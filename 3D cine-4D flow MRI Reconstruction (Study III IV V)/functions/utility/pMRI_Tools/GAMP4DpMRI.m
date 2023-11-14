classdef GAMP4DpMRI < hgsetget
    %GAMP4DPMRI2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        opts;
        xhatGAMP;
    end
    
    methods
        function obj = GAMP4DpMRI()
            obj.opts = GAMP4DpMRIOpts();
            obj.xhatGAMP = [];
        end
        
        % Method to crop data in the specified dimension "dim"
        function cropData(obj,cropVal,dim)
            
            if dim ==1
                obj.opts.data = ifft3_shift(obj.opts.data);
                obj.opts.data = obj.opts.data(cropVal+1:end-cropVal,:,:,:,:);
                obj.opts.data = fft3_shift(obj.opts.data);
                obj.opts.samp = obj.opts.samp(cropVal+1:end-cropVal,:,:,:);
                if ~isempty(obj.opts.maps)
                    obj.opts.maps = obj.opts.maps(cropVal+1:end-cropVal,:,:,:);
                end
                
            elseif dim ==2
                obj.opts.data = ifft3_shift(obj.opts.data);
                obj.opts.data = obj.opts.data(:,cropVal+1:end-cropVal,:,:,:);
                obj.opts.data = fft3_shift(obj.opts.data);
                obj.opts.samp = obj.opts.samp(:,cropVal+1:end-cropVal,:,:);
                if ~isempty(obj.opts.maps)
                    obj.opts.maps = obj.opts.maps(:,cropVal+1:end-cropVal,:,:);
                end
            else
                error('Use dim = 1 or 2');
            end
        end
        
        % Method that estimates sensitivity maps and combines coils
        function estimSensMaps(obj)
            
            % Coil combine
            obj.opts.data = coilCombine(obj.opts.data,13,'3dt');
            
            % Create time averaged image
            avg_image = sum(obj.opts.data,5);
            avg_pattern = sum(obj.opts.samp,4);
            avg_pattern(avg_pattern==0) = inf;
            avg_image = bsxfun(@rdivide,avg_image,avg_pattern);
%             avg_image = fftshift(fftshift(fftshift(avg_image,1),2),3);
%             avg_image = sqrt(size(avg_image,1)*size(avg_image,2))*ifft(ifft(ifft(avg_image,[],1),[],2),[],3);
%             avg_image = ifftshift(ifftshift(ifftshift(avg_image,1),2),3);
            avg_image = ifft3_shift(avg_image);
            sos = sqrt(sum(abs(avg_image).^2,4));
            sos = repmat(sos,[1,1,1,size(obj.opts.samp,4)]);
            sos = fftshift(fftshift(fftshift(sos,1),2),3);
            obj.opts.GAMPopt.xhat0 = sos(:)+1j*eps;

            
            p.mthd = 3; % '1' for espirit, '2' time-average spirit, '3' for walsh
            p.reEst = 0; % Res-estimating sensitivities
            p.fil = 3;
            [obj.opts.maps,~] = WalshCoilCombine3D(avg_image,p);

        end
        
        % Use single precision 
        function useSingle(obj)
            obj.opts.data = single(obj.opts.data);
            obj.opts.samp = single(obj.opts.samp);
            obj.opts.maps = single(obj.opts.maps);
            obj.opts.precision = 'single';
        end
        
        % Use GPU reconstruction
        function useGPU(obj)
            obj.opts.data = gpuArray(obj.opts.data);
            obj.opts.samp = gpuArray(obj.opts.samp);
            obj.opts.maps = gpuArray(obj.opts.maps);
            obj.opts.compute = 'gpu';
        end
        
        % Reconstruct Data
        function recon(obj)
            %GAMP4DPMRI Reconstructs 4D pMRI data using GAMP and 4D wavelets

            % fftshifts
            obj.opts.data = fftshift(fftshift(fftshift(obj.opts.data,1),2),3);
            obj.opts.samp = fftshift(fftshift(fftshift(obj.opts.samp,1),2),3);
            obj.opts.maps = fftshift(fftshift(fftshift(obj.opts.maps,1),2),3);
                      
            % Downsample the opts.data
            obj.opts.data = downsample_data(obj.opts.data,obj.opts.samp);

            % Normalize The columns of A to be unit norm
            R = numel(obj.opts.samp)/length(find(obj.opts.samp ==1));
            obj.opts.maps = obj.opts.maps*sqrt(R);
            obj.opts.data = obj.opts.data*sqrt(R);
            obj.opts.wvar = obj.opts.wvar*R;
            fprintf(sprintf('\nR = %0.3f\n',R))

            % Create Opterators
            pMRIOp = pMRI_Op_3D_t(obj.opts.maps,obj.opts.samp,'uniform_var',obj.opts.uniformVar,...
                     'compute',obj.opts.compute,'precision',obj.opts.precision);
            EstimIn = NullEstimIn(0,1);
            SparseTrans.wname = {'db1','db1','db1','db1'};
            SparseTrans.level = 1;
            W = ndDWTLinTrans4D(SparseTrans,size(obj.opts.samp),'uniform_var',...
                obj.opts.uniformVar,'compute',obj.opts.compute,...
                'precision',obj.opts.precision);

            % Concatonate All The Linear Transform Operator Together
            Op = LinTransConcat({pMRIOp;W},[1,1],obj.opts.precision,obj.opts.compute);  

            % Outputs esimators
            MeasEstimOut = CAwgnEstimOut(obj.opts.data(:),obj.opts.wvar,1);
            obj.opts.lambda = setLambda(size(obj.opts.samp),obj.opts.lambda);
            AnaEstimOut1 = CplxLaplaceEstimOut(obj.opts.lambda);
            EstimOut = EstimOutConcat({MeasEstimOut;AnaEstimOut1},[pMRIOp.M,W.M],...
                obj.opts.precision,obj.opts.compute);
            
            % Clean up memory
            clear MeasEstimOut AnaEstimOut1 pMRIOp avg_image avg_pattern sos
            
            % Reconstruct
            tic;
            obj.xhatGAMP = gampEst(EstimIn,EstimOut,Op,obj.opts.GAMPopt);
            display(sprintf('Reconstruction Completed in %s',num2str(toc)))
            obj.xhatGAMP = reshape(obj.xhatGAMP,size(obj.opts.samp));
            obj.xhatGAMP = ifftshift(ifftshift(ifftshift(obj.xhatGAMP,1),2),3);

            
            % fftshifts
            obj.opts.data = ifftshift(ifftshift(ifftshift(obj.opts.data,1),2),3);
            obj.opts.samp = ifftshift(ifftshift(ifftshift(obj.opts.samp,1),2),3);
            obj.opts.maps = ifftshift(ifftshift(ifftshift(obj.opts.maps,1),2),3);

        end
    end
    
end

