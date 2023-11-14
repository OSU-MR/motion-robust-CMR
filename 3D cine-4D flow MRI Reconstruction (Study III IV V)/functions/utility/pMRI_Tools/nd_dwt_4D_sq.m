%Calculates the nondecimated Wavelet transform using fast convlolution and
%   periodic boundry conditions.  The Filters are pre-computed and stored
%   in f_dec and f_rec to reduce computation time when performing Multiple
%   forward/backward transforms
%
%Methods:
%   ud_dwt_4D:  Constructor
%               Inputs: wname - Wavelet Filters to Use i.e. db1,db2,etc.
%                        Either a string 'db1' or a cell  of length 4 e.g. 
%                        {'db4','db1','db8,'db3'} where the first element in 
%                        the cell is filter for the first dimension, the 
%                        second element the seconde dimension, etc
%
%                       sizes - size of the 4D object [n1,n2,n3,n4]
%
%               Optional Inputs:
%                        pres_l2_norm -If set TRUE, the l2 norm in the 
%                        wavelet domain will be equal to the l2 norm in the
%                        signal domain.  Default is FALSE
% 
%                        compute - A string that sets the method of
%                        computation.  'mat' computes the wavelet tranform
%                        in Matlab.  'mex' computes the wavelet transfrom
%                        using a c mex file.  'gpu' computes the wavelets
%                        using matlab on a the GPU. When using 'gpu' the 
%                        input is expected to be a gpuArray and the output 
%                        also a gpuArray.  'gpu_off' computes the wavelets 
%                        using the gpu in matlab.  The calculation is offloaded
%                        to the GPU for calculation and brought back to 
%                        system memory afterwords. The input and output in 
%                        this case are both 'double' or 'single' depending
%                        on the precision used (see below)
%
%                        precision - a string that specifies double or
%                        single precision computation.  'double' or
%                        'single' are valid inputs
%
%   dec:        Multilevel Decomposition
%               Inputs: x - Image domain signal for decomposition
%
%                       levels - Number of decomposition Levels
%
%               Outputs: y - Multilevel non-decimated DWT coefficients in a
%                        5D array where the data is arranged [n1,n2,n3,n4,bands]
%                        The bands are orginized as follows.  Let "s1s2s3s4"
%                        represent the Bands.  The coefficients are ordered as follows,
%                        "LLLL", "HLLL", "LHLL", "HHLL", "LLHL",...,"HHHH
%                        where "H" denotes the high frequency filter and "L" 
%                        represents the low frequency filter. Successive 
%                        levels of decomposition are stacked such that the 
%                        highest "LLLL" is in [n1,n2,n3,n4,1]
%
%   rec:        Multilevel Reconstruction
%               Inputs: x - Wavelet coefficients in a 5D array size 
%                       [n1,n2,n3,n4,bands].
%
%               Outputs: y - Reconstructed 4D array.
%
%**************************************************************************
% The Ohio State University
% Written by:   Adam Rich 
% Last update:  7/5/2015
%**************************************************************************

classdef nd_dwt_4D_sq
    properties
        f_dec;          % Decomposition Filters
        sizes;          % Size of the 3D Image
        f_size;         % Length of the filters
        wname;          % Wavelet used
        pres_l2_norm;   % Binary indicator to preserver l2 norm of coefficients
        compute;        % How to compute the wavelet transform
        precision;      % which precision to use
    end
    
    methods
        % Constructor Object
        function obj = nd_dwt_4D_sq(wname,sizes,varargin)
            % Set Image size
            if length(sizes) ~=4
                error('The sizes vector must be length 4');
            else
                obj.sizes = sizes;
            end
            
            if ischar(wname)
                obj.wname = {wname,wname,wname,wname};
            elseif iscell(wname)
                obj.wname = wname;
            end
            
            % Set Default Options
            obj.pres_l2_norm = 0;
            obj.precision = 'double';
            obj.compute = 'mat';
            
            % Copy any optional inputs
            if mod(length(varargin),2)
                error('Optional inputs must come in pairs')
            end
            for ind = 1:2:length(varargin)
                switch lower(varargin{ind})
                    case 'pres_l2_norm'
                        obj.pres_l2_norm = varargin{ind+1};
                    case 'compute'
                        obj.compute = varargin{ind+1};
                    case 'precision'
                        obj.precision = varargin{ind+1};
                    otherwise
                        warning(sprintf('Unknown optional input #%d ingoring!',ind))
                end
            end
            
            if strcmpi(obj.compute,'mex') && strcmpi(obj.precision,'single')
                error('Single precsision is not currently supported for mex computation');
            end 
            
            % Get the Filter Coefficients
            [obj.f_dec,obj.f_size] = obj.get_filters(obj.wname);
            
            % Typecast the filters as single for single precision
            if strcmpi(obj.precision,'single')
                obj.f_dec = single(obj.f_dec);
            end
            
            % Put the filters on the GPU if using GPU computing
            if strcmpi(obj.compute,'gpu')
                obj.f_dec = gpuArray(obj.f_dec);
            end
        end
        
        % Multilevel Undecimated Wavelet Decomposition
        function y = dec(obj,x,level)
            
            % put the input array on the gpu
            if strcmpi(obj.compute,'gpu_off')
                x = gpuArray(x);
            end
            
            % Check if input is real
            if isreal(x)
                x_real = 1;
            else
                x_real = 0;
            end
            
            % Fourier Transform of Signal
            x = fftn(x);
            
			% Use Mex computation
            if strcmpi(obj.compute,'mex')
                y = nd_dwt_mex(x,obj.f_dec,0,level,obj.pres_l2_norm);
            
            % Use Matlab 
            else
                % Preallocate
                if ( strcmpi(obj.compute,'gpu_off') || strcmpi(obj.compute,'gpu') ) && strcmpi(obj.precision,'single')
                    y = zeros([obj.sizes, 16+15*(level-1)],'single');
                    y = gpuArray(y);
                elseif ( strcmpi(obj.compute,'gpu_off') || strcmpi(obj.compute,'gpu') ) && ~strcmpi(obj.precision,'single')
                    y = zeros([obj.sizes, 16+15*(level-1)],'gpuArray');
                elseif strcmpi(obj.precision,'single')
                    y = zeros([obj.sizes, 16+15*(level-1)],'single');
                else
                    y = zeros([obj.sizes, 16+15*(level-1)]);
                end
                
                % Calculate Mutlilevel Wavelet decomposition
                for ind = 1:level
                    % First Level
                    if ind ==1
                        y = level_1_dec(obj,x);
                    % Succssive Levels
                    else
                        y = cat(5,level_1_dec(obj,fftn(squeeze(y(:,:,:,:,1)))), y(:,:,:,:,2:end));
                    end
                end
            end
            
            % Take the real part if the input was real
            if x_real
                y = real(y);
            end
            
            % put arrays back in cpu memory
            if strcmpi(obj.compute,'gpu_off')
                y = gather(y);
            end
        end
        
        % Multilevel Undecimated Wavelet Reconstruction
        function y = rec(obj,x)
            
            % put the input array on the gpu
            if strcmpi(obj.compute,'gpu_off')
                x = gpuArray(x);
            end

            % Check if input is real
            if isreal(x)
                x_real = 1;
            else
                x_real = 0;
            end
            
            % Find the decomposition level
            level = 1+(size(x,5)-16)/15;
            
            % Fourier Transform of Signal
            x = fft(fft(fft(fft(x,[],1),[],2),[],3),[],4);
			
            if strcmpi(obj.compute,'mex')
                y = nd_dwt_mex(x,obj.f_dec,1,level,obj.pres_l2_norm);
            else
                % Reconstruct from Multiple Levels
                for ind = 1:level
                    % First Level
                    if ind ==1
                        y = level_1_rec(obj,x);
                        if ~obj.pres_l2_norm
                            y = y/16;
                        end
                    % Succssive Levels
                    else
                        y = fftn(y);
                        y = level_1_rec(obj,cat(5,y,x(:,:,:,:,17+(ind-2)*15:31+(ind-2)*15)));
                        if ~obj.pres_l2_norm
                            y = y/16;
                        end
                    end
                end 
            end
            
            % Take the real part if the input was real
            if x_real
                y = real(y);
            end
                        % put arrays back in cpu memory
            if strcmpi(obj.compute,'gpu_off')
                y = gather(y);
            end
        end
    end
    
     %% Private Methods
    methods (Access = protected,Hidden = true)
        
        % Returns the Filters and 
        function [f_dec,f_size] = get_filters(obj,wname)
            % Get Filters for the first domain
            [LO_D,HI_D]   = wave_filters(wname{1});
            [LO_D2,HI_D2] = wave_filters(wname{2});
            [LO_D3,HI_D3] = wave_filters(wname{3});
            [LO_D4,HI_D4] = wave_filters(wname{4});
            
            % Find the filter size
            f_size.s1 = length(LO_D);
            f_size.s2 = length(LO_D2);
            f_size.s3 = length(LO_D3);
            f_size.s4 = length(LO_D4);
            
            % Dimension Check
            if f_size.s1 > obj.sizes(1)
                error(['First Dimension of Data is shorter than the wavelet'...
                    ,' filter being used']);
            elseif f_size.s2 > obj.sizes(2)
                error(['Second Dimension of Data is shorter than the wavelet'...
                    ,' filter being used']);
            elseif f_size.s3 > obj.sizes(3)
                error(['Third Dimension of Data is shorter than the wavelet'...
                    ,' filter being used']);
            elseif f_size.s4 > obj.sizes(4)
                error(['Fourth Dimension of Data is shorter than the wavelet'...
                    ,' filter being used']);
            end
            
            % Get the 2D Filters by taking outer products
            dec_LL = LO_D.'*LO_D2;
            dec_HL = HI_D.'*LO_D2;
            dec_LH = LO_D.'*HI_D2;
            dec_HH = HI_D.'*HI_D2;
            
            % Take the Outerproducts for the third dimension
            dec_LLL = zeros(length(LO_D),length(LO_D2),length(LO_D3));
            dec_HLL = zeros(length(LO_D),length(LO_D2),length(LO_D3));
            dec_LHL = zeros(length(LO_D),length(LO_D2),length(LO_D3));
            dec_HHL = zeros(length(LO_D),length(LO_D2),length(LO_D3));
            dec_LLH = zeros(length(LO_D),length(LO_D2),length(LO_D3));
            dec_HLH = zeros(length(LO_D),length(LO_D2),length(LO_D3));
            dec_LHH = zeros(length(LO_D),length(LO_D2),length(LO_D3));
            dec_HHH = zeros(length(LO_D),length(LO_D2),length(LO_D3));
            
            for ind = 1:size(dec_LL,2)
                dec_LLL(:,ind,:) = dec_LL(:,ind)*LO_D3;
                dec_HLL(:,ind,:) = dec_HL(:,ind)*LO_D3;
                dec_LHL(:,ind,:) = dec_LH(:,ind)*LO_D3;
                dec_HHL(:,ind,:) = dec_HH(:,ind)*LO_D3;
                dec_LLH(:,ind,:) = dec_LL(:,ind)*HI_D3;
                dec_HLH(:,ind,:) = dec_HL(:,ind)*HI_D3;
                dec_LHH(:,ind,:) = dec_LH(:,ind)*HI_D3;
                dec_HHH(:,ind,:) = dec_HH(:,ind)*HI_D3;
            end

            % Take the Outerproducts for the fourth dimension
            for kk = 1:size(dec_LLL,3)
                for ii = 1:size(dec_LLL,2)
                    f_dec1(:,ii,kk,:,1) = dec_LLL(:,ii,kk)*LO_D4;
                    f_dec1(:,ii,kk,:,2) = dec_HLL(:,ii,kk)*LO_D4;
                    f_dec1(:,ii,kk,:,3) = dec_LHL(:,ii,kk)*LO_D4;
                    f_dec1(:,ii,kk,:,4) = dec_HHL(:,ii,kk)*LO_D4;
                    f_dec1(:,ii,kk,:,5) = dec_LLH(:,ii,kk)*LO_D4;
                    f_dec1(:,ii,kk,:,6) = dec_HLH(:,ii,kk)*LO_D4;
                    f_dec1(:,ii,kk,:,7) = dec_LHH(:,ii,kk)*LO_D4;
                    f_dec1(:,ii,kk,:,8) = dec_HHH(:,ii,kk)*LO_D4;
                    f_dec1(:,ii,kk,:,9) = dec_LLL(:,ii,kk)*HI_D4;
                    f_dec1(:,ii,kk,:,10) = dec_HLL(:,ii,kk)*HI_D4;
                    f_dec1(:,ii,kk,:,11) = dec_LHL(:,ii,kk)*HI_D4;
                    f_dec1(:,ii,kk,:,12) = dec_HHL(:,ii,kk)*HI_D4;
                    f_dec1(:,ii,kk,:,13) = dec_LLH(:,ii,kk)*HI_D4;
                    f_dec1(:,ii,kk,:,14) = dec_HLH(:,ii,kk)*HI_D4;
                    f_dec1(:,ii,kk,:,15) = dec_LHH(:,ii,kk)*HI_D4;
                    f_dec1(:,ii,kk,:,16) = dec_HHH(:,ii,kk)*HI_D4;
                end
            end

            % Add a circularshift of half the filter length to the 
            % reconstruction filters by adding phase to them
            phase1 = exp(1j*2*pi*f_size.s1/2*linspace(0,1-1/obj.sizes(1),obj.sizes(1)));
            phase2 = exp(1j*2*pi*f_size.s2/2*linspace(0,1-1/obj.sizes(2),obj.sizes(2)));
            phase3 = exp(1j*2*pi*f_size.s3/2*linspace(0,1-1/obj.sizes(3),obj.sizes(3)));
            phase4 = exp(1j*2*pi*f_size.s4/2*linspace(0,1-1/obj.sizes(4),obj.sizes(4)));

            % 2D Phase
            shift_2 = phase1.'*phase2;
            
            % 3D Phase
            shift_3 = zeros([obj.sizes(1),obj.sizes(2),obj.sizes(3)]);
            for ind = 1:size(shift_2,2)
               shift_3(:,ind,:) = shift_2(:,ind)*phase3; 
            end
            
            % 4D Phase
            shift = zeros([obj.sizes(1),obj.sizes(2),obj.sizes(3),obj.sizes(4)]);
            for ii = 1:size(shift_3,2)
                for kk = 1:size(shift_3,3)
                    shift(:,ii,kk,:) = shift_3(:,ii,kk)*phase4;
                end
            end
            
            % square the filter coefficients
            f_dec1 = f_dec1.^2;
            
            % Take the Fourier Transform of the Kernels for Fast
            % Convolution
            if obj.pres_l2_norm 
                scale = 1/16;
            else
                scale = 1;
            end
            if strcmpi(obj.compute,'mex')
                scale2 = 1/prod(obj.sizes);
            else
                scale2 = 1;
            end
            f_dec(:,:,:,:,1) = scale2*scale*shift.*fftn(f_dec1(:,:,:,:,1),[obj.sizes(1),obj.sizes(2),obj.sizes(3),obj.sizes(4)]);
            f_dec(:,:,:,:,2) = scale2*scale*shift.*fftn(f_dec1(:,:,:,:,2),[obj.sizes(1),obj.sizes(2),obj.sizes(3),obj.sizes(4)]);
            f_dec(:,:,:,:,3) = scale2*scale*shift.*fftn(f_dec1(:,:,:,:,3),[obj.sizes(1),obj.sizes(2),obj.sizes(3),obj.sizes(4)]);
            f_dec(:,:,:,:,4) = scale2*scale*shift.*fftn(f_dec1(:,:,:,:,4),[obj.sizes(1),obj.sizes(2),obj.sizes(3),obj.sizes(4)]);
            f_dec(:,:,:,:,5) = scale2*scale*shift.*fftn(f_dec1(:,:,:,:,5),[obj.sizes(1),obj.sizes(2),obj.sizes(3),obj.sizes(4)]);
            f_dec(:,:,:,:,6) = scale2*scale*shift.*fftn(f_dec1(:,:,:,:,6),[obj.sizes(1),obj.sizes(2),obj.sizes(3),obj.sizes(4)]);
            f_dec(:,:,:,:,7) = scale2*scale*shift.*fftn(f_dec1(:,:,:,:,7),[obj.sizes(1),obj.sizes(2),obj.sizes(3),obj.sizes(4)]);
            f_dec(:,:,:,:,8) = scale2*scale*shift.*fftn(f_dec1(:,:,:,:,8),[obj.sizes(1),obj.sizes(2),obj.sizes(3),obj.sizes(4)]);
            f_dec(:,:,:,:,9) = scale2*scale*shift.*fftn(f_dec1(:,:,:,:,9),[obj.sizes(1),obj.sizes(2),obj.sizes(3),obj.sizes(4)]);
            f_dec(:,:,:,:,10) = scale2*scale*shift.*fftn(f_dec1(:,:,:,:,10),[obj.sizes(1),obj.sizes(2),obj.sizes(3),obj.sizes(4)]);
            f_dec(:,:,:,:,11) = scale2*scale*shift.*fftn(f_dec1(:,:,:,:,11),[obj.sizes(1),obj.sizes(2),obj.sizes(3),obj.sizes(4)]);
            f_dec(:,:,:,:,12) = scale2*scale*shift.*fftn(f_dec1(:,:,:,:,12),[obj.sizes(1),obj.sizes(2),obj.sizes(3),obj.sizes(4)]);
            f_dec(:,:,:,:,13) = scale2*scale*shift.*fftn(f_dec1(:,:,:,:,13),[obj.sizes(1),obj.sizes(2),obj.sizes(3),obj.sizes(4)]);
            f_dec(:,:,:,:,14) = scale2*scale*shift.*fftn(f_dec1(:,:,:,:,14),[obj.sizes(1),obj.sizes(2),obj.sizes(3),obj.sizes(4)]);
            f_dec(:,:,:,:,15) = scale2*scale*shift.*fftn(f_dec1(:,:,:,:,15),[obj.sizes(1),obj.sizes(2),obj.sizes(3),obj.sizes(4)]);
            f_dec(:,:,:,:,16) = scale2*scale*shift.*fftn(f_dec1(:,:,:,:,16),[obj.sizes(1),obj.sizes(2),obj.sizes(3),obj.sizes(4)]);
        end
        
        % Single Level Redundant Wavelet Decomposition
        function y = level_1_dec(obj,x_f)
            % Preallocate
            if ( strcmpi(obj.compute,'gpu_off') || strcmpi(obj.compute,'gpu') ) && strcmpi(obj.precision,'single')
                y = zeros([obj.sizes,16],'single');
                y = gpuArray(y);
            elseif ( strcmpi(obj.compute,'gpu_off') || strcmpi(obj.compute,'gpu') ) && ~strcmpi(obj.precision,'single')
                y = zeros([obj.sizes,16],'gpuArray');
            elseif strcmpi(obj.precision,'single')
                y = zeros([obj.sizes,16],'single');
            else
                y = zeros([obj.sizes,16]);
            end
            
            % Calculate Wavelet Coefficents Using Fast Convolution
            y(:,:,:,:,1) =  ifftn(x_f.*obj.f_dec(:,:,:,:,1));
            y(:,:,:,:,2) =  ifftn(x_f.*obj.f_dec(:,:,:,:,2));
            y(:,:,:,:,3) =  ifftn(x_f.*obj.f_dec(:,:,:,:,3));
            y(:,:,:,:,4) =  ifftn(x_f.*obj.f_dec(:,:,:,:,4));
            y(:,:,:,:,5) =  ifftn(x_f.*obj.f_dec(:,:,:,:,5));
            y(:,:,:,:,6) =  ifftn(x_f.*obj.f_dec(:,:,:,:,6));
            y(:,:,:,:,7) =  ifftn(x_f.*obj.f_dec(:,:,:,:,7));
            y(:,:,:,:,8) =  ifftn(x_f.*obj.f_dec(:,:,:,:,8));
            y(:,:,:,:,9) =  ifftn(x_f.*obj.f_dec(:,:,:,:,9));
            y(:,:,:,:,10) = ifftn(x_f.*obj.f_dec(:,:,:,:,10));
            y(:,:,:,:,11) = ifftn(x_f.*obj.f_dec(:,:,:,:,11));
            y(:,:,:,:,12) = ifftn(x_f.*obj.f_dec(:,:,:,:,12));
            y(:,:,:,:,13) = ifftn(x_f.*obj.f_dec(:,:,:,:,13));
            y(:,:,:,:,14) = ifftn(x_f.*obj.f_dec(:,:,:,:,14));
            y(:,:,:,:,15) = ifftn(x_f.*obj.f_dec(:,:,:,:,15));
            y(:,:,:,:,16) = ifftn(x_f.*obj.f_dec(:,:,:,:,16));
            
        end
        
        % Single Level Redundant Wavelet Reconstruction
        function y = level_1_rec(obj,x_f)
            
            % Reconstruct the 3D array using Fast Convolution
            y = ifftn(squeeze(x_f(:,:,:,:,1)).*conj(obj.f_dec(:,:,:,:,1)));
            y = y + ifftn(squeeze(x_f(:,:,:,:,2)).*conj(obj.f_dec(:,:,:,:,2)));
            y = y + ifftn(squeeze(x_f(:,:,:,:,3)).*conj(obj.f_dec(:,:,:,:,3)));
            y = y + ifftn(squeeze(x_f(:,:,:,:,4)).*conj(obj.f_dec(:,:,:,:,4)));
            y = y + ifftn(squeeze(x_f(:,:,:,:,5)).*conj(obj.f_dec(:,:,:,:,5)));
            y = y + ifftn(squeeze(x_f(:,:,:,:,6)).*conj(obj.f_dec(:,:,:,:,6)));
            y = y + ifftn(squeeze(x_f(:,:,:,:,7)).*conj(obj.f_dec(:,:,:,:,7)));
            y = y + ifftn(squeeze(x_f(:,:,:,:,8)).*conj(obj.f_dec(:,:,:,:,8)));
            y = y + ifftn(squeeze(x_f(:,:,:,:,9)).*conj(obj.f_dec(:,:,:,:,9)));
            y = y + ifftn(squeeze(x_f(:,:,:,:,10)).*conj(obj.f_dec(:,:,:,:,10)));
            y = y + ifftn(squeeze(x_f(:,:,:,:,11)).*conj(obj.f_dec(:,:,:,:,11)));
            y = y + ifftn(squeeze(x_f(:,:,:,:,12)).*conj(obj.f_dec(:,:,:,:,12)));
            y = y + ifftn(squeeze(x_f(:,:,:,:,13)).*conj(obj.f_dec(:,:,:,:,13)));
            y = y + ifftn(squeeze(x_f(:,:,:,:,14)).*conj(obj.f_dec(:,:,:,:,14)));
            y = y + ifftn(squeeze(x_f(:,:,:,:,15)).*conj(obj.f_dec(:,:,:,:,15)));
            y = y + ifftn(squeeze(x_f(:,:,:,:,16)).*conj(obj.f_dec(:,:,:,:,16)));

        end
    end
end
