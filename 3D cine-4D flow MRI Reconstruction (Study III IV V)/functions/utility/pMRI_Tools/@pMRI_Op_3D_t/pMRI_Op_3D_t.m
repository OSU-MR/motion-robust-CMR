classdef pMRI_Op_3D_t < LinTrans
    %PMRI_OP_3D_T Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        frame_size;     % The Size of the image e.g. [256,256]
        N;              % Input Dimension Size
        M;              % Output Dimension Size
        Q;              % The Number of Frames
        C;              % The coil Sensitivity maps
        CSq;            % The squared sensitivity maps
        CSqTr;
        Fro2;           % Squared Forbenius Norm
        uniform_var;    % Binary indicator for uniform variance
        compute;        % How to compute the wavelet transform
        precision;      % which precision to use
        mask_patterns;
        mask_weights;   % Weights for weighted least squares
    end
    
    properties(Access = protected, Hidden = true)
        
        nddwt;
        s;              % 2D wavelet signiture matrix 
    end
    
    methods
        % Object Constructor
        function obj = pMRI_Op_3D_t(C,pattern,varargin)
                       
            % Set Default Options
            obj.precision = 'double';
            obj.compute = 'mat';
            obj.uniform_var = 0;
            obj.mask_weights = [];
            weights = [];
            
            % Copy any optional inputs 
            if mod(length(varargin),2)
                error('Optional inputs must come in pairs')
            end
            for ind = 1:2:length(varargin)
                switch lower(varargin{ind})
                    case 'uniform_var'
                        obj.uniform_var = varargin{ind+1};
                    case 'compute'
                        obj.compute = varargin{ind+1};
                    case 'precision'
                        obj.precision = varargin{ind+1};
                    case 'weights'
                        weights = varargin{ind+1};
                    otherwise
                        warning(sprintf('Unknown optional input #%d ingoring!',ind))
                end
            end
            
            obj.C = C;
%             obj.frame_size = size(squeeze(C(:,:,:,1,1)));
            obj.frame_size(1) = size(C,1);
            obj.frame_size(2) = size(C,2);
            obj.frame_size(3) = size(C,3);
            obj.Q = size(pattern,4);
            
            % Find the size of the measured Data
            patterns = permute(pattern,[1,2,3,5,4]);
            patterns = repmat(patterns,[1,1,1,size(C,4),1]);
            obj.mask_patterns = find(patterns==1);
            obj.M = length(obj.mask_patterns);
            
            % Find the size of the measured Data
            if ~isempty(weights)
                weights = permute(weights,[1,2,3,5,4]);
                weights = repmat(weights,[1,1,1,size(C,4),1]);
%                warning('look at this again')
                obj.mask_weights = weights(patterns==1);
            end
            
            % Find Image size
            obj.N = obj.frame_size(1)*obj.frame_size(2)*obj.frame_size(3)*obj.Q; 
            
            % Get Frobenius Norm
            obj.Fro2 = 1/(obj.N*size(obj.C,4));
            
            if ~obj.uniform_var
                obj.CSq = abs(obj.C.^2)/(obj.N/obj.Q);
                obj.CSqTr = obj.CSq(:,:,:,:,1);
                obj.CSq = obj.CSq(:,:,:,:,1);
                obj.CSq = reshape(obj.CSq,[size(obj.CSq,1)*size(obj.CSq,2)*size(obj.CSq,3),size(obj.CSq,4)]);
                obj.CSq = obj.CSq.';
            else
                obj.CSq = [];
            end
            
%             % Check to make sure some inputs are the right size
%             if sum(size(pattern(:,:,1)) ~= obj.frame_size) >=1
%                 error('Image size and sampling pattern must be the same size')
%             elseif size(C,1) ~= obj.frame_size(1) || size(C,2) ~= obj.frame_size(2)
%                 error('Image size and Coils must be the same size')
%             end
           
            if strcmpi(obj.precision,'single')
                obj.C = single(obj.C);
                obj.CSq = single(obj.CSq);
                obj.CSqTr = single(obj.CSqTr);
%                 obj.uniform_var = single(obj.uniform_var);
%                 obj.mask_patterns = single(obj.mask_patterns);
            end
            
            if strcmpi(obj.compute,'gpu_off') || strcmpi(obj.compute,'gpu')
                obj.C = gpuArray(obj.C);
                obj.CSq = gpuArray(obj.CSq);
                obj.CSqTr = gpuArray(obj.CSqTr);
                obj.uniform_var = gpuArray(obj.uniform_var);
                obj.mask_patterns = gpuArray(obj.mask_patterns);
                obj.mask_weights = gpuArray(obj.mask_weights);
%                 obj.M = gpuArray(obj.M );
%                 obj.N = gpuArray(obj.N ); 
%                 obj.frame_size = gpuArray(obj.frame_size);
%                 obj.Q = gpuArray(obj.Q);
            end
        end

        % Return the Dimensions of the Matrix
        function [m,n] = size(obj)
            n = obj.N;	
            m = obj.M;
        end
        
        % Forward 2D t Operator
        y = mult(obj,x)

        % Hermitian 2D t Operator
        y = multTr(obj,x)

        % Squared-Matrix multiply 
        y = multSq(obj,x)
        
        % Squared-Hermitian-Transposed Matrix multiply 
        y = multSqTr(obj,x)
    end
    
end





