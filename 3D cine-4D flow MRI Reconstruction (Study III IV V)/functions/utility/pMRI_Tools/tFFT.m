classdef tFFT< handle
    %TPCA Applies the forward and hermitian temporal principle component
    %decomposition
    %   Detailed explanation goes here
    
    properties
        sizes;  % Size of the input image
        M;      % Output vector size
    end
    
    methods
        % Constuctor
        function obj = tFFT(sizes)
            obj.sizes = sizes;
            obj.M = prod(sizes);
        end
        
        % Forward Multiply
        function y = mult(obj,x)
            x = reshape(x,[obj.sizes(1),obj.sizes(2),obj.sizes(3)]);
            y = 1/sqrt(obj.sizes(3))*fft(x,[],3);
            y = y(:);
        end
        
        % Backword Multiply
        function y = multTr(obj,x)
            x = reshape(x,[obj.sizes(1),obj.sizes(2),obj.sizes(3)]);
            y = sqrt(obj.sizes(3))*ifft(x,[],3);
            y = y(:);
        end
        
        % Forward squared Multiply
        function y = multSq(obj,x)
            y = ones(size(x(:)),class(x))*sum(x(:));
        end
        
        % Backword squared Multiply
        function y = multSqTr(obj,x)
            y = ones(size(x(:)),class(x))*sum(x(:));
        end
    end
    
end

