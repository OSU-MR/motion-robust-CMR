classdef tPCA < handle
    %TPCA Applies the forward and hermitian temporal principle component
    %decomposition
    %   Detailed explanation goes here
    
    properties
        V;      % 
        sizes;  % Size of the input image
        M;      % Output vector size
    end
    
    methods
        % Constuctor
        function obj = tPCA(sizes)
            obj.V = [];
            obj.sizes = sizes;
            obj.M = prod(sizes);
        end
        
        % Forward Multiply
        function y = mult(obj,x)
            x = reshape(x,[obj.sizes(1)*obj.sizes(2),obj.sizes(3)]);
%             [~,~,obj.V] = svd(x);
            covariance=cov(x);
            [obj.V,~] = eig(covariance);
            obj.V = obj.V(:,size(obj.V,1):-1:1);
            y = x*obj.V;
            y = y(:);
        end
        
        % Backword Multiply
        function y = multTr(obj,x)
            x = reshape(x,[obj.sizes(1)*obj.sizes(2),obj.sizes(3)]);
            y = x*(obj.V)';
            y = y(:);
        end
        
        % Forward Multiply with 
        function y = multSq(obj,x)
            x = reshape(x,[obj.sizes(1)*obj.sizes(2),obj.sizes(3)]);
            if isempty(obj.V)
                x = reshape(x,[obj.sizes(1)*obj.sizes(2),obj.sizes(3)]);
                covariance=cov(x);
                [obj.V,~] = eig(covariance);
                obj.V = obj.V(:,size(obj.V,1):-1:1);
            end
            y = x*abs(obj.V).^2;
            y = y(:);
        end
        
        % Backword Multiply
        function y = multSqTr(obj,x)
            x = reshape(x,[obj.sizes(1)*obj.sizes(2),obj.sizes(3)]);
            if isempty(obj.V)
                x = reshape(x,[obj.sizes(1)*obj.sizes(2),obj.sizes(3)]);
                covariance=cov(x);
                [obj.V,~] = eig(covariance);
                obj.V = obj.V(:,size(obj.V,1):-1:1);
            end
            y = x*abs(obj.V').^2;
            y = y(:);
        end
    end
    
end

