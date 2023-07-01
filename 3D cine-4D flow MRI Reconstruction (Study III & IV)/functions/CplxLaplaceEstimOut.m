classdef CplxLaplaceEstimOut < EstimOut
    %CPLXLAPLACEESTIMOUT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        lambda;
    end
    
    methods
        % Constructor
        function obj = CplxLaplaceEstimOut(lambda)
            obj.lambda = lambda;
        end
        
        % Estim Function
        function [zmean, zvar] = estim(obj, zmean0, zvar0)
            
            if isa(zmean0,'gpuArray')
                % Calculate Threshold
            
                % Calculate mean
                thresh_t = (abs(zmean0)-zvar0.*obj.lambda) > 0;
                zmean = arrayfun(@(a,b,c,d)(a.*(abs(b)-c.*d).*sign(b)),thresh_t,zmean0,zvar0,obj.lambda);
%                 zmean = thresh_t.*(abs(zmean0)-zvar0.*obj.lambda).*sign(zmean0);
%                 thresh_t.*(abs(zmean0)-zvar0.*obj.lambda).*sign(zmean0);
                clear zmean0
            
                % Caluclate Variance
                zvar = zvar0.*thresh_t;
                
            else
                % Calculate Threshold
                thresh = zvar0.*obj.lambda;
                
                % Calculate mean
                z_thresh = (abs(zmean0)-thresh);
                thresh_t = (abs(zmean0)-thresh) > 0;
                zmean = thresh_t.*z_thresh.*sign(zmean0);
                
                % Caluclate Variance
                zvar = zvar0.*thresh_t;
            end
                
        end
        
        % Log Likelihood Method
        function ll = logLike(obj,zhat,zvar)
%             ll = -abs(obj.lambda.*zhat) + log(obj.lambda.^2/(2*pi));
            ll = -abs(obj.lambda.*zhat);
%             ll = arrayfun(@(a,b)(-abs(a.*b)),obj.lambda,zhat);
        end
        
        % Log Scale Method. Same as loglike for max sum
        function ll = logScale(obj,zhat,zvar,var,pha)
%             ll = -abs(obj.lambda.*zhat) + log(obj.lambda.^2/(2*pi));
            ll = -abs(obj.lambda.*zhat);
        end
    end
    
end

