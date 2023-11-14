classdef sparseTransConcat < LinTrans
    %SPARSETRANSCONCAT Concatinate multiple sparsifying transforms together
    %   Detailed explanation goes here
    
    properties
        handles;        % A cell array of handles
        redu_factor;    % An array of redundancy factors
        im_sizes;       % The size of the input image
        N;
        M;
    end
    
    methods
        function obj = sparseTransConcat(handles,redu_factor,im_sizes)
            obj.handles = handles;
            obj.redu_factor = redu_factor;
            obj.im_sizes = im_sizes;
            obj.N = prod(im_sizes);
            obj.M = obj.N*sum(redu_factor);
        end
        
        % Return the Dimensions of the Matrix
        function [m,n] = size(obj)
            n = obj.N;	
            m = obj.M;
        end
    
        function y = mult(obj,x)
            y = [];
            for ind = 1:length(obj.handles)
                y = [y;obj.handles{ind}.mult(x)];
            end
            y = y/sqrt(sum(obj.redu_factor));
        end
        
        function y = multTr(obj,x)
            x = reshape(x,[obj.im_sizes,sum(obj.redu_factor)]);
            y = zeros(prod(obj.im_sizes),1,class(x));
            start_ind = 1;
            handle_ind = 1;
            for redu = [obj.redu_factor]
                y = y + redu*obj.handles{handle_ind}.multTr(x(:,:,:,start_ind:start_ind+redu-1));
                start_ind = start_ind+redu;
                handle_ind = handle_ind+1;
            end
%             y = y/sqrt(length(obj.redu_factor));
            y = y/sqrt(sum(obj.redu_factor));
        end
        
        function y = multSq(obj,x)
            y = [];
            for ind = 1:length(obj.handles)
                y = [y;obj.handles{ind}.multSq(x)];
            end
            y = y/(sum(obj.redu_factor));
        end
        
        function y = multSqTr(obj,x)
            x = reshape(x,[obj.im_sizes,sum(obj.redu_factor)]);
            y = zeros(prod(obj.im_sizes),1,class(x));
            start_ind = 1;
            handle_ind = 1;
            for redu = [obj.redu_factor]
                y = y + redu*obj.handles{handle_ind}.multSqTr(x(:,:,:,start_ind:start_ind+redu-1));
                start_ind = start_ind+redu;
                handle_ind = handle_ind+1;
            end
            y = y/(sum(obj.redu_factor));
        end
    end
    
end

