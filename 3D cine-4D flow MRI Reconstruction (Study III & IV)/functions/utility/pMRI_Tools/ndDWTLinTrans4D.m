classdef ndDWTLinTrans4D < LinTrans
    %NDDWTLINTRANS Linear transform class where the operator is a
    %non-decimated wavelet transform
    %   Detailed explanation goes here
    
%**************************************************************************
% The Ohio State University
% Written by:   Adam Rich 
% Email:        rich.178@osu.edu
% Last update:  4/13/2014
%**************************************************************************
    
% Public Properties
properties
    w_prop;
    sizes;
    Fro2;
    N;
    M;
    uniform_var;
    compute;
    precision;
end

% Private Properties
properties (Access = protected,Hidden = true)
    W;
    W_sq
end

methods
    % Constructor
    function obj = ndDWTLinTrans4D(w_prop,sizes,varargin)
        obj = obj@LinTrans(prod(sizes),prod(sizes)*(15*(w_prop.level-1) + 16));
        
        % Copy any optional inputs 
        if mod(length(varargin),2)
            error('Optional inputs must come in pairs')
        end

        % Set Default Options
        obj.precision = 'double';
        obj.compute = 'mat';
        obj.uniform_var = 0;
        
        % Set options based on inputs
        for ind = 1:2:length(varargin)
            switch lower(varargin{ind})
                case 'compute'
                    obj.compute = varargin{ind+1};
                case 'precision'
                    obj.precision = varargin{ind+1};
                case 'uniform_var'
                    obj.uniform_var = varargin{ind+1};
                otherwise
                    warning(sprintf('Unknown optional input #%d ignoring!',ind));
            end
        end
        
        % Set Class Properties
        obj.sizes = sizes;
        obj.w_prop = w_prop;
        obj.N = prod(sizes);
        obj.M = obj.N*(15*(w_prop.level-1) + 16);
        
        % Create nd-dwt operator
        obj.W = nd_dwt_4D(w_prop.wname,sizes,'pres_l2_norm',1,'compute',obj.compute,'precision',obj.precision);
        if obj.uniform_var ==0
            obj.W_sq = nd_dwt_4D_sq(w_prop.wname,sizes,'pres_l2_norm',1,'compute',obj.compute,'precision',obj.precision);
        end
        % Calculate the frobenious norm
        obj.Fro2 = 1/obj.M;
    end

    % Return the Dimensions of the Matrix
    function [m,n] = size(obj)
        n = obj.N;	
        m = obj.M;
    end

    % Multiply with W
    function y = mult(obj,x)
        x = reshape(x,obj.sizes);
        y = obj.W.dec(x,obj.w_prop.level);
% 		y = squeeze(permute(y,[2,3,4,1]));
        y = y(:);
    end

    % Multiply with W^T
    function y = multTr(obj,x)
		x = reshape(x,[obj.sizes(1),obj.sizes(2),obj.sizes(3),obj.sizes(4),16+15*(obj.w_prop.level-1)]);
% 		x = permute(x,[4,1,2,3]);
        x = obj.W.rec(x);
        y = x(:);
    end

    % Squared-Matrix multiply 
    function y = multSq(obj,x)
        if obj.uniform_var ==1
            y = ones(obj.M,1)*(obj.Fro2*sum(x(:),1));
        else
            x = reshape(x,obj.sizes);
            y = obj.W_sq.dec(x,obj.w_prop.level);
            y = y(:);
        end
    end

    % Squared-Hermitian-Transposed Matrix multiply 
    function y = multSqTr(obj,x)
        if obj.uniform_var ==1
            y = ones(obj.N,1)*(obj.Fro2*sum(x(:),1));
        else
            x = reshape(x,[obj.sizes(1),obj.sizes(2),obj.sizes(3),obj.sizes(4),16+15*(obj.w_prop.level-1)]);
            y = obj.W_sq.rec(x);
            y = y(:);
        end
    end
end
    
end

