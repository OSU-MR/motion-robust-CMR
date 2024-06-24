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

classdef harr_nddwt_4D
    properties
        sizes;          % Size of the 3D Image
        f_size;         % Length of the filters
        wname;          % Wavelet used
        pres_l2_norm;   % Binary indicator to preserver l2 norm of coefficients
        compute;        % How to compute the wavelet transform
        precision;      % which precision to use
        scale;
    end
    
    methods
        % Constructor Object
        function obj = harr_nddwt_4D(wname,sizes,varargin)
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
            
            if obj.pres_l2_norm
                obj.scale =1/2;
            else
                obj.scale = 1/sqrt(2);
            end
            
            if strcmpi(obj.compute,'mex') && strcmpi(obj.precision,'single')
                error('Single precsision is not currently supported for mex computation');
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
            
            % Calculate Wavelet Coefficents Using Convolution
            tmp = zeros(obj.sizes,'like',x_f);
            tmp1 = zeros(obj.sizes,'like',x_f);

            % =============================================================
            % app dim 1,2,3,4 
            % app dim1
            tmp(1:end-1,:,:,:) = obj.scale*(x_f(1:end-1,:,:,:)+x_f(2:end,:,:,:));
            tmp(end,:,:,:) = obj.scale*(x_f(1,:,:,:)+x_f(end,:,:,:));
  
            % app dim 2
            tmp1(:,1:end-1,:,:,:) = obj.scale*(tmp(:,1:end-1,:,:,:)+tmp(:,2:end,:,:,:));
            tmp1(:,end,:,:,:) = obj.scale*(tmp(:,1,:,:,:)+tmp(:,end,:,:,:));

            % app dim 3
            tmp(:,:,1:end-1,:) = obj.scale*(tmp1(:,:,1:end-1,:)+tmp1(:,:,2:end,:));
            tmp(:,:,end,:) = obj.scale*(tmp1(:,:,1,:)+tmp1(:,:,end,:));
            
            % app dim 4
            y(:,:,:,1:end-1,1) = obj.scale*(tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
            y(:,:,:,end,1) = obj.scale*(tmp(:,:,:,1)+tmp(:,:,:,end));  
            
            % =============================================================
            % det dim 1, app dim 2,3,4 
            % det dim1
            tmp(1:end-1,:,:,:) = obj.scale*(x_f(1:end-1,:,:,:)+-x_f(2:end,:,:,:));
            tmp(end,:,:,:) = obj.scale*(-x_f(1,:,:,:)+x_f(end,:,:,:));
  
            % app dim 2
            tmp1(:,1:end-1,:,:,:) = obj.scale*(tmp(:,1:end-1,:,:,:)+tmp(:,2:end,:,:,:));
            tmp1(:,end,:,:,:) = obj.scale*(tmp(:,1,:,:,:)+tmp(:,end,:,:,:));

            % app dim 3
            tmp(:,:,1:end-1,:) = obj.scale*(tmp1(:,:,1:end-1,:)+tmp1(:,:,2:end,:));
            tmp(:,:,end,:) = obj.scale*(tmp1(:,:,1,:)+tmp1(:,:,end,:));
            
            % app dim 4
            y(:,:,:,1:end-1,2) = obj.scale*(tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
            y(:,:,:,end,2) = obj.scale*(tmp(:,:,:,1)+tmp(:,:,:,end));
            
            % =============================================================
            % det dim 2, app dim 1,3,4 
            % app dim 1
            tmp(1:end-1,:,:,:) = obj.scale*(x_f(1:end-1,:,:,:)+x_f(2:end,:,:,:));
            tmp(end,:,:,:) = obj.scale*(x_f(1,:,:,:)+x_f(end,:,:,:));
  
            % det dim 2
            tmp1(:,1:end-1,:,:,:) = obj.scale*(tmp(:,1:end-1,:,:,:)-tmp(:,2:end,:,:,:));
            tmp1(:,end,:,:,:) = obj.scale*(-tmp(:,1,:,:,:)+tmp(:,end,:,:,:));

            % app dim 3
            tmp(:,:,1:end-1,:) = obj.scale*(tmp1(:,:,1:end-1,:)+tmp1(:,:,2:end,:));
            tmp(:,:,end,:) = obj.scale*(tmp1(:,:,1,:)+tmp1(:,:,end,:));
            
            % app dim 4
            y(:,:,:,1:end-1,3) = obj.scale*(tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
            y(:,:,:,end,3) = obj.scale*(tmp(:,:,:,1)+tmp(:,:,:,end));
            
            % =============================================================
            % det dim 1,2 app dim 3,4 
            % app dim1
            tmp(1:end-1,:,:,:) = obj.scale*(x_f(1:end-1,:,:,:)-x_f(2:end,:,:,:));
            tmp(end,:,:,:) = obj.scale*(-x_f(1,:,:,:)+x_f(end,:,:,:));
  
            % app dim 2
            tmp1(:,1:end-1,:,:,:) = obj.scale*(tmp(:,1:end-1,:,:,:)-tmp(:,2:end,:,:,:));
            tmp1(:,end,:,:,:) = obj.scale*(-tmp(:,1,:,:,:)+tmp(:,end,:,:,:));

            % app dim 3
            tmp(:,:,1:end-1,:) = obj.scale*(tmp1(:,:,1:end-1,:)+tmp1(:,:,2:end,:));
            tmp(:,:,end,:) = obj.scale*(tmp1(:,:,1,:)+tmp1(:,:,end,:));
            
            % app dim 4
            y(:,:,:,1:end-1,4) = obj.scale*(tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
            y(:,:,:,end,4) = obj.scale*(tmp(:,:,:,1)+tmp(:,:,:,end));  
            
            % =============================================================
            % det dim 3 app dim 1,2,4  
            % app dim1
            tmp(1:end-1,:,:,:) = obj.scale*(x_f(1:end-1,:,:,:)+x_f(2:end,:,:,:));
            tmp(end,:,:,:) = obj.scale*(x_f(1,:,:,:)+x_f(end,:,:,:));
  
            % app dim 2
            tmp1(:,1:end-1,:,:,:) = obj.scale*(tmp(:,1:end-1,:,:,:)+tmp(:,2:end,:,:,:));
            tmp1(:,end,:,:,:) = obj.scale*(tmp(:,1,:,:,:)+tmp(:,end,:,:,:));

            % app dim 3
            tmp(:,:,1:end-1,:) = obj.scale*(tmp1(:,:,1:end-1,:)-tmp1(:,:,2:end,:));
            tmp(:,:,end,:) = obj.scale*(-tmp1(:,:,1,:)+tmp1(:,:,end,:));
            
            % app dim 4
            y(:,:,:,1:end-1,5) = obj.scale*(tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
            y(:,:,:,end,5) = obj.scale*(tmp(:,:,:,1)+tmp(:,:,:,end)); 
            
            % =============================================================
            % det dim 1,3 app dim 2,4
            % app dim1
            tmp(1:end-1,:,:,:) = obj.scale*(x_f(1:end-1,:,:,:)-x_f(2:end,:,:,:));
            tmp(end,:,:,:) = obj.scale*(-x_f(1,:,:,:)+x_f(end,:,:,:));
  
            % app dim 2
            tmp1(:,1:end-1,:,:,:) = obj.scale*(tmp(:,1:end-1,:,:,:)+tmp(:,2:end,:,:,:));
            tmp1(:,end,:,:,:) = obj.scale*(tmp(:,1,:,:,:)+tmp(:,end,:,:,:));

            % app dim 3
            tmp(:,:,1:end-1,:) = obj.scale*(tmp1(:,:,1:end-1,:)-tmp1(:,:,2:end,:));
            tmp(:,:,end,:) = obj.scale*(-tmp1(:,:,1,:)+tmp1(:,:,end,:));
            
            % app dim 4
            y(:,:,:,1:end-1,6) = obj.scale*(tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
            y(:,:,:,end,6) = obj.scale*(tmp(:,:,:,1)+tmp(:,:,:,end));  
            
            % =============================================================
            % det dim 2,3 app dim 1, 4
            % app dim1
            tmp(1:end-1,:,:,:) = obj.scale*(x_f(1:end-1,:,:,:)+x_f(2:end,:,:,:));
            tmp(end,:,:,:) = obj.scale*(x_f(1,:,:,:)+x_f(end,:,:,:));
  
            % app dim 2
            tmp1(:,1:end-1,:,:,:) = obj.scale*(tmp(:,1:end-1,:,:,:)-tmp(:,2:end,:,:,:));
            tmp1(:,end,:,:,:) = obj.scale*(-tmp(:,1,:,:,:)+tmp(:,end,:,:,:));

            % app dim 3
            tmp(:,:,1:end-1,:) = obj.scale*(tmp1(:,:,1:end-1,:)-tmp1(:,:,2:end,:));
            tmp(:,:,end,:) = obj.scale*(-tmp1(:,:,1,:)+tmp1(:,:,end,:));
            
            % app dim 4
            y(:,:,:,1:end-1,7) = obj.scale*(tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
            y(:,:,:,end,7) = obj.scale*(tmp(:,:,:,1)+tmp(:,:,:,end));  
            
            % =============================================================
            % det dim 1,2,3, app dim 4
            % app dim1
            tmp(1:end-1,:,:,:) = obj.scale*(x_f(1:end-1,:,:,:)-x_f(2:end,:,:,:));
            tmp(end,:,:,:) = obj.scale*(-x_f(1,:,:,:)+x_f(end,:,:,:));
  
            % app dim 2
            tmp1(:,1:end-1,:,:,:) = obj.scale*(tmp(:,1:end-1,:,:,:)-tmp(:,2:end,:,:,:));
            tmp1(:,end,:,:,:) = obj.scale*(-tmp(:,1,:,:,:)+tmp(:,end,:,:,:));

            % app dim 3
            tmp(:,:,1:end-1,:) = obj.scale*(tmp1(:,:,1:end-1,:)-tmp1(:,:,2:end,:));
            tmp(:,:,end,:) = obj.scale*(-tmp1(:,:,1,:)+tmp1(:,:,end,:));
            
            % app dim 4
            y(:,:,:,1:end-1,8) = obj.scale*(tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
            y(:,:,:,end,8) = obj.scale*(tmp(:,:,:,1)+tmp(:,:,:,end));  
            
            % =============================================================
            % det dim 4 app dim 1,2,3
            % app dim1
            tmp(1:end-1,:,:,:) = obj.scale*(x_f(1:end-1,:,:,:)+x_f(2:end,:,:,:));
            tmp(end,:,:,:) = obj.scale*(x_f(1,:,:,:)+x_f(end,:,:,:));
  
            % app dim 2
            tmp1(:,1:end-1,:,:,:) = obj.scale*(tmp(:,1:end-1,:,:,:)+tmp(:,2:end,:,:,:));
            tmp1(:,end,:,:,:) = obj.scale*(tmp(:,1,:,:,:)+tmp(:,end,:,:,:));

            % app dim 3
            tmp(:,:,1:end-1,:) = obj.scale*(tmp1(:,:,1:end-1,:)+tmp1(:,:,2:end,:));
            tmp(:,:,end,:) = obj.scale*(tmp1(:,:,1,:)+tmp1(:,:,end,:));
            
            % app dim 4
            y(:,:,:,1:end-1,9) = obj.scale*(tmp(:,:,:,1:end-1)-tmp(:,:,:,2:end));
            y(:,:,:,end,9) = obj.scale*(-tmp(:,:,:,1)+tmp(:,:,:,end)); 
            
            % =============================================================
            % det dim 1,4, app dim 2,3
            % app dim1
            tmp(1:end-1,:,:,:) = obj.scale*(x_f(1:end-1,:,:,:)-x_f(2:end,:,:,:));
            tmp(end,:,:,:) = obj.scale*(-x_f(1,:,:,:)+x_f(end,:,:,:));
  
            % app dim 2
            tmp1(:,1:end-1,:,:,:) = obj.scale*(tmp(:,1:end-1,:,:,:)+tmp(:,2:end,:,:,:));
            tmp1(:,end,:,:,:) = obj.scale*(tmp(:,1,:,:,:)+tmp(:,end,:,:,:));

            % app dim 3
            tmp(:,:,1:end-1,:) = obj.scale*(tmp1(:,:,1:end-1,:)+tmp1(:,:,2:end,:));
            tmp(:,:,end,:) = obj.scale*(tmp1(:,:,1,:)+tmp1(:,:,end,:));
            
            % app dim 4
            y(:,:,:,1:end-1,10) = obj.scale*(tmp(:,:,:,1:end-1)-tmp(:,:,:,2:end));
            y(:,:,:,end,10) = obj.scale*(-tmp(:,:,:,1)+tmp(:,:,:,end));  
            
            % =============================================================
            % det dim 2,4, app dim 1,3
            % app dim1
            tmp(1:end-1,:,:,:) = obj.scale*(x_f(1:end-1,:,:,:)+x_f(2:end,:,:,:));
            tmp(end,:,:,:) = obj.scale*(x_f(1,:,:,:)+x_f(end,:,:,:));
  
            % app dim 2
            tmp1(:,1:end-1,:,:,:) = obj.scale*(tmp(:,1:end-1,:,:,:)-tmp(:,2:end,:,:,:));
            tmp1(:,end,:,:,:) = obj.scale*(-tmp(:,1,:,:,:)+tmp(:,end,:,:,:));

            % app dim 3
            tmp(:,:,1:end-1,:) = obj.scale*(tmp1(:,:,1:end-1,:)+tmp1(:,:,2:end,:));
            tmp(:,:,end,:) = obj.scale*(tmp1(:,:,1,:)+tmp1(:,:,end,:));
            
            % app dim 4
            y(:,:,:,1:end-1,11) = obj.scale*(tmp(:,:,:,1:end-1)-tmp(:,:,:,2:end));
            y(:,:,:,end,11) = obj.scale*(-tmp(:,:,:,1)+tmp(:,:,:,end));  
            
            % =============================================================
            % Det dim 1,2,4 app dim 3
            % app dim1
            tmp(1:end-1,:,:,:) = obj.scale*(x_f(1:end-1,:,:,:)-x_f(2:end,:,:,:));
            tmp(end,:,:,:) = obj.scale*(-x_f(1,:,:,:)+x_f(end,:,:,:));
  
            % app dim 2
            tmp1(:,1:end-1,:,:,:) = obj.scale*(tmp(:,1:end-1,:,:,:)-tmp(:,2:end,:,:,:));
            tmp1(:,end,:,:,:) = obj.scale*(-tmp(:,1,:,:,:)+tmp(:,end,:,:,:));

            % app dim 3
            tmp(:,:,1:end-1,:) = obj.scale*(tmp1(:,:,1:end-1,:)+tmp1(:,:,2:end,:));
            tmp(:,:,end,:) = obj.scale*(tmp1(:,:,1,:)+tmp1(:,:,end,:));
            
            % app dim 4
            y(:,:,:,1:end-1,12) = obj.scale*(tmp(:,:,:,1:end-1)-tmp(:,:,:,2:end));
            y(:,:,:,end,12) = obj.scale*(-tmp(:,:,:,1)+tmp(:,:,:,end));
            
            % =============================================================
            % det dim 3,4 app dim 1,2
            % app dim1
            tmp(1:end-1,:,:,:) = obj.scale*(x_f(1:end-1,:,:,:)+x_f(2:end,:,:,:));
            tmp(end,:,:,:) = obj.scale*(x_f(1,:,:,:)+x_f(end,:,:,:));
  
            % app dim 2
            tmp1(:,1:end-1,:,:,:) = obj.scale*(tmp(:,1:end-1,:,:,:)+tmp(:,2:end,:,:,:));
            tmp1(:,end,:,:,:) = obj.scale*(tmp(:,1,:,:,:)+tmp(:,end,:,:,:));

            % app dim 3
            tmp(:,:,1:end-1,:) = obj.scale*(tmp1(:,:,1:end-1,:)-tmp1(:,:,2:end,:));
            tmp(:,:,end,:) = obj.scale*(-tmp1(:,:,1,:)+tmp1(:,:,end,:));
            
            % app dim 4
            y(:,:,:,1:end-1,13) = obj.scale*(tmp(:,:,:,1:end-1)-tmp(:,:,:,2:end));
            y(:,:,:,end,13) = obj.scale*(-tmp(:,:,:,1)+tmp(:,:,:,end));  
            
            % =============================================================
            % det dim 1,3,4 app dim 2
            % app dim1
            tmp(1:end-1,:,:,:) = obj.scale*(x_f(1:end-1,:,:,:)-x_f(2:end,:,:,:));
            tmp(end,:,:,:) = obj.scale*(-x_f(1,:,:,:)+x_f(end,:,:,:));
  
            % app dim 2
            tmp1(:,1:end-1,:,:,:) = obj.scale*(tmp(:,1:end-1,:,:,:)+tmp(:,2:end,:,:,:));
            tmp1(:,end,:,:,:) = obj.scale*(tmp(:,1,:,:,:)+tmp(:,end,:,:,:));

            % app dim 3
            tmp(:,:,1:end-1,:) = obj.scale*(tmp1(:,:,1:end-1,:)-tmp1(:,:,2:end,:));
            tmp(:,:,end,:) = obj.scale*(-tmp1(:,:,1,:)+tmp1(:,:,end,:));
            
            % app dim 4
            y(:,:,:,1:end-1,14) = obj.scale*(tmp(:,:,:,1:end-1)-tmp(:,:,:,2:end));
            y(:,:,:,end,14) = obj.scale*(-tmp(:,:,:,1)+tmp(:,:,:,end));  
            
            % =============================================================
            % det dim 2, 3, 4 app dim 1
            % app dim1
            tmp(1:end-1,:,:,:) = obj.scale*(x_f(1:end-1,:,:,:)+x_f(2:end,:,:,:));
            tmp(end,:,:,:) = obj.scale*(x_f(1,:,:,:)+x_f(end,:,:,:));
  
            % app dim 2
            tmp1(:,1:end-1,:,:,:) = obj.scale*(tmp(:,1:end-1,:,:,:)-tmp(:,2:end,:,:,:));
            tmp1(:,end,:,:,:) = obj.scale*(-tmp(:,1,:,:,:)+tmp(:,end,:,:,:));

            % app dim 3
            tmp(:,:,1:end-1,:) = obj.scale*(tmp1(:,:,1:end-1,:)-tmp1(:,:,2:end,:));
            tmp(:,:,end,:) = obj.scale*(-tmp1(:,:,1,:)+tmp1(:,:,end,:));
            
            % app dim 4
            y(:,:,:,1:end-1,15) = obj.scale*(tmp(:,:,:,1:end-1)-tmp(:,:,:,2:end));
            y(:,:,:,end,15) = obj.scale*(-tmp(:,:,:,1)+tmp(:,:,:,end));  
            
            % =============================================================
            % det dim 1,2, 3, 4 app dim 1
            % app dim1
            tmp(1:end-1,:,:,:) = obj.scale*(x_f(1:end-1,:,:,:)-x_f(2:end,:,:,:));
            tmp(end,:,:,:) = obj.scale*(-x_f(1,:,:,:)+x_f(end,:,:,:));
  
            % app dim 2
            tmp1(:,1:end-1,:,:,:) = obj.scale*(tmp(:,1:end-1,:,:,:)-tmp(:,2:end,:,:,:));
            tmp1(:,end,:,:,:) = obj.scale*(-tmp(:,1,:,:,:)+tmp(:,end,:,:,:));

            % app dim 3
            tmp(:,:,1:end-1,:) = obj.scale*(tmp1(:,:,1:end-1,:)-tmp1(:,:,2:end,:));
            tmp(:,:,end,:) = obj.scale*(-tmp1(:,:,1,:)+tmp1(:,:,end,:));
            
            % app dim 4
            y(:,:,:,1:end-1,16) = obj.scale*(tmp(:,:,:,1:end-1)-tmp(:,:,:,2:end));
            y(:,:,:,end,16) = obj.scale*(-tmp(:,:,:,1)+tmp(:,:,:,end));
        end
        
        % Single Level Redundant Wavelet Reconstruction
        function y = level_1_rec(obj,x_f)
            
%             y = zeros(size(x_f(:,:,:,:,1)),'like',x_f);
            y = zeros(size(x_f),'like',x_f);
            tmp = zeros(obj.sizes,'like',x_f);
            tmp1 = zeros(obj.sizes,'like',x_f);
             % =============================================================
            % app dim 1,2,3,4 
            % app dim1
            tmp(2:end,:,:,:) = obj.scale*(x_f(1:end-1,:,:,:,1)+x_f(2:end,:,:,:,1));
            tmp(1,:,:,:) = obj.scale*(x_f(1,:,:,:,1)+x_f(end,:,:,:,1));
  
            % app dim 2
            tmp1(:,2:end,:,:,:) = obj.scale*(tmp(:,1:end-1,:,:,:)+tmp(:,2:end,:,:,:));
            tmp1(:,1,:,:,:) = obj.scale*(tmp(:,1,:,:,:)+tmp(:,end,:,:,:));

            % app dim 3
            tmp(:,:,2:end,:) = obj.scale*(tmp1(:,:,1:end-1,:)+tmp1(:,:,2:end,:));
            tmp(:,:,1,:) = obj.scale*(tmp1(:,:,1,:)+tmp1(:,:,end,:));
            
            % app dim 4
%             y(:,:,:,2:end) = obj.scale*(tmp(:,:,:,2:end)+tmp(:,:,:,2:end));
%             y(:,:,:,end) = obj.scale*(tmp(:,:,:,1)+tmp(:,:,:,end));  
            y(:,:,:,2:end,1) = obj.scale*(tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
            y(:,:,:,1,1) = obj.scale*(tmp(:,:,:,1)+tmp(:,:,:,end)); 
            
            % =============================================================
            % det dim 1, app dim 2,3,4 
            % det dim1
            tmp(2:end,:,:,:) = obj.scale*(-x_f(1:end-1,:,:,:,2)+x_f(2:end,:,:,:,2));
            tmp(1,:,:,:) = obj.scale*(x_f(1,:,:,:,2)-x_f(end,:,:,:,2));
  
            % app dim 2
            tmp1(:,2:end,:,:,:) = obj.scale*(tmp(:,1:end-1,:,:,:)+tmp(:,2:end,:,:,:));
            tmp1(:,1,:,:,:) = obj.scale*(tmp(:,1,:,:,:)+tmp(:,end,:,:,:));

            % app dim 3
            tmp(:,:,2:end,:) = obj.scale*(tmp1(:,:,1:end-1,:)+tmp1(:,:,2:end,:));
            tmp(:,:,1,:) = obj.scale*(tmp1(:,:,1,:)+tmp1(:,:,end,:));
            
            % app dim 4
%             y(:,:,:,2:end) = y(:,:,:,2:end) + obj.scale*(tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
%             y(:,:,:,1) = y(:,:,:,1) + obj.scale*(tmp(:,:,:,1)+tmp(:,:,:,end));
            y(:,:,:,2:end,2) = obj.scale*(tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
            y(:,:,:,1,2) = obj.scale*(tmp(:,:,:,1)+tmp(:,:,:,end)); 
            
            % =============================================================
            % det dim 2, app dim 1,3,4 
            % app dim 1
            tmp(2:end,:,:,:) = obj.scale*(x_f(1:end-1,:,:,:,3)+x_f(2:end,:,:,:,3));
            tmp(1,:,:,:) = obj.scale*(x_f(1,:,:,:,3)+x_f(end,:,:,:,3));
  
            % det dim 2
            tmp1(:,2:end,:,:,:) = obj.scale*(-tmp(:,1:end-1,:,:,:)+tmp(:,2:end,:,:,:));
            tmp1(:,1,:,:,:) = obj.scale*(tmp(:,1,:,:,:)-tmp(:,end,:,:,:));

            % app dim 3
            tmp(:,:,2:end,:) = obj.scale*(tmp1(:,:,1:end-1,:)+tmp1(:,:,2:end,:));
            tmp(:,:,1,:) = obj.scale*(tmp1(:,:,1,:)+tmp1(:,:,end,:));
            
            % app dim 4
%             y(:,:,:,2:end) = y(:,:,:,2:end) + obj.scale*(tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
%             y(:,:,:,1) = y(:,:,:,1) + obj.scale*(tmp(:,:,:,1)+tmp(:,:,:,end));
            y(:,:,:,2:end,3) = obj.scale*(tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
            y(:,:,:,1,3) = obj.scale*(tmp(:,:,:,1)+tmp(:,:,:,end)); 
            
            % =============================================================
            % det dim 1,2 app dim 3,4 
            % app dim1
            tmp(2:end,:,:,:) = obj.scale*(-x_f(1:end-1,:,:,:,4)+x_f(2:end,:,:,:,4));
            tmp(1,:,:,:) = obj.scale*(x_f(1,:,:,:,4)-x_f(end,:,:,:,4));
  
            % app dim 2
            tmp1(:,2:end,:,:,:) = obj.scale*(-tmp(:,1:end-1,:,:,:)+tmp(:,2:end,:,:,:));
            tmp1(:,1,:,:,:) = obj.scale*(tmp(:,1,:,:,:)-tmp(:,end,:,:,:));

            % app dim 3
            tmp(:,:,2:end,:) = obj.scale*(tmp1(:,:,1:end-1,:)+tmp1(:,:,2:end,:));
            tmp(:,:,1,:) = obj.scale*(tmp1(:,:,1,:)+tmp1(:,:,end,:));
            
            % app dim 4
%             y(:,:,:,2:end) = y(:,:,:,2:end)  + obj.scale*(tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
%             y(:,:,:,1) = y(:,:,:,1) + obj.scale*(tmp(:,:,:,1)+tmp(:,:,:,end));  
            y(:,:,:,2:end,4) = obj.scale*(tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
            y(:,:,:,1,4) = obj.scale*(tmp(:,:,:,1)+tmp(:,:,:,end)); 
            
            % =============================================================
            % det dim 3 app dim 1,2,4  
            % app dim1
            tmp(2:end,:,:,:) = obj.scale*(x_f(1:end-1,:,:,:,5)+x_f(2:end,:,:,:,5));
            tmp(1,:,:,:) = obj.scale*(x_f(1,:,:,:,5)+x_f(end,:,:,:,5));
  
            % app dim 2
            tmp1(:,2:end,:,:,:) = obj.scale*(tmp(:,1:end-1,:,:,:)+tmp(:,2:end,:,:,:));
            tmp1(:,1,:,:,:) = obj.scale*(tmp(:,1,:,:,:)+tmp(:,end,:,:,:));

            % app dim 3
            tmp(:,:,2:end,:) = obj.scale*(-tmp1(:,:,1:end-1,:)+tmp1(:,:,2:end,:));
            tmp(:,:,1,:) = obj.scale*(tmp1(:,:,1,:)-tmp1(:,:,end,:));
            
%             app dim 4
%             y(:,:,:,2:end) = y(:,:,:,2:end)  + obj.scale*(tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
%             y(:,:,:,1) = y(:,:,:,1) + obj.scale*(tmp(:,:,:,1)+tmp(:,:,:,end));  
            y(:,:,:,2:end,5) = obj.scale*(tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
            y(:,:,:,1,5) = obj.scale*(tmp(:,:,:,1)+tmp(:,:,:,end)); 
            
            % =============================================================
            % det dim 1,3 app dim 2,4
            % app dim1
            tmp(2:end,:,:,:) = obj.scale*(-x_f(1:end-1,:,:,:,6)+x_f(2:end,:,:,:,6));
            tmp(1,:,:,:) = obj.scale*(x_f(1,:,:,:,6)-x_f(end,:,:,:,6));
  
            % app dim 2
            tmp1(:,2:end,:,:,:) = obj.scale*(tmp(:,1:end-1,:,:,:)+tmp(:,2:end,:,:,:));
            tmp1(:,1,:,:,:) = obj.scale*(tmp(:,1,:,:,:)+tmp(:,end,:,:,:));

            % app dim 3
            tmp(:,:,2:end,:) = obj.scale*(-tmp1(:,:,1:end-1,:)+tmp1(:,:,2:end,:));
            tmp(:,:,1,:) = obj.scale*(tmp1(:,:,1,:)-tmp1(:,:,end,:));
            
            % app dim 4
%             y(:,:,:,2:end) = y(:,:,:,2:end)  + obj.scale*(tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
%             y(:,:,:,1) = y(:,:,:,1) + obj.scale*(tmp(:,:,:,1)+tmp(:,:,:,end)); 
            y(:,:,:,2:end,6) = obj.scale*(tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
            y(:,:,:,1,6) = obj.scale*(tmp(:,:,:,1)+tmp(:,:,:,end)); 
            
            % =============================================================
            % det dim 2,3 app dim 1, 4
            % app dim1
            tmp(2:end,:,:,:) = obj.scale*(x_f(1:end-1,:,:,:,7)+x_f(2:end,:,:,:,7));
            tmp(1,:,:,:) = obj.scale*(x_f(1,:,:,:,7)+x_f(end,:,:,:,7));
  
            % app dim 2
            tmp1(:,2:end,:,:,:) = obj.scale*(-tmp(:,1:end-1,:,:,:)+tmp(:,2:end,:,:,:));
            tmp1(:,1,:,:,:) = obj.scale*(tmp(:,1,:,:,:)-tmp(:,end,:,:,:));

            % app dim 3
            tmp(:,:,2:end,:) = obj.scale*(-tmp1(:,:,1:end-1,:)+tmp1(:,:,2:end,:));
            tmp(:,:,1,:) = obj.scale*(tmp1(:,:,1,:)-tmp1(:,:,end,:));
            
            % app dim 4
%             y(:,:,:,2:end) = y(:,:,:,2:end)  + obj.scale*(tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
%             y(:,:,:,1) = y(:,:,:,1) + obj.scale*(tmp(:,:,:,1)+tmp(:,:,:,end));   
            y(:,:,:,2:end,7) = obj.scale*(tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
            y(:,:,:,1,7) = obj.scale*(tmp(:,:,:,1)+tmp(:,:,:,end)); 
            
            % =============================================================
            % det dim 1,2,3, app dim 4
            % app dim1
            tmp(2:end,:,:,:) = obj.scale*(-x_f(1:end-1,:,:,:,8)+x_f(2:end,:,:,:,8));
            tmp(1,:,:,:) = obj.scale*(x_f(1,:,:,:,8)-x_f(end,:,:,:,8));
  
            % app dim 2
            tmp1(:,2:end,:,:,:) = obj.scale*(-tmp(:,1:end-1,:,:,:)+tmp(:,2:end,:,:,:));
            tmp1(:,1,:,:,:) = obj.scale*(tmp(:,1,:,:,:)-tmp(:,end,:,:,:));

            % app dim 3
            tmp(:,:,2:end,:) = obj.scale*(-tmp1(:,:,1:end-1,:)+tmp1(:,:,2:end,:));
            tmp(:,:,1,:) = obj.scale*(tmp1(:,:,1,:)-tmp1(:,:,end,:));
            
            % app dim 4
%             y(:,:,:,2:end) = y(:,:,:,2:end)  + obj.scale*(tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
%             y(:,:,:,1) = y(:,:,:,1) + obj.scale*(tmp(:,:,:,1)+tmp(:,:,:,end));
            y(:,:,:,2:end,8) = obj.scale*(tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
            y(:,:,:,1,8) = obj.scale*(tmp(:,:,:,1)+tmp(:,:,:,end)); 
            
            % =============================================================
            % det dim 4 app dim 1,2,3
            % app dim1
            tmp(2:end,:,:,:) = obj.scale*(x_f(1:end-1,:,:,:,9)+x_f(2:end,:,:,:,9));
            tmp(1,:,:,:) = obj.scale*(x_f(1,:,:,:,9)+x_f(end,:,:,:,9));
  
            % app dim 2
            tmp1(:,2:end,:,:,:) = obj.scale*(tmp(:,1:end-1,:,:,:)+tmp(:,2:end,:,:,:));
            tmp1(:,1,:,:,:) = obj.scale*(tmp(:,1,:,:,:)+tmp(:,end,:,:,:));

            % app dim 3
            tmp(:,:,2:end,:) = obj.scale*(tmp1(:,:,1:end-1,:)+tmp1(:,:,2:end,:));
            tmp(:,:,1,:) = obj.scale*(tmp1(:,:,1,:)+tmp1(:,:,end,:));
            
            % app dim 4
%             y(:,:,:,2:end) = y(:,:,:,2:end)  + obj.scale*(-tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
%             y(:,:,:,1) = y(:,:,:,1) + obj.scale*(tmp(:,:,:,1)-tmp(:,:,:,end));
            y(:,:,:,2:end,9) = obj.scale*(-tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
            y(:,:,:,1,9) = obj.scale*(tmp(:,:,:,1)-tmp(:,:,:,end));
            
            % =============================================================
            % det dim 1,4, app dim 2,3
            % app dim1
            tmp(2:end,:,:,:) = obj.scale*(-x_f(1:end-1,:,:,:,10)+x_f(2:end,:,:,:,10));
            tmp(1,:,:,:) = obj.scale*(x_f(1,:,:,:,10)-x_f(end,:,:,:,10));
  
            % app dim 2
            tmp1(:,2:end,:,:,:) = obj.scale*(tmp(:,1:end-1,:,:,:)+tmp(:,2:end,:,:,:));
            tmp1(:,1,:,:,:) = obj.scale*(tmp(:,1,:,:,:)+tmp(:,end,:,:,:));

            % app dim 3
            tmp(:,:,2:end,:) = obj.scale*(tmp1(:,:,1:end-1,:)+tmp1(:,:,2:end,:));
            tmp(:,:,1,:) = obj.scale*(tmp1(:,:,1,:)+tmp1(:,:,end,:));
            
            % app dim 4
%             y(:,:,:,2:end) = y(:,:,:,2:end) + obj.scale*(-tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
%             y(:,:,:,1) = y(:,:,:,1) + obj.scale*(tmp(:,:,:,1)-tmp(:,:,:,end)); 
            y(:,:,:,2:end,10) = obj.scale*(-tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
            y(:,:,:,1,10) = obj.scale*(tmp(:,:,:,1)-tmp(:,:,:,end));
            
            % =============================================================
            % det dim 2,4, app dim 1,3
            % app dim1
            tmp(2:end,:,:,:) = obj.scale*(x_f(1:end-1,:,:,:,11)+x_f(2:end,:,:,:,11));
            tmp(1,:,:,:) = obj.scale*(x_f(1,:,:,:,11)+x_f(end,:,:,:,11));
  
            % app dim 2
            tmp1(:,2:end,:,:,:) = obj.scale*(-tmp(:,1:end-1,:,:,:)+tmp(:,2:end,:,:,:));
            tmp1(:,1,:,:,:) = obj.scale*(tmp(:,1,:,:,:)-tmp(:,end,:,:,:));

            % app dim 3
            tmp(:,:,2:end,:) = obj.scale*(tmp1(:,:,1:end-1,:)+tmp1(:,:,2:end,:));
            tmp(:,:,1,:) = obj.scale*(tmp1(:,:,1,:)+tmp1(:,:,end,:));
            
            % app dim 4
%             y(:,:,:,2:end) = y(:,:,:,2:end)  + obj.scale*(-tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
%             y(:,:,:,1) = y(:,:,:,1) + obj.scale*(tmp(:,:,:,1)-tmp(:,:,:,end));
            y(:,:,:,2:end,11) = obj.scale*(-tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
            y(:,:,:,1,11) = obj.scale*(tmp(:,:,:,1)-tmp(:,:,:,end));
            
            % =============================================================
            % Det dim 1,2,4 app dim 3
            % app dim1
            tmp(2:end,:,:,:) = obj.scale*(-x_f(1:end-1,:,:,:,12)+x_f(2:end,:,:,:,12));
            tmp(1,:,:,:) = obj.scale*(x_f(1,:,:,:,12)-x_f(end,:,:,:,12));
  
            % app dim 2
            tmp1(:,2:end,:,:,:) = obj.scale*(-tmp(:,1:end-1,:,:,:)+tmp(:,2:end,:,:,:));
            tmp1(:,1,:,:,:) = obj.scale*(tmp(:,1,:,:,:)-tmp(:,end,:,:,:));

            % app dim 3
            tmp(:,:,2:end,:) = obj.scale*(tmp1(:,:,1:end-1,:)+tmp1(:,:,2:end,:));
            tmp(:,:,1,:) = obj.scale*(tmp1(:,:,1,:)+tmp1(:,:,end,:));
            
            % app dim 4
%             y(:,:,:,2:end) = y(:,:,:,2:end)  + obj.scale*(-tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
%             y(:,:,:,1) = y(:,:,:,1) + obj.scale*(tmp(:,:,:,1)-tmp(:,:,:,end));
            y(:,:,:,2:end,12) = obj.scale*(-tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
            y(:,:,:,1,12) = obj.scale*(tmp(:,:,:,1)-tmp(:,:,:,end));
            
            % =============================================================
            % det dim 3,4 app dim 1,2
            % app dim1
            tmp(2:end,:,:,:) = obj.scale*(x_f(1:end-1,:,:,:,13)+x_f(2:end,:,:,:,13));
            tmp(1,:,:,:) = obj.scale*(x_f(1,:,:,:,13)+x_f(end,:,:,:,13));
  
            % app dim 2
            tmp1(:,2:end,:,:,:) = obj.scale*(tmp(:,1:end-1,:,:,:)+tmp(:,2:end,:,:,:));
            tmp1(:,1,:,:,:) = obj.scale*(tmp(:,1,:,:,:)+tmp(:,end,:,:,:));

            % app dim 3
            tmp(:,:,2:end,:) = obj.scale*(-tmp1(:,:,1:end-1,:)+tmp1(:,:,2:end,:));
            tmp(:,:,1,:) = obj.scale*(tmp1(:,:,1,:)-tmp1(:,:,end,:));
            
            % app dim 4
%             y(:,:,:,2:end) = y(:,:,:,2:end)  + obj.scale*(-tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
%             y(:,:,:,1) = y(:,:,:,1) + obj.scale*(tmp(:,:,:,1)-tmp(:,:,:,end));
            y(:,:,:,2:end,13) = obj.scale*(-tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
            y(:,:,:,1,13) = obj.scale*(tmp(:,:,:,1)-tmp(:,:,:,end));
            
            % =============================================================
            % det dim 1,3,4 app dim 2
            % app dim1
            tmp(2:end,:,:,:) = obj.scale*(-x_f(1:end-1,:,:,:,14)+x_f(2:end,:,:,:,14));
            tmp(1,:,:,:) = obj.scale*(x_f(1,:,:,:,14)-x_f(end,:,:,:,14));
  
            % app dim 2
            tmp1(:,2:end,:,:,:) = obj.scale*(tmp(:,1:end-1,:,:,:)+tmp(:,2:end,:,:,:));
            tmp1(:,1,:,:,:) = obj.scale*(tmp(:,1,:,:,:)+tmp(:,end,:,:,:));

            % app dim 3
            tmp(:,:,2:end,:) = obj.scale*(-tmp1(:,:,1:end-1,:)+tmp1(:,:,2:end,:));
            tmp(:,:,1,:) = obj.scale*(tmp1(:,:,1,:)-tmp1(:,:,end,:));
            
            % app dim 4
%             y(:,:,:,2:end) = y(:,:,:,2:end)  + obj.scale*(-tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
%             y(:,:,:,1) = y(:,:,:,1) + obj.scale*(tmp(:,:,:,1)-tmp(:,:,:,end)); 
            y(:,:,:,2:end,14) = obj.scale*(-tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
            y(:,:,:,1,14) = obj.scale*(tmp(:,:,:,1)-tmp(:,:,:,end));
            
            % =============================================================
            % det dim 2, 3, 4 app dim 1
            % app dim1
            tmp(2:end,:,:,:) = obj.scale*(x_f(1:end-1,:,:,:,15)+x_f(2:end,:,:,:,15));
            tmp(1,:,:,:) = obj.scale*(x_f(1,:,:,:,15)+x_f(end,:,:,:,15));
  
            % app dim 2
            tmp1(:,2:end,:,:,:) = obj.scale*(-tmp(:,1:end-1,:,:,:)+tmp(:,2:end,:,:,:));
            tmp1(:,1,:,:,:) = obj.scale*(tmp(:,1,:,:,:)-tmp(:,end,:,:,:));

            % app dim 3
            tmp(:,:,2:end,:) = obj.scale*(-tmp1(:,:,1:end-1,:)+tmp1(:,:,2:end,:));
            tmp(:,:,1,:) = obj.scale*(+tmp1(:,:,1,:)-tmp1(:,:,end,:));
            
            % app dim 4
%             y(:,:,:,2:end) = y(:,:,:,2:end)  + obj.scale*(-tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
%             y(:,:,:,1) = y(:,:,:,1) + obj.scale*(tmp(:,:,:,1)-tmp(:,:,:,end));
            y(:,:,:,2:end,15) = obj.scale*(-tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
            y(:,:,:,1,15) = obj.scale*(tmp(:,:,:,1)-tmp(:,:,:,end));
            
            % =============================================================
            % det dim 1,2, 3, 4 app dim 1
            % app dim1
            tmp(2:end,:,:,:) = obj.scale*(-x_f(1:end-1,:,:,:,16)+x_f(2:end,:,:,:,16));
            tmp(1,:,:,:) = obj.scale*(x_f(1,:,:,:,16)-x_f(end,:,:,:,16));
  
            % app dim 2
            tmp1(:,2:end,:,:,:) = obj.scale*(-tmp(:,1:end-1,:,:,:)+tmp(:,2:end,:,:,:));
            tmp1(:,1,:,:,:) = obj.scale*(tmp(:,1,:,:,:)-tmp(:,end,:,:,:));

            % app dim 3
            tmp(:,:,2:end,:) = obj.scale*(-tmp1(:,:,1:end-1,:)+tmp1(:,:,2:end,:));
            tmp(:,:,1,:) = obj.scale*(+tmp1(:,:,1,:)-tmp1(:,:,end,:));
            
            % app dim 4
%             y(:,:,:,2:end) = y(:,:,:,2:end)  + obj.scale*(-tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
%             y(:,:,:,1) = y(:,:,:,1) + obj.scale*(tmp(:,:,:,1)-tmp(:,:,:,end));
            y(:,:,:,2:end,16) = obj.scale*(-tmp(:,:,:,1:end-1)+tmp(:,:,:,2:end));
            y(:,:,:,1,16) = obj.scale*(tmp(:,:,:,1)-tmp(:,:,:,end));
            
            y = sum(y,5);
        end
    end
end

