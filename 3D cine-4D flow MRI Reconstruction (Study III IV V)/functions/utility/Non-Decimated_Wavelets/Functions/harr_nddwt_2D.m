%Calculates the nondecimated Wavelet transform using fast convlolution and
%   periodic boundry conditions.  The Filters are pre-computed and stored
%   in f_dec and f_rec to reduce computation time when performing Multiple
%   forward/backward transforms
%
%
%Methods:
%   ud_dwt3D:   Constructor
%               Inputs: wname - Wavelet Filters to Use i.e. db1,db2,etc.
%                        Either a string 'db1' or a cell {'db4','db1'} where
%                        the first element in the string is filter for the 
%                        spatial domain and the second is the filter for the 
%                        time domain
%
%                       sizes - size of the 2D object [n1,n2]
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
%                       levels - Number of decomposition Levels
%               Outputs: y - Multilevel non-decimated DWT coefficients in a
%                        4D array where the data is arranged [n1,n2,bands]
%                        The bands are orginized as follows.  Let "HVT"
%                        represent the Horizontal, and Vertical Bands.  The
%                        coefficients are ordered as follows, "LL", "HL", 
%                        "LH", "HH", where "H" denotes the high frequency 
%                        filter and L represents the low frequency filter. 
%                        Successive levels of decomposition are stacked such 
%                        that the highest "LL" is in [n1,n2,1]
%
%   rec:        Multilevel Reconstruction
%               Inputs: x - Wavelet coefficients in a 3D array size 
%                       [n1,n2,bands].
%
%               Outputs: y - Reconstructed 2D array.%   
%
%**************************************************************************
% The Ohio State University
% Written by:   Adam Rich 
% Last update:  2/5/2015
%**************************************************************************

classdef harr_nddwt_2D
    properties
        sizes;          % Size of the 3D Image
        f_size;         % Length of the filters
        wname;          % Wavelet used
        pres_l2_norm;   % Binary indicator to preserver l2 norm of coefficients
        compute;        % How to compute the wavelet transform
        precision;      % which precision to use
        scale;
    end
    
    %% Public Methods
    methods
        % Constructor Object
        function obj = harr_nddwt_2D(wname,sizes,varargin)
            % Set Image size
            if length(sizes) ~=2
                error('The sizes vector must be length 2');
            else
                obj.sizes = sizes;
            end
            
            % Check wname input
            if mod(length(varargin),2)
                error('Optional inputs must come in pairs')
            end
            if ischar(wname)
                obj.wname = {wname,wname};
            elseif iscell(wname)
                if length(wname) ==2
                    obj.wname = wname;
                else
                    error(['You must specify two filter names in a cell array'...
                            ,'of length 2, or a single string for the same'...
                            ,' filter to be used in all dimensions']);
                end
            end
            
            % Set Default Options
            obj.pres_l2_norm = 0;
            obj.precision = 'double';
            obj.compute = 'mat';
            
            % Copy any optional inputs 
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
            
            % set kernel scale
            if obj.pres_l2_norm 
                obj.scale = 1/2;
            else
                obj.scale = 1/sqrt(2);
            end
            
            if strcmpi(obj.compute,'mex') && strcmpi(obj.precision,'single')
                error('Single precsision is not currently supported for mex computation');
            end            
        end
        
        % Multilevel Undecimated Wavelet Decomposition
        function y = dec(obj,x,level)
            
            if level ~=1
                error('Only single level decomposition supported for Harr')
            end
            
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
                    y = zeros([obj.sizes, 4+3*(level-1)],'single');
                    y = gpuArray(y);
                elseif ( strcmpi(obj.compute,'gpu_off') || strcmpi(obj.compute,'gpu') ) && ~strcmpi(obj.precision,'single')
                    y = zeros([obj.sizes, 4+3*(level-1)],'gpuArray');
                elseif strcmpi(obj.precision,'single')
                    y = zeros([obj.sizes, 4+3*(level-1)],'single');
                else
                    y = zeros([obj.sizes, 4+3*(level-1)]);
                end

                % Calculate Mutlilevel Wavelet decomposition
                for ind = 1:level
                    % First Level
                    if ind ==1
                        y = level_1_dec(obj,x);
                    % Succssive Levels
                    else
                        y = cat(3,level_1_dec(obj,fftn(squeeze(y(:,:,1)))), y(:,:,2:end));
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
            level = 1+(size(x,3)-4)/3;
            
            % Use c mex version if chosen
            if strcmpi(obj.compute,'mex')
                y = nd_dwt_mex(x,obj.f_dec,1,level,obj.pres_l2_norm);
            else
                
                % Reconstruct from Multiple Levels
                for ind = 1:level
                    % First Level
                    if ind ==1
                        y = level_1_rec(obj,x);
                        if ~obj.pres_l2_norm
                            y = y/4;
                        end
                    % Succssive Levels
                    else
                        y = fftn(y);
                        y = level_1_rec(obj,cat(3,y,x(:,:,5+(ind-2)*3:7+(ind-2)*3)));
                        if ~obj.pres_l2_norm
                            y = y/4;
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
                y = zeros([obj.sizes,4],'single');
                y = gpuArray(y);
            elseif ( strcmpi(obj.compute,'gpu_off') || strcmpi(obj.compute,'gpu') ) && ~strcmpi(obj.precision,'single')
                y = zeros([obj.sizes,4],'gpuArray');
            elseif strcmpi(obj.precision,'single')
                y = zeros([obj.sizes,4],'single');
            else
                y = zeros([obj.sizes,4]);
            end
            
            % Calculate Wavelet Coefficents Using Convolution
            ap_1 = zeros(obj.sizes,'like',x_f);
            det_1 = zeros(obj.sizes,'like',x_f);
            ap_1(1:end-1,:) = obj.scale*(x_f(1:end-1,:)+x_f(2:end,:));
            ap_1(end,:) = obj.scale*(x_f(1,:)+x_f(end,:));
            
            det_1(1:end-1,:) = obj.scale*(x_f(1:end-1,:)-x_f(2:end,:));
            det_1(end,:) = obj.scale*(-x_f(1,:)+x_f(end,:));
            
            % approximate dim 1 and 2
            y(:,1:end-1,1) = obj.scale*(ap_1(:,1:end-1)+ap_1(:,2:end));
            y(:,end,1) = obj.scale*(ap_1(:,1)+ap_1(:,end));
            
            % details dim1 approximate dim 2
            y(:,1:end-1,2) = obj.scale*(det_1(:,1:end-1)+det_1(:,2:end));
            y(:,end,2) = obj.scale*(det_1(:,1)+det_1(:,end));
            
            % approximate dim 1 details dim 2
            y(:,1:end-1,3) = obj.scale*(ap_1(:,1:end-1)-ap_1(:,2:end));
            y(:,end,3) = obj.scale*(-ap_1(:,1)+ap_1(:,end));
            
            % approximate dim 1 details dim 2
            y(:,1:end-1,4) = obj.scale*(det_1(:,1:end-1)-det_1(:,2:end));
            y(:,end,4) = obj.scale*(-det_1(:,1)+det_1(:,end));          
        end
        
        % Single Level Redundant Wavelet Reconstruction
        function y = level_1_rec(obj,x_f)
            
            ap_1 = zeros(obj.sizes,'like',x_f);
   
            % appriximate dim 1 and 2
            y = zeros(size(x_f(:,:,1)),'like',x_f);

            ap_1(2:end,:) = obj.scale*(x_f(1:end-1,:,1)+x_f(2:end,:,1));
            ap_1(1,:) = obj.scale*(x_f(1,:,1)+x_f(end,:,1));
            
            y(:,2:end) = obj.scale*(ap_1(:,1:end-1,1)+ap_1(:,2:end,1));
            y(:,1) = obj.scale*(ap_1(:,1,1)+ap_1(:,end,1));
            
            % Details dim 1 approximate dim 2
            ap_1(2:end,:) = obj.scale*(-x_f(1:end-1,:,2)+x_f(2:end,:,2));
            ap_1(1,:) = obj.scale*(x_f(1,:,2)-x_f(end,:,2));
            
            y(:,2:end) = y(:,2:end) + obj.scale*(ap_1(:,1:end-1,1)+ap_1(:,2:end,1));
            y(:,1) = y(:,1) + obj.scale*(ap_1(:,1,1)+ap_1(:,end,1));
            
            % approximate dim 1 details dim 2
            ap_1(2:end,:) = obj.scale*(x_f(1:end-1,:,3)+x_f(2:end,:,3));
            ap_1(1,:) = obj.scale*(x_f(1,:,3)+x_f(end,:,3));
            
            y(:,2:end) = y(:,2:end) + obj.scale*(-ap_1(:,1:end-1)+ap_1(:,2:end));
            y(:,1) = y(:,1) + obj.scale*(ap_1(:,1)-ap_1(:,end));
            
            % Details dim 1 and dim 2
            ap_1(2:end,:) = obj.scale*(-x_f(1:end-1,:,4)+x_f(2:end,:,4));
            ap_1(1,:) = obj.scale*(x_f(1,:,4)-x_f(end,:,4));
            
            y(:,2:end) = y(:,2:end) + obj.scale*(-ap_1(:,1:end-1)+ap_1(:,2:end));
            y(:,1) = y(:,1) + obj.scale*(ap_1(:,1)-ap_1(:,end));
        end
    end
end

