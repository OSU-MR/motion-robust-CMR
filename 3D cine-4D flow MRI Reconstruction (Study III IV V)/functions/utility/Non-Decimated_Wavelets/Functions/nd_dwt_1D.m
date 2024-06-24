%Calculates the 1 Dimensional nondecimated Wavelet transform using fast 
%    convlolution and periodic boundry conditions.  The Filters are 
%    pre-computed and stored in f_dec and f_rec to reduce computation time 
%    when performing Multiple forward/backward transforms
%
%Methods:
%   ud_dwt3D:   Constructor
%               Inputs: wname - Wavelet Filters to Use i.e. db1,db2,etc.
%                        Either a string 'db1' or a cell {'db4','db1'} where
%                        the first element in the string is filter for the 
%                        spatial domain and the second is the filter for the 
%                        time domain
%
%                       sizes - length of the 1D signal
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
%                        2D array where the data is arranged [n1,bands]
%                        The bands are orginized as follows.The coefficients 
%                        are ordered as follows,"L", "H" where"H" denotes 
%                        the high frequency filter and L represents the low 
%                        frequency filter. Successive levels of decomposition 
%                        are stacked such that the highest "L" is in [n1,1]
%
%   rec:        Multilevel Reconstruction
%               Inputs: x - Wavelet coefficients in a 2D array size 
%                       [n1,bands].
%
%               Outputs: y - Reconstructed 1D array.%   
%
%**************************************************************************
% The Ohio State University
% Written by:   Adam Rich 
% Last update:  2/5/2015
%**************************************************************************

classdef nd_dwt_1D
    %ND_DWT_1D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        f_dec;          % Decomposition Filters
        sizes;          % Size of the Signal
        f_size;         % Length of the filters
        wname;          % Wavelet used
        pres_l2_norm;   % Binary indicator to preserver l2 norm of coefficients
        compute;        % How to compute the wavelet transform
        precision;      % which precision to use
    end
    
    methods
        % Constructor Object
        function obj = nd_dwt_1D(wname,sizes,varargin)
            % Check Inputs
            if ischar(wname)
                obj.wname = {wname,wname};
            elseif iscell(wname)
                error('Wavelet Name Must be a string');
            end
            
            if length(sizes) ~= 1
                error('1D array length must be a scalar');
            else
                obj.sizes = sizes;
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
            if size(x,1) ==1
                x = x.';
            end
            x = fft(x);
            
            % Use Mex computation
			if strcmpi(obj.compute,'mex')
				y = nd_dwt_mex(x,obj.f_dec,0,level,obj.pres_l2_norm);
                
            % Use Matlab 
			else
                % Preallocate
                if ( strcmpi(obj.compute,'gpu_off') || strcmpi(obj.compute,'gpu') ) && strcmpi(obj.precision,'single')
                    y = zeros([obj.sizes, 2+(level-1)],'single');
                    y = gpuArray(y);
                elseif ( strcmpi(obj.compute,'gpu_off') || strcmpi(obj.compute,'gpu') ) && ~strcmpi(obj.precision,'single')
                    y = zeros([obj.sizes, 2+(level-1)],'gpuArray');
                elseif strcmpi(obj.precision,'single')
                    y = zeros([obj.sizes, 2+(level-1)],'single');
                else
                    y = zeros([obj.sizes, 2+(level-1)]);
                end

            	% Calculate Mutlilevel Wavelet decomposition
           	 	for ind = 1:level
              	  	% First Level
               	 	if ind ==1
                 	   y = level_1_dec(obj,x);
               		% Succssive Levels
               	 	else
                 	   y = cat(2,level_1_dec(obj,fftn(squeeze(y(:,1)))), y(:,2:end));
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
            level = ceil(size(x,2)-1);
            
            % Fourier Transform of Signal
            x = fft(x,[],1);
            
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
                        	y = y/2;
                    	end
                		% Succssive Levels
                	else
                    	y = fft(y);
                    	y = level_1_rec(obj,cat(2,y,x(:,3+(ind-2))));
                    	if ~obj.pres_l2_norm
                        	y = y/2;
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
        % Decomposition Filters
        
            % Get Filters for Spatial Domain
            [LO_D,HI_D] = wave_filters(wname{1});
            
            % Find the filter size
            f_size.s1 = length(LO_D);
            
            % Dimension Check
            if f_size.s1 > obj.sizes(1)
                error(['First Dimension of Data is shorter than the wavelet'...
                    ,' filter being used']);
            end
            
            % Add a circularshift of half the filter length to the 
            % reconstruction filters by adding phase to them
            shift = exp(1j*2*pi*f_size.s1/2*linspace(0,1-1/obj.sizes,obj.sizes));
            
            % Take the Fourier Transform of the Kernels for Fast
            % Convolution
            if obj.pres_l2_norm 
                scale = 1/sqrt(2);
            else
                scale = 1;
            end
            if strcmpi(obj.compute,'mex')
                scale2 = 1/prod(obj.sizes);
            else
                scale2 = 1;
            end
            f_dec(:,1) = scale2*scale*(shift.*fft(LO_D,obj.sizes)).';
            f_dec(:,2) = scale2*scale*(shift.*fft(HI_D,obj.sizes)).';
        end
        
        % Single Level Redundant Wavelet Decomposition
        function y = level_1_dec(obj,x_f)
            % Preallocate
            if ( strcmpi(obj.compute,'gpu_off') || strcmpi(obj.compute,'gpu') ) && strcmpi(obj.precision,'single')
                y = zeros([obj.sizes,2],'single');
                y = gpuArray(y);
            elseif ( strcmpi(obj.compute,'gpu_off') || strcmpi(obj.compute,'gpu') ) && ~strcmpi(obj.precision,'single')
                y = zeros([obj.sizes,2],'gpuArray');
            elseif strcmpi(obj.precision,'single')
                y = zeros([obj.sizes,2],'single');
            else
                y = zeros([obj.sizes,2]);
            end
            
            % Calculate Wavelet Coefficents Using Fast Convolution
            y(:,1) = ifft(x_f.*obj.f_dec(:,1));
            y(:,2) = ifft(x_f.*obj.f_dec(:,2));
            
        end
        
        % Single Level Redundant Wavelet Reconstruction
        function y = level_1_rec(obj,x_f)
            
            % Reconstruct the 3D array using Fast Convolution
            y = ifft(squeeze(x_f(:,1)).*conj(obj.f_dec(:,1)) + squeeze(x_f(:,2)).*conj(obj.f_dec(:,2)));
            
        end
    end
    
end

