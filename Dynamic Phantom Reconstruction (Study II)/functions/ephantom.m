function data = ephantom(tfun, param)

% Generates a moving 2D virtual phantom
%
% data: k-space data.  Assumes FOV = 1 in all dimensions
% tfun: a vector of length nt with numbers in the range -1 to 1
%                  that indecates the variation of the parameters (par +
%                  tfun*delta_par).  If not present, defaults to a sinusoide.
%
% Pablo Irarrazaval
% 25-11-03 Creation
% 27-11-03 Adds defined time function. Modifies calling prototype.
% 14-07-04 Fix small bug in limits of kx and ky
  
% check parameters
% ----------------

% check input arguments
if (nargin<1)||(nargin>2)
  error('Incorrect number of input arguments');
end

dim = [param.nx, param.ny, param.nt];

% prepare internal variables
% --------------------------
nx_min = -floor((dim(1)-1)/2)-1;
nx_max = ceil((dim(1)-1)/2)-1;
ny_min = -floor((dim(2)-1)/2)-1;
ny_max = ceil((dim(2)-1)/2)-1;
kx = (nx_min:nx_max)'*ones([1 dim(2)]);     %'
ky = ones([dim(1) 1])*(ny_min:ny_max);

% generates the k-space data for each frame
% -----------------------------------------
data = zeros(dim);
for time = 1:dim(3)
  frame = zeros(dim(1:2));
  % First I'll progrma the static part     %'
  frame = frame + 0.55*ephantom_elipse(kx,ky,...   % One on the bottom-left
				      [20 -10 24 18]/100,... % [center x, center y, size x, size y]
				      [0 0 0 0]/100,...
				      tfun(time),param);
  frame = frame + 0.5*ephantom_elipse(kx,ky,... % The rectangle
					  [-34 -3.5 4 5]/100,...
					  [0 0 0 0]/100,...
					  tfun(time),param);
  frame = frame + 0.5*ephantom_elipse(kx,ky,... % The rectangle
					  [-34 3.5 4 5]/100,...
					  [0 0 0 0]/100,...
					  tfun(time),param);
  frame = frame + 0.4*ephantom_elipse(kx,ky,... % The rectangle
					  [-29 0 4 5]/100,...
					  [0 0 0 0]/100,...
					  tfun(time),param);
%   frame = frame + 0.6*ephantom_elipse(kx,ky,... % The rectangle
% 					  [-32 0 8 12]/100,...
% 					  [0 0 0 0]/100,...
% 					  tfun(time),param);
 
% Now the moving bits
  frame = frame + 0.9*ephantom_elipse(kx,ky,...    % The outer ring
				   [0 0 90 86]/100,... % [center x, center y, size x, size y]
				   [0 2 0 2]/100,...   % [move x, move y, stretch x, stretch y]
 				   tfun(time),param);
  frame = frame - 0.7*ephantom_elipse(kx,ky,...  % The inner background
				   [0 0 80 76]/100,...
				   [0 2 0 2]/100,...
 				   tfun(time),param);
  frame = frame + 0.4*ephantom_elipse(kx,ky,...    % The big ellipse that carries the other two
				   [-10 -3.25 28 38]/100,...
   				   [0 3.5 0 3]/100,... % [0 2.5 0 7]/100,...
				   tfun(time),param);
  frame = frame + 0.24*ephantom_elipse(kx,ky,...
 				   [-10 -12.05 12 12]/100,...
 				   [0 3.5 0 3]/100,...
 				   tfun(time),param);
  frame = frame + 0.30*ephantom_elipse(kx,ky,...
 				   [-10 6 15 9]/100,...
 				   [0 1.85 0 2]/100,...
 				   tfun(time),param);
  frame = frame + 0.35*ephantom_elipse(kx,ky,...
 				   [20.5 10 20 12]/100,...
 				   [-1 0 -5 0]/100,...
 				   tfun(time),param);
 data(:,:,time) = frame;
end

return

