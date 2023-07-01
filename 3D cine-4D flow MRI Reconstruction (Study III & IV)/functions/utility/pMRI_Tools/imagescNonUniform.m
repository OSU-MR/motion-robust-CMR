function imagescNonUniform( image, x_dim_size, y_dim_size,clim)
% imagescNonUniform( image, x_dim_size, y_dim_size,clim)
%
% imagescNonUniform behaves like imagesc for non-uniformly sized pixels
%
%
% Inputs:
%   image:      The image to plot
%
%   x_dim_size: The size of the pixel in the x (horizontal) dimension
%
%   y_dim_size: The size of the pixel in the y (vertical) dimension
%
%   clim:       An optional input.  A size two vector of values to display 
%               in the image. It behaves exactly like CLIM in imagesc.  See
%               "help imagesc"
%
%
%**************************************************************************
% The Ohio State University
% Written by:   Adam Rich 
% Email:        rich.178@osu.edu
% Last update:  9/01/2015
%**************************************************************************

% Check if clim is provided
if nargin < 4
    clim = [min(image(:)),max(image(:))];
elseif length(clim)~=2
    error('clim must be length 2');
end

image(image(:)<clim(1)) = clim(1);
image(image(:)>clim(2)) = clim(2);

% Form the x and y vectors
y = linspace(0,round(size(image,1)*y_dim_size),size(image,1));
x = linspace(0,round(size(image,2)*x_dim_size),size(image,2));
z = 0;

% Form the mesh grid
[x_grid,y_grid,z_grid] = meshgrid(x,y,z);

% Plot
surface(x_grid,y_grid,z_grid,flipud(image),'FaceColor','texturemap','EdgeColor','none')
axis image

end

