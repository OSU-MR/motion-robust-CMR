clear all;
close all;

% Load data and set parameters
sizes = [32,32,16,16];              % Set the size of the 3D object
level = 2;                          % Set the Level of decomposition
wnames = {'db1','db3','db3','db5'}; % Set the wavelet used 
perserve_l2_norm = true;            % Choose wether to preserve the l2 norm or not (optional)    

% Make a random 4D complex signal
x = randn(sizes) + 1j*randn(sizes);

% Initialize the nd_dwt_4D class
nddwt = nd_dwt_4D(wnames,size(x),'pres_l2_norm',perserve_l2_norm);

% Perform a multilevel wavelet decomposition
x_trans = nddwt.dec(x,level);

% Perform a multilevel wavelet reconstruction
x_recon = nddwt.rec(x_trans);

% Print some stats to the screen
fprintf('Energy in signal domain = %s \t Energy in Wavelet Domain = %s\n',norm(x(:)), norm((x_trans(:))))
fprintf('Absolute Max Reconstruction Error = %s\n',max(abs((x_recon(:))-x(:))))
