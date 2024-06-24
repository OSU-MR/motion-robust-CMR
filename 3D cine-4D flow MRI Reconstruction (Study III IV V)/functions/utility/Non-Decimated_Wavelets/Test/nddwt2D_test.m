clear;
close all;
clc
% Load data and set parameters
sizes = [264,264];             % Set the size of the 2D object
level = 1;                      % Set the Level of decomposition
wnames = {'db1','db3'};   % Set the wavelet used 
perserve_l2_norm = true;        % Choose wether to preserve the l2 norm or not (optional)    

% Make a random 2D complex signal
x = randn(sizes) + 1j*randn(sizes);

% Initialize the nd_dwt_2D class
nddwt = nd_dwt_2D(wnames,size(x),'pres_l2_norm',perserve_l2_norm,'compute','mat','precision','double');

% Perform a multilevel wavelet decomposition
tic;
x_trans = nddwt.dec(x,level);

% Perform a multilevel wavelet reconstruction
x_recon = nddwt.rec(x_trans);
t1 = toc;

% Print some stats to the screen
fprintf('Matlab time = %s, output Class %s %s \n',num2str(t1),class(x_trans),class(x_recon))
fprintf('Energy in signal domain = %s \t Energy in Wavelet Domain = %s\n',norm(x(:)), norm((x_trans(:))))
fprintf('Absolute Max Reconstruction Error = %s\n\n',max(abs((x_recon(:))-x(:))))

%% Mex Double Computatoin
nddwt = nd_dwt_2D(wnames,size(x),'pres_l2_norm',perserve_l2_norm,'compute','mex','precision','double');

% Perform a multilevel wavelet decomposition
tic;
x_trans = nddwt.dec(x,level);

% Perform a multilevel wavelet reconstruction
x_recon = nddwt.rec(x_trans);
t1 = toc;

% Print some stats to the screen
fprintf('Mex time = %s, output Class %s %s \n',num2str(t1),class(x_trans),class(x_recon))
fprintf('Energy in signal domain = %s \t Energy in Wavelet Domain = %s\n',norm(x(:)), norm((x_trans(:))))
fprintf('Absolute Max Reconstruction Error = %s\n\n',max(abs((x_recon(:))-x(:))))

%% GPU Off Computatoin
nddwt = nd_dwt_2D(wnames,size(x),'pres_l2_norm',perserve_l2_norm,'compute','gpu_off','precision','double');

% Perform a multilevel wavelet decomposition
tic;
x_trans = nddwt.dec(x,level);

% Perform a multilevel wavelet reconstruction
x_recon = nddwt.rec(x_trans);
t1 = toc;

% Print some stats to the screen
fprintf('GPU_off time = %s, output Class %s %s \n',num2str(t1),class(x_trans),class(x_recon))
fprintf('Energy in signal domain = %s \t Energy in Wavelet Domain = %s\n',norm(x(:)), norm((x_trans(:))))
fprintf('Absolute Max Reconstruction Error = %s\n\n',max(abs((x_recon(:))-x(:))))

%% GPU computation
nddwt = nd_dwt_2D(wnames,size(x),'pres_l2_norm',perserve_l2_norm,'compute','gpu','precision','double');
% xGPU = gpuArray(x);
xGPU = x;
% Perform a multilevel wavelet decomposition
tic;
x_trans = nddwt.dec(xGPU,level);

% Perform a multilevel wavelet reconstruction
x_recon = nddwt.rec(x_trans);
t1 = toc;

% Print some stats to the screen
fprintf('GPU time = %s, output Class %s %s %s %s \n',num2str(t1),class(x_trans), class(gather(x_trans)), class(x_recon), class(gather(x_recon)))
fprintf('Energy in signal domain = %s \t Energy in Wavelet Domain = %s\n',norm(x(:)), norm((x_trans(:))))
fprintf('Absolute Max Reconstruction Error = %s\n\n',max(abs((x_recon(:))-x(:))))

%% Matlab Single
nddwt = nd_dwt_2D(wnames,size(x),'pres_l2_norm',perserve_l2_norm,'compute','mat','precision','single');
x = single(x);
% Perform a multilevel wavelet decomposition
tic;
x_trans = nddwt.dec(x,level);

% Perform a multilevel wavelet reconstruction
x_recon = nddwt.rec(x_trans);
t1 = toc;

% Print some stats to the screen
fprintf('Matlab Single time = %s, output Class %s %s \n',num2str(t1),class(x_trans),class(x_recon))
fprintf('Energy in signal domain = %s \t Energy in Wavelet Domain = %s\n',norm(x(:)), norm((x_trans(:))))
fprintf('Absolute Max Reconstruction Error = %s\n\n',max(abs((x_recon(:))-x(:))))

%% GPU_off Single
nddwt = nd_dwt_2D(wnames,size(x),'pres_l2_norm',perserve_l2_norm,'compute','gpu_off','precision','single');
 
% Perform a multilevel wavelet decomposition
tic;
x_trans = nddwt.dec(x,level);

% Perform a multilevel wavelet reconstruction
x_recon = nddwt.rec(x_trans);
t1 = toc;

% Print some stats to the screen
fprintf('GPU_off single time = %s, output Class %s %s \n',num2str(t1),class(x_trans),class(x_recon))
fprintf('Energy in signal domain = %s \t Energy in Wavelet Domain = %s\n',norm(x(:)), norm((x_trans(:))))
fprintf('Absolute Max Reconstruction Error = %s\n\n',max(abs((x_recon(:))-x(:))))

%% Initialize the nd_dwt_2D class
nddwt = nd_dwt_2D(wnames,size(x),'pres_l2_norm',perserve_l2_norm,'compute','gpu','precision','single');
xGPU = gpuArray(x);

% Perform a multilevel wavelet decomposition
tic;
x_trans = nddwt.dec(xGPU,level);

% Perform a multilevel wavelet reconstruction
x_recon = nddwt.rec(x_trans);
t1 = toc;

% Print some stats to the screen
fprintf('GPU single time = %s, output Class %s %s %s %s \n',num2str(t1),class(x_trans), class(gather(x_trans)), class(x_recon), class(gather(x_recon)))
fprintf('Energy in signal domain = %s \t Energy in Wavelet Domain = %s\n',norm(x(:)), norm((x_trans(:))))
fprintf('Absolute Max Reconstruction Error = %s\n\n',max(abs((x_recon(:))-x(:))))

