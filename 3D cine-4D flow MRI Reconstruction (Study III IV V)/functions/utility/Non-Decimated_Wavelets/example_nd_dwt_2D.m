clear all;
close all;
addpath('./Functions');
% Load data and set parameters
x = double(imread('cameraman.tif'));   % Load a noisy test singal
level = 1;                            % Set the Level of decomposition
wnames = {'db1','db4'};               % Set the wavelet used 
perserve_l2_norm = true;              % Choose wether to preserve the l2 norm or not (optional) 

% Make the signal complex
x = x + 1j*ones(256,256);

% Initialize the nd_dwt_2D class
nddwt = nd_dwt_2D(wnames,size(x),'pres_l2_norm',perserve_l2_norm);

% Perform a multilevel wavelet decomposition
x_trans = nddwt.dec(x,level);

% Perform a multilevel wavelet reconstruction
x_recon = nddwt.rec(x_trans);

% Plot a single level of wavelet decomposition
figure(1)
clf
coeffs = [x_trans(:,:,1)/max(max(x_trans(:,:,1))),...
    x_trans(:,:,2)/max(max(x_trans(:,:,2)));...
    x_trans(:,:,3)/max(max(x_trans(:,:,3))),...
    x_trans(:,:,4)/max(max(x_trans(:,:,4)))];
imagesc(abs(coeffs))
colormap(gray)
axis image
title('One Level of the Wavelet Decomposition')

% Plot the original and reconstructed signal
figure(2)
clf
subplot(211)
imagesc(abs(x))
title('Original Signal')
colormap(gray)
axis image

subplot(212)
imagesc(abs(x_recon))
title('Reconstructed Signal')
colormap(gray)
axis image

% Print some stats to the screen
fprintf('Energy in signal domain = %s \t Energy in Wavelet Domain = %s\n',norm(x(:)), norm((x_trans(:))))
fprintf('Absolute Max Reconstruction Error = %s\n',max(abs((x_recon(:))-x(:))))