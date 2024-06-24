clear all;
close all;

% Load data and set parameters
load noisbloc;           % Load a noisy test singal
level = 3;               % Set the Level of decomposition
wname = 'db3';           % Set the wavelet used 
perserve_l2_norm = true; % Choose wether to preserve the l2 norm or not (optional)    

% Make the signal complex
x = noisbloc.*(ones(size(noisbloc)) + 1j*ones(size(noisbloc)));

% Initialize the nd_dwt_1D class
nddwt = nd_dwt_1D(wname,length(x),'perserve_l2_norm',perserve_l2_norm);

% Perform a multilevel wavelet decomposition
x_trans = nddwt.dec(x,level);

% Perform a multilevel wavelet reconstruction
x_recon = nddwt.rec(x_trans);

% Plot the original and reconstructed signal
figure(1)
clf
subplot(211)
plot(abs(x))
title('Original Signal')

subplot(212)
plot(abs(x_recon))
title('Reconstructed Signal')

% Plot the wavelet coefficients
figure(2)
subplot(level+1,1,1)
plot(abs(x_trans(:,1)))
title('Approximate Coefficients')
for ind = 1:level
    subplot(level+1,1,ind+1)
    plot(abs(x_trans(:,ind+1)))
    title(sprintf('Detail Coefficients level=%d',level-ind+1))
end

% Print some stats to the screen
fprintf('Energy in signal domain = %s \t Energy in Wavelet Domain = %s\n',norm(x), norm((x_trans(:))))
fprintf('Absolute Max Reconstruction Error = %s\n',max(abs((x_recon(:))-x(:))))