function [ lambda_vec ] = setLambdaWavFour( sizes,lambda )
% [ lambda_vec ] = setLambda( sizes,lambda )
%
% setLambda creates a vector of regularization parameters the size of the
% 3D non-decimated wavelet transform with a different values for each
% subband
%
%
% Inputs:
%   sizes:  The size of the 3 Dimensional Image
%
%   lambda: The base regularization parameter
%
% Output
%   lambda_vec: A Nx1 vector the size of the 3D nodecimated wavelet
%               transform with different regularization for each subband
%
%**************************************************************************
% The Ohio State University
% Written by:   Adam Rich 
% Email:        rich.178@osu.edu
% Last update:  10/09/2014
%**************************************************************************

if length(sizes)==3
    k = 1.2;
    lambda_band = zeros(8,1);
    lambda_band(1) = 5*lambda;
    lambda_band(2) = 1*lambda;
    lambda_band(3) = 1*lambda;
    lambda_band(4) = 2*lambda;
    lambda_band(5) = 1*k*lambda;
    lambda_band(6) = 2*k*lambda;
    lambda_band(7) = 2*k*lambda;
    lambda_band(8) = 3*k*lambda;

    % High Temporal Weights
%     k = 1.2;
%     lambda_band = zeros(8,1);
%     lambda_band(1) = 0.01*lambda;
%     lambda_band(2) = 1*lambda;
%     lambda_band(3) = 1*lambda;
%     lambda_band(4) = 4*lambda;
%     lambda_band(5) = 1*k*lambda;
%     lambda_band(6) = 4*k*lambda;
%     lambda_band(7) = 4*k*lambda;
%     lambda_band(8) = 9*k*lambda;

    lambda_vec = zeros([sizes,8]);
    for ind = 1:8
       lambda_vec(:,:,:,ind) = lambda_band(ind);
    end
    lambda_vec(:,:,1,1) = 0.001*lambda_vec(:,:,1,1);
    
elseif length(sizes)==4

%     warning('lower t lambda')
    % High Temporal Weights
    k = 10;
    lambda_band = zeros(8,1);
    lambda_band(1) = 0.001*lambda;
    lambda_band(2) = 1*lambda;
    lambda_band(3) = 1*lambda;
    lambda_band(4) = 1*lambda;
    lambda_band(5) = 1*lambda;
    lambda_band(6) = 1*lambda;
    lambda_band(7) = 1*lambda;
    lambda_band(8) = 2*lambda;
    lambda_band(9) = k*lambda;
    lambda_band(10) = k*lambda;
    lambda_band(11) = k*lambda;
    lambda_band(12) = k*lambda;
    lambda_band(13) = k*lambda;
    lambda_band(14) = k*lambda;
    lambda_band(15) = k*lambda;
    lambda_band(16) = 1.5*k*lambda;

    lambda_vec = zeros([sizes,16]);
    for ind = 1:16
       lambda_vec(:,:,:,:,ind) = lambda_band(ind);
    end
else 
    error('error uknow dimension size')
end
lambda_vec = lambda_vec(:);
end

