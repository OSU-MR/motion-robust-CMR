function [ lambda_vec ] = setLambdatPCA( sizes,lambda )
% [ lambda_vec ] = setLambda( sizes,lambda )
%
% setLambda creates a vector of regularization parameters the size of the
%   temporal PCA decomposition
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
% Last update:  3/14/2018
%**************************************************************************


lambda_vec = ones(sizes(1)*sizes(2),sizes(3));
lambda_vec(:,end) = 0.1*lambda;
lambda_vec(:,1:end-1) = lambda;
lambda_vec = lambda_vec(:);




