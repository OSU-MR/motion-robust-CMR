function [ path ] = getPath( fname )
%GETPATH Summary of this function goes here
%   Detailed explanation goes here

path = which(fname);
path = path(1:end-length(fname));
end

