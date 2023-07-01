classdef ReVEALData
% The ReVEALData Class is used to store PC-MRI during ReVEAL reconstruction
%
% Methods:
%
% Properties
% 	Yb:		Background encoded data
%
% 	Yx:		x direction velocity encoded data
%
% 	Yy:		y direction velocity encoded data
%
% 	Yz:		z direction velocity encoded data
%
% 	sampB:	b direction velocity encoded sampling pattern
%
% 	sampX:	x direction velocity encoded sampling pattern
%
% 	sampY: 	y direction velocity encoded sampling pattern
%
% 	sampZ:	z direction velocity encoded sampling pattern
%
%	maxwellCorrX	Maxwell Correction for x velocity encoding
%
%	maxwellCorrY:	Maxwell Correction for y velocity encoding
%
%	maxwellCorrZ:	Maxwell Correction for z velocity encoding
%
%   sensMaps;       Estimated Sensitivity Maps
%
%   x0;             Time averaged image for initialization
%
%	scanParam:		MRI Scan Parameters Tr, Te, etc
%
%**************************************************************************
% The Ohio State University
% Written by:   Adam Rich 
% Written on:   2/28/2015
% Last update:  3/2/2015
%**************************************************************************

	properties
		Yb;				% Background encoded data
		Yx;				% x direction velocity encoded data
		Yy;				% y direction velocity encoded data
		Yz;				% z direction velocity encoded data
		sampB;			% b direction velocity encoded sampling pattern
		sampX;			% x direction velocity encoded sampling pattern
		sampY; 			% y direction velocity encoded sampling pattern
		sampZ;			% z direction velocity encoded sampling pattern
		maxwellCorrX;	% Maxwell Correction for x velocity encoding
		maxwellCorrY;	% Maxwell Correction for y velocity encoding
		maxwellCorrZ;	% Maxwell Correction for z velocity encoding
        sensMaps;       % Estimated Sensitivity Maps
        x0;             % Time averaged image for initialization
		scanParam;		% MRI Scan Parameters Tr, Te, etc
        noiseChannel;   % The last channel after coil compression
        weightsB;       % Weights for nonlinear least squares
        weightsX;       % Weights for nonlinear least squares
        weightsY;       % Weights for nonlinear least squares
        weightsZ;       % Weights for nonlinear least squares
        chunk;
        chunkInds;
        chunkOutputs;
	end
	
	methods
		% Constructor Method
		function obj = ReVEALData()
			obj.Yb = [];
			obj.Yx = [];
			obj.Yy = [];
			obj.Yz = [];
			obj.sampB = [];
			obj.sampX = [];
			obj.sampY = [];
			obj.sampZ = [];
			obj.maxwellCorrX = [];
			obj.maxwellCorrY = [];
			obj.maxwellCorrZ = [];
            obj.sensMaps = [];
			obj.scanParam = [];
            obj.noiseChannel = [];
            obj.weightsB = [];
            obj.weightsX = [];
            obj.weightsY = [];
            obj.weightsZ = [];
            obj.chunk = [];
		end
	end
end