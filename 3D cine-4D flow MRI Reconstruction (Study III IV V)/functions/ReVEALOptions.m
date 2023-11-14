classdef ReVEALOptions
% The ReVEAL Class to store algorithm options
%
% Methods:
%
% Properties
%   data
%
%   samplingPatterns
%
%   maxwellCorrections
%
%   options
%
%   parameters
%
%   scan_parameters
%
%**************************************************************************
% The Ohio State University
% Written by:   Adam Rich
% Written on:   3/2/2015
% Last update:  3/2/2015
%**************************************************************************
    
    properties
        sensParamFil;   % Sensitivity map blur parameter
        is1Dir;         % Flag for 1 Directional velocity measurements
        isPlanar;       % Flag for planar imaging
        ReVEALOpts;
        GAMPOpt;
        GAMPOptB;
        GAMPOptX;
        GAMPOptY;
        GAMPOptZ;
        SparseTrans;
    end
    
    methods
        % Constructor
        function obj = ReVEALOptions
            obj.is1Dir = [];
            obj.isPlanar = [];
            obj.sensParamFil = 3;
            
            % ReVEAL options
            obj.ReVEALOpts.nit = 8;         % Number of iterations
            obj.ReVEALOpts.tol = 0.0001;    % Stopping Tolerance
            obj.ReVEALOpts.L1 = 1;          % Use L1 output estimator
            obj.ReVEALOpts.gamma = 0.20;    % Probability of Velocity
            obj.ReVEALOpts.lambda0 = 2;     % Regularization on L t-bands
            obj.ReVEALOpts.sigmaSq = [];    % Variance connecting the graphs
            obj.ReVEALOpts.wvar = 1e-10;    % Noise Variance
            obj.ReVEALOpts.uniform_var = 0; % Use Uniform Variance Assumption
            obj.ReVEALOpts.MAP = 0;         % Use MAP or MMSE
            obj.ReVEALOpts.num_coils = 12;  % Number of coils
            obj.ReVEALOpts.compute = 'mat'; % Use GPU acceleration
            obj.ReVEALOpts.precision = 'single'; % Use double or single precision
            obj.ReVEALOpts.data_size = [];
            obj.ReVEALOpts.bg_corr = 0;     % background phase correction: (0) no, (1) yes
            
            obj.GAMPOpt = GampOpt();        % generate default options
            obj.GAMPOpt.verbose = 0;        % Report step size etc
            obj.GAMPOpt.nit = 10;           % Number of Iterations
            obj.GAMPOpt.stepIncr = 1.1;     % Step size Increase
            obj.GAMPOpt.step = 0.1;         % Initial Step size
            obj.GAMPOpt.stepTol = 0.005;    % Min Step size befor stopping
            obj.GAMPOpt.stepWindow = 0;     % Stepsize Window
            obj.GAMPOpt.tol = 1e-4;         % Stopping Tolerance
            obj.GAMPOpt.legacyOut = 0;      % Legacy output format
            
            % GAMP Options
            obj.GAMPOptB = GampOpt();        % generate default options
            obj.GAMPOptB.verbose = 0;        % Report step size etc
            obj.GAMPOptB.nit = 10;           % Number of Iterations
            obj.GAMPOptB.stepIncr = 1.1;     % Step size Increase
            obj.GAMPOptB.step = 0.1;         % Initial Step size
            obj.GAMPOptB.stepTol = 0.005;    % Min Step size befor stopping
            obj.GAMPOptB.stepWindow = 0;     % Stepsize Window
            obj.GAMPOptB.tol = 1e-4;         % Stopping Tolerance
            obj.GAMPOptB.legacyOut = 0;      % Legacy output format
            
            % GAMP Options
            obj.GAMPOptX = GampOpt();        % generate default options
            obj.GAMPOptX.verbose = 0;        % Report step size etc
            obj.GAMPOptX.nit = 10;           % Number of Iterations
            obj.GAMPOptX.stepIncr = 1.1;     % Step size Increase
            obj.GAMPOptX.step = 0.1;         % Initial Step size
            obj.GAMPOptX.stepTol = 0.005;    % Min Step size befor stopping
            obj.GAMPOptX.stepWindow = 0;     % Stepsize Window
            obj.GAMPOptX.tol = 1e-4;         % Stopping Tolerance
            obj.GAMPOptX.legacyOut = 0;      % Legacy output format
            
            % GAMP Options
            obj.GAMPOptY = GampOpt();        % generate default options
            obj.GAMPOptY.verbose = 0;        % Report step size etc
            obj.GAMPOptY.nit = 10;           % Number of Iterations
            obj.GAMPOptY.stepIncr = 1.1;     % Step size Increase
            obj.GAMPOptY.step = 0.1;         % Initial Step size
            obj.GAMPOptY.stepTol = 0.005;    % Min Step size befor stopping
            obj.GAMPOptY.stepWindow = 0;     % Stepsize Window
            obj.GAMPOptY.tol = 1e-4;         % Stopping Tolerance
            obj.GAMPOptY.legacyOut = 0;      % Legacy output format
            
            % GAMP Options
            obj.GAMPOptZ = GampOpt();        % generate default options
            obj.GAMPOptZ.verbose = 0;        % Report step size etc
            obj.GAMPOptZ.nit = 10;           % Number of Iterations
            obj.GAMPOptZ.stepIncr = 1.1;     % Step size Increase
            obj.GAMPOptZ.step = 0.1;         % Initial Step size
            obj.GAMPOptZ.stepTol = 0.005;    % Min Step size befor stopping
            obj.GAMPOptZ.stepWindow = 0;     % Stepsize Window
            obj.GAMPOptZ.tol = 1e-4;         % Stopping Tolerance
            obj.GAMPOptZ.legacyOut = 0;      % Legacy output format
            
            % Sparse Trans Options
            obj.SparseTrans.type = 'nd-dwt';        % Discrete Wavelet Transform
            obj.SparseTrans.level = 1;              % DWT levels
            obj.SparseTrans.wname = 'db1';
%            obj.SparseTrans.wname = {'db1','db1','db1'};  % DWT wavelet

        end
        
        % Make a copy of a handle object.
        function new = copy(obj,handle)
            % Instantiate new object of the same class.
            new = feval(class(handle));

            % Copy all non-hidden properties.
            p = properties(handle);
            for i = 1:length(p)
                new.(p{i}) = handle.(p{i});
            end
            
            new.GAMPOptB = obj.GAMPOptZ.copy(obj.GAMPOptB);
            new.GAMPOptX = obj.GAMPOptZ.copy(obj.GAMPOptX);
            new.GAMPOptY = obj.GAMPOptZ.copy(obj.GAMPOptY);
            new.GAMPOptZ = obj.GAMPOptZ.copy(obj.GAMPOptZ);
        end
    end
    
end

