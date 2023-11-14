classdef GAMP4DpMRIOpts
    %GAMP4DPMRIOPTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        data;
        samp;
        maps;
        compute;
        uniformVar;
        GAMPopt;
        wvar;
        lambda;
        precision;
    end
    
    methods
        % Constructor MEthod
        function obj = GAMP4DpMRIOpts()
            obj.data = [];
            obj.samp = [];
            obj.maps = [];
            obj.compute = 'mat';
            obj.uniformVar = 1;
            obj.precision = 'double';
            obj.compute = 'mat';
            obj.GAMPopt = GampOpt();
            obj.GAMPopt.nit = 40;
            obj.GAMPopt.stepWindow = 0;
            obj.GAMPopt.step = 0.1;
            obj.GAMPopt.verbose = 1;
            obj.GAMPopt.stepIncr = 1.1;     % Step size Increase
            obj.GAMPopt.stepTol = 0.005;    % Min Step size befor stopping
            obj.GAMPopt.tol = 1e-4;         % Stopping Tolerance
            
            obj.wvar = [];
            obj.lambda = [];
        end
        
    end
end

