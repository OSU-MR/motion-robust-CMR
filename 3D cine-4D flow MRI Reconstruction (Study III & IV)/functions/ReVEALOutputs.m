classdef ReVEALOutputs
    %REVEALOUTPUTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        xHatb;
        xHatx;
        xHaty;
        xHatz;
        thetaX;
        thetaY;
        thetaZ;
        thetaMAPX;
        thetaMAPY;
        thetaMAPZ;
        vX;
        vY;
        vZ;
        sensMaps;
        scanParam
        raw_header
    end
    
    methods
        function obj = ReVEALOutputs()
            obj.xHatb = [];
            obj.xHatx = [];
            obj.xHaty = [];
            obj.xHatz = [];
            obj.thetaX = [];
            obj.thetaY = [];
            obj.thetaZ = [];
            obj.thetaMAPX = [];
            obj.thetaMAPY = [];
            obj.thetaMAPZ = [];
            obj.vX = [];
            obj.vY = [];
            obj.vZ = [];
            obj.sensMaps = [];
            obj.scanParam = [];
            obj.raw_header = [];
        end
    end
    
end

