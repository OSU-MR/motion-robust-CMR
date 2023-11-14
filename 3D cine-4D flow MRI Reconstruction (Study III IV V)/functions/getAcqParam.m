function [acqParam] = getAcqParam(acqParam_mat) 

% FOV
acqParam.FOV = char(acqParam_mat(1).Private_0051_100c);

% Spatial resolution
acqParam.SpatialResolutionX = str2double(acqParam.FOV(5:7))/double(acqParam_mat(1).Rows); % spatial resolution
acqParam.SpatialResolutionY = str2double(acqParam.FOV(9:11))/double(acqParam_mat(1).Columns); 

% Venc
VENC = char(acqParam_mat(1).Private_0051_1014);
acqParam.VENC = str2double(VENC(2:4));
if isnan(acqParam.VENC) % if venc is less than 100
    acqParam.VENC = str2double(VENC(2:3));
end
                
% In plane phase encoding direction                    
acqParam.InPlanePhaseEncodingDirection = char(acqParam_mat(1).InPlanePhaseEncodingDirection);
                    
% Trigger Time                        
for j = 1 : length(acqParam_mat)
    
    acqParam.TT(j) = acqParam_mat(j).TriggerTime;
    
end                        
end

