function [NetFlow, valveRegurgFrac] = flow_quant(phase_cine, contour, opt)

contourMask = zeros(size(phase_cine));

for i = 1 : size(phase_cine, 3)    
    contourMask(:,:,i) = poly2mask(contour(:, 2, i), contour(:, 1, i), size(phase_cine, 1), size(phase_cine, 2));
end

for k = 1 : size(phase_cine, 3)

    phase_cine_temp = phase_cine(:,:,k); 
    
    contourMask_temp = contourMask(:,:,k);
    
    contourMask_temp = contourMask_temp(:);
    
    vel(k,1) = mean(phase_cine_temp(find(contourMask_temp ==1)));
    
    vel(k,1) = vel(k,1)*opt.acqParam.VENC; 
 
    area(k,1) = sum(contourMask_temp ==1); 
    
end


area = area .* ((opt.acqParam.SpatialResolutionX)*(opt.acqParam.SpatialResolutionY))/100; 

delta_TT = opt.acqParam.TT(1,2:size(opt.acqParam.TT, 2))-opt.acqParam.TT(1,1:size(opt.acqParam.TT, 2)-1); 

delta_TT(1, size(delta_TT, 2)+1) = delta_TT(1, size(delta_TT, 2));

flow = vel .* area ;

NetFlow = sum(flow.*delta_TT')/1000; 

% Find valve regurgitation fraction

% Check sign of flow
if mean(flow) >= 0 
    indRegurg = find(flow < 0);
elseif mean(flow) < 0
    indRegurg = find(flow > 0);
end


valveRegurgFrac = -100*(sum(flow(indRegurg).*delta_TT(indRegurg)')/1000) / NetFlow;

end

