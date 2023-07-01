function [S0] = getInitialMask(magMask, opt )

% Initialize exclusion operator, S
S0 = zeros(size(magMask));

switch opt.initialMaskLoc
    case "mid"
        
        if size(S0, 3) > 1
            mid1 = floor(size(magMask,2)/2)+1; % Center
            mid2 = floor(size(magMask,3)/2)+1;
            extent1 = floor((size(magMask,2)*opt.midFOVFrac)/2);
            extent2 = floor((size(magMask,3)*0.65)/2);
            S0(:,mid1-extent1:mid1+extent1,mid2-extent2:mid2+extent2) = 1;
            S0 = logical(S0);
        else
            switch opt.acqParam.InPlanePhaseEncodingDirection
                case 'ROW'
                    mid = floor(size(magMask,2)/2)+1; % Center
                    extent = floor((size(magMask,2)*opt.midFOVFrac)/2);
                    S0(:,mid-extent:mid+extent) = 1;
                    S0 = logical(S0);
                case 'COL'
                    mid = floor(size(magMask,1)/2)+1; % Center
                    extent = floor((size(magMask,1)*opt.midFOVFrac)/2);
                    S0(mid-extent:mid+extent,:) = 1;
                    S0 = logical(S0);
            end
        end
                   
        
        
        
        
        
        
        
        
        
        
    case "sideLo"
        switch opt.acqParam.InPlanePhaseEncodingDirection
            case 'ROW' % Phase encode across rows in image
                extent = floor((size(magMask,2)*opt.midFOVFrac));
                S0(:,1:extent) = 1;
                S0 = logical(S0);
            case 'COL' % Phase encode across columns in image
                extent = floor((size(magMask,1)*opt.midFOVFrac));
                S0(1:extent,:) = 1;
                S0 = logical(S0);
        end
    case "sideHi"
        switch opt.acqParam.InPlanePhaseEncodingDirection
            case 'ROW' % Phase encode across rows in image
                extent = floor((size(magMask,2)*opt.midFOVFrac));
                S0(:,end-extent:end) = 1;
                S0 = logical(S0);
            case 'COL' % Phase encode across columns in image
                extent = floor((size(magMask,1)*opt.midFOVFrac));
                S0(end-extent:end,:) = 1;
                S0 = logical(S0);
        end
end

