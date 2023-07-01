function [ mask ] = readContour_segment( cName, imsize )
% readContour_segment.m [FUNCTION]
% Reads contours generated using Segment by Medviso

% Load the segment data
cIn = load(cName);

% use the first ROI that is complete
for ind = size(cIn.setstruct(1).Roi,2)
    if ~sum(sum(sum(isnan(cIn.setstruct(1).Roi(ind).X))))
        contourx = cIn.setstruct(1).Roi(ind).X;
        contoury = cIn.setstruct(1).Roi(ind).Y;
    end
end

C = zeros(size(contourx,1),2,size(contourx,2));
for i=1:size(contourx,1)
    for k = 1:size(contourx,2)
        C(i,1,k) = contourx(i,k);
        C(i,2,k) = contoury(i,k);
    end
end

% create a mask fromt the contour
mask = zeros([imsize,size(C,3)]);
for ind = 1:size(C,3)
    tmp = C(:,2,ind);
    if sum(isnan(tmp(:)))
        keyboard
    end
    mask(:,:,ind) = poly2mask(C(:,2,ind),C(:,1,ind),imsize(1),imsize(2));
end

end
