function [ contour ] = readContourSegment( cName )
% readContourSegment
% Reads contours generated using Segment by Medviso

cIn = load(cName);
contourx = cIn.setstruct(1).Roi.X;
contoury = cIn.setstruct(1).Roi.Y;
for i=1:size(contourx,1)
    for k = 1:size(contourx,2)
        contour(i,1,k) = contourx(i,k);
        contour(i,2,k) = contoury(i,k);
    end
end
end

