
function create_GIF( x,file_name, frame_rate, clip_factor )
% create_GIF( x,file_name, frame_rate, clip_factor )
% Creates a GIF from a 3D Array
% Inputs:
% x:            Real Valued 3D array
%
% file_name:    File name for the movie
%
% frame_rate:   The frame rate of the movie in frames per second
%
% clip_factor:  The amount to clip on the high end of the image.  e.g. a
%               clip factor of 2 will plot from 0-1/2*max(x(:))
%
%**************************************************************************
% The Ohio State University
% Written by:   Adam Rich 
% Email:        rich.178@osu.edu
% Last update:  4/26/2014
%**************************************************************************

file_name = [file_name,'.gif'];
frame_rate = 1/frame_rate; % Delay between GIF frames, it controls the speed of the movie
 
s = size(x);
x = x/max(x(:));
out = uint8(clip_factor*255*x);
% out = x;

% 3D Image
if ndims(x) ==3
    imwrite(out(:,:,1),file_name, 'DelayTime', frame_rate, 'LoopCount', Inf)
    for i = 2:s(3)-1
        imwrite(out(:,:,i), file_name,'WriteMode','append','DelayTime', frame_rate)  
    end
    imwrite(out(:,:,s(3)), file_name,'WriteMode','append','DelayTime', frame_rate);
    
% 4D image
elseif ndims(x) ==4
    
    n_rows = ceil(size(out,3)/4);
    ind = 1;
    outFlat = [];
    for j = 1:n_rows
        rowTMP = [];
        for k = 1:4
            if ind <=size(out,3)
                rowTMP = cat(2,rowTMP,out(:,:,ind,:));
            else
                rowTMP = cat(2,rowTMP,zeros(size(out(:,:,1,:))) );
            end
            ind = ind+1;
        end
        outFlat = cat(1,outFlat,rowTMP);
    end
    outFlat = squeeze(outFlat);
    imwrite(outFlat(:,:,1),file_name, 'DelayTime', frame_rate, 'LoopCount', Inf)
    for i = 2:s(4)-1
        imwrite(outFlat(:,:,i), file_name,'WriteMode','append','DelayTime', frame_rate)  
    end
    imwrite(outFlat(:,:,s(4)), file_name,'WriteMode','append','DelayTime', frame_rate);
end

end
    

