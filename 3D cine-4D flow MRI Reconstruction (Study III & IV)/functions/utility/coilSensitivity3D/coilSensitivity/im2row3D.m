function res = im2row3D(im, winSize)
%res = im2row(im, winSize)
[sx,sy,sz,nc] = size(im);

res = zeros((sx-winSize(1)+1)*(sy-winSize(2)+1)*(sz-winSize(3)+1),prod(winSize),nc);
count=0;
for z=1:winSize(3)
for y=1:winSize(2)
    for x=1:winSize(1)
        count = count+1;
        res(:,count,:) = reshape(im(x:sx-winSize(1)+x,y:sy-winSize(2)+y,z:sz-winSize(3)+z,:),...
            (sx-winSize(1)+1)*(sy-winSize(2)+1)*(sz-winSize(3)+1),1,nc);
    end
end
end
