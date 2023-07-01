function [S, cI] = coilSen(kdata, samp, param)

mthd = param.mthd;
fft2Scal = 1/sqrt(size(kdata,1)*size(kdata,2));

if ndims(kdata) == 4
    kdata_avg = sum(kdata,4)./(repmat(sum(samp,3),[1,1,size(kdata,3)])+eps);
elseif ndims(kdata) == 3
    kdata_avg = kdata;
end

if mthd == 1 % espirit
    S = espirit_sens(kdata,samp,3,0.01);
%     S = flipdim(flipdim(circshift(S,[-1,-1,0]),2),1);
        
elseif mthd ==2 % espirit after time-averaging the data
    S = espirit_sens(kdata_avg,logical(sum(samp,3)),3,0.01);
%     S = flipdim(flipdim(circshift(S,[-1,-1,0]),2),1);

elseif mthd == 3 % Walsh
    I = kdata2im(kdata_avg,logical(sum(samp,3)), fft2Scal); % Find low-res images
%     figure; imagesc(abs(I(:,:,2))); axis('image');
    [cI,S] = WalshCoilCombine(I, param);
%     S = circshift(flipdim(flipdim(S,2),1),[1,1,0]);
else
    error('Undefined method for coil sensisitive estimation');
end


% Find low-res image from the data 'd'
function img = kdata2im(d, s, fft2Scal)
[CalibSize, ~] = getCalibSize(s);
fil = g2d(size(d(:,:,1)), CalibSize);
d = d.* repmat(fil, [1,1,size(d,3)]);
img = 1/fft2Scal * fftshift(fftshift(ifft2(ifftshift(ifftshift(d,1),2)),2),1);



function F = g2d(S1, S2)
[x,y] = ndgrid(1:S1(1), 1:S1(2));
cntr=floor(S1/2)+1;
sigx = S2(1)/2;
sigy = S2(2)/2;
F = 1/(2*pi*sigx*sigy)*exp(-((x-cntr(1)).^2)./(2*sigx^2) - ((y-cntr(2)).^2)./(2*sigy^2));
F = F/max(F(:));
% figure; imagesc(F); axis('image');




% function [S] = coilSen(kdata, samp, param)



% if param.mthd == 1 % espirit
%     [S,~] = espirit_recon(kdata,samp);
% %     S = flipdim(flipdim(circshift(S,[-1,-1,0]),2),1);
%         
% elseif param.mthd ==2 % espirit after time-averaging the data
%     [S,~] = espirit_recon(sum(kdata,4)./(repmat(sum(samp,3),[1,1,size(kdata,3)])+eps),logical(sum(samp,3)));
% %     S = flipdim(flipdim(circshift(S,[-1,-1,0]),2),1);
% 
% elseif param.mthd == 3 % Walsh
%     I = kdata2im(sum(kdata,4)./(repmat(sum(samp,3),[1,1,size(kdata,3)])+eps),logical(sum(samp,3))); % Find low-res images
% %     figure; imagesc(abs(I(:,:,1))); axis('image');
%     [~,S] = WalshCoilCombine(I, param);
%     S = circshift(flipdim(flipdim(S,2),1),[1,1,0]);
% else
%     error('Undefined method for coil sensisitive estimation');
% end
% 
% 
% % Find low-res image from the data 'd'
% function img = kdata2im(d, s)
% [CalibSize, ~] = getCalibSize(s);
% fil = g2d(size(d(:,:,1)), CalibSize);
% d = d.* repmat(fil, [1,1,size(d,3)]);
% img = fftshift(ifft2(ifftshift(d)));
% 
% 
% 
% function F = g2d(S1, S2)
% [x,y] = ndgrid(1:S1(1), 1:S1(2));
% cntr=floor(S1/2)+1;
% sigx = S2(1)/4;
% sigy = S2(2)/4;
% F = 1/(2*pi*sigx*sigy)*exp(-((x-cntr(1)).^2)./(2*sigx^2) - ((y-cntr(2)).^2)./(2*sigy^2));
% F = F/max(F(:));
% % figure; imagesc(F); axis('image');

