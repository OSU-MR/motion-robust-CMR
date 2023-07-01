function [S, cI] = coilSen_walsh(kdata, samp,param)

fft2Scal = 1/sqrt(size(kdata,1)*size(kdata,2));

if ndims(kdata) == 4
    kdata_avg = sum(kdata,4)./(repmat(sum(samp,3),[1,1,size(kdata,3)])+eps);
%     kdata_avg1 = sum(kdata(:,:,:,1:end/2),4)    ./(repmat(sum(samp(:,:,1:end/2),3),    [1,1,size(kdata,3)])+eps);
%     kdata_avg2 = sum(kdata(:,:,:,end/2+1:end),4)./(repmat(sum(samp(:,:,end/2+1:end),3),[1,1,size(kdata,3)])+eps);
%     msk = repmat(logical(sum(samp(:,:,1:end/2),3)),    [1,1,size(kdata,3)]) .* ...
%           repmat(logical(sum(samp(:,:,end/2+1:end),3)),[1,1,size(kdata,3)]);
%     figure; imagesc(msk(:,:,1));
%     kdata_avg = (kdata_avg1 + kdata_avg2)/2.*msk;
%     warning('temporarily changed, comment lines 8 - 13 and uncomment line 7 instead');
elseif ndims(kdata) == 2 || ndims(kdata) == 3
    kdata_avg = kdata;
end

    I = kdata2im(kdata_avg,logical(sum(samp,3)), fft2Scal, param.ACSco); % Find low-res images
    [cI,S] = WalshCoilCombine(I, param);
    [cI,S] = sensCorrect(cI, S, param.avgPhs); % Remove phase from time average



% Find low-res image from the data 'd'
function img = kdata2im(d, s, fft2Scal, co)
[CalibSize, ~] = getCalibSize(s);
% CalibSize = [40,40];
fil = g2d(size(d(:,:,1)), CalibSize, co);
d = d.* repmat(fil, [1,1,size(d,3)]);
img = 1/fft2Scal * fftshift(fftshift(ifft2(ifftshift(ifftshift(d,1),2)),2),1);
% cctest = 1;



function F = g2d(S1, S2, co)
[x,y] = ndgrid(1:S1(1), 1:S1(2));
cntr=floor(S1/2)+1;
sigx = S2(1)*co(1);
sigy = S2(2)*co(2);
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

