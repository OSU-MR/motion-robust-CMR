%%
%Input:
%kdata: data in k-space
%       kx, ky, kz, coil, time 
%samp: sampling pattern, same for all coils, so no coil information
%      it is fully sampled in RO direction always.
%      two dimensions less than kdata
%      kx, ky, kz, time
%      type: boolean, true indicating the point was acquired
%param:
%     .mthd  : method
%     .fil   : order of smoothing matrix, default 9
%     .opt   : type of phase correction

%Output:
%S: sensitivity map
%cI: combined Image

% an example of how to use coilSen3D
% param.mthd = 3; % Walsh method
% [S, cI] = coilSen3D(kdata, samp, param);

%%

function [S, cI] = coilSen3D(kdata, samp, param)

mthd = param.mthd;
fft3Scal = 1/sqrt(size(kdata,1)*size(kdata,2)*size(kdata,3));

if ndims(kdata) == 5 % kx, ky, kz, coil, time
    % average over time for each coil
    kdata_avg = sum(kdata,5)./(repmat(sum(samp,4),[1,1,1,size(kdata,4)])+eps);
elseif ndims(kdata) == 4
    kdata_avg = kdata;
end

if mthd == 1 % espirit
    disp('not implemented for 3D.');
    [cI, S] = espirit_sens3D(kdata,samp,6,param);
%     S = flipdim(flipdim(circshift(S,[-1,-1,0]),2),1);
        
elseif mthd ==2 % espirit after time-averaging the data
    disp('not implemented for 3D.');
    % Perform espirit on low-res image to save on memory
    N = [size(kdata_avg,1), size(kdata_avg,2),size(kdata_avg,3)];
    factor = [2,2,2];
    cN = [floor(N(1)/factor(1)),floor(N(2)/factor(2)),floor(N(3)/factor(3))];
    center = [floor(N(1)/2)+1, floor(N(2)/2)+1,floor(N(3)/2)+1];
    kdata_avg_crop = kdata_avg(center(1)-ceil(cN(1)/2):center(1)+ceil(cN(1)/2)-1,center(2)-ceil(cN(2)/2):center(2)+ceil(cN(2)/2)-1,center(3)-ceil(cN(3)/2):center(3)+ceil(cN(3)/2)-1,:);
    samp_crop = samp(center(1)-ceil(cN(1)/2):center(1)+ceil(cN(1)/2)-1,center(2)-ceil(cN(2)/2):center(2)+ceil(cN(2)/2)-1,center(3)-ceil(cN(3)/2):center(3)+ceil(cN(3)/2)-1,:);

    tic;
    [cI, S] = espirit_sens3D(kdata_avg_crop,logical(sum(samp_crop,4)),param.ACSsz,param);
    toc
    
    cN = size(kdata_avg_crop); cN = cN(1:3);
    
    % apply window
    ham1(:,1,1) = hamming(cN(1)); 
    ham2(1,:,1) = hamming(cN(2));
    ham3(1,1,:) = hamming(cN(3));
    
    window = bsxfun(@times, bsxfun(@times,ham1,ham2),ham3);
    S_kspace = fft3_shift(S);
    S_kspace = bsxfun(@times,S_kspace,window);
    
    padSize = ceil((N-cN)/2);
    S_kspace_pad = padarray(S_kspace,padSize,'both');
    
    S_new = ifft3_shift(S_kspace_pad);
    
    S = S_new;
    
    placeholder = [];
%     S = flipdim(flipdim(circshift(S,[-1,-1,0]),2),1);

elseif mthd == 3 % Walsh
    I = kdata2im(kdata_avg,logical(sum(samp,4)), fft3Scal); % Find low-res images
%     figure; imagesc(abs(I(:,:,2))); axis('image');
    [S,cI] = WalshCoilCombine3D(I, param);
%     S = circshift(flipdim(flipdim(S,2),1),[1,1,0]);
else
    error('Undefined method for coil sensisitive estimation');
end


% Find low-res image from the data 'd'
function img = kdata2im(d, s, fft3Scal)
[CalibSize, ~] = getCalibSize3D(s);
fil = g3d(size(d(:,:,:,1)), CalibSize); % size(d(:, :, :, 1) is the acquisition matrix
% sz = size(d);
% calib_d = zeros(sz);
% cx = sz(1)/2;
% cy = ceil(sz(2)/2);
% cz = sz(3)/2;
% xrgn = 1:64;%ceil(cx-CalibSize(1)/2+1):ceil(cx+CalibSize(1)/2);
% yrgn = ceil(cy-CalibSize(2)/2+1):ceil(cy+CalibSize(2)/2);
% zrgn = ceil(cz-CalibSize(3)/2+1):ceil(cz+CalibSize(3)/2);
% calib_d(xrgn, yrgn, zrgn, :) = d(xrgn, yrgn, zrgn, :);
% apply Gaussian filter to data from each coil
% calib_d = calib_d.* repmat(fil, [1,1,1,size(calib_d,4)]);
d = d.* repmat(fil, [1,1,1,size(d,4)]);
img = zeros(size(d));
for c = 1:size(d,4)
    img(:,:,:,c) = 1/fft3Scal * fftshift(fftshift(fftshift(ifftn(ifftshift(ifftshift(ifftshift(d(:,:,:,c),1),2),3)),3),2),1);
end



function F = g3d(S1, S2)
[x,y,z] = ndgrid(1:S1(1), 1:S1(2), 1:S1(3));
cntr=floor(S1/2)+1; % centre
sigx = S2(1)/2; % sigma x
sigy = S2(2)/2; % sigma y
sigz = S2(3)/2; % sigma z
F = 1/(2*pi*sigx*sigy*sigz)*exp(-((x-cntr(1)).^2)./(2*sigx^2) - ((y-cntr(2)).^2)./(2*sigy^2) - ((z-cntr(3)).^2)./(2*sigz^2));
F = F/max(F(:));
% figure; imagesc(F); axis('image');


% function [CalibSize, dummy] = getCalibSize(s)
% CalibSize = size(s(:,:,:));
% dummy = CalibSize;

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

