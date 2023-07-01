function [ image_seq_unwrap ] = phaseUnwrap( image_seq, mag )
%PHASEUNWRAP Summary of this function goes here
%   Detailed explanation goes here
figure(1)
clf
imagesc(image_seq(:,:,1))
ell = imellipse;
display('Selecte an area to unwrap and press ENTER')
pause
mask = ell.createMask;

if nargin<2
    mag = [];
    debug = 0;
end
debug = 1;
%%
[x_ind,y_ind] = find(mask==1);

image_seq_unwrap = image_seq;
for ii = 1:length(x_ind)
    t_series = squeeze(image_seq(x_ind(ii),y_ind(ii),:));
    t_unwrap = unwrap(t_series);
    image_seq_unwrap(x_ind(ii),y_ind(ii),:) = t_unwrap;
    
    if debug
        mag_t_series = squeeze(mag(x_ind(ii),y_ind(ii),:));
        figure(1)
        clf
        plot(t_series)
        hold on
        plot(t_unwrap)
        plot(mag_t_series)
        legend('original','unwrapped','mag')
        x_ind(ii)
        y_ind(ii)
        pause
    end
end

end

