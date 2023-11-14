function [scale,time_av] = normImage(kdata,samp)

    weights = repmat(sum(samp,3),[1,1,size(kdata,3)]);
    weights(find(weights==0)) = Inf;

    % Normalize k_data by weights and take an inverse Fourier Transform
    time_av = ifft2_shift(sum(kdata,4)./weights);
    [time_av, ~] = WalshCoilCombine(time_av,1);
    
%     figure
%     imagesc(abs(time_av)/max(abs(time_av(:))))
    
    x_50 = prctile(abs(time_av(:)),60);
    tmp = time_av(abs(time_av)>=x_50);
    scale = mean(abs(tmp));
%     im = time_av;
%     im(abs(time_av)<x_50) = 0;
    
%     level = graythresh(abs(time_av)/max(abs(time_av(:))));
% %     mask = imbinarize(abs(time_av));
%     mask = im2bw(abs(time_av)/max(abs(time_av(:))),level);
    
%     scale = mean(time_av(mask==1));
%     mean(abs(time_av(:)))
end