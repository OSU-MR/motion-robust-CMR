function [ rand_samp, R ] = varRandSamp( sizes ,R_desired,sigma)
%VARRANDSAMP create a variabler density random sampling pattern
%   Detailed explanation goes here

envelope1 = normpdf(linspace(-1,1,sizes(1)),0,sigma);
envelope2 = normpdf(linspace(-1,1,sizes(2)),0,sigma);
envelope = envelope1'*envelope2;
% envelope1 = kaiser(sizes(1));
% envelope2 = kaiser(sizes(2));
% envelope1 = triang(sizes(1));
% envelope2 = triang(sizes(2));
% envelope = envelope1*envelope2';

% envelope = envelope/max(envelope(:));
rand_samp = rand(sizes).*envelope;

thresh = max(0.9*rand_samp(:));
keep_searching = 1;
search_ind =0;
thresh_step = 0.1;
while keep_searching && search_ind <50
    rand_samp_t = rand_samp;
    rand_samp_t(rand_samp<thresh) = 0;
    rand_samp_t(rand_samp>=thresh) = 1;
    
    figure(1)
    clf
    imagesc(rand_samp_t);
    R = numel(rand_samp_t)/length(find(rand_samp_t ==1))
    
   
    if R >= R_desired
        thresh = thresh-thresh_step;
    else
        thresh = thresh +thresh_step;
        thresh_step = 0.5*thresh_step;
        
    end
    
    if abs(R_desired - R) < 0.2
        keep_searching = 0;
    end
    search_ind = search_ind + 1;
end

rand_samp = rand_samp_t;
R = numel(rand_samp)/length(find(rand_samp ==1));
end

