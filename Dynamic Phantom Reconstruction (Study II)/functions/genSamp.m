function [samp, sampInd] =  genSamp(p)
% Generate sampling pattern
    sd = p.sd;
    R  = p.R;
    N  = p.N;
    acs= p.acs;

    rng(sd);
    samp = zeros(N);
    samp(:, randsample(N(2), round(N(2)/R))) = 1;
    samp(:, end/2-acs/2+1 : end/2+acs/2) = 1;
    sampInd = find(samp~=0); % sampling indices
end
