function [ lambda_vec, theta1, theta2, sigma1_vec, sigma2_vec,gm] = estBGParams( x )
%ESTBGPARAMS Summary of this function goes here
%   Detailed explanation goes here

% nddwt = nd_dwt_2D('db1',[size(x,1),size(x,2)]);
% 
% w = nddwt.dec(x(:,:,1),1);
% size(w)
% tmp = w(:,:,1);
% 
% % [ lambda, theta, sigma, omega,gm] = GM_init_BG( real(tmp(:)), 2, 1 )
% gm = fitgmdist(abs(tmp(:)), 2);

sizes = size(x);
nddwt = nd_dwt_3D('db1',size(x));
w = nddwt.dec(gather(x),1);

for ind = 1:8
    tmp = squeeze(real(w(:,:,1,ind)));
    [ lambda(ind), theta1(ind), theta2(ind), sigma1(ind), sigma2(ind), gm] = GM_init(tmp(:), 2, 1 );
% gm = fitgmdist(abs(tmp(:)), 2);
end
warning('remove lines')
lambda = repmat(lambda.',[2,1]);
size(lambda)
% lambda(5:end) = 0.5*lambda(1:4);

sigma1 = repmat(sigma1.',[2,1]);
% sigma1(5:end) = 0.5*sigma1(1:4);
% sigma1(1:4) = 0.3*sigma1(1:4);

sigma2 = repmat(sigma2.',[2,1]);
% sigma2(5:end) = 0.5*sigma2(1:4);
% sigma1(1:4) = 0.3*sigma1(1:4);

lambda_vec = zeros([sizes,8]);
for ind = 1:8
   lambda_vec(:,:,:,ind) = lambda(ind);
end
lambda_vec = lambda_vec(:);

sigma1_vec = zeros([sizes,8]);
for ind = 1:8
   sigma1_vec(:,:,:,ind) = sigma1(ind);
end
sigma1_vec = sigma1_vec(:);


sigma2_vec = zeros([sizes,8]);
for ind = 1:8
   sigma2_vec(:,:,:,ind) = sigma2(ind);
end
sigma2_vec = sigma2_vec(:);

end

