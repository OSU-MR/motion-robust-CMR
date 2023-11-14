function [ phi1,phi2,eta] = CGMUpdate_local(zhat,zvar,sizes,phi1,phi2,eta)
%CGMUPDATE Summary of this function goes here
%   Detailed explanation goes here

% Define Complex Normal Distribution
CN = @(mu,sig) 1./(pi*sig).*exp( -1./sig.*abs(mu).^2);

% Truncate zhat and zvar to Wavelet Coefficients
zhat = zhat(end-8*prod(sizes)+1:end);
zvar = zvar(end-8*prod(sizes)+1:end);

S1 = (1-eta).*CN(zhat,phi1+zvar)+eps;
S2 = eta.*CN(zhat,phi2+zvar)+eps;

b1 = S1./(S1+S2);
b2 = S2./(S1+S2);

gamma1 = zhat.*phi1./(zvar + phi1);
gamma2 = zhat.*phi2./(zvar + phi2);

v1 = (zvar.*phi1)./(zvar + phi1);
v2 = (zvar.*phi2)./(zvar + phi2);


% phi1 = sum( b1.*(abs(gamma1).^2 + v1) )/sum(b1);
% phi2 = sum( b2.*(abs(gamma2).^2 + v2) )/sum(b2);
% eta = mean(S2./(S1+S2));


b1 = reshape(b1,[sizes,8]);
b2 = reshape(b2,[sizes,8]);
gamma1 = reshape(gamma1,[sizes,8]);
gamma2 = reshape(gamma2,[sizes,8]);
v1 = reshape(v1,[sizes,8]);
v2 = reshape(v2,[sizes,8]);
phi1 = reshape(b1.*(abs(gamma1).^2 + v1),[sizes,8]);
phi2 = reshape(b2.*(abs(gamma2).^2 + v2)./b2,[sizes,8]);
eta  = reshape(S2./(S1+S2),[sizes,8]);

% Creat A Gaussian Filter
filt_size = 3;
filt_size_t = 3;
filt = normpdf(linspace(-1,1,filt_size),0,1);
filt_t = normpdf(linspace(-1,1,filt_size_t),0,1);
filter_2 = filt.'*filt;
filter_3 = zeros(filt_size,filt_size,filt_size_t);
for ind = 1:size(filter_3,2)
    filter_3(:,ind,:) = filter_2(:,ind)*filt_t;
end
filter_3 = filter_3/sum(filter_3(:));
filter_3 = fftn(filter_3,size(b1(:,:,:,1)));
x_conv = @(x) real(circshift(ifftn(fftn(x).*filter_3),[-1,-1,-1]));

for ind = 1:8
	phi1(:,:,:,ind) = x_conv(phi1(:,:,:,ind))./x_conv(b1(:,:,:,1));
	phi2(:,:,:,ind) = x_conv(phi2(:,:,:,ind))./x_conv(b2(:,:,:,1));
	eta(:,:,:,ind)  = x_conv(eta(:,:,:,ind));
end

phi1 = phi1(:);
phi2 = phi2(:);
eta = abs(eta(:))/(1.05*max(abs(eta(:))));
end

