function [ phi1,phi2,eta] = CGMUpdate(zhat,zvar,sizes,phi1,phi2,eta)
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

for ind = 1:8
	phi1(:,:,:,ind) = sum(sum(sum(phi1(:,:,:,ind))))./sum(sum(sum(b1(:,:,:,ind))));
	phi2(:,:,:,ind) = sum(sum(sum(phi2(:,:,:,ind))))./sum(sum(sum(b2(:,:,:,ind))));
	eta(:,:,:,ind)  = mean(mean(mean(eta(:,:,:,ind))));
end

% warning('scaling low band')
% phi1(:,:,:,1) = 50*phi1(:,:,:,1);
% phi2(:,:,:,1) = 50*phi2(:,:,:,1);
% warning('remove these lines')
% phi1(:,:,:,5:end) = 5*phi1(:,:,:,5:end);
% phi2(:,:,:,5:end) = 5*phi2(:,:,:,5:end);

phi1 = phi1(:);
phi2 = phi2(:);
eta = eta(:);
end

