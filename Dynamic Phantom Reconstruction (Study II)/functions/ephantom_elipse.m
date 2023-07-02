function data = ephantom_elipse(kx,ky,par,dpar,factor,param)
%
% generates a 2D k-space elipse (called by ephantom2)
%
% Usage: data = ephantom2_elipse(kx,ky,par,dpar,factor);
%
% It does not check parameters! You shoud do that in ephantom2
%
% kx,ky: kspace position of samples
% par: [center_x center_y semiaisx_x semiaxis_yy]
% dpar: [dcenter_x dcenter_y dsemiaisx_x dsemiaxis_yy]
% factor: portion of dpar to add to par
% sm: width of apodization

% Pablo Irarrazaval
% 25-11-03 Creation
% 27-11-03 Add dpar and factor

sig  = param.sig;
row  = param.row;
ctr  = param.ctr;
nor  = param.nor;

ipar = par+factor*dpar;
arg = sqrt(ipar(3)^2*kx.^2 + ipar(4)^2*ky.^2);
[x y] = find(abs(arg)<1e-8); arg(x,y) = 1e-8; % avoid zeros

data = prod(size(kx))*ipar(3)*ipar(4)*...
       0.5*besselj(1,pi*arg)./arg;
data = exp(-i*2*pi*(ipar(1)*kx+ipar(2)*ky)).*data;
% sig  = min(size(kx))/2.4;
data = data .* gauss2d([size(kx,1),size(kx,2)], ctr, sig, row, nor);
