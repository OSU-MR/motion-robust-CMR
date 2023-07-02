function data = ephantom_rectangle(kx,ky,par,dpar,factor,param)
%
% generates a 2D k-space rectangle (called by ephantom2)
%
% Usage: data = ephantom2_rectangle(kx,ky,par,dpar,factor);
%
% It does not check parameters! You shoud do that in ephantom2
%
% kx,ky: kspace position of samples
% par: [centerx centery sidex sidey]
% dpar: [dcenterx dcentery dsidex dsidey]
% factor: portion of dpar to add to par
% sm: width of apodization
%
% Pablo Irarrazaval
% 25-11-03
% 27-11-03 Add dpar and factor

sig  = param.sig;
row  = param.row;
ctr  = param.ctr;
nor  = param.nor;

ipar = par+factor*dpar;
data = prod(size(kx))*ipar(3)*ipar(4)*sinc(ipar(3)*kx).*sinc(ipar(4)*ky);
data = exp(-i*2*pi*(ipar(1)*kx+ipar(2)*ky)).*data;
data = data .* gauss2d([size(kx,1),size(kx,2)], ctr, sig, row, nor);
