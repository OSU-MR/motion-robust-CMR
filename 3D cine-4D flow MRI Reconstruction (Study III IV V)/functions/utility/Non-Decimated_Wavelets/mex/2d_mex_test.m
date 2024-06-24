clear all;
%%
x = randn(128,128) + 1j*randn(128,128);

nddwt = nd_dwt_2D('db1',size(x),1,0);
nddwt_mex = nd_dwt_2D('db1',size(x),1,1);

tic;
y_mat = nddwt.dec(x,4);
t1 =toc;


% tic;
% for ind = 1:1;
% 	y_mex = nddwt_mex.dec(x,4);
% end
% t2 = toc;

% display(sprintf('%f',100*(t2-t1)/t1))
% display(sprintf('%f',t1/t2))
% %%
% max(abs(y_mat(:)-y_mex(:)))
