clear all;
%%
x = randn(129,131) + 1j*randn(129,131);
level = 1;
nddwt = nd_dwt_2D('db1',size(x),1,0);
nddwt_mex = nd_dwt_2D('db1',size(x),1,1);
pause(0.5)
tic;
y_mat = nddwt.dec(x,level);
t1 =toc;

tic;
% for ind = 1:1;
	y_mex = nddwt_mex.dec(x,level);
% end
t2 = toc;

display(sprintf('Forward test %f percent faster, %f as fast, error %s',100*(t2-t1)/t1,t1/t2,num2str(max(abs(y_mat(:)-y_mex(:))))))

tic;
x_mat = nddwt.rec(y_mat);
t1 =toc;

tic;
% for ind = 1:1;
x_mex = nddwt_mex.rec(y_mex);
% end
t2 = toc;
fprintf('\n');
display(sprintf('Backward test %f percent faster, %f as fast, error %s',100*(t2-t1)/t1,t1/t2,num2str(max(abs(x_mat(:)-x_mex(:))))))
