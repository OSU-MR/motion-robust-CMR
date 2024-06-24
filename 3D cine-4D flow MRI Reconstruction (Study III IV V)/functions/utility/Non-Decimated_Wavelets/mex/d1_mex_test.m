clear all;
%%
n = 10000;
x = randn(n,1) + 1j*randn(n,1);
level = 3;
preserve_l2 = 0;
nddwt = nd_dwt_1D('db1',n,preserve_l2,0);
nddwt_mex = nd_dwt_1D('db1',n,preserve_l2,1);

num_test = 10;
tic;
for ind = 1:num_test;
	y_mat = nddwt.dec(x,level);
end
t1 =toc;


tic;
for ind = 1:num_test;
	y_mex = nddwt_mex.dec(x,level);
end
t2 = toc;

display(sprintf('Forward 1D test %f percent faster, %f as fast, error %s',100*(t2-t1)/t1,t1/t2,num2str(max(abs(y_mat(:)-y_mex(:))))))

tic;
for ind = 1:num_test;
	x_mat = nddwt.rec(y_mat);
end
t3 =toc;

tic;
for ind = 1:num_test
	x_mex = nddwt_mex.rec(y_mex);
end
t4 = toc;
fprintf('\n');
display(sprintf('Backward 1D test %f percent faster, %f as fast, error %s',100*(t4-t3)/t3,t3/t4,num2str(max(abs(x_mat(:)-x_mex(:))))))
