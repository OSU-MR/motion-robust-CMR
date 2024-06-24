clear all;
%%
x = randn(128,68,12,20) + 1j*randn(128,68,12,20);
level = 3;
l2 = 1;
nddwt = nd_dwt_4D('db3',size(x),l2,0);
nddwt_mex = nd_dwt_4D('db3',size(x),l2,1);
num_test = 1;
pause(0.5)
tic;
for ind = 1:num_test
    y_mat = nddwt.dec(x,level);
end
t1 =toc/num_test;

tic;
for ind = 1:num_test;
	y_mex = nddwt_mex.dec(x,level);
end
t2 = toc/num_test;

display(sprintf('Forward test %f percent faster, %f as fast, error %s',100*(t2-t1)/t1,t1/t2,num2str(max(abs(y_mat(:)-y_mex(:))))))

tic;
for ind = 1:num_test
    x_mat = nddwt.rec(y_mat);
end
t1 =toc/num_test;

tic;
for ind = 1:num_test
x_mex = nddwt_mex.rec(y_mex);
end
t2 =toc/num_test;

fprintf('\n');
display(sprintf('Backward test %f percent faster, %f as fast, error %s',100*(t2-t1)/t1,t1/t2,num2str(max(abs(x_mat(:)-x_mex(:))))))
