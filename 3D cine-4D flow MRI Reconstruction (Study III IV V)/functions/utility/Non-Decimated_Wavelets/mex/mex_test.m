clear
clc;

level = 1;
l2 = 0;
num_test = 10;

%% ========================================================================
% 1D Test
% =========================================================================
n = 10000;
x = randn(n,1) + 1j*randn(n,1);
nddwt = nd_dwt_1D('db1',n,l2,0);
nddwt_mex = nd_dwt_1D('db1',n,l2,1);

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
display(sprintf('Backward 1D test %f percent faster, %f as fast, error %s',100*(t4-t3)/t3,t3/t4,num2str(max(abs(x_mat(:)-x_mex(:))))))
fprintf('\n');
fprintf('\n');

%% ========================================================================
% 2D Test
% =========================================================================
x = randn(129,131) + 1j*randn(129,131);
nddwt = nd_dwt_2D('db1',size(x),l2,0);
nddwt_mex = nd_dwt_2D('db1',size(x),l2,1);

tic;
for ind = 1:num_test
    y_mat = nddwt.dec(x,level);
end
t1 =toc;

tic;
for ind = 1:num_test
	y_mex = nddwt_mex.dec(x,level);
end
t2 = toc;

display(sprintf('Forward 2D test %f percent faster, %f as fast, error %s',100*(t2-t1)/t1,t1/t2,num2str(max(abs(y_mat(:)-y_mex(:))))))

tic;
for ind = 1:num_test
    x_mat = nddwt.rec(y_mat);
end
t1 =toc;

tic;
for ind = 1:num_test
    x_mex = nddwt_mex.rec(y_mex);
end
t2 = toc;
display(sprintf('Backward 2D test %f percent faster, %f as fast, error %s',100*(t2-t1)/t1,t1/t2,num2str(max(abs(x_mat(:)-x_mex(:))))))
fprintf('\n');
fprintf('\n');

%% ========================================================================
% 3D Test
% =========================================================================
x = randn(131,128,30) + 1j*randn(131,128,30);
nddwt = nd_dwt_3D('db1',size(x),l2,'mat');
nddwt_mex = nd_dwt_3D('db1',size(x),l2,'mex');

tic;
for ind = 1:num_test;
	y_mat = nddwt.dec(x,level);
end
t1 = toc/num_test;

tic;
for ind = 1:num_test;
	y_mex = nddwt_mex.dec(x,level);
end
t2 = toc/num_test;

display(sprintf('Forward 3D test %f percent faster, %f as fast, error %s',100*(t2-t1)/t1,t1/t2,num2str(max(abs(y_mat(:)-y_mex(:))))))

tic;
for ind = 1:num_test;
	x_mat = nddwt.rec(y_mat);
end
t3 =toc/num_test;

tic;
for ind = 1:num_test
	x_mex = nddwt_mex.rec(y_mex);
end
t4 = toc/num_test;
display(sprintf('Backward 3D test %f percent faster, %f as fast, error %s',100*(t4-t3)/t3,t3/t4,num2str(max(abs(x_mat(:)-x_mex(:))))))
fprintf('\n');
fprintf('\n');

%% ========================================================================
% 4D Test
% =========================================================================
x = randn(128,68,8,8) + 1j*randn(128,68,8,8);
nddwt = nd_dwt_4D('db3',size(x),l2,0);
nddwt_mex = nd_dwt_4D('db3',size(x),l2,1);

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

display(sprintf('Forward 4D test %f percent faster, %f as fast, error %s',100*(t2-t1)/t1,t1/t2,num2str(max(abs(y_mat(:)-y_mex(:))))))

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

display(sprintf('Backward 4D test %f percent faster, %f as fast, error %s',100*(t2-t1)/t1,t1/t2,num2str(max(abs(x_mat(:)-x_mex(:))))))
