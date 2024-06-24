clear all

mex -v nd_dwt_mex.c nddwt.c -I"/usr/local/include/" -L"/usr/local/lib" -lfftw3 -lm -lfftw3_threads

