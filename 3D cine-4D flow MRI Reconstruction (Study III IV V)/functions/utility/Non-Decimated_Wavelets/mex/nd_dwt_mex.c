#include "mex.h"
#include "matrix.h"
#include <stdlib.h>
#include <fftw3.h>
#include "nddwt.h"

/*  the gateway routine.  */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    /* */
    int num_dims,ind, numel,num_dims_kernel,level,l2_pres;
    const int *dims_mat,*dims_kernel;
    int *dims_c;
    double *imageR, *imageI, *kernelR, *kernelI,*outR,*outI;
    int *dims_out;
	
    /* Input Checks*/
    if (nrhs < 5) {
        mexErrMsgIdAndTxt("MATLAB:FFT2mx:invalidNumInputs",
                          "Four Inputs Required");
    }
    if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1])) {
        mexErrMsgIdAndTxt( "MATLAB:FFT2mx:invalidNumInputs",
            "Arrays must be double");
    }
    if (!mxIsComplex(prhs[0]) || !mxIsComplex(prhs[1]) ) {
        mexErrMsgIdAndTxt( "MATLAB:FFT2mx:invalidNumInputs",
            "Arrays must be complex");
    }
    
    /* Forward Transform */
    if (mxGetScalar(prhs[2]) ==0) {
        dims_mat = mxGetDimensions(prhs[0]);            /* array of image dimension sizes */
		if (mxGetNumberOfDimensions(prhs[0]) == 2){
			if (dims_mat[1]==1 && mxGetNumberOfElements(prhs[1]) != mxGetNumberOfElements(prhs[0])*2 ){
				mexErrMsgIdAndTxt( "MATLAB:FFT2mx:invalidNumInputs",
	                          "FIlter size and image size not consistant");
							  
			}
			else if (dims_mat[1]!=1 && mxGetNumberOfElements(prhs[1]) != mxGetNumberOfElements(prhs[0])*(1<<(mxGetNumberOfDimensions(prhs[0])))) {
				mexErrMsgIdAndTxt( "MATLAB:FFT2mx:invalidNumInputs",
	                          "FIlter size and image size not consistant");
							  
			}
		}
       	else if(mxGetNumberOfElements(prhs[1]) !=
            mxGetNumberOfElements(prhs[0])*(1<<(mxGetNumberOfDimensions(prhs[0]))) ) {
	            mexErrMsgIdAndTxt( "MATLAB:FFT2mx:invalidNumInputs",
	                          "FIlter size and image size not consistant");
        }
        

        /* Get pointers to input array */
        imageR =  (double *) mxGetPr(prhs[0]);
        imageI =  (double *) mxGetPi(prhs[0]);
        kernelR = (double *) mxGetPr(prhs[1]);
        kernelI = (double *) mxGetPi(prhs[1]);
    
        /* Get input dimesions */
        num_dims = mxGetNumberOfDimensions(prhs[0]);    /* number of image dimensions */
        dims_mat = mxGetDimensions(prhs[0]);            /* array of image dimension sizes */
        numel = mxGetNumberOfElements(prhs[0]);         /* number of elements in the array */
        dims_kernel = mxGetDimensions(prhs[1]);         /* array of kernel dimension sizes */
        num_dims_kernel = mxGetNumberOfDimensions(prhs[1]);/* number of kernel dimensions */
        level = mxGetScalar(prhs[3]);                   /* level of decomposition */
    
		if (dims_mat[1]==1){
			num_dims = 1;
		}
		
        /* C uses row major array storage, reverse dimension order to use column major*/
        dims_c = (int *) malloc(sizeof(int)*num_dims);
        for (ind =0;ind<num_dims;ind++){
            dims_c[ind] = dims_mat[ind];
        }
        
        /* set the output dims vector */
        dims_out = (int *) malloc(sizeof(int)*num_dims_kernel);
        for (ind =0;ind<num_dims_kernel;ind++){
            dims_out[ind] = dims_kernel[ind];
        }
        dims_out[num_dims_kernel-1] = ((1<<num_dims) + ((1<<num_dims)-1)*(level-1));
        
        /* Create an output array */
        plhs[0] = mxCreateNumericArray(num_dims_kernel, dims_out, mxDOUBLE_CLASS, mxCOMPLEX);
        outR = (double *) mxGetPr(plhs[0]);
        outI = (double *) mxGetPi(plhs[0]);
    
        if (level==1) {
            /* Take the wavelet transform */
            nd_dwt_dec_1level(outR, outI, imageR, imageI, kernelR, kernelI,num_dims,
                              dims_c);
        }
        else{
            nd_dwt_dec(outR, outI, imageR, imageI, kernelR, kernelI,num_dims,
                              dims_c,level);
        }
    
        /* free memory */
        free(dims_c);
        free(dims_out);
    }
    
    /* Inverse Transform */
    else{
		        
        /* Get pointers to input array */
        imageR =  (double *) mxGetPr(prhs[0]);
        imageI =  (double *) mxGetPi(prhs[0]);
        kernelR = (double *) mxGetPr(prhs[1]);
        kernelI = (double *) mxGetPi(prhs[1]);
        
        /* Get input dimesions */
        num_dims = mxGetNumberOfDimensions(prhs[0])-1;
        dims_mat = mxGetDimensions(prhs[0]);
        numel = mxGetNumberOfElements(prhs[0]);
        dims_kernel = mxGetDimensions(prhs[1]);
        num_dims_kernel = mxGetNumberOfDimensions(prhs[1])-1;
        level = mxGetScalar(prhs[3]);                   /* level of decomposition */
		l2_pres = mxGetScalar(prhs[4]);
				
        /* Input Check */
        if (mxGetNumberOfElements(prhs[1]) != (mxGetNumberOfElements(prhs[0]) - (level-1)*(mxGetNumberOfElements(prhs[1]) -mxGetNumberOfElements(prhs[1])/(1<<num_dims))) ) {
            mexErrMsgIdAndTxt( "MATLAB:FFT2mx:invalidNumInputs",
                              "FIlter size and image size not consistant");
        }
        
        /* C uses row major array storage, reverse dimension order to use column major*/
        dims_c = (int *) malloc(sizeof(int)*num_dims);
        for (ind =0;ind<num_dims;ind++){
            dims_c[ind] = dims_mat[ind];
        }
        
        /* Create an output array */
        plhs[0] = mxCreateNumericArray(num_dims_kernel, dims_kernel, mxDOUBLE_CLASS, mxCOMPLEX);
        outR = (double *) mxGetPr(plhs[0]);
        outI = (double *) mxGetPi(plhs[0]);
		
        /* Take the Inverse wavelet transform */
        if (level ==1) {
            nd_dwt_rec_1level(outR, outI, imageR, imageI, kernelR, kernelI,num_dims,
                              dims_c,l2_pres);
        }
        else {
            nd_dwt_rec(outR, outI, imageR, imageI, kernelR, kernelI,num_dims,
                              dims_c,level,l2_pres);
        }
        
        /* free dims_c*/
        free(dims_c);
    }
}

