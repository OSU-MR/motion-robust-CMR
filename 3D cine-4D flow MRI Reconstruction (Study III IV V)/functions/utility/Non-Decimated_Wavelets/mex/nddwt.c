//
//  nddwt.c
//  wavelets_3D
//
//  Created by Adam V. Rich on 8/14/15.
//  Copyright (c) 2015 Adam. All rights reserved.
//

//#include "nddwt.h"
#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>

/* =================================================================================================*/
fftw_plan init_fftw_plan(double *inR, double *inI, double *outR, double *outI, int *dims, int numDims, int howmany_trans){
    /*
     inR:        pointer to Real portion of input array
     inI:        pointer to Imaginary portion of input array
     outR:       pointer to Real portion of output array
     outI:       pointer to Imaginary portion of output array
     dims:       pointer to array of data dimensions
     numDims:    integer number of dimensions in the array
     */
    
    /* Initialize fftw plan */
    fftw_plan fftwplan;
    fftw_iodim dims_guru[numDims],howmany_dims[1];
    unsigned int numel = 1; /* number of elements in the array to be transformed */
    
    /* Set dims for guru interface */
    for (int ind = 0; ind<numDims; ind++) {
        dims_guru[ind].n = dims[ind];
        if (ind ==0) {
            dims_guru[ind].is = 1;
            dims_guru[ind].os = 1;
        }
        else{
            dims_guru[ind].is = dims_guru[ind-1].is*dims_guru[ind-1].n;
            dims_guru[ind].os = dims_guru[ind-1].is*dims_guru[ind-1].n;
//            printf("%d\n",dims_guru[ind].is);
        }
        numel *= dims[ind];
    }
    
    howmany_dims[0].n = howmany_trans;
    howmany_dims[0].is = numel;
    howmany_dims[0].os = numel;
    
    //    FILE *fp;
    //    fp=fopen("/Users/adam/Documents/c_code/wavelets/wavelets/wisdom", "w");
    //     fftw_import_wisdom_from_filename("/Users/adam/Documents/c_code/wavelets/wavelets/wisdom");
    //    fftw_import_wisdom_from_file(fp);
    
    /* Initialize Plan */
    fftwplan = fftw_plan_guru_split_dft(numDims, dims_guru,1,howmany_dims,inR,inI,outR,outI,FFTW_ESTIMATE);
    //    fftw_export_wisdom_to_file(fp);
    //    fclose(fp);
    
    /* Return an fftwplan */
    return fftwplan;
}

/* =================================================================================================*/
void pointByPoint(double *in1R, double *in1I,double *in2R,double *in2I,double *outR,double *outI, int size,int conj){
    /*  Point by point multiply two complex vectors in1 and in2 and store them in out
        if conj ==1, then the second input will be conjugate before multiplying */
    //    R[(a+ib)(c+id)]	=	ac-bd
    //    I[(a+ib)(c+id)]	=	(a+b)(c+d)-ac-bd. */
    register double ac,bd,a_b,c_d;
    
    /* normal point by point */
    if (conj ==0) {
        #pragma omp parallel for private(ac,bd,a_b,c_d)
        for (int ind = 0; ind < size; ind++) {
            ac = in1R[ind]*in2R[ind];
            bd = in1I[ind]*in2I[ind];
            a_b = in1R[ind] + in1I[ind];
            c_d = in2R[ind] + in2I[ind];
            outR[ind] = ac-bd;
            outI[ind] = a_b*c_d-ac-bd;
        }
    }
    /* conjugate the secon inpute then point by point */
    else{
       #pragma omp parallel for private(ac,bd,a_b,c_d)
        for (int ind = 0; ind < size; ind++) {
            ac = in1R[ind]*in2R[ind];
            bd = in1I[ind]*(-in2I[ind]);
            a_b = in1R[ind] +in1I[ind];
            c_d = in2R[ind] -in2I[ind];
            outR[ind] = ac-bd;
            outI[ind] = a_b*c_d-ac-bd;
        }
    }
}

/* =================================================================================================*/
void nd_dwt_dec_1level(double *outR, double *outI, double *imageR, double *imageI,
                       double *kernelR, double *kernelI, int num_dims,int *dims){
    /* fftw plan init */
    int int_threads = fftw_init_threads();
//    printf("%d\n",int_threads);
    fftw_plan_with_nthreads(8);
    fftw_plan ifftw_plan_in_place,fftw_plan_in_place;
    double numel, dims_pow;
    numel= 1;
    dims_pow = (1<<num_dims);
    
    /* Make fftw planes */
    fftw_plan_in_place = init_fftw_plan(imageR, imageI, imageR, imageI, dims, num_dims,1);
    ifftw_plan_in_place = init_fftw_plan(outI, outR, outI, outR, dims, num_dims,dims_pow);
    
    /* find the number of elements in the image */
    for (int ind =0; ind<num_dims; ind++) {
        numel *=dims[ind];
    }

    /* Take the fft of the input image */
    // fftw_execute_split_dft(fftw_plan_in_place, imageR, imageI, imageR, imageI);
 
    /* Loop over each subband. multiply be the kernel */
//    #pragma omp parallel for
    for (int ind = 0; ind<numel*dims_pow; ind+=numel) {
        pointByPoint(imageR, imageI, &kernelR[ind], &kernelI[ind], &outR[ind], &outI[ind], numel,0);
    }
    
    /* Take ifft of the result */
    fftw_execute_split_dft(ifftw_plan_in_place, outI, outR, outI, outR);
    
//     /* Normalize the DFT */
//     for (int ind = 0; ind<numel*8; ind++) {
//         outR[ind] = outR[ind]/numel;
//         outI[ind] = outI[ind]/numel;
//     }
    
    /* Destory fftw_plans*/
    fftw_destroy_plan(fftw_plan_in_place);
    fftw_destroy_plan(ifftw_plan_in_place);
}

/* =================================================================================================*/
void nd_dwt_rec_1level(double *outR, double *outI, double *imageR, double *imageI,
                       double *kernelR, double *kernelI, int num_dims,int *dims,int l2_pres){
    /* fftw plan init */
    int int_threads = fftw_init_threads();
    fftw_plan_with_nthreads(8);
    fftw_plan ifftw_plan_in_place,fftw_plan_in_place;
    unsigned int numel,dims_pow;
    numel= 1;
    dims_pow = (1<<num_dims);
    
    /* find the number of elements in the image */
    for (int ind =0; ind<num_dims; ind++) {
        numel *=dims[ind];
    }
    
    /* Make fftw planes */
    fftw_plan_in_place = init_fftw_plan(imageR, imageI, imageR, imageI, dims, num_dims,dims_pow);
    ifftw_plan_in_place = init_fftw_plan(imageI, imageR, imageI, imageR, dims, num_dims,dims_pow);
    
    /* Take the fft of the input image */
    // fftw_execute_split_dft(fftw_plan_in_place, imageR, imageI, imageR, imageI); /* */
    pointByPoint(imageR, imageI, kernelR, kernelI, imageR, imageI, numel*dims_pow,1);
    fftw_execute_split_dft(ifftw_plan_in_place, imageI, imageR, imageI, imageR);
    
    /* Loop over each subband and add them up */
//    #pragma omp parallel for
    for (int ind = 0; ind<numel*dims_pow; ind = ind + numel) {
        for (int k = 0; k<numel; k++) {
            outR[k] += imageR[k+ind];
            outI[k] += imageI[k+ind];
        }
    }
    
    /* Normalize */
	if (l2_pres==0){
		register double dims_inv = (double) 1/dims_pow;
		for (int ind = 0; ind<numel; ind++) {
			outR[ind] = outR[ind]*dims_inv;
			outI[ind] = outI[ind]*dims_inv;
		}
	}
    /* Destory fftw_plans*/
    fftw_destroy_plan(fftw_plan_in_place);
    fftw_destroy_plan(ifftw_plan_in_place);
}

/* =================================================================================================*/
void nd_dwt_dec(double *outR, double *outI, double *imageR, double *imageI,
                double *kernelR, double *kernelI, int num_dims,int *dims,int level){
    
    double *approxR, *approxI;
    unsigned int size_trans,level_start,numel;
    fftw_plan fftw_plan_in_place;
		
    numel = 1;
    /* find the number of elements in the image */
    for (int ind =0; ind<num_dims; ind++) {
        numel *=dims[ind];
    }
	
	/* allocate memory for old approximate coefficients */
    approxR = (double *) malloc(sizeof(double)*numel);
    approxI = (double *) malloc(sizeof(double)*numel);
	
	/* create plan for old approximate coefficients */
	fftw_plan_in_place = init_fftw_plan(approxR, approxI, approxR, approxI, dims, num_dims,1);
    
    size_trans = (1<<num_dims) + ((1<<num_dims)-1)*(level-1);
    level_start = ((1<<num_dims)-1)*(level-1);
    /* printf("size trans =%d  levels start = %d\n",size_trans*numel,level_start*numel); */
    
	/* Take fft and store the old aproximate coefficient in approx */
    nd_dwt_dec_1level(&outR[level_start*numel], &outI[level_start*numel], imageR, imageI,
                      kernelR, kernelI, num_dims, dims);
    for (int ind =0; ind<numel; ind++) {
        approxR[ind]= outR[level_start*numel + ind];
        approxI[ind]= outI[level_start*numel + ind];
    }
	
	/* take fft of old approximate coefficients, this could be avoided! */
	fftw_execute_split_dft(fftw_plan_in_place, approxR, approxI, approxR, approxI);
		
	/* loop over remaining levels */
    for (int level_ind =level-1; level_ind >0; level_ind--) {
        level_start = ((1<<num_dims)-1)*(level_ind-1);
        nd_dwt_dec_1level(&outR[level_start*numel], &outI[level_start*numel], approxR, approxI,
                          kernelR, kernelI, num_dims, dims);
        for (int ind =0; ind<numel; ind++) {
            approxR[ind]= outR[level_start*numel + ind];
            approxI[ind]= outI[level_start*numel + ind];
        }
		fftw_execute_split_dft(fftw_plan_in_place, approxR, approxI, approxR, approxI);
    }
    
    /* Free Memory */
    free(approxR);
    free(approxI);
}

/* =================================================================================================*/
void nd_dwt_rec(double *outR, double *outI, double *imageR, double *imageI,
                double *kernelR, double *kernelI, int num_dims,int *dims, int level,int l2_pres){
    unsigned int level_start,numel, dims_pow;
    numel = 1;
	fftw_plan fftw_plan_in_place;
    // dims_pow = (1<<num_dims);
	
    /* find the number of elements in the image */
    for (int ind =0; ind<num_dims; ind++) {
        numel *=dims[ind];
    }
	fftw_plan_in_place = init_fftw_plan(outR, outI, outR, outI, dims, num_dims,1);
	
    /* reconstruct lowest levels */
    nd_dwt_rec_1level(outR, outI, imageR, imageI, kernelR, kernelI, num_dims,dims,l2_pres);
    
    /* Take the fft of the input image */
    fftw_execute_split_dft(fftw_plan_in_place, outR, outI, outR, outI); 
	
    /* copy previous reconstructed level to approximate coefficent locations */
	level_start = ((1<<num_dims)-1);
    for (int ind = 0;ind<numel; ind++) {
        imageR[level_start*numel + ind] = outR[ind];
        imageI[level_start*numel + ind] = outI[ind];
        outR[ind] = 0.0;
        outI[ind] = 0.0;
    }
	
    /* loop through sucessive levels */
	for (int level_ind =2; level_ind <=level; level_ind++) {
        level_start = ((1<<num_dims)-1)*(level_ind-1);
              
        /* reconstruct  next lowest levels */
        nd_dwt_rec_1level(outR, outI, &imageR[level_start*numel], &imageI[level_start*numel], kernelR, kernelI, num_dims,dims,l2_pres);
    
		if (level_ind != level){
		    /* Take the fft of the input image */
		    fftw_execute_split_dft(fftw_plan_in_place, outR, outI, outR, outI); 
	
		    /* copy previous reconstructed level to approximate coefficent locations */
			level_start = ((1<<num_dims)-1)*(level_ind);
		    for (int ind = 0;ind<numel; ind++) {
		        imageR[level_start*numel + ind] = outR[ind];
		        imageI[level_start*numel + ind] = outI[ind];
		        outR[ind] = 0.0;
		        outI[ind] = 0.0;
		    }
		}
	}

}


