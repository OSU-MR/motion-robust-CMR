//
//  nddwt.h
//  wavelets_3D
//
//  Created by Adam V. Rich on 8/14/15.
//  Copyright (c) 2015 Adam. All rights reserved.
//

#ifndef nddwt_H
#define nddwt_H

/* Create a plan for a multidimensional fft */
fftw_plan init_fftw_plan(double *inR, double *inI, double *outR, double *outI, int *dims, int numDims, int howmany_trans);

/* complex multiply a*b = c */
void pointByPoint(double *aR,double *aI,double *bR,double *bI,double *cR,double *cI, int size,int conj);

/* single level decomposition */
void nd_dwt_dec_1level(double *outR, double *outI, double *imageR, double *imageI, double *kernelR, double *kernelI, int num_dims,
                       int *dims);

/* single level reconstruction */
void nd_dwt_rec_1level(double *outR, double *outI, double *imageR, double *imageI, double *kernelR, double *kernelI, int num_dims,
                       int *dims,int l2_pres);

/* multilevel decomposition */
void nd_dwt_dec(double *outR, double *outI, double *imageR, double *imageI, double *kernelR, double *kernelI, int num_dims,
                       int *dims,int level);

/* mulitleve level reconstruction */
void nd_dwt_rec(double *outR, double *outI, double *imageR, double *imageI, double *kernelR, double *kernelI, int num_dims,
                       int *dims, int level, int l2_pres);

#endif /* defined(__wavelets_3D__nddwt__) */
