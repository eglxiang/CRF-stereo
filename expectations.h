#ifndef EXPECTATIONS_H
#define EXPECTATIONS_H

// expectations.h
// $Id: expectations.h,v 1.10 2007/12/11 13:28:48 weinman Exp $
//
// code to compute the empirical and model expectations for the gradient update
// Daniel Scharstein

#include "MRFEnergy.h"
#include "SyncMeanField.h"
#include <vector>

typedef std::vector<float> fvec;


void computeMaskedImage(CByteImage disp, CByteImage trueimage, int ignoreVal=0);

float computeMaskedLikelihood(CByteImage disp, SyncMeanField *mrf);

// Using occlusion model
void computeEmpirDistU(CByteImage disp, fvec &distU, 
                       uchar occlVal=0, int ignoreVal=-1);
void computeModelDistU(CByteImage disp, fvec &distU, MRFEnergy* mrf,
                       int ignoreVal=-1);

// No occlusion model
void computeEmpirDistU(CByteImage truedisp, CByteImage disp, fvec &distU, 
                       int ignoreVal=-1);
void computeModelDistU(CByteImage truedisp, CByteImage disp, fvec &distU, 
                       MRFEnergy* mrf, int ignoreVal=-1);

// Gradient-modulated Occlusion model
void computeEmpirDistVPI(CByteImage validdisp, CByteImage im1grad, 
                         fvec &distV, fvec &distP, 
                         uchar occlVal=0, int ignoreVal=-1);
void computeModelDistVPI(CByteImage disp, CByteImage im1grad, 
                         fvec &distV, fvec &distP, MRFEnergy* mrf, 
                         int ignoreVal=-1);

// Simple occlusion model
void computeEmpirDistV(CByteImage disp, CByteImage im1grad, 
                       fvec &distV, fvec &distP, 
                       uchar occlVal=0, int ignoreVal=-1, int bs =1);
void computeModelDistV(CByteImage disp, CByteImage im1grad, 
                       fvec &distV, fvec &distP, 
                       MRFEnergy* mrf, int ignoreVal=-1, int bs = 1);

// No occlusion model
void computeEmpirDistV(CByteImage truedisp, CByteImage disp, CByteImage im1grad,
                       fvec &distV, int ignoreVal=-1, int bs=1);
void computeModelDistV(CByteImage truedisp, CByteImage disp, CByteImage im1grad,
                       fvec &distV, MRFEnergy* mrf, int ignoreVal=-1, int bs =1);




// vector element-wise product: prod = a .* b
void vecEltWiseProd(fvec a, fvec b, fvec &prod);

// vector element-wise quotient: quot = a ./ b
void vecEltWiseQuot(fvec a, fvec b, fvec &quot);

// vector difference: diff = a - b
void vecDiff(fvec a, fvec b, fvec &diff);

// vector copy: dst = src
void vecCopy(fvec src, fvec &dst);

// vector sum: sum = a + b
void vecSum(fvec a, fvec b, fvec &sum);

// scale vector: t = s * t
void vecScale(fvec &t, float s);

// vector L2 norm
float vecNorm(fvec a);

// limit vector elements to specified range
void vecBound(fvec &a, float minval, float maxval);

// print vector to stderr
void printfVec(std::vector<int> a, char *name);
void printfVec(fvec a, char *name);

// concatenate vectors
void concatVec(fvec a, fvec b, fvec &ab);
void concatVec(fvec a, fvec b, fvec c, fvec &abc);

// split vectors
void splitVec(fvec ab, fvec &a, fvec &b);
void splitVec(fvec ab, fvec &a, fvec &b, fvec &c);

#endif