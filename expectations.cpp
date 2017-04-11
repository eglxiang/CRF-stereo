// expectations.cpp
// $Id: expectations.cpp,v 1.11 2007/12/11 13:28:48 weinman Exp $
//
// code to compute the empirical and model expections for the gradient update
// Daniel Scharstein

#include "crfstereo.h"
#include "expectations.h"
#include <assert.h>
#include <math.h>
#include <iostream>

using std::cout;


#define ABS(x) ((x)>0? (x) : (-(x)))


// CPAL - compute a better valid disp image
void computeMaskedImage(CByteImage disp, CByteImage truedisp, int ignoreVal)
{

  if (ignoreVal<0)
    return;

  CShape sh = disp.Shape();
  int width = sh.width, height = sh.height;

  for (int y = 0; y < height; y++) {
	uchar *di = &disp.Pixel(0, y, 0);
	uchar *ti = &truedisp.Pixel(0, y, 0);

	for (int x = 0; x < width; x++) {
      if (ti[x]==ignoreVal) { // only evaluate non-ignored pixels
		di[x]=ignoreVal;
      } 
	}
  }
}

float computeMaskedLikelihood(CByteImage disp, SyncMeanField *mrf)
{
  
  CShape sh = disp.Shape();
  int width = sh.width, height = sh.height;

  register float loglik = 0;

  for (int y = 0; y < height; y++) {

	uchar *ti = &disp.Pixel(0, y, 0);

	for (int x = 0; x < width; x++) {

      if (ti[x]!=0)
        loglik += mrf->getBelief(mrf->pixelIndex(y,x),ti[x]);
    }
  }

  return loglik;
}

// Occlusion model
// return data cost at occluded pixels
void computeEmpirDistU(CByteImage disp, fvec &distU, uchar occlVal, int ignoreVal)
{
  CShape sh = disp.Shape();

  int width = sh.width, height = sh.height;

  int nU = (int)distU.size();

  if (nU == 0)
	return; // nothing to do


  float *costs;
  int nPixels, nStates;

  getCostArray(costs, nPixels, nStates);

  assert(nPixels == width * height);

  float totCost = 0;
  float totOccl = 0;

  for (int y = 0; y < height; y++) {

	uchar *di = &disp.Pixel(0, y, 0);

	for (int x = 0; x < width; x++) {

      if ((int)(di[x])!=ignoreVal) {
        
        if (di[x]!=occlVal){
          totCost += costs[di[x]];
        }
        else{
          // cost of occluded (delta function)
          ++totOccl;
          
        }
      }

      // move costs pointer to next pixel
      costs += nStates;
	}
  }

  distU[0] = totOccl;

  if(nU == 2)
	distU[1] = totCost;

}

// No occlusion model
// return data cost at non-occluded pixels
void computeEmpirDistU(CByteImage validdisp, CByteImage disp, fvec &distU, int ignoreVal)
{
  CShape sh = validdisp.Shape();

  int width = sh.width, height = sh.height;

  int nU = (int)distU.size();

  if (nU == 0)
	return; // nothing to do


  float *costs;
  int nPixels, nStates;

  getCostArray(costs, nPixels, nStates);

  assert(nPixels == width * height);

  float totCost = 0;

  for (int y = 0; y < height; y++) {

	uchar *di = &disp.Pixel(0, y, 0);
	uchar *tdi = &validdisp.Pixel(0, y, 0);


	for (int x = 0; x < width; x++) {

      if ((int)(tdi[x])!=ignoreVal && tdi[x]!=0){
        totCost += costs[di[x]];
      }
		
      // move costs pointer to next pixel
      costs += nStates;
	}
  }
  distU[0] = totCost;
}

// Occlusion model
// return expected data cost
void computeModelDistU(CByteImage disp, fvec &distU, MRFEnergy* mrf, int ignoreVal)
{

  CShape sh = disp.Shape();


  int width = sh.width, height = sh.height;

  int nU = (int)distU.size();

  if (nU == 0)
	return; // nothing to do
  

  float *costs;
  int nPixels, nStates;
  
  getCostArray(costs, nPixels, nStates);

  assert(nPixels == width * height);

  float totCost = 0;
  float occCost = 0;

  float *belief = new float[nStates];

  for (int y = 0; y < height; y++) {

	uchar *di = &disp.Pixel(0, y, 0);

	for (int x = 0; x < width; x++) {

      if ((int)(di[x])!=ignoreVal) {

        mrf->getBelief(mrf->pixelIndex(y,x),belief);
        
        // cost for non-occluded pixel
        for (int b=0 ; b<nStates -1 ; b++)
          totCost += exp(belief[b])*costs[b];
        
        // cost for occluded pixel
        occCost += exp(belief[nStates-1]);
        
      }
      
      // move costs pointer to next pixel
      costs += nStates;
	}
  }

  delete[] belief;

  distU[0] = occCost; 

  if(nU == 2)
    distU[1] = totCost;

}

// No occlusion model
void computeModelDistU(CByteImage validdisp, CByteImage disp, fvec &distU, MRFEnergy* mrf, int ignoreVal)
{

  CShape sh = validdisp.Shape();


  int width = sh.width, height = sh.height;

  int nU = (int)distU.size();

  if (nU == 0)
	return; // nothing to do
  

  float *costs;
  int nPixels, nStates;
  
  getCostArray(costs, nPixels, nStates);

  assert(nPixels == width * height);

  float totCost = 0;
  float occCost = 0;

  float *belief = new float[nStates];

  for (int y = 0; y < height; y++) {

    uchar *tdi = &validdisp.Pixel(0, y, 0);

	for (int x = 0; x < width; x++) {

      if ((int)(tdi[x])!=ignoreVal && tdi[x]!=0) { // skip occluded pixels

        
        mrf->getBelief(mrf->pixelIndex(y,x),belief);
        
        // cost for non-occluded pixel
        for (int b=0 ; b<nStates -1 ; b++)
          totCost += exp(belief[b])*costs[b];
        
        // cost for occluded pixel
        occCost += exp(belief[nStates-1]);
      }

      // move costs pointer to next pixel
      costs += nStates;
	}
  }

  delete[] belief;

  distU[0] = occCost; 

  if(nU == 2)
    distU[1] = totCost;

}


// return smoothness cost component for each thetaV parameter at nonoccluded pixels
void computeEmpirDistVPI(CByteImage disp, CByteImage im1grad, fvec &distV, fvec &distP, uchar occlVal, int ignoreVal)
{
	CShape sh = disp.Shape();
	int width = sh.width, height = sh.height;
	int nV = (int)distV.size();
	int nP = (int)distP.size();

	if (nV == 0 && nP == 0)
		return; // nothing to do

	for (int k=0; k < nV; k++)
		distV[k] = 0;

	for (int k=0; k < nP; k++)
		distP[k] = 0;

	for (int y = 0; y < height-1; y++) {

		uchar *di = &disp.Pixel(0, y, 0);
		uchar *di2 = &disp.Pixel(0, y+1, 0);
		uchar *gr = &im1grad.Pixel(0, y, 0);

		for (int x = 0; x < width-1; x++) {

			int g = gr[x]; // quantized gradient at the pixel
          
            if ((int)(di[x])==ignoreVal)
              continue;

			// Horizontal Edge

			if (di[x] == occlVal && di[x+1] == occlVal){ // both are occluded
              distP[g]++;
			}
			else if (di[x] == occlVal || di[x+1] == occlVal){ // one is occluded
              distP[2*nV + g]++;
			}
			else if ((int)(di[x+1])!=ignoreVal) { // neither is occluded

					int dh = di[x] - di[x+1]; // horizontal disparity
					dh = ABS(dh); // disparity magnitude
					if (dh > 0){   // count if STRICTLY positive
						distV[g]++;

					}	
			}
          
			// Vertical Edge

			if (di[x] == occlVal && di2[x] == occlVal){ // both are occluded
              distP[nV + g]++;
			}
			else if(di[x] == occlVal || di2[x] == occlVal){ // one is occluded
              distP[3*nV + g]++;
			}
			else if ((int)(di2[x])!=ignoreVal) { // neither is occluded

					int dv = di[x] - di2[x]; // vertical disparity
					dv = ABS(dv); // disparity magnitude
					if (dv > 0){ // count if STRICTLY positive         
						distV[g]++;
					}
			}
		} // end for x
    } // for y	
}


// Simple occlusion model
// return smoothness cost component for each thetaV parameter at nonoccluded pixels
void computeEmpirDistV(CByteImage disp, CByteImage im1grad, fvec &distV, fvec &distP, uchar occlVal, int ignoreVal, int bs)
{
	CShape sh = disp.Shape();
	int width = sh.width/bs, height = sh.height/bs;
	int nV = (int)distV.size();
	int nP = (int)distP.size();

	if (nV == 0 && nP==0)
		return; // nothing to do

	for (int k=0; k < nV; k++)
		distV[k] = 0;

	for (int k=0; k < nP; k++)
		distP[k] = 0;

    int cnt = 0;
	for (int y = 0; y < height-1; y++) {
		uchar *di = &disp.Pixel(0, y, 0);
		uchar *di2 = &disp.Pixel(0, y+1, 0);
		uchar *gr = &im1grad.Pixel(0, y, 0);

		for (int x = 0; x < width-1; x++) {

          if ((int)(di[x])==ignoreVal)
            continue;

			int g = gr[x]; // quantized gradient at the pixel
          
			// Horizontal Edge
			if (di[x] == occlVal && di[x+1] == occlVal){ // both are occluded
				distP[0] ++;
			}
			else if (di[x] == occlVal || di[x+1] == occlVal){ // one is occluded
				distP[1] ++;
			}
			else if ((int)(di[x+1])!=ignoreVal) { // neither is occluded
					int dh = di[x] - di[x+1]; // horizontal disparity
					dh = ABS(dh); // disparity magnitude
					if (dh > 0){   // count if STRICTLY positive
						distV[g]++;
					}
                    else
                      cnt++;
			}
          
			// Vertical Edge

			if (di[x] == occlVal && di2[x] == occlVal){ // both are occluded
		 		distP[0] ++;
			}
			else if(di[x] == occlVal || di2[x] == occlVal){ // one is occluded
				distP[1] ++;
			}
			else if ((int)(di2[x])!=ignoreVal){
					int dv = di[x] - di2[x]; // vertical disparity
					dv = ABS(dv); // disparity magnitude
					if (dv > 0){ // count if STRICTLY positive         
						distV[g]++;
					}
                    else
                      cnt++;

			}
		} // end for x
    } // for y	
    std::cout<<"count: " << cnt << "\n";
}

// No occlusion model
// return smoothness cost component for each thetaV parameter at nonoccluded pixels
void computeEmpirDistV(CByteImage validdisp, CByteImage disp, CByteImage im1grad, fvec &distV, int ignoreVal, int bs)
{
  CShape sh = disp.Shape();

  int width = sh.width/bs, height = sh.height/bs;

  int nV = (int)distV.size();

  if (nV == 0)
	return; // nothing to do

  for (int k=0; k < nV; k++)
	distV[k] = 0;


  for (int y = 0; y < height-1; y++) 
    {
      uchar *di = &disp.Pixel(0, y, 0);
      uchar *di2 = &disp.Pixel(0, y+1, 0);
      uchar *tdi = &validdisp.Pixel(0, y, 0);
      uchar *tdi2 = &validdisp.Pixel(0, y+1, 0);

      uchar *gr = &im1grad.Pixel(0, y, 0);

      for (int x = 0; x < width-1; x++) 
        {
          
          if ((int)(tdi[x])==ignoreVal || tdi[x] == 0) // skip occluded pixels
            continue;

          int g = gr[x]; // quantized gradient at the pixel
          
          //
          // Horizontal Edge
          //

          // neighbor right is not occluded
          if ((int)(tdi[x+1])!=ignoreVal && tdi[x+1] != 0)
            {

              int dh = di[x] - di[x+1]; // horizontal disparity

              dh = ABS(dh); // disparity magnitude

              if (dh > 0)   // count if STRICTLY positive
                {
                  distV[g]++;
                }
            } 

          //
          // Vertical Edge
          //

          // neighbor down is not occluded 
          if ((int)(tdi2[x])!=ignoreVal && tdi2[x] != 0) 
            {
              int dv = di[x] - di2[x]; // vertical disparity
              dv = ABS(dv); // disparity magnitude

              if (dv > 0) // count if STRICTLY positive
                {
                  distV[g]++;
                }
            } 
        }
    }
}

// Gradient-modulated occlusion model
// return expected smoothness cost component for each thetaV parameter 
void computeModelDistVPI(CByteImage disp, CByteImage im1grad, fvec &distV, fvec &distP, MRFEnergy* mrf, int ignoreVal, int bs)
{
  CShape sh = disp.Shape();

  int width = sh.width, height = sh.height;

  int nV = (int)distV.size();
  int nP = (int)distP.size();

  if (nV == 0 & nP == 0)
	return; // nothing to do

  for (int k=0; k < nV; k++)
	distV[k] = 0;

  for (int k=0; k < nP; k++)
    distP[k] = 0;

  int nStates = mrf->getNumLabels();

  float totBelEqDisp,totBel1Occ,belBothOcc, belNotEq;

  float *belEqDisp;

  float *belief = new float[nStates*nStates];

  for (int y = 0; y < height-1; y++) 
    {
      uchar *gr = &im1grad.Pixel(0, y, 0);
      uchar *di = &disp.Pixel(0, y, 0);
      uchar *di2 = &disp.Pixel(0,y+1,0);

      for (int x = 0; x < width-1; x++) 
        {
          
          if ((int)(di[x])==ignoreVal)
            continue;
          
          int g = gr[x]; // quantized gradient at the pixel
          
          //
          // Horizontal Edge
          //

          if ((int)(di[x+1])!=ignoreVal) {

            mrf->getEdgeBelief( mrf->pixelIndex(y,x), mrf->pixelIndex(y,x+1), 
                                belief);

            totBelEqDisp = 0;
            totBel1Occ = 0;
            belBothOcc = 0;
            
            belEqDisp = belief;
            
            // Sum Probability along diagonal to get Pr[d_p==d_q]
            // for non-occluded states
            for (int b=0 ; b<nStates-1 ; b++, belEqDisp+=(nStates+1))
            {
              totBelEqDisp += exp(*belEqDisp);
            }
          
            // Probability of both occluded
            belBothOcc = exp(*belEqDisp);
            
            // Sum probability of either pixel being occluded
            // JW: could speed with pointer arithmetic. i'm lazy and error prone.
            for (int b=0 ; b<nStates-1 ; b++)
            {
              totBel1Occ += exp(belief[b*(nStates)+nStates-1]); // [b][nDisps]
              totBel1Occ += exp(belief[(nStates-1)*(nStates)+b]); //[nDisps][b]
            }
            
            belNotEq = 1-totBelEqDisp-belBothOcc-totBel1Occ;
            
            
            if (belBothOcc>0)
              distP[g] += (belBothOcc>1 ? 1 : belBothOcc);
            
            
            if (totBel1Occ>0)
              distP[2*nV + g] += (totBel1Occ>1 ? 1 : totBel1Occ);
            
          
            if (belNotEq>0)
              distV[g] += (belNotEq>1 ? 1 : belNotEq);
          }

          
          //
          // Vertical Edge
          //
          
          if ((int)(di2[x])!=ignoreVal) {

            mrf->getEdgeBelief( mrf->pixelIndex(y,x), mrf->pixelIndex(y+1,x), 
                                belief);
            
            totBelEqDisp = 0;
            totBel1Occ = 0;
            belBothOcc = 0;
            
            belEqDisp = belief;
            
            // Sum Probability along diagonal to get Pr[d_p==d_q]
            // for non-occluded states
            for (int b=0 ; b<nStates-1 ; b++, belEqDisp+=(nStates+1))
            {
              totBelEqDisp += exp(*belEqDisp);
            }
          
            // Probability of both occluded
            belBothOcc = exp(*belEqDisp);
            
            // Sum probability of either pixel being occluded
            // JW: could speed with pointer arithmetic. i'm lazy and error prone.
            for (int b=0 ; b<nStates-1 ; b++)
            {
              totBel1Occ += exp(belief[b*(nStates)+nStates-1]); // [b][nDisps]
              totBel1Occ += exp(belief[(nStates-1)*(nStates)+b]); //[nDisps][b]
            }
            
            belNotEq = 1-totBelEqDisp-belBothOcc-totBel1Occ;
            
            
          if (belBothOcc>0)
            distP[nV + g] += (belBothOcc>1 ? 1 : belBothOcc);
          
          if (totBel1Occ>0)
            distP[3*nV + g] += (totBel1Occ>1 ? 1 : totBel1Occ);
          
          if (belNotEq>0)
            distV[g] += (belNotEq>1 ? 1 : belNotEq);
          }
        }
    }
  delete[] belief;
}

// Simple occlusion model
// return expected smoothness cost component for each thetaV parameter 

void computeModelDistV(CByteImage disp, CByteImage im1grad, fvec &distV, fvec &distP, MRFEnergy* mrf, int ignoreVal, int bs)
{
  CShape sh = disp.Shape();

  int width = sh.width, height = sh.height;

  int nV = (int)distV.size();
  int nP = (int)distP.size();

  if (nV == 0 & nP == 0)
	return; // nothing to do

  for (int k=0; k < nV; k++)
	distV[k] = 0;

  for (int k=0; k < nP; k++)
    distP[k] = 0;

  int nStates = mrf->getNumLabels();

  float totBelEqDisp,totBel1Occ,belBothOcc, belNotEq;

  float *belEqDisp;

  float *belief = new float[nStates*nStates];

  for (int y = 0; y < height-1; y++) 
    {
      uchar *gr = &im1grad.Pixel(0, y, 0);
      uchar *di = &disp.Pixel(0,y,0);
      uchar *di2 = &disp.Pixel(0,y+1,0);

      for (int x = 0; x < width-1; x++) 
        {

          if ((int)(di[x])==ignoreVal)
            continue;

          int g = gr[x]; // quantized gradient at the pixel
          
          //
          // Horizontal Edge
          //
          
          if ((int)(di[x+1])!=ignoreVal) {

            mrf->getEdgeBelief( mrf->pixelIndex(y,x), mrf->pixelIndex(y,x+1), 
                                belief);
            
            totBelEqDisp = 0;
            totBel1Occ = 0;
            belBothOcc = 0;
            
            belEqDisp = belief;
            
            // Sum Probability along diagonal to get Pr[d_p==d_q]
            // for non-occluded states
            for (int b=0 ; b<nStates-1 ; b++, belEqDisp+=(nStates+1))
            {
              totBelEqDisp += exp(*belEqDisp);
            }
          
            // Probability of both occluded
            belBothOcc = exp(*belEqDisp);
            
            // Sum probability of either pixel being occluded
            // JW: could speed with pointer arithmetic. i'm lazy and error prone.
            for (int b=0 ; b<nStates-1 ; b++)
            {
              totBel1Occ += exp(belief[b*(nStates)+nStates-1]); // [b][nDisps]
              totBel1Occ += exp(belief[(nStates-1)*(nStates)+b]); //[nDisps][b]
            }

            belNotEq = 1-totBelEqDisp-belBothOcc-totBel1Occ;
            
            
            if (belBothOcc>0)
              distP[0] += (belBothOcc>1 ? 1 : belBothOcc);
            
            if (totBel1Occ>0)
            distP[1] += (totBel1Occ>1 ? 1 : totBel1Occ);

            if (belNotEq>0)
              distV[g] += (belNotEq>1 ? 1 : belNotEq);
          }
          
          //
          // Vertical Edge
          //

          if ((int)(di2[x])!=ignoreVal) {

            mrf->getEdgeBelief( mrf->pixelIndex(y,x), mrf->pixelIndex(y+1,x), 
                                belief);
            
            totBelEqDisp = 0;
            totBel1Occ = 0;
            belBothOcc = 0;
            
            belEqDisp = belief;
            
            // Sum Probability along diagonal to get Pr[d_p==d_q]
          // for non-occluded states
            for (int b=0 ; b<nStates-1 ; b++, belEqDisp+=(nStates+1))
            {
              totBelEqDisp += exp(*belEqDisp);
            }
          
            // Probability of both occluded
            belBothOcc = exp(*belEqDisp);
            
            // Sum probability of either pixel being occluded
            // JW: could speed with pointer arithmetic. i'm lazy and error prone.
            for (int b=0 ; b<nStates-1 ; b++)
            {
              totBel1Occ += exp(belief[b*(nStates)+nStates-1]); // [b][nDisps]
              totBel1Occ += exp(belief[(nStates-1)*(nStates)+b]); //[nDisps][b]
            }
            
            belNotEq = 1-totBelEqDisp-belBothOcc-totBel1Occ;
            
            if (belBothOcc>0)
              distP[0] += (belBothOcc>1 ? 1 : belBothOcc);
            
            if (totBel1Occ>0)
              distP[1] += (totBel1Occ>1 ? 1 : totBel1Occ);
            
            if (belNotEq>0)
              distV[g] += (belNotEq>1 ? 1 : belNotEq);
          }      
          }
        }
      
  delete[] belief;
}

// No occlusion model
// return expected smoothness cost component for each thetaV parameter 
// at nonoccluded pixels
void computeModelDistV(CByteImage validdisp, CByteImage disp, CByteImage im1grad, fvec &distV, MRFEnergy* mrf, int ignoreVal)
{
  CShape sh = validdisp.Shape();

  int width = sh.width, height = sh.height;

  int nV = (int)distV.size();

  if (nV == 0)
	return; // nothing to do

  for (int k=0; k < nV; k++)
	distV[k] = 0;

  int nStates = mrf->getNumLabels();

  float totBelEqDisp,belNotEq;

  float *belEqDisp;

  float *belief = new float[nStates*nStates];

  for (int y = 0; y < height-1; y++) 
    {
      uchar *gr = &im1grad.Pixel(0, y, 0);

      uchar *tdi = &validdisp.Pixel(0, y, 0);
      uchar *tdi2 = &validdisp.Pixel(0, y+1, 0);

      for (int x = 0; x < width-1; x++) 
        {
          
          int g = gr[x]; // quantized gradient at the pixel

          if ((int)(tdi[x])==ignoreVal || tdi[x]==0) // skip occluded pixels
            continue;

          //
          // Horizontal Edge
          //
          
          // neighbor right is not occluded
          if ((int)(tdi[x+1])!=ignoreVal && tdi[x+1]!=0) { 
            
            mrf->getEdgeBelief( mrf->pixelIndex(y,x), mrf->pixelIndex(y,x+1), belief);
            
            totBelEqDisp = 0;
            
            belEqDisp = belief;
          
            // Sum Probability along diagonal to get Pr[d_p==d_q]
            for (int b=0 ; b<nStates ; b++, belEqDisp+=(nStates+1))
              {
                totBelEqDisp += exp(*belEqDisp);
              }
            
            belNotEq = 1-totBelEqDisp;

            if (belNotEq>0)
              distV[g] += (belNotEq>1 ? 1 : belNotEq);
          }
          
          //
          // Vertical Edge
          //

          // neighbor down is not occluded
          if ((int)(tdi2[x])!=ignoreVal && tdi2[x]!=0) {

            mrf->getEdgeBelief( mrf->pixelIndex(y,x), mrf->pixelIndex(y+1,x), belief);
            
            totBelEqDisp = 0;

            belEqDisp = belief;
          
            // Sum Probability along diagonal to get Pr[d_p==d_q]
            for (int b=0 ; b<nStates-1 ; b++, belEqDisp+=(nStates+1))
            {
              totBelEqDisp += exp(*belEqDisp);
            }
            

            belNotEq = 1-totBelEqDisp;

            if (belNotEq>0)
              distV[g] += (belNotEq>1 ? 1 : belNotEq);

          }
        }
    }

  delete[] belief;
}



// vector element-wise product: prod = a .* b
void vecEltWiseProd(fvec a, fvec b, fvec &prod)
{
  int nK = (int)a.size();
  assert((int)b.size() == nK);
  prod.resize(nK);
  for (int k = 0; k < nK; k++)
	prod[k] = a[k] * b[k];
}

// vector element-wise quotient: quot = a ./ b
void vecEltWiseQuot(fvec a, fvec b, fvec &quot)
{
  int nK = (int)a.size();
  assert((int)b.size() == nK);
  quot.resize(nK);
  for (int k = 0; k < nK; k++)
	quot[k] = (a[k] == 0) ? 0 // don't divide by b if a is already 0 (avoids 0/0 case)
      : a[k] / b[k];
}

// vector difference: diff = a - b
void vecDiff(fvec a, fvec b, fvec &diff)
{
  int nK = (int)a.size();
  assert((int)b.size() == nK);
  diff.resize(nK);
  for (int k = 0; k < nK; k++)
	diff[k] = a[k] - b[k];
}

// vector sum: sum = a + b
void vecSum(fvec a, fvec b, fvec &sum)
{
  int nK = (int)a.size();
  assert((int)b.size() == nK);
  sum.resize(nK);
  for (int k = 0; k < nK; k++)
	sum[k] = a[k] + b[k];
}

// vector copy: dst = src
void vecCopy(fvec src, fvec &dst)
{
  int nK = (int)src.size();
  dst.resize(nK);
  for (int k = 0; k < nK; k++)
	dst[k] = (float)src[k];
}

// scale vector: t = s* t
void vecScale(fvec &t, float s)
{
  int nK = (int)t.size();
  for (int k = 0; k < nK; k++)
	t[k] *= s;
}

// vector L2 norm
float vecNorm(fvec a)
{
  float ssum = 0;
  int nK = (int)a.size();
  for (int k = 0; k < nK; k++)
	ssum += a[k] * a[k];

  return sqrt(ssum);
}

// limit vector elements to specified range
void vecBound(fvec &a, float minval, float maxval)
{
    int nK = (int)a.size();
    bool flagmin = false;
    bool flagmax = false;
    for (int k = 0; k < nK; k++) {
	if (a[k] < minval) {
	    a[k] = minval;
	    flagmin = true;
	}
	if (a[k] > maxval) {
	    a[k] = maxval;
	    flagmax = true;
	}
    }
    if (flagmin)
	fprintf(stderr, "vecBound: set to minval %g\n", minval);
    if (flagmax)
	fprintf(stderr, "vecBound: set to maxval %g\n", maxval);
}


// print vector to stderr
void printfVec(std::vector<int> a, char *name)
{
  int nK = (int)a.size();
  DEBUG_OUT2(verbose, debugfile, "%s: (%s", name, nK==0? ")\n" : "");
  for (int k = 0; k < nK; k++)
	DEBUG_OUT2(verbose, debugfile, "%d%s", a[k], k==nK-1? ")\n" : ",");
}

void printfVec(fvec a, char *name)
{
  int nK = (int)a.size();
  DEBUG_OUT2(verbose, debugfile, "%s: (%s", name, nK==0? ")\n" : "");
  for (int k = 0; k < nK; k++)
	DEBUG_OUT2(verbose, debugfile, "%g%s", a[k], k==nK-1? ")\n" : ",");
}

// concatenate vectors
void concatVec(fvec a, fvec b, fvec c, fvec &abc)
{
  int nA = (int)a.size();
  int nB = (int)b.size();
  int nC = (int)c.size();
  abc.resize(nA + nB + nC);
  int k, i = 0;
  for (k = 0; k < nA; k++)
	abc[i++] = a[k];
  for (k = 0; k < nB; k++)
	abc[i++] = b[k];
  for (k = 0; k < nC; k++)
	abc[i++] = c[k];
}

// split vectors
void splitVec(fvec abc, fvec &a, fvec &b, fvec &c)
{
  int nA = (int)a.size();
  int nB = (int)b.size();
  int nC = (int)c.size();

  assert((int)abc.size() == nA + nB + nC);
  int k, i = 0;
  for (k = 0; k < nA; k++)
	a[k] = abc[i++];
  for (k = 0; k < nB; k++)
	b[k] = abc[i++];
  for( k = 0; k < nC; k++)
    c[k] = abc[i++];
}

// concatenate vectors
void concatVec(fvec a, fvec b, fvec &ab)
{
  int nA = (int)a.size();
  int nB = (int)b.size();

  ab.resize(nA + nB );
  int k, i = 0;
  for (k = 0; k < nA; k++)
	ab[i++] = a[k];
  for (k = 0; k < nB; k++)
	ab[i++] = b[k];
}

// split vectors
void splitVec(fvec ab, fvec &a, fvec &b)
{
  int nA = (int)a.size();
  int nB = (int)b.size();

  assert((int)ab.size() == nA + nB);
  int k, i = 0;
  for (k = 0; k < nA; k++)
	a[k] = ab[i++];
  for (k = 0; k < nB; k++)
	b[k] = ab[i++];
}


// vector difference: diff = a - b
