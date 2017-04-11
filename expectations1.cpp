// expectations.cpp
// $Id: expectations1.cpp,v 1.1.1.1 2007/12/01 01:01:14 cpal Exp $
//
// code to compute the empirical and model expections for the gradient update
// Daniel Scharstein

#include "crfstereo.h"
#include "expectations.h"
#include <assert.h>
#include <math.h>
#include <iostream>


#define ABS(x) ((x)>0? (x) : (-(x)))

// CPAL - compute a better valid disp image
void computeValidDistImage(CByteImage disp)
{
  CShape sh = disp.Shape();
  int width = sh.width, height = sh.height;

  for (int y = 0; y < height; y++) {
	uchar *di = &disp.Pixel(0, y, 0);

	for (int x = 0; x < width; x++) {
      if (x < 50 || x > width - 50) { // only evaluate non-occluded pixels
		di[x]=0;
      } 
	}
  }
}

// CPAL - compute a better valid disp image
void computeMaskedImage(CByteImage disp, CByteImage truedisp)
{
  CShape sh = disp.Shape();
  int width = sh.width, height = sh.height;

  for (int y = 0; y < height; y++) {
	uchar *di = &disp.Pixel(0, y, 0);
	uchar *ti = &truedisp.Pixel(0, y, 0);

	for (int x = 0; x < width; x++) {
      if (ti[x]==0) { // only evaluate non-occluded pixels
		di[x]=0;
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


void computeModelDistU(CByteImage disp, CByteImage validdisp, fvec &distU)
{
  CShape sh = disp.Shape();
  if (sh != validdisp.Shape())
	throw CError("image shapes don't match");
  int width = sh.width, height = sh.height;

  int nU = (int)distU.size();
  if (nU == 0)
	return; // nothing to do

  assert((int)distU.size() == 1); // for now just single param

  float *costs;
  int nPixels, nDisps;
  getCostArray(costs, nPixels, nDisps);
  assert(nPixels == width * height);

  float totcost = 0;
  float occludedTotal = 0;
  float max = 0;;

  for (int y = 0; y < height; y++) {
	uchar *di = &disp.Pixel(0, y, 0);
	uchar *tdi = &validdisp.Pixel(0, y, 0);
	for (int x = 0; x < width; x++) {
      // CPAL - edges as well!!!
      // (int x = 50; x < width-50; x++) {
      if (di[x] < 80) { // only evaluate non-occluded pixels
		totcost += costs[di[x]];
      }
	  else{
		occludedTotal++;
	  }
      // move costs pointer to next pixel
      costs += nDisps;
	}
  }
  distU[0] = occludedTotal; // 100 - hmmm, need to adjust the scale of thetaU wrt. thetaV
  // CPAL removed the /100
}
// return data cost at nonoccluded pixels
void computeEmpirDistU(CByteImage disp, CByteImage validdisp, fvec &distU, float occ)
{
  CShape sh = disp.Shape();
  if (sh != validdisp.Shape())
	throw CError("image shapes don't match");
  int width = sh.width, height = sh.height;

  int nU = (int)distU.size();
  if (nU == 0)
	return; // nothing to do

  assert((int)distU.size() == 1); // for now just single param

  float *costs;
  int nPixels, nDisps;
  getCostArray(costs, nPixels, nDisps);
  assert(nPixels == width * height);

  float totcost = 0;
  float occludedTotal = 0;
  float max = 0;;

  for (int y = 0; y < height; y++) {
	uchar *di = &disp.Pixel(0, y, 0);
	uchar *tdi = &validdisp.Pixel(0, y, 0);
	for (int x = 0; x < width; x++) {
      // CPAL - edges as well!!!
      // (int x = 50; x < width-50; x++) {
      if (di[x] > occ) { // only evaluate non-occluded pixels
		totcost += costs[di[x]];
      }
	  else{
		occludedTotal++;
	  }
      // move costs pointer to next pixel
      costs += nDisps;
	}
  }
  distU[0] = occludedTotal; // 100 - hmmm, need to adjust the scale of thetaU wrt. thetaV
  // CPAL removed the /100
}

// return data cost at nonoccluded pixels
void computeEmpirDistU(CByteImage disp, CByteImage validdisp, fvec &distU)
{
  CShape sh = disp.Shape();
  if (sh != validdisp.Shape())
	throw CError("image shapes don't match");
  int width = sh.width, height = sh.height;

  int nU = (int)distU.size();
  if (nU == 0)
	return; // nothing to do

  assert((int)distU.size() == 1); // for now just single param

  float *costs;
  int nPixels, nDisps;
  getCostArray(costs, nPixels, nDisps);
  assert(nPixels == width * height);

  float totcost = 0;
  float occludedTotal = 0;
  float max = 0;;

  for (int y = 0; y < height; y++) {
	uchar *di = &disp.Pixel(0, y, 0);
	uchar *tdi = &validdisp.Pixel(0, y, 0);
	for (int x = 0; x < width; x++) {
      // CPAL - edges as well!!!
      // (int x = 50; x < width-50; x++) {
      if (di[x] > 0) { // only evaluate non-occluded pixels
		totcost += costs[di[x]];
      }
	  else{
		occludedTotal++;
	  }
      // move costs pointer to next pixel
      costs += nDisps;
	}
  }
  distU[0] = occludedTotal; // 100 - hmmm, need to adjust the scale of thetaU wrt. thetaV
  // CPAL removed the /100
}

// return smoothness cost component for each thetaV parameter at nonoccluded pixels
void computeEmpirDistV(CByteImage disp, CByteImage validdisp, CByteImage im1grad, fvec &distV)
{
  CShape sh = disp.Shape();

  if (sh != validdisp.Shape())
	throw CError("image shapes don't match");

  int width = sh.width, height = sh.height;

  int nV = (int)distV.size();

  if (nV == 0)
	return; // nothing to do

  for (int k=0; k < nV; k++)
	distV[k] = 0;

  int cnt = 0;

  for (int y = 0; y < height-1; y++) 
    {
      uchar *di = &disp.Pixel(0, y, 0);
      uchar *tdi = &validdisp.Pixel(0, y, 0);
      uchar *di2 = &disp.Pixel(0, y+1, 0);
      uchar *tdi2 = &validdisp.Pixel(0, y+1, 0);
      uchar *gr = &im1grad.Pixel(0, y, 0);

      for (int x = 0; x < width-1; x++) 
        {
          
          if (di[x] == 0) // skip occluded pixels
			  //conintue;
		  {
			  // Lam 
			  if( di[x+1] == 0){
				distV[nV-2] ++;
			  }
			  else{
				distV[nV-1] ++;
			  }
		  }
		  if( di[x+1] == 0){
				distV[nV-1] ++;
		  }

		  if (di2[x] == 0) // skip occluded pixels
		  {
			  // Lam 
			  if( di2[x+1] == 0){
				distV[nV-2] ++;
			  }
			  else{
				distV[nV-1] ++;
			  }
		  }
		  if( di2[x+1] == 0){
				distV[nV-1] ++;
		  }
		  if(di[x+1] == 0 || di[x] == 0 || di2[x+1] == 0 || di2[x] == 0)
			  continue;
		  // End Lam

          int g = gr[x]; // quantized gradient at the pixel
          
          //
          // Horizontal Edge
          //

          if (di[x+1] != 0) // neighbor right is not occluded
            {

              int dh = di[x] - di[x+1]; // horizontal disparity

              dh = ABS(dh); // disparity magnitude
              if (dh > 0)   // count if STRICTLY positive
                {
                  //distV[dh > 1]++;
                  distV[g]++;
                  cnt++;
                }
            } 
          else 
            { // count valid pixels next to occluded areas as two disp jumps of >= 2
              //distV[1] += 2;
              distV[g] += 1;
            }
          
          //
          // Vertical Edge
          //

          if (di2[x] != 0) // neighbor down is not occluded 
            {
              int dv = di[x] - di2[x]; // vertical disparity
              dv = ABS(dv); // disparity magnitude
              if (dv > 0) // count if STRICTLY positive
                {
                //distV[dv > 1]++;
                distV[g]++;
                cnt++;
                }
            } 
          else 
            { // count valid pixels next to occluded areas as two disp jumps of >= 2
              //distV[1] += 2;
              distV[g] += 1;
            }
        }
    }
}



// return smoothness cost component for each thetaV parameter at nonoccluded pixels
void computeModelDistV(CByteImage disp, CByteImage validdisp, CByteImage im1grad, fvec &distV)
{
  CShape sh = disp.Shape();

  if (sh != validdisp.Shape())
	throw CError("image shapes don't match");

  int width = sh.width, height = sh.height;

  int nV = (int)distV.size();

  if (nV == 0)
	return; // nothing to do

  for (int k=0; k < nV; k++)
	distV[k] = 0;

  int cnt = 0;
  int max = 0;

  for (int y = 0; y < height-1; y++) 
    {
      uchar *di = &disp.Pixel(0, y, 0);
      uchar *tdi = &validdisp.Pixel(0, y, 0);
      uchar *di2 = &disp.Pixel(0, y+1, 0);
      uchar *tdi2 = &validdisp.Pixel(0, y+1, 0);
      uchar *gr = &im1grad.Pixel(0, y, 0);

      for (int x = 0; x < width-1; x++) 
        {
          if(di[x] > max)
			max = di[x];
          if (di[x] == 80) // skip occluded pixels
			  //conintue;
		  {
			  // Lam 
			  if( di[x+1] == 80){
				distV[nV-2] ++;
			  }
			  else{
				distV[nV-1] ++;
			  }
		  }
		  if( di[x+1] == 80){
				distV[nV-1] ++;
		  }

		  if (di2[x] == 80) // skip occluded pixels
		  {
			  // Lam 
			  if( di2[x+1] == 80){
				distV[nV-2] ++;
			  }
			  else{
				distV[nV-1] ++;
			  }
		  }
		  if( di2[x+1] == 80){
				distV[nV-1] ++;
		  }
		  if(di[x+1] == 80 || di[x] == 80 || di2[x+1] == 80 || di2[x] == 80)
			  continue;
		  // End Lam

          int g = gr[x]; // quantized gradient at the pixel
          
          //
          // Horizontal Edge
          //

          if (di[x+1] != 0) // neighbor right is not occluded
            {

              int dh = di[x] - di[x+1]; // horizontal disparity

              dh = ABS(dh); // disparity magnitude
              if (dh > 0)   // count if STRICTLY positive
                {
                  //distV[dh > 1]++;
                  distV[g]++;
                  cnt++;
                }
            } 
          else 
            { // count valid pixels next to occluded areas as two disp jumps of >= 2
              //distV[1] += 2;
              distV[g] += 2;
            }
          
          //
          // Vertical Edge
          //

          if (di2[x] != 0) // neighbor down is not occluded 
            {
              int dv = di[x] - di2[x]; // vertical disparity
              dv = ABS(dv); // disparity magnitude
              if (dv > 0) // count if STRICTLY positive
                {
                //distV[dv > 1]++;
                distV[g]++;
                cnt++;
                }
            } 
          else 
            { // count valid pixels next to occluded areas as two disp jumps of >= 2
              //distV[1] += 2;
              distV[g] += 2;
            }
        }
    }
  std::cout << max << "\n";
}

// return expected data cost at nonoccluded pixels
void computeModelDistU(CByteImage disp, CByteImage validdisp, fvec &distU,
                       MRFEnergy* mrf)
{

  CShape sh = disp.Shape();

  if (sh != validdisp.Shape())
	throw CError("image shapes don't match");

  int width = sh.width, height = sh.height;

  int nU = (int)distU.size();

  if (nU == 0)
	return; // nothing to do
  
  assert((int)distU.size() == 1); // for now just single param

  float *costs;
  int nPixels, nDisps;
  
  getCostArray(costs, nPixels, nDisps);
  
  assert(nPixels == width * height);

  float totcost = 0;
  float occCost = 0;

  float *belief = new float[nDisps];

  for (int y = 0; y < height; y++) {

	uchar *tdi = &validdisp.Pixel(0, y, 0);

	for (int x = 0; x < width; x++) {
      // CPAL - edges as well!!!
      // (int x = 50; x < width-50; x++) {

      mrf->getBelief(mrf->pixelIndex(y,x),belief);

      // cost for non-occluded pixel
      //for (int b=0 ; b<nDisps ; b++)
      //totcost += exp(belief[b])*costs[b];

      // cost for occluded pixel
      occCost += exp(belief[nDisps]);

      // move costs pointer to next pixel
      costs += nDisps;
	}
  }
  delete[] belief;
  distU[0] = occCost; // 100 - hmmm, need to adjust the scale of thetaU wrt. thetaV
  // CPAL removed the /100
}

// return expected smoothness cost component for each thetaV parameter 
// at nonoccluded pixels
void computeModelDistV(CByteImage disp, CByteImage validdisp, CByteImage im1grad,
                       fvec &distV, MRFEnergy* mrf)
{
  CShape sh = disp.Shape();

  if (sh != validdisp.Shape())
	throw CError("image shapes don't match");

  int width = sh.width, height = sh.height;

  int nV = (int)distV.size();

  if (nV == 0)
	return; // nothing to do

  for (int k=0; k < nV; k++)
	distV[k] = 0;

  int nDisps = mrf->getNumLabels(); // is there another way to get this? -JJW

  float totBelEqDisp;
  float *belEqDisp;
  
  float totOCC = 0.0;

  //float *belief = new float[nDisps*nDisps];
  float *belief = new float[nDisps*nDisps];

  for (int y = 0; y < height-1; y++) {
	  
	  uchar *di = &disp.Pixel(0, y, 0);
      uchar *tdi = &validdisp.Pixel(0, y, 0);
      uchar *di2 = &disp.Pixel(0, y+1, 0);
      uchar *tdi2 = &validdisp.Pixel(0, y+1, 0);
      uchar *gr = &im1grad.Pixel(0, y, 0);


      for (int x = 0; x < width-1; x++) 
        {
			mrf->getEdgeBelief( mrf->pixelIndex(y,x), mrf->pixelIndex(y,x+1), belief);
			totBelEqDisp = 0;
			totOCC = 0;
          
			// Horizontal Edge
			if (di[x] == 0){
				// skip occluded pixels
				//continue;
				// Lam
				if(di[x+1] == 0){ // both are occluded
					belEqDisp = belief;
					distV[nV-2] += exp(*(belEqDisp+(nDisps*nDisps-1)));
				}
				else{
					belEqDisp = belief;
					for (int b=(nDisps - 1)*(nDisps) ; b< (nDisps-1) ; b++, belEqDisp++)
						totOCC += exp(*belEqDisp);
					distV[nV-1] += totOCC;
				}
				continue;
				// Lam
			}

			if (di[x+1] == 0){

				belEqDisp = belief;
				for (int b=nDisps - 1 ; b< (nDisps-1) ; b++, belEqDisp+=nDisps)
					totOCC += exp(*belEqDisp);
				distV[nV-1] += totOCC;
				continue;
			}


          int g = gr[x]; // quantized gradient at the pixel
          //
          //
          // Horizontal Edge
          //

          

          belEqDisp = belief;
          
          // Sum Probability along diagonal to get Pr[d_p==d_q]
          for (int b=0 ; b<nDisps ; b++, belEqDisp+=(nDisps+1))
            totBelEqDisp += exp(*belEqDisp);

#ifndef NDEBUG
          float *bel = belief;
          if (0)//x==100 && y==100)
            {
               printf("(%d,%d) horiz ",x,y);
//               for (int b=0 ; b<nDisps*nDisps ; b++, bel++)
//                 printf("%f\t",*bel);
//               printf("\n");
              printf("sum diag: %f\n",totBelEqDisp);
            }
#endif

		  // needed?
          if (totBelEqDisp<1)
          // count valid pixels next to occluded areas as two disp jumps of >= 2
            distV[g] += (1-totBelEqDisp) * (tdi[x+1] != 0 ? 1 : 2);
          
          //
          // Vertical Edge
          //

          mrf->getEdgeBelief( mrf->pixelIndex(y,x), mrf->pixelIndex(y+1,x), belief);
          totBelEqDisp = 0;
          belEqDisp = belief;
          
          // Sum Probability along diagonal to get Pr[d_p==d_q]
          for (int b=0 ; b<nDisps ; b++, belEqDisp+=(nDisps+1))
            totBelEqDisp += exp(*belEqDisp);

		  /*for (int b=nDisps - 1 ; b< (nDisps-1) ; b++, belEqDisp+=nDisps)
            totOCC += exp(*belEqDisp);

		  belEqDisp = belief;
		  for (int b=(nDisps - 1)*(nDisps) ; b< (nDisps-1) ; b++, belEqDisp++)
            totOCC += exp(*belEqDisp);*/


#ifndef NDEBUG
          bel = belief;
          if (x==171 && y==25)
            {
              printf("(%d,%d) vert ",x,y);
               for (int b=0 ; b<nDisps*nDisps ; b++, bel++)
                 printf("%f\t",*bel);
               printf("\n");
              printf("sum diag: %f\n",totBelEqDisp);
            }
#endif
          
          // count valid pixels next to occluded areas as two disp jumps of >= 2
          if (totBelEqDisp<1)
            distV[g] += (1-totBelEqDisp) * (tdi2[x] != 0 ? 1 : 2);



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
void concatVec(fvec a, fvec b, fvec &ab)
{
  int nA = (int)a.size();
  int nB = (int)b.size();
  ab.resize(nA + nB);
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
