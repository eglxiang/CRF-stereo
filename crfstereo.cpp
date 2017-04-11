/* crfstereo.cpp
 * $Id: crfstereo.cpp,v 1.10 2007/12/08 14:11:48 weinman Exp $
 *
 * MRF stereo matcher for CRF stereo paper
 * derived from mrfstereo.cpp
 * uses MRF library (GC/expansion)
 *
 * Daniel Scharstein
 * 11/13/2006 initial version
 */



// maximum number of iterations for graph cut / BP 
// (can be big, since termination is controlled via closeEnoughPercent parameter)
#define MAXITER 1

#include "crfstereo.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "GCoptimization.h"
#include "SyncSumProd.h"
#include "SparseSyncSumProd.h"
#include "MeanField.h"
#include "SparseMeanField.h"
#include "TRW-S.h"
#include <vector>
#include "convert.h"
#include "CImg.h"

using namespace cimg_library;

using std::vector;

int verbose;
FILE *debugfile;


// global variables to keep track of what needs to be freed later
MRF::CostVal *dsiArray = NULL;
vector<MRF::CostVal *> dsiArray2;
vector<int> globalNpixels;
int globalNdisps = 0;
bool randomv = false;

unsigned int cueIndex = 0;
vector<MRF::CostVal *> hCueArrayV;
vector<MRF::CostVal *> vCueArrayV;

vector<MRF::CostVal *> hBOcclArrayV;
vector<MRF::CostVal *> vBOcclArrayV;
vector<MRF::CostVal *> hSOcclArrayV;
vector<MRF::CostVal *> vSOcclArrayV;

vector<MRF::CostVal *> dsiArrayV;

bool pairwise = false;
bool localOccl = false;
bool pairwiseInteraction = false;

vector<int> globalwidth; // for debugging
//float globallambda1 = 1;
//float globallambda2 = 1;

int inferencer = USE_BP;

MRF::CostVal occluded1 = 0.0;
MRF::CostVal occluded2 = 0.0;

void setGlobalNpixels(vector<int>& w, vector<int>& h){
	for(unsigned int i = 0; i < w.size(); ++i)
		globalNpixels.push_back(w[i]*h[i]);
}
void setGlobalWidth(vector<int>& v){
	globalwidth = v;
}
// make sure MRF library is compiled with float costs
void checkForFloatCosts() {
  if (typeid(MRF::EnergyVal) != typeid(float) ||
      typeid(MRF::CostVal) != typeid(float)) {
	fprintf(stderr, "Need to define MRF::EnergyVal and MRF::CostVal as float in mrf.h\n");
	exit(1);
  }
}

// compute global data cost dsiArray - matching cost volume (absolute color diffs)
// also compute best data cost index per pixel (i.e. WTA disparities)
void initializeDataCostVector(vector<CByteImage> im1,       // reference (left) images
                        vector<CByteImage> im2,       // target (right) images
                        int nD,               // number of disparities
                        bool birchfield,      // use Birchfield/Tomasi costs (absolute differences by default)
                        bool squaredDiffs,    // use squared differences (absolute differences by default)
                        int truncDiffs,       // truncate diffs (before squaring), use 255 for no truncation
                        vector<CByteImage> &wta,   // winner-take-all disparities
						int bs, // block size
						float pn) //p-norm)    
{
  checkForFloatCosts();
  if(localOccl == true)
	  globalNdisps = nD +1;		
  else
	  globalNdisps = nD;

  for(unsigned int i = 0; i < im2.size(); ++i){	
	  CShape sh = im1[i].Shape();	
	  sh.width /= bs; sh.height /= bs;
	  int width = sh.width, height = sh.height, nB = sh.nBands;
	  int tempglobalNpixels = width * height;
	  int nColors = __min(3, nB);
	  dsiArrayV.push_back(new MRF::CostVal[tempglobalNpixels * globalNdisps]);

	  DEBUG_OUT4(verbose, debugfile, "Image size: %d x %d (%d colors), nDisp = %d\n", width, height, nColors, nD);
	  DEBUG_OUT1(verbose, debugfile, "Data term: %s differences ", squaredDiffs ? "squared" : "absolute");

	  if (truncDiffs < 255)
		  DEBUG_OUT1(verbose, debugfile, "truncated at %d", truncDiffs);
	  DEBUG_OUT1(verbose, debugfile, "%s\n", birchfield ? " (Birchfield/Tomasi)" : " ");
  
	  sh.nBands = 1;
	  wta[i].ReAllocate(sh);


	  // worst value for sumdiff below 
	  MRF::CostVal worst_match = (MRF::CostVal)nColors * pow(255,pn); //(squaredDiffs ? 255 * 255 : 255);
	  // truncation threshold - NOTE: if squared, don't multiply by nColors (Eucl. dist.)
	  MRF::CostVal maxsumdiff = (MRF::CostVal)nColors*pow(abs(truncDiffs),pn);//(squaredDiffs ? truncDiffs * truncDiffs : nColors * abs(truncDiffs));
	  // value for out-of-bounds matches
	  MRF::CostVal badcost = __min(worst_match, maxsumdiff);
	
	  bool inBound = true;
	  int dsiIndex = 0;
	  for (int y = 0; y < height; y++) {
		  uchar *WTArow = &wta[i].Pixel(0, y, 0);
		  for (int x = 0; x < width; x++) {
		//	  uchar *pix1 = &im1[i].Pixel(x, y, 0);
			  int bestd = 0;
			  MRF::CostVal bestval = badcost;
			  for (int d = 0; d < nD; d++) {
				  // blockmatching
				float sumdiff = 0;
				inBound = (x-d >= 0 ? true : false); // if the boundary is out of bound than the whole block is out of bound.
				MRF::CostVal dsiValue=0;
				for(int bsx = 0; bsx < bs; ++bsx){
					for(int bsy = 0; bsy < bs; ++bsy){
						uchar *pix1 = &im1[i].Pixel(x*bs+bsx, y*bs+bsy, 0);
						int x2 = (x-d)*bs+bsx;
						
						if (inBound && d < nD) { // in bounds
							uchar *pix2 = &im2[i].Pixel(x2, y*bs+bsy, 0);
					  		
					  		for (int b = 0; b < nColors; b++) {
						  		int diff = 0;
						  		if (birchfield) {
						  	  		// Birchfield/Tomasi cost
							  		int im1c = pix1[b];
							  		int im1l = x == 0?   im1c : (im1c + pix1[b - nB]) / 2;
							  		int im1r = x == width-1? im1c : (im1c + pix1[b + nB]) / 2;
							  		int im2c = pix2[b];
							  		int im2l = x2 == 0?   im2c : (im2c + pix2[b - nB]) / 2;
							  		int im2r = x2 == width-1? im2c : (im2c + pix2[b + nB]) / 2;
							  		int min1 = __min(im1c, __min(im1l, im1r));
							  		int max1 = __max(im1c, __max(im1l, im1r));
							  		int min2 = __min(im2c, __min(im2l, im2r));
							  		int max2 = __max(im2c, __max(im2l, im2r));
							  		int di1 = __max(0, __max(im1c - max2, min2 - im1c));
							  		int di2 = __max(0, __max(im2c - max1, min1 - im2c));
							  		diff = __min(di1, di2);
						  		} else {
							  		// simple absolute difference
							  		int di = pix1[b] - pix2[b];
							  		diff = abs(di);
						  		}
						  		// square diffs if requested (Birchfield too...)
						  		// sumdiff += (squaredDiffs ? diff * diff : diff);
								sumdiff += pow(diff,pn);
							} // finish one pixel for the block.

					  		// truncate diffs
					  		dsiValue += ((float)__min(sumdiff, maxsumdiff));
						} 
				  		else { // out of bounds: use maximum truncated cost
					  		dsiValue += badcost;
				  		}
					}
				} // finish calculate the cost for the block.
				  		
				// The cost of pixel p and label l is stored at dsiArray[p*nLabels+l]
				(dsiArrayV[i])[dsiIndex++] = pow((float)dsiValue,1/pn);
				//costs.push_back((float)dsiValue);

				if (dsiValue < bestval) {
					bestval = dsiValue;  
					bestd = d;  
				}
			  } // end block matching 
			  // Lam
			  if(localOccl == true) 
				  dsiIndex++; // Skip initializing occluded cost, will be assigned later.
			  WTArow[x] = bestd;
		  }
	  }	
  }
}


// compute global data cost dsiArray - matching cost volume (absolute color diffs)
// also compute best data cost index per pixel (i.e. WTA disparities)
void initializeDataCost(CByteImage im1,       // reference (left) image
                        CByteImage im2,       // target (right) image
                        int nD,               // number of disparities
                        bool birchfield,      // use Birchfield/Tomasi costs (absolute differences by default)
                        bool squaredDiffs,    // use squared differences (absolute differences by default)
                        int truncDiffs,       // truncate diffs (before squaring), use 255 for no truncation
                        CByteImage &wta)      // winner-take-all disparities
{
  checkForFloatCosts();

  CShape sh = im1.Shape();
  int width = sh.width, height = sh.height, nB = sh.nBands;
  int nColors = __min(3, nB);

  DEBUG_OUT4(verbose, debugfile, "Image size: %d x %d (%d colors), nDisp = %d\n", width, height, nColors, nD);
  DEBUG_OUT1(verbose, debugfile, "Data term: %s differences ", squaredDiffs ? "squared block" : "absolute block");
  if (truncDiffs < 255)
	DEBUG_OUT1(verbose, debugfile, "truncated at %d", truncDiffs);
  DEBUG_OUT1(verbose, debugfile, "%s\n", birchfield ? " (Birchfield/Tomasi)" : " ");

  int tempglobalNpixels = width * height;

  if(localOccl == true)
	  globalNdisps = nD +1;
  else
	  globalNdisps = nD;
  
  dsiArray = new MRF::CostVal[tempglobalNpixels * globalNdisps];
  sh.nBands = 1;
  wta.ReAllocate(sh);

  // worst value for sumdiff below
  MRF::CostVal worst_match = (MRF::CostVal)nColors * (squaredDiffs ? 255 * 255 : 255);
  // truncation threshold - NOTE: if squared, don't multiply by nColors (Eucl. dist.)
  MRF::CostVal maxsumdiff = (MRF::CostVal)(squaredDiffs ? truncDiffs * truncDiffs : nColors * abs(truncDiffs));
  // value for out-of-bounds matches
  MRF::CostVal badcost = __min(worst_match, maxsumdiff);
  float pn = 0.5;
  int bs = 3;
  int dsiIndex = 0;
  for (int y = 0; y < height-bs+1; y++) {
	uchar *WTArow = &wta.Pixel(0, y, 0);
		for (int x = 0; x < width-bs+1; x++) {
			int bestd = 0;
			MRF::CostVal bestval = badcost;
			for (int d = 0; d < nD; d++) {
				// blockmatching
				int sumdiff = 0;
				MRF::CostVal dsiValue=0;
				for(int bsx = 0; bsx < bs; ++bsx){
					for(int bsy = 0; bsy < bs; ++bsy){
						uchar *pix1 = &im1.Pixel(x+bsx, y+bsy, 0);
						int x2 = x-d+bsx;
						if (x2 >= 0 && d < nD) { // in bounds
							uchar *pix2 = &im2.Pixel(x2, y+bsy, 0);
							for (int b = 0; b < nColors; b++) {
								int diff = 0;
								if (birchfield) {
									// Birchfield/Tomasi cost
									int im1c = pix1[b];
									int im1l = x == 0?   im1c : (im1c + pix1[b - nB]) / 2;
									int im1r = x == width-1? im1c : (im1c + pix1[b + nB]) / 2;
									int im2c = pix2[b];
									int im2l = x2 == 0?   im2c : (im2c + pix2[b - nB]) / 2;
									int im2r = x2 == width-1? im2c : (im2c + pix2[b + nB]) / 2;
									int min1 = __min(im1c, __min(im1l, im1r));
									int max1 = __max(im1c, __max(im1l, im1r));
									int min2 = __min(im2c, __min(im2l, im2r));
									int max2 = __max(im2c, __max(im2l, im2r));
									int di1 = __max(0, __max(im1c - max2, min2 - im1c));
									int di2 = __max(0, __max(im2c - max1, min1 - im2c));
									diff = __min(di1, di2);
								} else {
									// simple absolute difference
									int di = pix1[b] - pix2[b];
									diff = abs(di);
								}
								// square diffs if requested (Birchfield too...)
								// sumdiff += (squaredDiffs ? diff * diff : diff);
							}
							// truncate diffs
							dsiValue += pow(((float) __min(sumdiff, maxsumdiff)),pn);
						} 
						else { // out of bounds: use maximum truncated cost
							dsiValue = badcost;
						}
					}
				}
				// The cost of pixel p and label l is stored at dsiArray[p*nLabels+l]	
				dsiArray[dsiIndex++] = pow((float)dsiValue,1/pn);
				//costs.push_back((float)dsiValue);
						
				if (dsiValue < bestval) {
					bestval = dsiValue;
					bestd = d;
				}
				
				// Lam
					if(localOccl == true) 
						dsiIndex++; // Skip initializing occluded cost, will be assigned later.
					WTArow[x] = bestd;
				}
		}
  }
}

// create data cost from global cost array dsiArray
DataCost *computeDataCost(int index)
{
  return new DataCost(computeDataCostArray(index));
}

MRF::CostVal* computeDataCostArray(int )
{
  if (dsiArray == NULL)
	throw CError("call initializeDataCost first");

  return dsiArray;
}

// create data cost using parameters thetaU
DataCost *computeDataCost(std::vector<float> thetaU, int index)
{
  return new DataCost ( computeDataCostArray (thetaU, index));
}

MRF::CostVal* computeDataCostArray(std::vector<float> thetaU, int index)
{
	if (dsiArrayV.empty() == true)
		throw CError("call initializeDataCost first");
	int nU = (int)thetaU.size();

	if (nU == 0) {
		return dsiArrayV[index];
	} else if (nU == 1) {
		// pair wise with no learning scaler
		float local = thetaU[0];

		int arrSize = globalNpixels[index] * globalNdisps;

		// allocate second array when called the first time
		if (dsiArray2.size() == index)
			dsiArray2.push_back(new MRF::CostVal[arrSize]);

		// make scaled copy
		for (int i = 0; i < arrSize; ++i) {
			if( ( (i+1) % globalNdisps == 0))
				dsiArray2[index][i] = local;
			else
				dsiArray2[index][i] = dsiArrayV[index][i];
		}
		return dsiArray2[index];
	} else if (nU == 2) {
		float scale = thetaU[1];
		// scale cost values by single factor
	
		float local = thetaU[0];

		int arrSize = globalNpixels[index] * globalNdisps;

		// allocate second array when called the first time
		if (dsiArray2.size() == index)
			dsiArray2.push_back(new MRF::CostVal[arrSize]);

		// make scaled copy
		for (int i = 0; i < arrSize; ++i) {
			if( ( (i+1) % globalNdisps == 0))
				dsiArray2[index][i] = local;
			else
				dsiArray2[index][i] = scale * dsiArrayV[index][i];
		}
		return dsiArray2[index];
	}else { // nU > 1
		throw CError("more than one thetaU parameter not yet supported");
	}
}

// 
void getCostArray(float* &costs, int &nPixels, int &nDisps)
{
	costs = (dsiArray2.empty() == true)? dsiArray : dsiArray2[cueIndex];
  nPixels = globalNpixels[cueIndex];
  nDisps = globalNdisps;
}


// compute quantized absolute color gradients of input image
void computeQuantizedGradients(CByteImage im1, CByteImage &im1grad, std::vector<int> gradThresh, int bs)
  /*
    use nK = gradThresh.size() levels of quantization.
    quantization is done as follows:
    compute absolute gradient (RMS difference of color bands) in x and y
    given gradient g, assign largest k such that gradThresh[k] <= g
    the last entry (gradThresh[nK-1]) is ignored and is assumed to be infinity
    in im1grad, use band 0 for x gradient, and band 1 for y gradients
  */
{
  CShape sh = im1.Shape();
  sh.width /= bs;
  sh.height /= bs;
  CByteImage temp = im1;

  ScaleAndOffset(im1, temp, (float)1.0/bs, 0); // scale to integer disps
  int width = sh.width, height = sh.height, nB = sh.nBands;
  int nColors = __min(3, nB);

  CImg<unsigned char> a;
  a = im1;
  a.resize(width,height,1,nColors,5);
  assign(im1, a);
  a = im1;

  gradThresh.push_back(256 * 256 * nColors); // make sure last threshold is big enough

  int nK = (int)gradThresh.size();
  int k;

  DEBUG_OUT1(verbose, debugfile, "%d Gradient bins:\t0-", nK);
  for (k = 0; k < nK-1; k++)
	DEBUG_OUT2(verbose, debugfile, "%d\t%d-", gradThresh[k]-1, gradThresh[k]);
  DEBUG_OUT0(verbose, debugfile, "inf\n");

  for (k = 0; k < nK-1; k++) // compute sum of squared color differences below,
	gradThresh[k] *= gradThresh[k] * nColors;  //  need to adjust thresholds to get RMS distance

  sh.nBands = 2;  // band 0: x gradient, band 1: y gradient
  im1grad.ReAllocate(sh);
  for (int y = 0; y < height; y++) {
	//printf("y=%d\n",y);
	for (int x = 0; x < width; x++) {

		int sx = 0, sy = 0;
		uchar *grad   = &im1grad.Pixel(x, y, 0);
		//for(int bsx = 0; bsx < bs; ++bsx){
		//	for(int bsy = 0; bsy < bs; ++bsy){
								  
				//uchar *pix   = &im1.Pixel(x*bs+bsx, y*bs+bsy, 0);
    //  			uchar *pix1x = &im1.Pixel(bs*(x + (x < width-1)) +bsx, y*bs+bsy, 0);
    //  			uchar *pix1y = &im1.Pixel(x*bs+bsx, (y + (y < height-1))*bs+bsy, 0);

				uchar *pix   = &im1.Pixel(x, y, 0);
      			uchar *pix1x = &im1.Pixel((x + (x < width-1)) , y, 0);
      			uchar *pix1y = &im1.Pixel(x, (y + (y < height-1)), 0);
				
      			for (int b = 0; b < nColors; b++) {
					int dx = pix[b] - pix1x[b];
					int dy = pix[b] - pix1y[b];
					sx += dx * dx;
					sy += dy * dy;
      			}


			//}
	  //}
	  //sx = sx/(bs*bs*);
	  //sy = sy/(bs*bs);

      for (k = 0; gradThresh[k] <= sx; k++); // find lowest bin for x gradient
      assert(k < nK);
      grad[0] = k;
      for (k = 0; gradThresh[k] <= sy; k++)
		; // find lowest bin for y gradient
      assert(k < nK);
      grad[1] = k;
	}
  }
  if (1) {
	CByteImage g;
	BandSelect(im1grad, g, 0, 0);
	ScaleAndOffset(g, g, (float)255.0/nK, 0);
	WriteImageVerb(g, "grad0.png", true);
	BandSelect(im1grad, g, 1, 0);
	ScaleAndOffset(g, g, (float)255.0/nK, 0);
	WriteImageVerb(g, "grad1.png", true);
	//exit(1);
  }
}



// function-based smoothness cost computation
MRF::CostVal fnCost(int pix1, int pix2, int di, int dj)
{
	if (pairwise){
		if(di == globalNdisps-1 && dj == globalNdisps-1)
			return occluded2;
		else if(di == globalNdisps-1 || dj == globalNdisps-1)
			return occluded1;
	}
	else if (pairwiseInteraction){
		if(di == globalNdisps-1 && dj == globalNdisps-1){		
			if (pix1 < pix2) {
				if (pix2 - pix1 == 1)
					return hBOcclArrayV[cueIndex][pix1];
				else {
					assert(pix2 - pix1 == globalwidth[cueIndex]);
					return vBOcclArrayV[cueIndex][pix1];
				}
			} else { // (pix2 < pix1)
				if (pix1 - pix2 == 1)
					return hBOcclArrayV[cueIndex][pix2];
				else {
					assert(pix1 - pix2 == globalwidth[cueIndex]);
					return vBOcclArrayV[cueIndex][pix2];
				}
			}
		}
		else if(di == globalNdisps-1 || dj == globalNdisps-1){
			if (pix1 < pix2) {
				if (pix2 - pix1 == 1)
					return hSOcclArrayV[cueIndex][pix1];
				else {
					assert(pix2 - pix1 == globalwidth[cueIndex]);
					return vSOcclArrayV[cueIndex][pix1];
				}
			} else { // (pix2 < pix1)
				if (pix1 - pix2 == 1)
					return hSOcclArrayV[cueIndex][pix2];
				else {
					assert(pix1 - pix2 == globalwidth[cueIndex]);
					return vSOcclArrayV[cueIndex][pix2];
				}
			}
		}
	}

  if (di == dj)
	return 0;

  //int deltad = di < dj ? dj-di : di-dj;
  //if (deltad == 1)
  //return globallambda1;
  //return globallambda2;
  // later, multiuply this by gradient cost:

  // else deltad > 0

  MRF::CostVal gradcost = 0;
  if (pix1 < pix2) {
	if (pix2 - pix1 == 1)
      gradcost = hCueArrayV[cueIndex][pix1];
	else {
      assert(pix2 - pix1 == globalwidth[cueIndex]);
      gradcost = vCueArrayV[cueIndex][pix1];
	}
  } else { // (pix2 < pix1)
	if (pix1 - pix2 == 1)
      gradcost = hCueArrayV[cueIndex][pix2];
	else {
      assert(pix1 - pix2 == globalwidth[cueIndex]);
      gradcost = vCueArrayV[cueIndex][pix2];
	}
  }
  return gradcost;
}

SmoothnessCost *computeSmoothnessCost(CByteImage im1grad, std::vector<float> theta, int i)
{
  std::vector<float> thetaP;

  return computeSmoothnessCost(im1grad,theta,thetaP, i);
}

SmoothnessCost *computeSmoothnessCost(CByteImage im1grad, std::vector<float> theta, std::vector<float> thetaP, int i)
{
 
  return new SmoothnessCost( computeSmoothnessCostFunction (im1grad, theta, thetaP, i) );
}

MRF::SmoothCostGeneralFn computeSmoothnessCostFunction(CByteImage im1grad, 
													   std::vector<float> theta, int i)
{
  std::vector<float> thetaP;
  return computeSmoothnessCostFunction(im1grad,theta,thetaP, i);
}

MRF::SmoothCostGeneralFn computeSmoothnessCostFunction(CByteImage im1grad, 
													   std::vector<float> theta, 
                                                       std::vector<float> thetaP, int i)
{
  checkForFloatCosts();
  int nK = (int)theta.size();
  int nP = (int)thetaP.size();

  CShape sh = im1grad.Shape();
  int width = sh.width, height = sh.height;

   if( nP == 2){
	occluded1 = thetaP[1]; // Only one of the pair is occluded.
	occluded2 = thetaP[0]; // Both Pair are occluded.
  }
   else if( nP == nK*4){
	   if (hBOcclArrayV.size() <= i)
		   hBOcclArrayV.push_back(new MRF::CostVal[width * height]);
	   if (vBOcclArrayV.size() <= i)
		   vBOcclArrayV.push_back(new MRF::CostVal[width * height]);
	   if (hSOcclArrayV.size() <= i)
		   hSOcclArrayV.push_back(new MRF::CostVal[width * height]);
	   if (vSOcclArrayV.size() <= i)
		   vSOcclArrayV.push_back(new MRF::CostVal[width * height]);
   }
  
  if (nK == 0) // run without smoothness params
	theta.push_back(1.0f);
  
  DEBUG_OUT1(verbose, debugfile, "%d Gradient weights:\t", nK);
  for (int k = 0; k < nK; k++)
	DEBUG_OUT1(verbose, debugfile, "%g\t", theta[k]);
  DEBUG_OUT0(verbose, debugfile, "\n");
  
  
  //DEBUG_OUT(verbose, debugfile, "Smoothness term: L%d norm, truncated at %d\n", smoothexp, smoothmax);
  
  // allocate the arrays when called the first time
  if (hCueArrayV.size() <= i)
	  hCueArrayV.push_back(new MRF::CostVal[width * height]);
  if (vCueArrayV.size() <= i)
	  vCueArrayV.push_back(new MRF::CostVal[width * height]);

  //globallambda1 = theta[0];
  //globallambda2 = theta[1];
  //DEBUG_OUT(verbose, debugfile, "Smoothness term: lambda1 = %g, lambda2 = %g\n", globallambda1, globallambda2);
  
  int n = 0;
  for (int y = 0; y < height; y++) {
	for (int x = 0; x < width; x++) {
      uchar *grad = &im1grad.Pixel(x, y, 0);
      hCueArrayV[i][n] = theta[grad[0]];
      vCueArrayV[i][n] = theta[grad[1]];
      n++;
	}
  }


  if( nP == 4*nK){
	n = 0;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			uchar *grad = &im1grad.Pixel(x, y, 0);

			hBOcclArrayV[i][n] = thetaP[grad[0]];
			vBOcclArrayV[i][n] = thetaP[nK + grad[1]];

			hSOcclArrayV[i][n] = thetaP[2*nK + grad[0]];
			vSOcclArrayV[i][n] = thetaP[3*nK + grad[1]];
			n++;
		}
	}
  }

  // FUNCTION WRAPPING BREAKS THE ABILITY FOR THIS RETURN TYPE - JJW
//   if (0) { // define smoothnest cost using truncated linear
// 	MRF::CostVal lambda = 1;
// 	int smoothexp = 1;
// 	MRF::CostVal smoothmax = 1;
// 	// Potts model
// 	return new SmoothnessCost(smoothexp, smoothmax, lambda, hCueArray, vCueArray);
//   }
//   else { // define smoothness cost using general function
	return fnCost;
    //  }
}


// delete any arrays allocated for data and smoothness costs
void deleteGlobalArrays() {
  delete [] dsiArray;


  for(vector<MRF::CostVal *>::reverse_iterator b = dsiArray2.rbegin(); b != dsiArray2.rend(); ++b){
	  free(*b);
  }

  for(vector<MRF::CostVal *>::reverse_iterator b = hCueArrayV.rbegin(); b != hCueArrayV.rend(); ++b){
	  free(*b);
  }

  for(vector<MRF::CostVal *>::reverse_iterator b = vCueArrayV.rbegin(); b != vCueArrayV.rend(); ++b){
	  free(*b);
  }

  for(vector<MRF::CostVal *>::reverse_iterator b = hBOcclArrayV.rbegin(); b != hBOcclArrayV.rend(); ++b){
	  free(*b);
  }

  for(vector<MRF::CostVal *>::reverse_iterator b = vBOcclArrayV.rbegin(); b != vBOcclArrayV.rend(); ++b){
	  free(*b);
  }

  for(vector<MRF::CostVal *>::reverse_iterator b = hSOcclArrayV.rbegin(); b != hSOcclArrayV.rend(); ++b){
	  free(*b);
  }

  for(vector<MRF::CostVal *>::reverse_iterator b = vSOcclArrayV.rbegin(); b != vSOcclArrayV.rend(); ++b){
	  free(*b);
  }

  for(vector<MRF::CostVal *>::reverse_iterator b = dsiArrayV.rbegin(); b != dsiArrayV.rend(); ++b){
	  free(*b);
  }
  /*delete [][] hCueArray;
  delete [][] vCueArray;
  delete [][] hBOcclArray;
  delete [][] vBOcclArray;
  delete [][] hSOcclArray;
  delete [][] vSOcclArray;*/
}

// main stereo matcher
// if WTAdisp is a valid image, it is used to initialize the labels

void crfstereo(int width, int height, int nD, DataCost *dcost, SmoothnessCost *scost,
               CByteImage WTAdisp, CByteImage &disp, float closeEnoughPercent)
{
  crfstereo (width, height, nD, dcost, scost, WTAdisp, disp, closeEnoughPercent, 
             USE_GC,1);
}

void crfstereo(int width, int height, int nD, DataCost *dcost, SmoothnessCost *scost,
               CByteImage WTAdisp, CByteImage &disp, float closeEnoughPercent, 
               int inferencer, int innerIter)
{

  if (innerIter == -1)
    innerIter = 5*width*height;

  EnergyFunction *energy = new EnergyFunction(dcost, scost);

  MRF* mrf;

  switch (inferencer){
  case USE_TRWS:
    mrf = new TRWS(width, height, nD, energy);
	DEBUG_OUT0(verbose, debugfile, "Creating TRWS\n");
	break;
  case USE_MF:
    mrf = new MeanField(width, height, nD, energy);
    DEBUG_OUT0(verbose, debugfile, "Creating MeanField\n");
    break;
  case USE_SMF:
    mrf = new SparseMeanField(width, height, nD, energy);
    DEBUG_OUT0(verbose, debugfile, "Creating SparseMeanField\n");
    break;
  case USE_BP:
    mrf = new SyncSumProd(width, height, nD, energy);
    DEBUG_OUT0(verbose, debugfile, "Creating SyncSumProd\n");
    break;
  case USE_SBP:
    {
      mrf = new SparseSyncSumProd(width, height, nD, energy);
      FloatType damper = 0.9;
      mrf->setParameters(PARAM_DAMPER, &damper);
      
      FloatType msgTol = 0.0001 * 2 * (width*(height-1) + (width-1)*height);
      mrf->setParameters(PARAM_MSGDTOL, &msgTol);
      
      FloatType divTol = 0.001;
      mrf->setParameters(PARAM_DIVTOL, &divTol);
      
      DEBUG_OUT0(verbose, debugfile, "Creating SparseSyncSumProd\n");
      break;
    }
  case USE_GC:
    mrf = new Expansion(width, height, nD, energy);
    DEBUG_OUT0(verbose, debugfile, "Creating Graph Cuts / Expansion\n");
    break;
  default:
    fprintf(stderr, "Unknown inferencer type");
    exit(1);
  }
  
  if(inferencer != USE_TRWS)
	mrf->setParameters(1,&randomv);
  mrf->initialize();

  if (WTAdisp.Shape().width == width)
    {
      // initialize labels to WTA disps
      int n = 0;

      for (int y = 0; y < height; y++) 
        {
          uchar *row = &WTAdisp.Pixel(0, y, 0);

          for (int x = 0; x < width; x++) 
            {
              mrf->setLabel(n++, row[x]);
            }
        }
    } 
  else 
    {
      mrf->clearAnswer();
    }

  delete energy;

  crfstereo(width, height, disp, mrf, closeEnoughPercent, innerIter );

  delete mrf;
}

void crfstereo(int width, int height, CByteImage &disp, MRF* mrf, 
               float closeEnoughPercent, int innerIter, int maxIter )
{
  MRF::EnergyVal E, Ed, Es, Eold;
  Ed = mrf->dataEnergy();	
  Es = mrf->smoothnessEnergy(); 
  E = Ed + Es; // mrf->totalEnergy();

  DEBUG_OUT3(verbose, debugfile, "Energy = %f (Ed=%g, Es=%g) at start\n", E, Ed, Es);

 

  Eold = E;

  float t, tot_t = 0;
  
  for (int iter = 0; iter < maxIter; iter++)
    {
      mrf->optimize(innerIter, t);
      tot_t += t ;
      
      Ed = mrf->dataEnergy(); 
      Es = mrf->smoothnessEnergy(); 
      E = Ed + Es; // mrf->totalEnergy();

      float percentDecrease = (float)(100.0*(1.0 - E/Eold));

      DEBUG_OUT5(verbose, debugfile, "Energy = %f (Ed=%g, Es=%g), %.3f secs, %.3f%% decrease\n", E, Ed, Es, tot_t, percentDecrease);

      if (percentDecrease < closeEnoughPercent)
	    break;
      
      if (E > Eold)
	    fprintf(stderr, "Warning: energy is increasing!\n");

      Eold = E;
    }
  
  // get disparities
  CShape sh(width, height, 1);
  disp.ReAllocate(sh);

  int n = 0;
  for (int y = 0; y < height; y++) 
    {
      uchar *row = &disp.Pixel(0, y, 0);
      for (int x = 0; x < width; x++) 
        {
          row[x] = mrf->getLabel(n++);
        }
    }
}

void writeDisparities(CByteImage disp, int outscale, char *outstem)
{
  CByteImage disp2;
  DEBUG_OUT1(verbose, debugfile, "scaling disparities by %d\n", outscale);
  ScaleAndOffset(disp, disp2, (float)outscale, 0);

  char dispname[1000];
  sprintf(dispname, "%s.png", outstem);
  DEBUG_OUT1(verbose, debugfile, "Writing image %s\n", dispname);
  WriteImage(disp2, dispname);
}

void setPairwise(bool v){
	pairwise = v;
}

void setPairwiseInteraction(bool v){
	pairwiseInteraction = v;
}

void setLocal(bool v){
	localOccl = v;
}

bool getLocal(){
	return localOccl;
}


void updateThetaP(std::vector<float> theta){
	int nP = theta.size();
	if(nP == 2){
		occluded1 = theta[1]; // Only one of the pair is occluded.
		occluded2 = theta[0]; // Both Pair are occluded.
	}
}

int getNumStates(){
	return globalNdisps;
}

float rangeNormailize(float rawMin, float rawMax, float min, float max, float real){
	if(real > rawMax){
		return max;
	}
	else if(real < rawMin){
		return min;
	}
	else{
		return (max - min) * (real - rawMin )/ ( rawMax - rawMin) + min;
	}
}

void updateCueIndex(unsigned int i){
	cueIndex = i;
}

void setRandom(bool rv){
	randomv = rv;
}
