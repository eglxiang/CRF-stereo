#ifndef CRFSTEREO_H
#define CRFSTEREO_H

// crfstereo.h
// $Id: crfstereo.h,v 1.8 2007/12/07 17:23:22 weinman Exp $

#include "mrf.h"
#include "imageLib.h"
#include <vector>

#define USE_GC 0
#define USE_BP 1
#define USE_SBP 2
#define USE_MF 3
#define USE_SMF 4
#define USE_TRWS 5

class ImageWriter;

// global variables for debugging output
extern int verbose;
extern FILE *debugfile;

#ifdef _MSC_VER        // Visual C++
// unfortunately, variable argument macros not supported until VCC 2005...
#define DEBUG_OUT0(verbose, debugfile, fmt) { \
  if (verbose) fprintf(stderr, fmt); \
  if (debugfile) fprintf(debugfile, fmt); }
#define DEBUG_OUT1(verbose, debugfile, fmt, arg1) { \
  if (verbose) fprintf(stderr, fmt, arg1); \
  if (debugfile) fprintf(debugfile, fmt, arg1); }
#define DEBUG_OUT2(verbose, debugfile, fmt, arg1, arg2) { \
  if (verbose) fprintf(stderr, fmt, arg1, arg2); \
  if (debugfile) fprintf(debugfile, fmt, arg1, arg2); }
#define DEBUG_OUT3(verbose, debugfile, fmt, arg1, arg2, arg3) { \
  if (verbose) fprintf(stderr, fmt, arg1, arg2, arg3); \
  if (debugfile) fprintf(debugfile, fmt, arg1, arg2, arg3); }
#define DEBUG_OUT4(verbose, debugfile, fmt, arg1, arg2, arg3, arg4) { \
  if (verbose) fprintf(stderr, fmt, arg1, arg2, arg3, arg4); \
  if (debugfile) fprintf(debugfile, fmt, arg1, arg2, arg3, arg4); }
#define DEBUG_OUT5(verbose, debugfile, fmt, arg1, arg2, arg3, arg4, arg5) { \
  if (verbose) fprintf(stderr, fmt, arg1, arg2, arg3, arg4, arg5); \
  if (debugfile) fprintf(debugfile, fmt, arg1, arg2, arg3, arg4, arg5); }
#else
// easier in g++:
#define DEBUG_OUT(verbose, debugfile, args...) { \
  if (verbose) fprintf(stderr, args); \
  if (debugfile) fprintf(debugfile, args); }
#define DEBUG_OUT0 DEBUG_OUT
#define DEBUG_OUT1 DEBUG_OUT
#define DEBUG_OUT2 DEBUG_OUT
#define DEBUG_OUT3 DEBUG_OUT
#define DEBUG_OUT4 DEBUG_OUT
#define DEBUG_OUT5 DEBUG_OUT
#endif

// prototypes


void setGlobalNpixels(vector<int>& w, vector<int>& h);
void setGlobalWidth(vector<int>& v);

// make sure MRF library is compiled with float costs
void checkForFloatCosts();

// compute global data cost array - matching cost volume (absolute color diffs)
// also compute best data cost index per pixel (i.e. WTA disparities)
void initializeDataCost(CByteImage im1, CByteImage im2, int nD,  bool birchfield,  bool squaredDiffs, int truncDiffs,  CByteImage &wta); 

void initializeDataCostVector(vector<CByteImage>, vector<CByteImage>, int, bool, bool, int, vector<CByteImage> &, int , float);
// create data cost from global cost array
DataCost *computeDataCost(int index);
MRF::CostVal *computeDataCostArray(int index);

// create data cost using parameters thetaU
DataCost *computeDataCost(std::vector<float> thetaU, int index);
MRF::CostVal *computeDataCostArray(std::vector<float> thetaU, int index);

// get pointer to cost array and its dimensions
void getCostArray(float* &costs, int &nPixels, int &nDisps);

// compute quantized absolute color gradients of input image
void computeQuantizedGradients(CByteImage im1, CByteImage &im1Grad, std::vector<int> gradThresh, int bs);
/*
  use nK = gradThresh.size()+1 levels of quantization.
  quantization is done as follows:
  compute absolute gradient (Euclidean distance in color space) in x and y
  given gradient g, assign largest k such that g <= gradThresh[k]
  entry gradThresh[nK] is assumed to be infinity
  in im1Grad, use band 0 for x gradient, and band 1 for y gradients
*/

SmoothnessCost *computeSmoothnessCost(CByteImage im1grad, std::vector<float> theta);
MRF::SmoothCostGeneralFn computeSmoothnessCostFunction(CByteImage im1grad, std::vector<float> theta);
// I don't know how to use a default object argument, so I leave the above two functions in  - JJW
SmoothnessCost *computeSmoothnessCost(CByteImage im1grad, std::vector<float> theta, std::vector<float> thetaP, int i);
MRF::SmoothCostGeneralFn computeSmoothnessCostFunction(CByteImage im1grad, std::vector<float> theta, std::vector<float> thetaP, int i);

/*
the parameters that determine the smoothness cost are stored in a
vector so they can be handled uniformly during learning:

  std::vector<float> theta;

  int nK = theta.size();
  //theta[0]: cost for |dp-dq| == 1
  //theta[k], k=1..nK-1: cost for |dp-dq| > 1 and gradient_pq < gradThresh[k]

OOOPS, I just realized that the MRF library doesn't support this, except if you set
up the smoothness cost as a general function, which is much slower.  If you
use the hCue and vCue parameters, then the smoothness cost (which only depends
on the two labels) simply gets multiplied by constants...

I can modify the graph cut code later, but for now let's stick with the old model:
  //theta[k], k=0..nK-1: cost for |dp-dq| > 0 and gradient_pq < gradThresh[k]
*/

// delete any arrays allocated for data and smoothness costs
void deleteGlobalArrays();


// main stereo matcher
// if WTAdisp is a valid image, it is used to initialize the labels
void crfstereo(int width, int height, int nD, DataCost *dcost, SmoothnessCost *scost, 
	       CByteImage WTAdisp, CByteImage &disp, float closeEnoughPercent);

void crfstereo(int width, int height, int nD, DataCost *dcost, SmoothnessCost *scost, 
               CByteImage WTAdisp, CByteImage &disp, float closeEnoughPercent, 
               int inferencer, int innerIter = -1);

void crfstereo(int width, int height, CByteImage &disp, MRF *mrf, 
               float closeEnoughPercent, int innerIter = 1, 
               int maxIter = 50);


void writeDisparities(CByteImage disp, int outscale, char *outstem);

void setPairwise(bool v);
void setPairwiseInteraction(bool v);
void setLocal(bool v);
void updateThetaP(std::vector<float> theta);
int getNumStates();
float rangeNormailize(float rawMin, float rawMax, float min, float max, float real);

void updateCueIndex(unsigned int i);
bool getLocal();
void setRandom(bool rv);

#endif
