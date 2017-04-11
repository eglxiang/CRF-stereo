// gradientdescent.cpp
// $Id: gradientdescentMPE1.cpp,v 1.1.1.1 2007/12/01 01:01:14 cpal Exp $
//
// finds best parameters theta using gradient descent
//
// Daniel Scharstein

// if you want to have a thetaV parameter for equal neighboring disparities, 
// define HAVE_COST_FOR_DELTA_D_EQUAL_ZERO in crfstereo.h

// mask occluded areas in saved disparity maps
#define WRITE_MASKED_DISPS 0

static char usage[] = "\n\
 usage: %s [options] imL imR truedisp outstem \n\
\n\
  reads imL, imR, and true disparities (in png or pgm/ppm format)\n\
  runs MRF stereo\n\
  writes disparities corresponding imL to outstem.png and parameters to outstem.txt\n\
\n\
  options:\n\
    -n nD           disparity levels, by default 16 (i.e. disparites 0..15)\n\
    -b              use Birchfield/Tomasi costs\n\
    -s              use squared differences (absolute differences by default)\n\
    -t trunc        truncate differences to <= 'trunc'\n\
    -u u1           initial data cost param (thetaU) - can specify one or none\n\
	-c c1 -c c2 ... initial costs for occluded states\n\
    -v v1 -v v2 ... initial costs (thetaV's) for each bin (can specify nV >= 0)\n\
    -g t1 -g t2 ... gradient thresholds between bins (must specify nV-1) \n\
    -r rate         learning rate (for thetaVs) during gradient descent\n\
    -f factorU      factor for learning rate for thetaU\n\
    -z sigma        use Gaussian regularization with scale sigma\n\
    -i maxiter      maximum number of iterations for gradient descent\n\
    -x closeEnough  sufficiently small energy decrease percentage to terminate GC/BP\n\
    -o outscale     scale factor for disparities (full range by default)\n\
    -q              quiet (turn off debugging output)\n\
\n\
  Example usage:\n\
  gradientdescent -i10 -b -r.0001 -v10 -v10 -v10 -g4 -g8 -n20 -o8 scenes/venus/im* scenes/venus/truedisp.png outvenus\n\
\n";

// -l and -m options currently not used:
//    -l lambda1     weight of smoothness term for delta_d==1 (20 by default)\n
//    -m lambda2     weight of smoothness term for delta_d>=2  (20 by default)\n


#include <stdio.h>
#ifdef _MSC_VER        // Visual C++
# include "XGetopt.h"  // replacement for getopt()
#else
# include <getopt.h>
#endif
#include "crfstereo.h"
#include "expectations.h"
#include "evaldisps.h"
#include <iostream>

// global variables for debugging output
int verbose;
FILE *debugfile;

int main(int argc, char **argv)
{
    // make sure MRF library is compiled with float costs
    checkForFloatCosts();

    // parameters controlled via command-line options:
    int nD = 20;               // disparity levels (d = 0 .. nD-1)
    bool birchfield = false;   // use Birchfield/Tomasi costs
    bool squaredDiffs = false; // use squared differences (absolute differences by default)
    int truncDiffs = 255;      // truncated differences (before squaring), by default not
    //float lambda1 = 20;         // weight of smoothness term (20 by default)
    //float lambda2 = 20;         // weight of smoothness term (20 by default)
    //int gradThresh = 8;        // intensity gradient cue threshold (8 by default)
    //float gradPenalty = 2;     // if grad < gradThresh, multiply smoothness cost by this
    float rate = 0.001f;
    float gaussSigma = 0;      // Scale of gaussian  regularizer (zero for none)
    float factorU = 1;
    int maxiter = 100;          // maximum number of iterations for gradient descent
    float closeEnoughPercent = 2.0f; // when to stop GC or BP iterations
    int outscale = -1;         // scale factor for disparities; -1 means full range 255.0/(nD-1)
    verbose = 1;               // print messages to stderr

    std::vector<int> gradThreshVec;
    fvec thetaU; // data weights (at most one for now)
    fvec thetaV; // smoothness weights associated with gradient bins
    fvec theta;  // concatenation of the two
	fvec thetaO;  // Lam

    // parse command-line arguments using the "getopt" utility
    int o;
	while ((o = getopt(argc, argv, "n:bst:c:p:u:v:g:r:f:i:x:o:z:q")) != -1)
	switch (o) {
	case 'n': nD = atoi(optarg); break;
	case 'b': birchfield = 1; break;
	case 's': squaredDiffs = 1; break;
	case 't': truncDiffs = atoi(optarg); break;
	    //case 'l': lambda1 = atof(optarg); break;
	    //case 'm': lambda2 = atof(optarg); break;
	case 'u': thetaU.push_back((float)atof(optarg)); break;
	case 'p': thetaO.push_back((float)atof(optarg)); break;
	case 'v': thetaV.push_back((float)atof(optarg)); break;
	case 'g': gradThreshVec.push_back(atoi(optarg)); break;
	case 'r': rate = (float)atof(optarg); break;
	case 'f': factorU = (float)atof(optarg); break;
	case 'i': maxiter = atoi(optarg); break;
	case 'x': closeEnoughPercent = (float)atof(optarg); break;
	case 'o': outscale = atoi(optarg); break;
	case 'q': verbose = 0; break;
    case 'z': gaussSigma = (float)atof(optarg); break;
	default: 
	    fprintf(stderr, "Ignoring unrecognized option %s\n", argv[optind-1]);
	}

    if (optind != argc-4) {
	fprintf(stderr, usage, argv[0]);
	exit(1);
    }

    char *im1name = argv[optind++];
    char *im2name = argv[optind++];
    char *truedispname = argv[optind++];
    char *outstem = argv[optind++];

    if (outscale < 0)
	outscale = 255 / (nD-1);

    try {
	char debugname[1000];
	sprintf(debugname, "%s.txt", outstem);
	debugfile = fopen(debugname, "w");
	if (debugfile == NULL) 
	    throw CError("Cannot write to %s", debugname);
	if (verbose)
	    fprintf(stderr, "writing parameters to %s\n", debugname);
	for (int i=0; i<argc; i++)
	    fprintf(debugfile, "%s ", argv[i]);
	fprintf(debugfile, "\n");

	CByteImage im1, im2;      // input images (gray or color)
	CByteImage truedisp;       // quantized gradients of image 1
	CByteImage im1grad;       // quantized gradients of image 1
	CByteImage WTAdisp;       // WTA disparities (from matching cost only)
	CByteImage disp;          // output disparities
	DEBUG_OUT1(verbose, debugfile, "Reading image %s\n", im1name);
	ReadImage(im1, im1name);
	DEBUG_OUT1(verbose, debugfile, "Reading image %s\n", im2name);
	ReadImage(im2, im2name);
	DEBUG_OUT1(verbose, debugfile, "Reading image %s\n", truedispname);
	ReadImage(truedisp, truedispname);
	ScaleAndOffset(truedisp, truedisp, (float)1.0/outscale, 0); // scale to integer disps
  
	CShape sh = im1.Shape();
	if (sh != im2.Shape())
	    throw CError("image shapes don't match");
	int width = sh.width, height = sh.height;

	// ******** data cost

	// initialize data cost (also computes WTA disparities)
	initializeDataCost(im1, im2, nD, birchfield, squaredDiffs, truncDiffs, WTAdisp, thetaU[0]);

	printfVec(thetaU, "thetaU");
	int nU = (int)thetaU.size(); // number of data parameters (0 or 1 for now)

	if (nU > 1)
	    throw CError("Must give 1 or 0 thetaU values");

	fvec empirDistU, modelDistU;
	empirDistU.resize(nU);
	modelDistU.resize(nU);
	
	// ******** smoothness cost

	printfVec(gradThreshVec, "thresholds");
	printfVec(thetaV, "thetaV");
#ifdef HAVE_COST_FOR_DELTA_D_EQUAL_ZERO
	DEBUG_OUT0(verbose, debugfile, "last thetaV is for deltaD==0\n");
#endif
	int nT = (int)gradThreshVec.size(); // number of thresholds
	int nV = (int)thetaV.size(); // number of smoothness parameters

#ifdef HAVE_COST_FOR_DELTA_D_EQUAL_ZERO
	if (nV != nT+2)
	    throw CError("Must give exactly TWO more thetaV values than thresholds (last one for deltaD==0)");
#else
	if (nV > 0 && nV != nT+1)
	    throw CError("Must give exactly one more thetaV values than thresholds");
#endif

    DEBUG_OUT1(verbose,debugfile,"Gaussian Sigma = %f\n",gaussSigma);
    
	// compute quantized image gradients
	computeQuantizedGradients(im1, im1grad, gradThreshVec);

	// set up the smoothness cost
	//theta[0] = lambda1; // for delta_d == 1
	//theta[1] = lambda2; // for delta_d >= 2

	fvec empirDistV, modelDistV;
	

	if (nU + nV < 1)
	    throw CError("Must give at least on thetaV or thetaU");

	// Lam
	for( fvec::iterator b = thetaO.begin(); b != thetaO.end(); ++b)
		thetaV.push_back(*b);
	nV = (int)thetaV.size();
	// end Lam

	empirDistV.resize(nV);
	modelDistV.resize(nV);

	writeDisparities(truedisp, outscale, outstem);
	CByteImage tmpDisp;
	CByteImage previousDisp, maskeddisp;
	CopyPixels(WTAdisp, previousDisp);

	float oldnorm = -1;
	fvec oldTheta;

	concatVec(thetaU, thetaV, theta);
	printfVec(theta, "*********************** initial theta ");

	for (int iter=0; iter < maxiter; iter++) {

	    // split theta into U and V components
	    splitVec(theta, thetaU, thetaV);

	    // run the stereo matcher, starting from the previous solution
	    DataCost *dcost = computeDataCost(thetaU);
	    SmoothnessCost *scost = computeSmoothnessCost(im1grad, thetaV);
	    crfstereo(width, height, nD + 1, dcost, scost, previousDisp, disp, closeEnoughPercent);
	    delete scost;
	    delete dcost;

	    // evaluate matching errors
	    float bad1=0, rms=0;
	    CByteImage errormap;
	    //evaldisps(disp, truedisp, errormap, maskeddisp, bad1, rms, 1);
        evaldisps(disp, truedisp, errormap, bad1, rms, 1);
	    DEBUG_OUT2(verbose, debugfile, "bad1= %g   rms= %g\n", bad1, rms);
	    
	    if (maxiter > 1) {
		char outname[1000];
		sprintf(outname, "%s-iter%03d", outstem, iter);
		writeDisparities(disp, outscale, outname);
	    }

	    // compute empirical data distribution (of truedisp)
	    computeEmpirDistU(truedisp, truedisp, empirDistU);
	    // compute model data distribution (of computed disp)
	    computeModelDistU(disp, truedisp, modelDistU);

	    // compute empirical smoothness distribution (of truedisp)
	    computeEmpirDistV(truedisp, truedisp, im1grad, empirDistV);
	    // compute model smoothness distribution (of computed disp)
	    computeModelDistV(disp, truedisp, im1grad, modelDistV);

	    fvec empirDist, modelDist, diffDist, relDiff, regTheta,gradTheta;
	    concatVec(empirDistU, empirDistV, empirDist);
	    concatVec(modelDistU, modelDistV, modelDist);

	    vecDiff(modelDist, empirDist, diffDist);
	    printfVec(empirDist, "empirDist (GT)");
	    printfVec(modelDist, "modelDist     ");
	    printfVec(diffDist,  "diffDist      ");
        
        if (gaussSigma!=0)
          {
            vecCopy(theta,regTheta);
            vecScale(regTheta,-1/(gaussSigma*gaussSigma));
            vecSum(diffDist,regTheta,gradTheta);
            
            printfVec(regTheta,  "gaussReg      ");
            printfVec(gradTheta,  "gradient      ");
          }
        else
          vecCopy(diffDist,gradTheta);

	    // terminate if gradient is small enough percentage of modeldist
	    vecEltWiseQuot(gradTheta, modelDist, relDiff);

	    float norm = vecNorm(relDiff);
	    DEBUG_OUT1(verbose, debugfile, "norm of gradient = %.1f%%\n", 100 * norm);

	    // DON'T terminate to show overfitting!
	    if (100*norm < 0.5) // terminate when norm of diff is less than 1 percent
          break;

	    if (oldnorm < 0 || norm < 2*oldnorm) {
		// if current solution is better than previous one,
		// we're fine, keep going with bigger and bigger steps in same direction

		vecCopy(theta, oldTheta); // remember current point

		// update theta
		vecScale(gradTheta, rate);
		if (nU == 1) // param 0 is thetaU
		    gradTheta[0] *= factorU; // use different scale for thetaU
		printfVec(gradTheta, " scaled grad   ");
		// theta = theta + step * gradTheta
		vecSum(theta, gradTheta, theta);

		// make sure thetas don't go to or below 0
		vecBound(theta, 1e-3f, 1e10f);

		oldnorm = norm;
		rate *= 1.2f;
		DEBUG_OUT1(verbose, debugfile, "rate = %g\n", rate);

	    } else { // norm > 2*oldnorm
		// if norm gets much bigger, go back half-way and try again with much smaller step size
		vecSum(theta, oldTheta, theta);  //theta = (theta + oldTheta)/2
		vecScale(theta, 0.5);

		rate *= 0.25;
		DEBUG_OUT1(verbose, debugfile, "*** BACKING UP, rate = %g\n", rate);
	    }

	    printfVec(theta, "*********************** new theta ");

	    // start next iteration with current solution
	    CopyPixels(disp, previousDisp);

		if(maxiter > 1){
			char outname[1000];
			CopyPixels(disp, tmpDisp);
			sprintf(outname, "%s-masked-iter%03d", outstem, iter);
			computeMaskedImage(tmpDisp, truedisp);
			writeDisparities(tmpDisp, outscale, outname);

		}
	}

	if (WRITE_MASKED_DISPS)
	    writeDisparities(maskeddisp, outscale, outstem);
	else
	    writeDisparities(disp, outscale, outstem);

	fclose(debugfile);

	deleteGlobalArrays();
    }
    catch (CError &err) {
	fprintf(stderr, err.message);
	fprintf(stderr, "\n");
	exit(1);
    }
    catch (bad_alloc) {
	fprintf(stderr, "*** Error: not enough memory\n");
	exit(1);
    }

    return 0;
}
