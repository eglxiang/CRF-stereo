// crfmain.cpp
// console application to run crfstereo

static char usage[] = "\n crfstereo version 1.1 \n\
\n\
 usage: %s [options] imL imR outstem \n\
\n\
  reads imL and imR (in png or pgm/ppm format)\n\
  runs MRF stereo\n\
  writes disparities corresponding imL to outstem.png and parameters to outstem.txt\n\
\n\
  options:\n\
    -n nD           disparity levels, by default 16 (i.e. disparites 0..15)\n\
    -a              use graph cuts minimization\n\
    -f              use mean field inference\n\
    -c              initialize mean field at graph cuts\n\
    -p              use belief propagation inference\n\
    -k              use sparse message passing\n\
    -b              use Birchfield/Tomasi costs\n\
    -s              use squared differences (absolute differences by default)\n\
    -t trunc        truncate differences to <= 'trunc'\n\
    -u u1           data cost param (thetaU) - can specify one or none\n\
    -v v1 -v v2 ... costs (thetaV's) for each bin (can specify nV >= 0)\n\
    -g t1 -g t2 ... gradient thresholds between bins (must specify nV-1) \n\
    -x closeEnough  sufficiently small energy decrease percentage to terminate GC/BP\n\
    -o outscale     scale factor for disparities (full range by default)\n\
    -q              quiet (turn off debugging output)\n\
\n\
  Example usage:\n\
  crfstereo  -v10 -v10 -v10 -g4 -g8 -n20 -o8 scenes/venus/im* scenes/venus/truedisp.png outvenus\n\
\n";

#include <stdio.h>
#include <fstream>

#ifdef _MSC_VER        // Visual C++
# include "XGetopt.h"  // replacement for getopt()
#else
# include <getopt.h>
#endif

#include "crfstereo.h"
#include "expectations.h"
#include "evaldisps.h"

#include "GCoptimization.h"
#include "AsyncSumProd.h"
#include "MeanField.h"
#include "SparseAsyncSumProd.h"
#include "SparseMeanField.h"

#include "LabelImageWriter.h"

// global variables for debugging output
int verbose;
FILE *debugfile;
    
int main(int argc, char **argv)
{

    // make sure MRF library is compiled with float costs
    checkForFloatCosts();

    // parameters controlled via command-line options:
    int inferencer = USE_BP;
    bool useSparse = false;
    int nD = 20;               // disparity levels (d = 0 .. nD-1)
    bool birchfield = false;   // use Birchfield/Tomasi costs
    bool useGC = false;         // initialize MF at graph cuts
    bool squaredDiffs = false; // use squared differences (absolute differences by default)
    int truncDiffs = 255;      // truncated differences (before squaring), by default not
    //float lambda1 = 20;         // weight of smoothness term (20 by default)
    //float lambda2 = 20;         // weight of smoothness term (20 by default)
    //int gradThresh = 8;        // intensity gradient cue threshold (8 by default)
    //float gradPenalty = 2;     // if grad < gradThresh, multiply smoothness cost by this
    int maxInner = 1;          // maximum number of inner iterations 
    int maxOuter = 50;          // maximum number of outer iterations 
    float closeEnoughPercent = 2.0f; // when to stop GC or BP iterations
    int outscale = -1;         // scale factor for disparities; -1 means full range 255.0/(nD-1)
    verbose = 1;               // print messages to stderr

    std::vector<int> gradThreshVec;
    fvec thetaV; // smoothness weights associated with gradient bins
    fvec thetaU; // data weights associated with gradient bins

    // parse command-line arguments using the "getopt" utility
    int o;
    while ((o = getopt(argc, argv, "acn:fpkbst:u:v:g:r:i:x:o:q")) != -1)
	switch (o) {
	case 'n': nD = atoi(optarg); break;
    case 'a': inferencer = USE_GC; break;
    case 'f': inferencer = USE_MF; break;
    case 'c': useGC = true; break;
    case 'p': inferencer = USE_BP; break;
    case 'k': useSparse = true; break;
	case 'b': birchfield = 1; break;
	case 's': squaredDiffs = 1; break;
	case 't': truncDiffs = atoi(optarg); break;
 	case 'u': thetaU.push_back((float)atof(optarg)); break;
	case 'v': thetaV.push_back((float)atof(optarg)); break;
	case 'g': gradThreshVec.push_back(atoi(optarg)); break;
	case 'x': closeEnoughPercent = (float)atof(optarg); break;
	case 'o': outscale = atoi(optarg); break;
	case 'q': verbose = 0; break;
	default: 
	    fprintf(stderr, "Ignoring unrecognized option %s\n", argv[optind-1]);
	}

    if (useSparse == true)
        inferencer = ( inferencer== USE_BP ? USE_SBP : USE_SMF);
    
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
    char logname[1000];
	sprintf(debugname, "%s.txt", outstem);
    sprintf(logname,"%s-log.txt",outstem);

	debugfile = fopen(debugname, "w");

    ofstream logStream(logname);

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
	initializeDataCost(im1, im2, nD, birchfield, squaredDiffs, truncDiffs, WTAdisp);

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
	int nT = (int)gradThreshVec.size(); // number of thresholds
	int nV = (int)thetaV.size(); // number of smoothness parameters

	if (nV > 0 && nV != nT+1)
	    throw CError("Must give exactly one more thetaV values than thresholds");

	// compute quantized image gradients
	computeQuantizedGradients(im1, im1grad, gradThreshVec);

	// set up the smoothness cost
	//theta[0] = lambda1; // for delta_d == 1
	//theta[1] = lambda2; // for delta_d >= 2

	fvec empirDistV, modelDistV;
	empirDistV.resize(nV);
	modelDistV.resize(nV);

	if (nU + nV < 1)
	    throw CError("Must give at least on thetaV or thetaU");

	CByteImage previousDisp;
	CopyPixels(WTAdisp, previousDisp);


	// CPAL - added to filter edges
	//computeValidDistImage(truedisp);
	CByteImage tmpDisp;
	CopyPixels(WTAdisp, previousDisp);

    // Create MRF object
    MRF::CostVal* dcostArr = computeDataCostArray(thetaU);
    MRF::SmoothCostGeneralFn scostFn = computeSmoothnessCostFunction
      (im1grad, thetaV);

    EnergyFunction *energy = new EnergyFunction (new DataCost(dcostArr),
                                                 new SmoothnessCost (scostFn) );

    MRFEnergy *mrf;
    
    switch (inferencer)
      {
      case USE_GC:
        break;
      case USE_MF:
        {
          mrf = new MeanField(width, height, nD, energy);
          
          FloatType damper = 1;
          mrf->setParameters(PARAM_DAMPER, &damper);
          
          FloatType msgTol = 0.001;//0.001 * width * height;
          mrf->setParameters(PARAM_MSGDTOL, &msgTol);

          MeanField::DiffType msgDiffFun = MeanField::ENERGY;

          printf("msgDiffFun = %d\n",msgDiffFun);

          mrf->setParameters(PARAM_MSGDFUN, &msgDiffFun);

          maxInner = 10;
          maxOuter = 1;
          

          DEBUG_OUT0(verbose, debugfile, "Creating MeanField\n");
          DEBUG_OUT2(verbose, debugfile, "damper = %0.1f  tol = %g\n",
                     damper,msgTol);

          break;
        }
      case USE_SMF:
        {
          mrf = new SparseMeanField(width, height, nD, energy);

          FloatType damper = 1;
          mrf->setParameters(PARAM_DAMPER, &damper);
          
          FloatType msgTol = 0.001;//0.001 * width * height;
          mrf->setParameters(PARAM_MSGDTOL, &msgTol);
          
          MeanField::DiffType msgDiffFun = MeanField::ENERGY;
          mrf->setParameters(PARAM_MSGDFUN, &msgDiffFun);

          FloatType divTol = 0.001;
          mrf->setParameters(PARAM_DIVTOL, &divTol);

          maxInner = 500;
          maxOuter = 1;

          DEBUG_OUT0(verbose, debugfile, "Creating SparseMeanField\n");
          DEBUG_OUT3(verbose, debugfile, "damper = %0.1f  tol = %g  div = %g\n",
                     damper,msgTol,divTol);
          break;
        }
      case USE_BP:
        {
          mrf = new AsyncSumProd(width, height, nD, energy);
          
          FloatType damper = 1.0;
          mrf->setParameters(PARAM_DAMPER, &damper);
          
          FloatType msgTol = 1;//0.01 * 2 * (width*(height-1) + (width-1)*height);
          mrf->setParameters(PARAM_MSGDTOL, &msgTol);

          maxInner = 1;//5*width*height;
          maxOuter = 50;
          
          DEBUG_OUT0(verbose, debugfile, "Creating AsyncSumProd\n");
          DEBUG_OUT2(verbose, debugfile, "damper = %0.1f  tol = %g\n",damper,msgTol);
          break;
        }
      case USE_SBP:
        {
          mrf = new SparseAsyncSumProd(width, height, nD, energy);
          
          FloatType damper = 1.0;
          mrf->setParameters(PARAM_DAMPER, &damper);
          
          FloatType msgTol = 1;//0.0001 * 2 * (width*(height-1) + (width-1)*height);
          mrf->setParameters(PARAM_MSGDTOL, &msgTol);

          FloatType divTol = 1e-6;
          mrf->setParameters(PARAM_DIVTOL, &divTol);
          
          maxInner = 1;//5*width*height;
          maxOuter = 50;
          
          DEBUG_OUT0(verbose, debugfile, "Creating SparseAsyncSumProd\n");
          DEBUG_OUT3(verbose, debugfile, "damper = %0.1f  tol = %g  div = %g\n",
                     damper,msgTol,divTol);
          break;
        }
      default:
        fprintf(stderr, "Unknown inferencer type");
        return -1;
      }

    if (inferencer!=USE_GC)
      {
        mrf->initialize();
        mrf->setLog(&logStream);
      }

    delete energy;

          // Initialize mean field beliefs with graph-cuts
    if (inferencer==USE_GC || 
        (useGC && (inferencer==USE_MF || inferencer==USE_SMF)))
        {
          EnergyFunction *gcenergy = new EnergyFunction 
            (new DataCost(dcostArr), new SmoothnessCost (scostFn) );

          MRF *gcmrf = new Expansion(width,height,nD,gcenergy);
          
          bool useRand = false;
          gcmrf->setParameters(1,&useRand);

          gcmrf->initialize();

          crfstereo(width,height,disp,gcmrf,closeEnoughPercent,1,50);

          MRF::Label *gclabel = gcmrf->getAnswerPtr();

          if (inferencer!=USE_GC)
            {
              MRF::Label *mflabel = mrf->getAnswerPtr();
              
              int nPixels = width*height; 
              
              for (int i=0 ; i<nPixels ; i++, gclabel++, mflabel++)
                *mflabel = *gclabel;
            }

          delete energy;
          delete gcmrf;

          if (inferencer==USE_MF)
            {
              MeanField *mfmrf = dynamic_cast<MeanField*>(mrf);
              mfmrf->initializeAlg(MeanField::LABELPOINT);
            }
          else if (inferencer==USE_SMF)
            {
              SparseMeanField *mfmrf = 
                dynamic_cast<SparseMeanField*>(mrf);
              mfmrf->initializeAlg(MeanField::LABELPOINT);
            }
          

        } else if(0) {

          if (inferencer==USE_MF)
            {
              MeanField *mfmrf = dynamic_cast<MeanField*>(mrf);
              mfmrf->initializeAlg(MeanField::LOCAL);
            }
          else if (inferencer==USE_SMF)
            {
              SparseMeanField *mfmrf = 
                dynamic_cast<SparseMeanField*>(mrf);
              mfmrf->initializeAlg(MeanField::LOCAL);
            }
        }

    if (inferencer!=USE_GC)
      {
        char iterstem[1000];
        sprintf(iterstem,"%s-inneriter",outstem);
        LabelImageWriter imWriter(iterstem,width,height,outscale);
        mrf->setLabelImageWriter(&imWriter);
        
        // run the stereo matcher 
        crfstereo(width, height, disp, mrf, closeEnoughPercent, 
                  maxInner, maxOuter);
      }

    // evaluate matching errors
    float bad1=0, rms=0;
    CByteImage errormap;
    evaldisps(disp, truedisp, errormap, bad1, rms, 1);
    DEBUG_OUT2(verbose, debugfile, "bad1= %g   rms= %g\n", bad1, rms);
    
    
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
