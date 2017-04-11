// gradientdescent.cpp
// $Id: gradientdescent.cpp,v 1.13 2007/12/11 13:28:48 weinman Exp $
//
// finds best parameters theta using gradient descent
//
// Daniel Scharstein

static char usage[] = "\n\
 usage: %s [options] imL imR truedisp outstem \n\
\n\
  reads imL, imR, and true disparities (in png or pgm/ppm format)\n\
  runs MRF stereo\n\
  writes disparities corresponding imL to outstem.png and parameters to outstem.txt\n\
\n\
  options:\n\
    -f f1			scale factor \n\
    -n nD           disparity levels, by default 16 (i.e. disparites 0..15)\n\
    -d              value of true disparity to mask/ignore\n\
	-l 1/0			1 read images from directory, 0 read single image. \n\
    -m              use mean field inference\n\
    -c              initialize mean field at graph cuts\n\
    -e              don't use graph cuts for first iteration\n\
    -t              use belief propagation inference\n\
    -k              use sparse message passing\n\
    -b              use Birchfield/Tomasi costs\n\
    -s              use squared differences (absolute differences by default)\n\
    -t trunc        truncate differences to <= 'trunc'\n\
    -u u1           initial data cost param (thetaU) - can specify one or none\n\
    -v v1 -v v2 ... initial costs (thetaV's) for each bin (can specify nV >= 0)\n\
    -g t1 -g t2 ... gradient thresholds between bins (must specify nV-1) \n\
	-p p1 -p p2     pairwise cost (must specify 2) \n\
    -r rate         learning rate during gradient descent\n\
    -z sigma        use Gaussian regularization with scale sigma\n\
    -i maxiter      maximum number of iterations for gradient descent\n\
    -x closeEnough  sufficiently small energy decrease percentage to terminate GC/BP\n\
    -o outscale     scale factor for disparities (full range by default)\n\
    -q              quiet (turn off debugging output)\n\
	-h 1/0			1 to enable random and 0 to disable random \n\
	-j				run with additional features. \n\
					0. No occlusion \n\
					wxyz, w bit activate occlusion model, x bit activate local pairwise,\n\
					y bit activate gradient modulated, z bit activate learning rate \n\
\n\
  Example usage:\n\
  gradientdescent -i10 -b -r.0001 -v10 -v10 -v10 -g4 -g8 -n20 -o8 scenes/venus/im* scenes/venus/truedisp.png outvenus\n\
\n";

// -l and -m options currently not used:
//    -l lambda1     weight of smoothness term for delta_d==1 (20 by default)\n
//    -m lambda2     weight of smoothness term for delta_d>=2  (20 by default)\n


#include <stdio.h>
#include <fstream>
#include <limits>

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
#include <iostream>
#include "helperLib.h"

#include "LabelImageWriter.h"
#include <iostream>
#include <algorithm>
#include "helperLib.h"

// global variables for debugging output
extern int verbose;
extern FILE *debugfile;
    
int main(int argc, char **argv)
{
    // make sure MRF library is compiled with float costs
    checkForFloatCosts();

    // parameters controlled via command-line options:
    int inferencer = USE_BP;
    bool useSparse = false;
    int nD = 20;               // disparity levels (d = 0 .. nD-1)
	int nP = 0;
    int nF = 0;
    int ignoreVal = -1;        // Value of ground truth disp image to ignore
    bool gcExceptFirst = false; 
    bool birchfield = false;   // use Birchfield/Tomasi costs
    bool useGC = false;         // initialize MF at graph cuts
    bool squaredDiffs = false; // use squared differences (absolute differences by default)
    int truncDiffs = 255;      // truncated differences (before squaring), by default not
    //float lambda1 = 20;         // weight of smoothness term (20 by default)
    //float lambda2 = 20;         // weight of smoothness term (20 by default)
    //int gradThresh = 8;        // intensity gradient cue threshold (8 by default)
    //float gradPenalty = 2;     // if grad < gradThresh, multiply smoothness cost by this
    float rate = 0.001f;
    float gaussSigma = 0;      // Scale of gaussian  regularizer (zero for none)
    int maxiter = 100;          // maximum number of iterations for gradient descent
    float closeEnoughPercent = 2.0f; // when to stop GC or BP iterations
    int outscale = -1;         // scale factor for disparities; -1 means full range 255.0/(nD-1)
    verbose = 1;               // print messages to stderr
	fvec factorU;
	int methods = 0;		   // use to activate new added features.		
	int dirimage = 0;		   // 1 to read images from directory, 0 to read single image file.
	int random = 0;			   // use to enable random in graphcut
	bool useRand = false;
	bool testImage = false;

    std::vector<int> gradThreshVec;
    fvec thetaV; // smoothness weights associated with gradient bins
    fvec thetaU; // data weights associated with gradient bins
	fvec thetaP; // data weights associated with gradient bins

    // parse command-line arguments using the "getopt" utility
    int o;
	while ((o = getopt(argc, argv, "cd:n:emf:w:p:kbst:u:v:l:h:j:g:r:i:x:o:qz:")) != -1)
	switch (o) {
	case 'n': nD = atoi(optarg); break;
    case 'd': ignoreVal = atoi(optarg);
    case 'm': inferencer = USE_MF; break;
    case 'c': useGC = true; break;
    case 'e': gcExceptFirst = true; break;
    case 'a': inferencer = USE_BP; break;
    case 'k': useSparse = true; break;
	case 'b': birchfield = 1; break;
	case 's': squaredDiffs = 1; break;
	case 't': truncDiffs = atoi(optarg); break;
	    //case 'l': lambda1 = atof(optarg); break;
	    //case 'm': lambda2 = atof(optarg); break;
	case 'u': thetaU.push_back((float)atof(optarg)); break;
	case 'v': thetaV.push_back((float)atof(optarg)); break;
	case 'p': thetaP.push_back((float)atof(optarg)); break;
	case 'g': gradThreshVec.push_back(atoi(optarg)); break;
	case 'r': rate = (float)atof(optarg); break;
	case 'f': factorU.push_back((float)atof(optarg)); break;
	case 'i': maxiter = atoi(optarg); break;
	case 'x': closeEnoughPercent = (float)atof(optarg); break;
	case 'o': outscale = atoi(optarg); break;
	case 'q': verbose = 0; break;
    case 'z': gaussSigma = (float)atof(optarg); break;
	case 'j': methods = (int) atoi(optarg); break;
	case 'h': random = (int) atoi(optarg); break;
	case 'l': dirimage = atoi(optarg); break;
	case 'w': testImage = (int) atoi(optarg); break;
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
	FILE* logPlot;

    try {
	char debugname[1000];
    char logname[1000];
	char plotName[1000];


	sprintf(debugname, "%s.txt", outstem);
    sprintf(logname,"%s-log.txt",outstem);
	sprintf(plotName, "%sPlot.txt", outstem);

	debugfile = fopen(debugname, "w");
    ofstream logStream(logname);

	if (debugfile == NULL) 
	    throw CError("Cannot write to %s", debugname);
	if (verbose){
		fprintf(stderr, "writing parameters to %sPlot\n", logname);
	    fprintf(stderr, "writing parameters to %s\n", debugname);
	}

	for (int i=0; i<argc; i++)
	    fprintf(debugfile, "%s ", argv[i]);
	fprintf(debugfile, "\n");

	if (logPlot == NULL) 
	    throw CError("Cannot write to %s", plotName);

	logPlot = fopen(plotName, "w");
	DEBUG_OUT0(verbose, logPlot, "iter \t training rms \t var rms\t test rms \n");
	fclose(logPlot);

	vector<CByteImage> set1Image;
	vector<CByteImage> set2Image;
	vector<CByteImage> setTrueDisp;

	if(testImage == 1 && dirimage == 0){
		throw CError("Please activate the read directory image option and test.png exists");
	}

	if(dirimage == 1){
		readImages(im1name, set1Image, testImage);
		readImages(im2name, set2Image, testImage);
		readImages(truedispname, setTrueDisp, testImage);

		if( (set1Image.size() != set2Image.size()) || (set2Image.size() != setTrueDisp.size()) )
			 throw CError("Not all directory contains the same number of images.\n");
	}
	else{
		CByteImage t1, t2, d1;
		ReadImage(t1, im1name);
		ReadImage(t2, im2name);
		ReadImage(d1, truedispname);

		set1Image.push_back(t1);
		set2Image.push_back(t2);
		setTrueDisp.push_back(d1);
	}

	nP = (int)thetaP.size();
	int nU = (int)thetaU.size(); // number of data parameters (0 or 1 for now)
	int nV = (int)thetaV.size(); // number of smoothness parameters

	// activiate occlusion method.
	if( methods > 1111 || methods < 0){
		throw CError("Incorrect model has been selected\n");
	}
	else{
		if((methods/1000) == 1){
			// 1. Activiate occlusion model
			if(nU > 0)
				setLocal(true);
			else 
				throw CError("Please specfic a localcost (example: -u 12) \n");
		}
		else{
			if(nU != 0)
				throw CError("Please remove -u from command line, or activiate local occlusion with -m\n");
		}

		if(methods >= 1 &&  methods < 1000){
			throw CError("Please activiate localcost before pairwise or gradient features \n");
		}

		if((methods/100)%10 == 1){
			// 2.	Simple pairwise occlusion model
			if(nP == 2)
				setPairwise(true);
			else 
				throw CError("Please specfic a pairwise cost (example: -p 12 -p 12) \n");
		}
		if((methods/10)%10 == 1){
			// 3.	Gradient modulated occlusion model
			if(nP == 4*nV)
				setPairwiseInteraction(true);
			else
				throw CError("Please specfic a 12 pairwise costs (example: -p 12 -p 12 .. -p 12) \n");
		}
		
		if(methods%10 == 1){
			// 4.	Learning scalar occlusion model
		}
	}

	// keep total for all training images.
	fvec toatlEmpirDistU, totalModelDistU;
	toatlEmpirDistU.resize(nU);
	totalModelDistU.resize(nU);

	fvec totalEmpirDistP, totalModelDistP;
	totalEmpirDistP.resize(nP);
	totalModelDistP.resize(nP);

	fvec totalEmpirDistV, totalModelDistV;
	totalEmpirDistV.resize(nV);
	totalModelDistV.resize(nV);

	vector<CByteImage> setIm1grad;
	setIm1grad.resize(set1Image.size());

	vector<CByteImage> setWTAdisp;  
	setWTAdisp.resize(set1Image.size());

	vector<CByteImage> setDisp;  
	setDisp.resize(set1Image.size());

	for(unsigned int i = 0; i < setTrueDisp.size(); ++i){
		ScaleAndOffset(setTrueDisp[i], setTrueDisp[i], (float)1.0/outscale, 0); // scale to integer disps
	}
	
	CShape sh;
	vector<int> width, height;
	for(unsigned int i = 0; i < set1Image.size(); ++i){
		sh = set1Image[i].Shape();
		width.push_back(sh.width);
		height.push_back(sh.height);
		if (sh != set2Image[i].Shape())
			throw CError("image shapes don't match");
	}
	

	setGlobalNpixels(width,height);
	setGlobalWidth(width);



	nF = (int) factorU.size();
	int nT = (int)gradThreshVec.size(); // number of thresholds

	if(nF != nU){
		std::cout << "Scaling vector must have the same size of thetaU\n";
		factorU.resize(nU);
		for(fvec::iterator b = factorU.begin(); b != factorU.end(); ++b){
			*b = 1.0f;
		}	
		nF = (int) factorU.size();
	}

	if (nU + nV < 1)
	    throw CError("Must give at least on thetaV or thetaU");

	if (nU > 2)
	    throw CError("Must give 1 or 0 thetaU values");

	printfVec(thetaP, "thetaP");
	printfVec(factorU, "Scaling Factor");

	// initialize data cost (also computes WTA disparities)
	initializeDataCostVector(set1Image, set2Image, nD, birchfield, squaredDiffs, truncDiffs, setWTAdisp);

	printfVec(thetaU, "thetaU");
	std::cout << "number " << nU  <<  "\n" << std::endl;

	fvec empirDistU, modelDistU;
	empirDistU.resize(nU);
	modelDistU.resize(nU);

	fvec empirDistP, modelDistP;
	empirDistP.resize(nP);
	modelDistP.resize(nP);

	for(vector<CByteImage>::size_type i = 0; i < set1Image.size(); ++i){
		// compute quantized image gradients
		computeQuantizedGradients(set1Image[i], setIm1grad[i], gradThreshVec);
	}

	fvec empirDistV, modelDistV;
	empirDistV.resize(nV);
	modelDistV.resize(nV);

	vector<CByteImage> tmpPreviousDisp, previousDisp, maskeddisp;
	previousDisp.resize(setWTAdisp.size());
	for(unsigned int i = 0; i < setWTAdisp.size(); ++i)
		CopyPixels(setWTAdisp[i], previousDisp[i]);

	printfVec(factorU, "Scaling Factor");
	printfVec(thetaU, "thetaU");
	printfVec(thetaP, "thetaP");

	if (nU > 2)
	    throw CError("Must give 2, 1 or 0 thetaU values");
	
	// ******** smoothness cost

	printfVec(gradThreshVec, "thresholds");
	printfVec(thetaV, "thetaV");

	if (nV > 0 && nV != nT+1)
	    throw CError("Must give exactly one more thetaV values than thresholds");

    DEBUG_OUT1(verbose,debugfile,"Gaussian Sigma = %f\n",gaussSigma);

	if (nU + nV < 1)
	    throw CError("Must give at least one thetaV or thetaU");

	vector<MRF::CostVal*> dcostArr;
	vector<MRF::SmoothCostGeneralFn> scostFn;
	vector<MRFEnergy *> mrf;
	int innerIter, outerIter;

	vector<FloatType> damper;
	vector<FloatType> msgTol;
	vector<MeanField::DiffType> msgDiffFun;
	vector<FloatType> divTol;

	for(unsigned int i = 0; i < set1Image.size(); ++i){
		// Create MRF object
		updateCueIndex(i);
		dcostArr.push_back(computeDataCostArray(thetaU,i));
		scostFn.push_back(computeSmoothnessCostFunction(setIm1grad[i], thetaV, thetaP,i));

		EnergyFunction *energy = new EnergyFunction (new DataCost(dcostArr[i]),
                                             new SmoothnessCost (scostFn[i]) );

	
		switch (inferencer){
			case USE_MF:{
				printf("USE_MF \n");  
				if(nU > 0)
					mrf.push_back(new MeanField(width[i], height[i], nD+1, energy));
				else
					mrf.push_back(new MeanField(width[i], height[i], nD, energy));
          
				damper.push_back(1);
				mrf[i]->setParameters(PARAM_DAMPER, &(damper[i]));
          
				msgTol.push_back(0.001); //0.001 * width * height;
				mrf[i]->setParameters(PARAM_MSGDTOL, &(msgTol[i]));

				msgDiffFun.push_back(MeanField::ENERGY);
				mrf[i]->setParameters(PARAM_MSGDFUN, &(msgDiffFun[i]));

				innerIter = 500;
				outerIter = 1;

				DEBUG_OUT0(verbose, debugfile, "Creating MeanField\n");
				DEBUG_OUT2(verbose, debugfile, "damper = %0.1f  tol = %g\n",
					 damper[i],msgTol[i]);

				break;
			}
			case USE_SMF:
			{
				printf("USE_SMF  \n");
				if(nU > 0 )
					mrf.push_back(new SparseMeanField(width[i], height[i], nD+1, energy));
				 else
					 mrf.push_back(new SparseMeanField(width[i], height[i], nD, energy));;

				damper.push_back(1);
				mrf[i]->setParameters(PARAM_DAMPER, &(damper[i]));
	      
				msgTol.push_back(0.01); //0.001 * width * height;
				mrf[i]->setParameters(PARAM_MSGDTOL, &(msgTol[i]));
	      
				msgDiffFun.push_back(MeanField::ENERGY);
				mrf[i]->setParameters(PARAM_MSGDFUN, &(msgDiffFun[i]));

				divTol.push_back(0.001);
				mrf[i]->setParameters(PARAM_DIVTOL, &(divTol[i]));

				innerIter = 500;
				outerIter = 1;
	      
				DEBUG_OUT0(verbose, debugfile, "Creating SparseMeanField\n");
				DEBUG_OUT2(verbose, debugfile, "damper = %0.1f  tol = %g\n",
					 damper[i],msgTol[i]);
				break;
			}
			case USE_BP:
			{
				if(nU > 0)
					mrf.push_back(new AsyncSumProd(width[i], height[i], nD+1, energy));
				else
					mrf.push_back(new AsyncSumProd(width[i], height[i], nD, energy));
	      
				damper.push_back(0.9);
				mrf[i]->setParameters(PARAM_DAMPER, &(damper[i]));
	      
				msgTol.push_back(1); //0.01 * 2 * (width*(height-1) + (width-1)*height);
				mrf[i]->setParameters(PARAM_MSGDTOL, &(msgTol[i]));
	      
				innerIter = 50;
				outerIter = 10;
	      
				DEBUG_OUT0(verbose, debugfile, "Creating AsyncSumProd\n");
				DEBUG_OUT2(verbose, debugfile, "damper = %0.1f  tol = %g\n",damper,msgTol);
				break;
			}
			case USE_SBP:
			{
				if(nU > 0)
					mrf.push_back(new SparseAsyncSumProd(width[i], height[i], nD+1 , energy));
				else
					mrf.push_back(new SparseAsyncSumProd(width[i], height[i], nD, energy));
	      
				damper.push_back(0.9);
				mrf[i]->setParameters(PARAM_DAMPER, &(damper[i]));
	      
				msgTol.push_back(0.01); //0.0001 * 2 * (width*(height-1) + (width-1)*height);
				mrf[i]->setParameters(PARAM_MSGDTOL, &(msgTol[i]));
	
				divTol.push_back(0.001);
				mrf[i]->setParameters(PARAM_DIVTOL, &(divTol[i]));

				innerIter = 50;
				outerIter = 10;
	      
				DEBUG_OUT0(verbose, debugfile, "Creating SparseAsyncSumProd\n");
				DEBUG_OUT2(verbose, debugfile, "damper = %0.1f  tol = %g\n",
					 damper[i],msgTol[i]);
				break;
			}
			default:
				fprintf(stderr, "Unknown inferencer type");
				return -1;
			}

			mrf[i]->initialize();
			mrf[i]->setLog(&logStream);
    
			delete energy;
		}

    float oldnorm = -1;

    fvec oldTheta;
	float testRms;
	float meanRms;
	float varRms;

	for (int iter=0; iter < maxiter; iter++) {
		for(int i = 0; i < toatlEmpirDistU.size(); ++i){
			toatlEmpirDistU[i] = 0;
			totalModelDistU[i] = 0;
		}

		for(int i = 0; i < totalEmpirDistP.size(); ++i){
			totalEmpirDistP[i] = 0;
			totalModelDistP[i] = 0;
		}

		for(int i = 0; i < totalEmpirDistV.size(); ++i){
			totalEmpirDistV[i] = 0;
			totalModelDistV[i] = 0;
		}
		vector<float> trainRms;

		for(unsigned int i = 0; i < set1Image.size(); ++i){
			if(testImage == 1 && i == 0)
				DEBUG_OUT0(verbose, debugfile, "****** Result of Test Image *****\n\n\n");

			updateCueIndex(i); // quick hack, so we don't have to change code in  graph cut files.
			// Initialize mean field beliefs with graph-cuts
			if (useGC && (inferencer==USE_MF || inferencer==USE_SMF) && 
				!(gcExceptFirst && iter==0))
			{
				EnergyFunction *gcenergy = new EnergyFunction 
					(new DataCost(dcostArr[i]), new SmoothnessCost (scostFn[i]) );

				MRF *gcmrf;

				if(nU > 0)
					gcmrf = new Expansion(width[i],height[i],nD+1,gcenergy);
				else
					gcmrf = new Expansion(width[i],height[i],nD  ,gcenergy);
          
				gcmrf->setParameters(1,&random);

				gcmrf->initialize();

				crfstereo(width[i],height[i],setDisp[i],gcmrf,closeEnoughPercent,1);

				MRF::Label *gclabel = gcmrf->getAnswerPtr();
				MRF::Label *mflabel = mrf[i]->getAnswerPtr();

				 int nPixels = width[i]*height[i]; 
          
				 for (int i=0 ; i<nPixels ; i++, gclabel++, mflabel++)
					*mflabel = *gclabel;
          
				// delete energy; // ?
				delete gcenergy;
				delete gcmrf;

		
				if (inferencer==USE_MF){
					MeanField *mfmrf = dynamic_cast<MeanField*>(mrf[i]);
					mfmrf->initializeAlg(MeanField::LABELPOINT);
				}
				else{
					SparseMeanField *mfmrf = dynamic_cast<SparseMeanField*>(mrf[i]);
					mfmrf->initializeAlg(MeanField::LABELPOINT);
				}
			} else if (0){

				if (inferencer==USE_MF){
					MeanField *mfmrf = dynamic_cast<MeanField*>(mrf[i]);
					mfmrf->initializeAlg(MeanField::LOCAL);	
				}
				else{
					SparseMeanField *mfmrf = dynamic_cast<SparseMeanField*>(mrf[i]);
					mfmrf->initializeAlg(MeanField::LOCAL);
				}
			}

//       char iterstem[1000];
//       sprintf(iterstem,"%s-outeriter-%d-inneriter",outstem,iter);
//       LabelImageWriter imWriter(iterstem,width,height,outscale);
//       mrf->setLabelImageWriter(&imWriter);

      // run the stereo matcher 
      crfstereo(width[i], height[i], setDisp[i], mrf[i], closeEnoughPercent, innerIter, outerIter);
	  computeMaskedImage(setDisp[i],setTrueDisp[i],ignoreVal);
      logStream.flush();

      
      // evaluate matching errors
      float bad1=0, rms=0;
      vector<CByteImage> errormap;
	  errormap.resize(setDisp.size());
	  evaldisps(setDisp[i], setTrueDisp[i], errormap[i], bad1, rms, 1);
	  confusionMatrix(setTrueDisp[i], setDisp[i]);
      DEBUG_OUT2(verbose, debugfile, "bad1= %g   rms= %g\n", bad1, rms);
	    
       if (maxiter > 1) {
			char outname[1000];
			if(testImage == 1 && i == 0)
				sprintf(outname, "%s-test-%d-%03d", outstem, i, iter);
			else
				sprintf(outname, "%s-iter-%d-%03d", outstem, i, iter);
			writeDisparities(setDisp[i], outscale, outname);
		}

	   if(testImage == 1 && i == 0){
			testRms = rms;
			DEBUG_OUT0(verbose, debugfile, "***End of Result of Test Image***** \n\n");
			continue;
		}
		else{
			trainRms.push_back(rms);
		}

       if (nU>0) { // occlusion model
         // compute empirical distribution (of truedisp)
         computeEmpirDistU(setTrueDisp[i], empirDistU);
         // compute model distribution
         computeModelDistU(setDisp[i], modelDistU, mrf[i]);
       }
       else { // no occlusion model
         computeEmpirDistU(setTrueDisp[i], setTrueDisp[i], empirDistU);
         computeModelDistU(setTrueDisp[i], setDisp[i], modelDistU, mrf[i]);
       }

	   // Sum of all empire and model dists for all training images.
		vecSum(toatlEmpirDistU, empirDistU, toatlEmpirDistU);
		vecSum(totalModelDistU, modelDistU, totalModelDistU);

       if (nP>0) { // occlusion model
         if (nP == 4*nV){ // gradient modulated terms
           computeEmpirDistVPI(setTrueDisp[i], setIm1grad[i], empirDistV, empirDistP);
           computeModelDistVPI(setDisp[i], setIm1grad[i], modelDistV, modelDistP, mrf[i]);
         }
         else{ // simple occlusion terms
           computeEmpirDistV(setTrueDisp[i], setIm1grad[i], empirDistV, empirDistP);
           computeModelDistV(setDisp[i], setIm1grad[i], modelDistV, modelDistP, mrf[i]);
         }
       }
       else { // no occlusion model
         computeEmpirDistV(setTrueDisp[i], setTrueDisp[i], setIm1grad[i], empirDistV);
         computeModelDistV(setTrueDisp[i], setDisp[i], setIm1grad[i], modelDistV, mrf[i]);
       }

	   // sum of all empire and model dists
		vecSum(totalEmpirDistV, empirDistV, totalEmpirDistV);
		vecSum(totalModelDistV, modelDistV, totalModelDistV);

		// sum of all empire and model dists
		vecSum(totalEmpirDistP, empirDistP, totalEmpirDistP);
		vecSum(totalModelDistP, modelDistP, totalModelDistP);
	}	// next image
    

      fvec empirDist, modelDist, diffDist, regTheta,relDiff, gradTheta, theta;

      concatVec(toatlEmpirDistU, totalEmpirDistV, totalEmpirDistP, empirDist);
      concatVec(totalModelDistU, totalModelDistV, totalModelDistP, modelDist);
      concatVec(thetaU, thetaV, thetaP, theta);

      vecDiff(modelDist, empirDist, diffDist);

      printfVec(empirDist, "empirDist (GT)");
      printfVec(modelDist, "modelDist     ");
      printfVec(diffDist,  "diffDist      ");

	  DEBUG_OUT0(verbose, debugfile, "Result of the set of images\n");
	  // lam , wrong
	  float bad2=0, rms2=0;
		vector<CByteImage> errormap2;
		errormap2.resize(setDisp.size());
		evaldispsSet(setDisp, setTrueDisp, errormap2, bad2, rms2, 1);
		confusionMatrix(setTrueDisp, setDisp);
		DEBUG_OUT2(verbose, debugfile, "Set training: bad1= %g   rms= %g\n", bad2, rms2);

		// calculate mean
		meanRms = 0;
		varRms = 0;
		for(vector<float>::iterator b = trainRms.begin(); b != trainRms.end(); ++b)
			meanRms += *b;

		// mean rms.
		meanRms /= trainRms.size();

		for(vector<float>::iterator b = trainRms.begin(); b != trainRms.end(); ++b)
			varRms += (*b-meanRms)*(*b-meanRms);

		// calculate variance.
		varRms /= trainRms.size();

		logPlot = fopen(plotName, "a");
		DEBUG_OUT2(verbose, logPlot, "%i\t%g\t%g\t%g\n", iter, meanRms, varRms, testRms);
		fclose(logPlot);

		

		// end calculation.
		// Lam


      
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
      vecEltWiseQuot(gradTheta, empirDist, relDiff);
      
      
      
      float norm = vecNorm(relDiff);
      DEBUG_OUT1(verbose, debugfile, "norm of gradient = %.1f%%   rate = %g\n", 100 * norm, rate);
      if (100*norm < 0.5) // terminate when norm of diff is less than 1 percent
		break;
      


      if (oldnorm<0 || norm<oldnorm) {
		// if current solution is better than previous one,
		// we're fine, keep going with bigger and bigger steps in same direction
		vecCopy(theta, oldTheta); // remember current point

        // Scale local terms
        for (int u=0 ; u<nU ; u++)
          gradTheta[u] *= factorU[u];

		// update theta
		vecScale(gradTheta, rate);

		printfVec(gradTheta, " scaled grad   ");

		// theta = theta + step * gradTheta
		vecSum(theta, gradTheta, theta);

		// make sure thetas don't go to or below 0
        //if (useGC)
         // vecBound(theta, 1e-3f, 1e10f);


        rate *= 1.05f;

		oldnorm = norm;

		DEBUG_OUT1(verbose, debugfile, "rate = %g\n", rate);

      } else { // norm > oldnorm
		// if norm gets bigger, go back half-way and try again with much smaller step size
		vecSum(theta, oldTheta, theta);  //theta = (theta + oldTheta)/2
		vecScale(theta, 0.5);

        if (norm > 2*oldnorm)
          rate *= 0.5;
        else
          rate *= 0.85;

		DEBUG_OUT1(verbose, debugfile, "*** BACKING UP, rate = %g\n", rate);
      }

      printfVec(theta, "*********************** new theta ");

	
	    
      // split back into components
      splitVec(theta, thetaU, thetaV, thetaP);

	   for(unsigned int i = 0; i < setDisp.size(); ++i){

		 // start next iteration with current solution
        CopyPixels(setDisp[i], previousDisp[i]);

		
       // CPAL 
       if (maxiter > 1) {
		   CByteImage tmpDisp;
         char outname[1000];
         CopyPixels(setDisp[i],tmpDisp);
         sprintf(outname, "%s-masked-iter%03d", outstem, iter);
         computeMaskedImage(tmpDisp,setTrueDisp[i]);
         writeDisparities(tmpDisp, outscale, outname);
       }

      
		// Update the energy function based on new thetas
		dcostArr[i] = computeDataCostArray(thetaU, i);
		scostFn[i] = computeSmoothnessCostFunction(setIm1grad[i], thetaV, thetaP, i);

	 
      
		mrf[i]->setData(dcostArr[i]);
		mrf[i]->setSmoothness(scostFn[i]);
	   

		if (inferencer==USE_SMF){
			SparseMeanField *mfmrf = 
				dynamic_cast<SparseMeanField*>(mrf[i]);
			mfmrf->initializeAlg(MeanField::RESET);
		}
		else
			mrf[i]->initializeAlg();
	}
	}// iteration loop

	for(unsigned int i = 0; i < setDisp.size(); ++i){
		char outname[1000];
		sprintf(outname, "%s-iter-%d-final", outstem, i);

		/*if (WRITE_MASKED_DISPS)
			writeDisparities(maskeddisp[i], outscale, outname);
		else*/
			writeDisparities(setDisp[i], outscale, outname);
	}

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
