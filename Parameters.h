#ifndef PARAMETERS_FILE
#define PARAMETERS_FILE

#include <stdio.h>
#include <vector>
#include <fstream>
#include <assert.h>
#include <iostream>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <string>

#ifdef _MSC_VER        // Visual C++
# include "XGetopt.h"  // replacement for getopt()
#else
# include <getopt.h>
#endif

using std::vector;
using std::ifstream;
using std::cout;
using std::ostream;
using std::endl;
using std::string;

typedef std::vector<float> fvec;


static char usage[] = "\n\
 usage: %s [options] imL imR truedisp outstem \n\
\n\
  reads imL, imR, and true disparities (in png or pgm/ppm format)\n\
  runs MRF stereo\n\
  writes disparities corresponding imL to outstem.png and parameters to outstem.txt\n\
\n\
  options:\n\
    -n nD           disparity levels, by default 16 (i.e. disparites 0..15)\n\
    -d              value of true disparity to mask/ignore\n\
    -b              use Birchfield/Tomasi costs\n\
    -s              use squared differences (absolute differences by default)\n\
    -t trunc        truncate differences to <= 'trunc'\n\
    -u u1           initial data cost param (thetaU) - can specify one or none\n\
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


class Parameters{
private:
	int getNumParameters(char *);
	void convertTo2DArray(char **, char*);
public:
	Parameters();
	Parameters(Parameters& data);

	void parseCommandLine(int, char **);
	void readParameterFile(char *);
	vector<string> getImageFileName();

	friend ostream& operator<<(ostream&, const Parameters&);

	int getDisparityLevel() const;
	int getNumberPairwise() const;
	int getScaleFactor() const;
	int getOutScale() const;

	unsigned int getThetaUCount() const;
	const fvec& getThetaU() const;

	
	string getImage1Name() const;
	string getImage2Name() const;
	string getTrueDispName() const;
	string getOutStem() const;


private:
	// parameters controlled via command-line options:
    int nD;					// disparity levels (d = 0 .. nD-1)
	int nP;					// Number of pairwise.
	int nF;					// Number of scaling factor.
	int nU;					// Number of local pair.
	int nV;					// Number of smoothness parameters
	int nT;					// Number of gradient thresholds.

    int ignoreVal;			// Value of ground truth disp image to ignore
    bool birchfield;		// use Birchfield/Tomasi costs
    bool squaredDiffs;		// use squared differences (absolute differences by default)
    int truncDiffs;			// truncated differences (before squaring), by default not
    float rate;
    float gaussSigma;		// Scale of gaussian  regularizer (zero for none)
    fvec factorU;           // Scaling factor for learning rate
    int maxiter;			// maximum number of iterations for gradient descent
    float closeEnoughPercent; // when to stop GC or BP iterations
    int outscale;			// scale factor for disparities; -1 means full range 255.0/(nD-1)
	int verbose;

    std::vector<int> gradThreshVec;
    fvec thetaU;			// data weights (at most one for now)
    fvec thetaV;			// smoothness weights associated with gradient bins
    fvec theta;				// concatenation of the two
	fvec thetaP;			// Pairwise Occluded Cost

	string im1name;
    string im2name;
    string truedispname;
    string outstem;

	friend class CRF;
	friend class Log;
};

#endif