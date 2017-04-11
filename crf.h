#include "Parameters.h"
#include "Data.h"
#include "expectations.h"
#include "evaldisps.h"

#include <iostream>

using std::cout;

class CRF{

private:
	void computeDistU();
	void computeDistV();
	bool update(const fvec&,  float&, CByteImage&, fvec&);

public:
	CRF(Parameters &, Data&, Log& l);
	void initialize();
	void checkInputs();
	void run();

private:
	Parameters *p;
	Data *d;
	Log *l;
	bool localOcclusion;		// Used to activate local occlusion model.
	bool pairwise;				// Used to activate local pairwise  model.
	bool pairwiseInteraction;	// Used to activate local pairwise interaction model.

	fvec empirDistV;	
	fvec empirDistP;
	fvec empirDistU;

	fvec modelDistP;
	fvec modelDistV;
	fvec modelDistU;
};
