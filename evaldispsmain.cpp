// evaldispsmain.cpp
// $Id: evaldispsmain.cpp,v 1.1.1.1 2007/12/01 01:01:14 cpal Exp $
//
// main program for running evaldisps from the commandline

static char usage[] = "\nUsage: %s disp.png truedisp.png dispscale [errormap.png]\n";

#include <stdio.h>
#include "imageLib.h"
#include "evaldisps.h"
    
int main(int argc, char **argv)
{
    try {
	if (argc < 4 || argc > 5)
	    throw CError(usage, argv[0]);

	int verbose = 0;
	CByteImage disp, truedisp, errormap;
	int argn=1;
	char *dispname = argv[argn++];
	char *truedispname = argv[argn++];
	int dispscale = atoi(argv[argn++]);
	char *errmapname = argn < argc ? argv[argn++] : NULL;

	ReadImageVerb(disp, dispname, verbose);
	ReadImageVerb(truedisp, truedispname, verbose);

	// convert to integer disparities
	ScaleAndOffset(disp, disp, 1.0/dispscale, 0);
	ScaleAndOffset(truedisp, truedisp, 1.0/dispscale, 0);

	float bad1=0, rms=0;
	evaldisps(disp, truedisp, errormap, bad1, rms, verbose);
	//printf("bad=%g  rms=%g\n", bad1, rms);

	if (errmapname)
	    WriteImageVerb(errormap, errmapname, verbose);
    }
    catch (CError &err) {
	fprintf(stderr, err.message);
	fprintf(stderr, "\n");
	exit(1);
    }
    return 0;
}
