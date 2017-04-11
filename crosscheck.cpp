// crosscheck.cpp
// $Id: crosscheck.cpp,v 1.1.1.1 2007/12/01 01:01:14 cpal Exp $
//
// read left and right disparity map, convert to integer disparities,
// crosscheck, and save disparity map with occluded pixels marked 0

static char usage[] = "\nUsage: %s disp0.png disp1.png dispscale outdisp0.png\n";

#include <stdio.h>
#include "imageLib.h"
    
int main(int argc, char **argv)
{
    try {
	if (argc != 5)
	    throw CError(usage, argv[0]);

	int verbose = 1;
	CByteImage disp0, disp1, outdisp;
	int argn=1;
	char *disp0name = argv[argn++];
	char *disp1name = argv[argn++];
	int dispscale = atoi(argv[argn++]);
	char *outdispname = argv[argn++];

	ReadImageVerb(disp0, disp0name, verbose);
	ReadImageVerb(disp1, disp1name, verbose);
	CShape sh = disp0.Shape();
	if (sh != disp1.Shape())
	    throw CError("image shapes don't match");
	int width = sh.width, height = sh.height;

	outdisp.ReAllocate(sh);

	if (0) { // convert to integer first
	    // round to nearest integer disparities
	    ScaleAndOffset(disp0, disp0, 1.0/dispscale, 0.5);
	    ScaleAndOffset(disp1, disp1, 1.0/dispscale, 0.5);

	    for (int y = 0; y < height; y++) {
		uchar *di0 = &disp0.Pixel(0, y, 0);
		uchar *di1 = &disp1.Pixel(0, y, 0);
		uchar *out = &outdisp.Pixel(0, y, 0);
		for (int x = 0; x < width; x++) {
		    int d0 = di0[x];
		    int x1 = x - d0;
		    int d1 = x1 >= 0? di1[x1] : 0;
		    //if (d0 == 0 || d1 == 0 || d0 != d1)
		    if (d0 == 0 || d1 == 0 || d0-d1 > 1 || d1-d0 > 1)
			out[x] = 0;
		    else
			out[x] = d0;
		}
	    }
	} else { // convert to integer last
	    float fs = (float)dispscale;
	    for (int y = 0; y < height; y++) {
		uchar *di0 = &disp0.Pixel(0, y, 0);
		uchar *di1 = &disp1.Pixel(0, y, 0);
		uchar *out = &outdisp.Pixel(0, y, 0);
		for (int x = 0; x < width; x++) {
		    float fd0 = di0[x]/fs;
		    int d0 = (int)(fd0 + .5);
		    int x1 = x - d0;
		    float fd1 = x1 >= 0? di1[x1]/fs : 0;
		    int d1 = (int)(fd1 + .5);
		    float fdd = fd0 - fd1;
		    if (fdd < 0)
			fdd = -fdd;
		    if (d0 == 0 || d1 == 0 || fdd > .5)
			out[x] = 0;
		    else
			out[x] = d0;
		}
	    }
	}

	// scale back up
	ScaleAndOffset(outdisp, outdisp, dispscale, 0);
	WriteImageVerb(outdisp, outdispname, verbose);


    }
    catch (CError &err) {
	fprintf(stderr, err.message);
	fprintf(stderr, "\n");
	exit(1);
    }
    return 0;
}
