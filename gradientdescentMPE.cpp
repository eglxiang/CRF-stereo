// gradientdescent.cpp
// $Id: gradientdescentMPE.cpp,v 1.14 2007/12/11 13:28:48 weinman Exp $
//
// finds best parameters theta using gradient descent
//
// Daniel Scharstein

// if you want to have a thetaV parameter for equal neighboring disparities, 
// define HAVE_COST_FOR_DELTA_D_EQUAL_ZERO in crfstereo.h

// mask occluded areas in saved disparity maps
#define WRITE_MASKED_DISPS 0
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <fstream>
#include "imageLib.h"
//#include "rightTrans.cpp"
using namespace std;
using std::fstream;


int main(int argc,char *argv[]) {
	// check input args
	if(argc != 5){
		cout << argv[0] << " " << argc<< " left_image right_image left_disp right_disp \n";
		exit(1);
	}
	CByteImage right_img, left_img, right_disp, left_disp;
	ReadImage(right_img, argv[1]);
	ReadImage(left_img, argv[2]);
	ReadImage(right_disp, argv[3]);
	ReadImage(left_disp, argv[4]);
}






