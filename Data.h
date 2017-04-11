#ifndef DATA_FILE
#define DATA_FILE

#include "Parameters.h"
#include "crfstereo.h"
//#include "expectations.h"
//#include "evaldisps.h"
#include  <iostream>
#include "Logger.h"
using std::endl;

class Data{
public:
	Data();
	void initialize(Parameters&, Log &);
	void checkInputs();

private:
	CByteImage im1, im2;      // input images (gray or color)
	CByteImage truedisp;       // quantized gradients of image 1
	CByteImage im1grad;       // quantized gradients of image 1
	CByteImage WTAdisp;       // WTA disparities (from matching cost only)
	CByteImage disp;          // output disparities

	CShape sh;
	int width, height;

	friend class CRF;
};
#endif