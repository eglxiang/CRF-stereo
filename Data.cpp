#include "Data.h"
Data::Data(){
	
}

void Data::initialize(Parameters& p, Log& log){
	log << "Reading image " << p.getImage1Name() << "\n";
	ReadImage(im1, p.getImage1Name().c_str());
	log << "Reading image " << p.getImage2Name() << "\n";
	ReadImage(im2, p.getImage2Name().c_str());
	log << "Reading image " << p.getTrueDispName() << "\n";
	ReadImage(truedisp, p.getTrueDispName().c_str());

	ScaleAndOffset(truedisp, truedisp, (float)1.0/p.getOutScale(), 0); // scale to integer disps

	sh = im1.Shape();
	if (sh != im2.Shape())
	    throw CError("image shapes don't match");
	width = sh.width;
	height = sh.height;
}