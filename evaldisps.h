// $Id: evaldisps.h,v 1.2 2007/12/07 17:23:22 weinman Exp $
// computes and prints disparity error statistics, creates errormap
// returns percentages of bad pixels (error > 1) and RMS disp error (both in nonoccluded areas)
// verbose==1: print lots, verbose==0: print single line, verbose<0: print nothing
void evaldisps(CByteImage disp, CByteImage truedisp, CByteImage &errormap,
	       float &bad1, float &rms, int verbose);

void confusionMatrix(CByteImage truedisp, CByteImage disp);

void evaldispsSet(vector<CByteImage> disp, vector<CByteImage> truedisp, vector<CByteImage> &errormap,
               float &bad1, float &rms, int verbose);

void confusionMatrix(vector<CByteImage> &truedisp, vector<CByteImage> &disp);