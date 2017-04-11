// evaldisps.cpp
// $Id: evaldisps.cpp,v 1.5 2007/12/07 17:23:22 weinman Exp $

#include <stdio.h>
#include <math.h>
#include <iostream>
#include "imageLib.h"
#include "crfstereo.h"


#define ERRLIMIT 4  /* highest absolute disparity error to count */
#define UNK 0     /* color for unknown true disparity and occlusion */
int errcol[ERRLIMIT+1] = {255, 200, 40, 20, 0};  // colors for err==0, 1, 2, 3, and >3

using std::cout;

// computes and prints disparity error statistics, creates errormap
// returns percentages of bad pixels (error > 1) and RMS disp error (both in nonoccluded areas)
// verbose==1: print lots, verbose==0: print single line, verbose<0: print nothing
void evaldisps(CByteImage disp, CByteImage truedisp, CByteImage &errormap,
               float &bad1, float &rms, int verbose)
{
  CShape sh = disp.Shape();
 // if (sh != truedisp.Shape())

	//throw CError("image shapes don't match");
  int width = sh.width, height = sh.height;

  errormap.ReAllocate(sh);

  int cnt[ERRLIMIT+1];
  int tot = 0;
  int k;
  int nStates =  getNumStates(); //I CAN'T GET THIS TO LINK! -JW
  for (k = 0; k <= ERRLIMIT; k++)
	cnt[k] = 0;
  float sserr = 0;

  for (int y = 0; y < height; y++) {
	uchar *di = &disp.Pixel(0, y, 0);
	uchar *tdi = &truedisp.Pixel(0, y, 0);
	uchar *errmap = &errormap.Pixel(0, y, 0);
	for (int x = 0; x < width; x++) {
      int d = di[x];
      int td = tdi[x];
      if (td == 0 || d == (nStates -1) ) { // unknown true disp
		errmap[x] = UNK;
      } else {
		int err = td - d;
		float ferr = (float)err;
		if (err < 0) err = -err;
		if (err > ERRLIMIT) err = ERRLIMIT;
		tot++;
		cnt[err]++;
		sserr += ferr * ferr;
		errmap[x] = errcol[err];
      }
	}
  }

  if (verbose > 0) {
	printf("err thresh:\t");
	for (k = 0; k < ERRLIMIT; k++)
      printf("%d\t", k);
	printf("\nbad disp %%:\t");
  }
  int ct = 0;
  for (k = 0; k < ERRLIMIT; k++) {
	ct += cnt[k];
	float errpercent = (float)(100.0 * (1.0 - (float)ct/(float)tot));
	if (k == 1)
      bad1 = errpercent;
	if (verbose > 0)
      printf("%.1f\t", errpercent);
	else if (verbose >= 0)
      printf("%d:%5.1f   ", k, errpercent);
  }

  rms = sqrt(sserr / tot);
  if (verbose > 0)
	printf("\nRMS disp error = %.2f\n", rms);
  else if (verbose >= 0)
	printf("rms: %.2f\n", rms);
}

void confusionMatrix(CByteImage truedisp, CByteImage disp)
{
	int d, td;
	uchar *di, *tdi;
	int model_emp = 0;
	int model_notEmp = 0;
	int notMod_emp = 0;
	int notMod_notEmp = 0;
	int nState = getNumStates(); //I CAN'T GET THIS TO LINK!
	if(getLocal() == false) nState++;

	CShape sh = disp.Shape();
	//if (sh != truedisp.Shape())
	//	throw CError("image shapes don't match");
	int width = sh.width, height = sh.height;

	 for (int y = 0; y < height; y++) {
			di = &disp.Pixel(0, y, 0);
			tdi = &truedisp.Pixel(0, y, 0);
	
			for (int x = 0; x < width; x++) {
				d = di[x];
				td = tdi[x];
				if(td == 0 & d == nState - 1)
					model_emp++;
				else if(td != 0 & d == nState - 1)
					model_notEmp++;
				else if(td == 0 & d != nState - 1)
					notMod_emp++;
				else
					notMod_notEmp++;
			}
	 }
	 DEBUG_OUT0(verbose, debugfile, "\n Confusion Matrix \n");
	 DEBUG_OUT0(verbose, debugfile, "Model \t O \t D \t\t O \t D \n");
	 DEBUG_OUT0(verbose, debugfile, "Truth \n");
	 DEBUG_OUT2(verbose, debugfile, "O \t %d \t %d \t\t", model_emp , notMod_emp );
	 DEBUG_OUT2(verbose, debugfile, "%f \t %f \n", model_emp/(float)(notMod_emp+model_emp) , notMod_emp/(float)(notMod_emp+model_emp) );
	 DEBUG_OUT2(verbose, debugfile, "D \t %d \t %d \t\t", model_notEmp , notMod_notEmp );
	 DEBUG_OUT2(verbose, debugfile, "%f \t %f \n", model_notEmp/(float)(model_notEmp+notMod_notEmp) , notMod_notEmp/(float)(model_notEmp+notMod_notEmp) );
}


// computes and prints disparity error statistics, creates errormap
// returns percentages of bad pixels (error > 1) and RMS disp error (both in nonoccluded areas)
// verbose==1: print lots, verbose==0: print single line, verbose<0: print nothing
void evaldispsSet(vector<CByteImage> disp, vector<CByteImage> truedisp, vector<CByteImage> &errormap,
               float &bad1, float &rms, int verbose)
{
	int cnt[ERRLIMIT+1];
	int tot = 0;
	int k;
	float sserr = 0;
	int nStates = getNumStates(); // I CAN'T GET THIS TO LINK! -JW
	for (k = 0; k <= ERRLIMIT; k++)
		cnt[k] = 0;


	for(unsigned int i = 0; i < disp.size(); ++i){ 
		CShape sh = disp[i].Shape();
		if (sh != truedisp[i].Shape())
			throw CError("image shapes don't match");
		int width = sh.width, height = sh.height;

		errormap[i].ReAllocate(sh);
		
		for (int y = 0; y < height; y++) {
			uchar *di = &(disp[i]).Pixel(0, y, 0);
			uchar *tdi = &(truedisp[i]).Pixel(0, y, 0);
			uchar *errmap = &(errormap[i]).Pixel(0, y, 0);
			for (int x = 0; x < width; x++) {
				int d = di[x];
				int td = tdi[x];				
				if (td == 0 || d == (nStates -1) ) { // unknown true disp
					errmap[x] = UNK;
				} else {
					int err = td - d;
					float ferr = (float)err;
					if (err < 0) err = -err;
					if (err > ERRLIMIT) err = ERRLIMIT;
					tot++;
					cnt[err]++;
					sserr += ferr * ferr;
					errmap[x] = errcol[err];
				}
			}
		}
	}

	if (verbose > 0) {
		printf("err thresh:\t");
		for (k = 0; k < ERRLIMIT; k++)
			printf("%d\t", k);
		printf("\nbad disp %%:\t");
	}
 
	int ct = 0; 
	for (k = 0; k < ERRLIMIT; k++) {
		ct += cnt[k];
		float errpercent = (float)(100.0 * (1.0 - (float)ct/(float)tot));
		if (k == 1)
			bad1 = errpercent;
		if (verbose > 0)
			printf("%.1f\t", errpercent);
		else if (verbose >= 0)
			printf("%d:%5.1f   ", k, errpercent);
	}

	rms = sqrt(sserr / tot);	
	if (verbose > 0)
		printf("\nRMS disp error = %.2f\n", rms);
	else if (verbose >= 0)
		printf("rms: %.2f\n", rms);
}

void confusionMatrix(vector<CByteImage> &truedisp, vector<CByteImage> &disp)
{
	int d, td;
	uchar *di, *tdi;
	int model_emp = 0;
	int model_notEmp = 0;
	int notMod_emp = 0;
	int notMod_notEmp = 0;
	int nState = getNumStates();

	if(getLocal() == false) nState++;


	for(unsigned int i = 0; i < disp.size(); ++i){
		CShape sh = disp[i].Shape();
		if (sh != truedisp[i].Shape())
			throw CError("image shapes don't match");
		int width = sh.width, height = sh.height;

		for (int y = 0; y < height; y++) {	
			di = &(disp[i]).Pixel(0, y, 0);
			tdi = &(truedisp[i]).Pixel(0, y, 0);
				
			for (int x = 0; x < width; x++) {
				d = di[x];
				td = tdi[x];
				if(td == 0 & d == nState - 1)
					model_emp++; // both occluded
				else if(td != 0 & d == nState - 1)
					model_notEmp++; // true is not, but model is 
				else if(td == 0 & d != nState - 1)
					notMod_emp++; // model is not occluded but true is occluded
				else
					notMod_notEmp++; // both are not occluded.
			}
		 }
	}
	 DEBUG_OUT0(verbose, debugfile, "\n Confusion Matrix \n");
	 DEBUG_OUT0(verbose, debugfile, "Model \t O \t D \t\t O \t D \n");
	 DEBUG_OUT0(verbose, debugfile, "Truth \n");
	 DEBUG_OUT2(verbose, debugfile, "O \t %d \t %d \t\t", model_emp , notMod_emp );
	 DEBUG_OUT2(verbose, debugfile, "%f \t %f \n", model_emp/(float)(notMod_emp+model_emp) , notMod_emp/(float)(notMod_emp+model_emp) );
	 DEBUG_OUT2(verbose, debugfile, "D \t %d \t %d \t\t", model_notEmp , notMod_notEmp );
	 DEBUG_OUT2(verbose, debugfile, "%f \t %f \n", model_notEmp/(float)(model_notEmp+notMod_notEmp) , notMod_notEmp/(float)(model_notEmp+notMod_notEmp) );
}
