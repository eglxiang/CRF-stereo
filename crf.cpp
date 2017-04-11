#include "crf.h"

CRF::CRF(Parameters &newp, Data& newd, Log& newl){
	p =&newp;
	d = &newd;
	l = &newl;
}

void CRF::checkInputs(){	
	// Check 
	// d->checkInputs();
	
}

void CRF::computeDistU(){
	if (localOcclusion == true) { // occlusion model
		computeEmpirDistU(d->truedisp, empirDistU,  0, p->ignoreVal);
		computeEmpirDistU(    d->disp, modelDistU, p->nD, p->ignoreVal);
	}
	else { // no occlusion model
		computeEmpirDistU(d->truedisp, d->truedisp, empirDistU,p->ignoreVal);
		computeEmpirDistU(d->truedisp,     d->disp, modelDistU,p->ignoreVal);
	}
}

void CRF::computeDistV(){
	if(localOcclusion == true){ // occlusion model
		if(p->nP == 4*p->nV){ // gradient modulated terms
			computeEmpirDistVPI(d->truedisp, d->im1grad, empirDistV, empirDistP,
                                0, p->ignoreVal);
			computeEmpirDistVPI(    d->disp, d->im1grad, modelDistV, modelDistP,
                                p->nD, p->ignoreVal);
		}	
		else{ // simple occlusion terms
			computeEmpirDistV(d->truedisp, d->im1grad, empirDistV, empirDistP, 
				               0, p->ignoreVal);
			computeEmpirDistV(    d->disp, d->im1grad, modelDistV, modelDistP, 
                              p->nD, p->ignoreVal);
		}
	}else{ // no occlusion model
		computeEmpirDistV(d->truedisp, d->truedisp, d->im1grad, empirDistV, p->ignoreVal);
		computeEmpirDistV(d->truedisp,     d->disp, d->im1grad, modelDistV, p->ignoreVal);
	}
}

bool CRF::update(const fvec& relDiff, float &oldnorm, CByteImage& previousDisp, fvec& gradTheta){
	float norm = vecNorm(relDiff);
	fvec oldTheta;
//	    DEBUG_OUT1(verbose, debugfile, "norm of gradient = %.1f%%\n", 100 * norm);

	if (100*norm < 0.5) // terminate when norm of diff is less than 1 percent
          return true;

	if (oldnorm < 0 || norm < oldnorm) {
		// if current solution is better than previous one,
		// we're fine, keep going with bigger and bigger steps in same direction
		vecCopy(p->theta, oldTheta); // remember current point

		// Scale local terms
		for (int u=0 ; u< p->nU ; u++)
			gradTheta[u] *= p->factorU[u];

		// update theta
		vecScale(gradTheta, p->rate);

		printfVec(gradTheta, " scaled grad   ");

		// theta = theta + step * gradTheta
		vecSum(p->theta, gradTheta, p->theta);		

		// make sure thetas don't go to or below 0
		vecBound(p->theta, 1e-3f, 1e10f);

		p->rate *= 1.05f;
			
		oldnorm = norm;

//			DEBUG_OUT1(verbose, debugfile, "rate = %g\n", rate);

		CopyPixels(d->disp, previousDisp);
   
	} else { // norm > 2*oldnorm
		// if norm gets much bigger, go back half-way and try again with much smaller step size
		vecSum(p->theta, oldTheta, p->theta);  //theta = (theta + oldTheta)/2
		vecScale(p->theta, 0.5);
	
		if (norm>2*oldnorm) 
			p->rate *= 0.5f;
		else
			p->rate *= 0.85f;

//			DEBUG_OUT1(verbose, debugfile, "*** BACKING UP, rate = %g\n", rate);	
	}
}

void CRF::run(){
	float oldnorm = -1;
	CByteImage tmpPreviousDisp, previousDisp, maskeddisp;
	CopyPixels(d->WTAdisp, previousDisp);

	concatVec(p->thetaU, p->thetaV, p->thetaP, p->theta);
	printfVec(p->theta, "*********************** initial theta ");

	for (int iter=0; iter < p->maxiter; iter++) {
		 // split theta into U, V,  and P components
	    splitVec(p->theta, p->thetaU, p->thetaV, p->thetaP);
	    // run the stereo matcher, starting from the previous solution
	    DataCost *dcost = computeDataCost(p->thetaU);
	    SmoothnessCost *scost = computeSmoothnessCost(d->im1grad, p->thetaV, p->thetaP);

		if (localOcclusion == true) // occlusion model ==> extra state			
			crfstereo(d->width, d->height, p->nD + 1, dcost, scost, previousDisp, d->disp, p->closeEnoughPercent);
		else	
			crfstereo(d->width, d->height, p->nD, dcost, scost, previousDisp, d->disp, p->closeEnoughPercent);

	    delete scost;
	    delete dcost;

        if (p->ignoreVal>= 0)
			computeMaskedImage(d->disp,d->truedisp,p->ignoreVal);

	    // evaluate matching errors
	    float bad1=0, rms=0;
	    CByteImage errormap;
	    //evaldisps(disp, truedisp, errormap, maskeddisp, bad1, rms, 1);
        evaldisps(d->disp, d->truedisp, errormap, bad1, rms, 1);
		confusionMatrix(d->truedisp, d->disp);
//	    DEBUG_OUT2(verbose, debugfile, "bad1= %g   rms= %g\n", bad1, rms);
	    
	    if (p->maxiter > 1) {
			char outname[1000];
			sprintf(outname, "%s-iter%03d", p->outstem, iter);
			writeDisparities(d->disp, p->outscale, outname);
	    }

		computeDistU();
		computeDistV();

		fvec empirDist, modelDist, diffDist, relDiff, regTheta,gradTheta;
	    concatVec(empirDistU, empirDistV, empirDistP, empirDist);
	    concatVec(modelDistU, modelDistV, modelDistP, modelDist);

	    vecDiff(modelDist, empirDist, diffDist);
	    printfVec(empirDist, "empirDist (GT)");
	    printfVec(modelDist, "modelDist     ");
	    printfVec(diffDist,  "diffDist      ");
        
        if (p->gaussSigma!=0){
            vecCopy(p->theta,regTheta);
            vecScale(regTheta,-1/(p->gaussSigma*p->gaussSigma));
            vecSum(diffDist,regTheta,gradTheta);
            
            printfVec(regTheta,  "gaussReg      ");
            printfVec(gradTheta,  "gradient      "); 
		}
        else
			vecCopy(diffDist,gradTheta);


		// terminate if gradient is small enough percentage of empirdist
	    vecEltWiseQuot(gradTheta, empirDist, relDiff);

		if(update(relDiff, oldnorm, previousDisp, gradTheta))
			break;

		printfVec(p->theta, "*********************** new theta ");
	    // start next iteration with current solution
	}

	//if (WRITE_MASKED_DISPS)
	//	writeDisparities(maskeddisp, p->outscale, p->outstem.c_str());
	//else
	//	writeDisparities(d->disp, p->outscale, p->outstem.c_str());

//	fclose(debugfile);
}


void CRF::initialize(){
	
	localOcclusion = p->nD > 0; // Active Local Occlusion, if theta vector is greater that 0


	//// initialize data cost (also computes WTA disparities)
	initializeDataCost(d->im1, d->im2, p->nD, p->birchfield, p->squaredDiffs, p->truncDiffs, d->WTAdisp);

	log << "number " << p->nD  <<  "\n" << std::endl;

	if (p->nU > 2)
	    throw CError("Must give 1 or 0 thetaU values");

	// compute quantized image gradients
	computeQuantizedGradients(d->im1, d->im1grad, p->gradThreshVec);

	empirDistU.resize(p->nU);
	modelDistU.resize(p->nU);

	empirDistP.resize(p->nP);
	modelDistP.resize(p->nP);

	empirDistV.resize(p->nV);
	modelDistV.resize(p->nV);

	if(p->nP == 2)
		pairwise = true; 
	else if(p->nP == 4*p->nV)
		pairwiseInteraction = true;
}