#include "Parameters.h"

// return the number of parameters in string.
int Parameters::getNumParameters(char *string){
	int count = 0;
	char *index = string;			// use to index string
	bool space = isspace(*index);	// the state of the being of string

	while( *index != '\0'){
		
		if(space != isspace(*index)){
			// The state has changed, this implies that either 
			// case 1: we were at a command and we encounter a whitespace (-n \t)
			// case 2: we were at a whitespace and encounter a command (\t -n)
			if(space == false)
				count++;	// case 2, we found another command.
			space = isspace(*index);
		}
		index++;
	}
	return count;
}

// convert an array of parameters into a multiple dimensional array
// that contains a parameter in each dimension.  This is required
// in order to use getOpt.
void Parameters::convertTo2DArray(char **des, char* source){
	char *sindex = source;	// use to index source
	char **dindex = des;	// use to index des
	int start = 0, end = 0;	// use to track the beginning and ending of the next command.
	bool space = isspace(*sindex);	// the state of the being of string

	while(*sindex != '\0'){
		if(space != isspace(*sindex)){
			// The state has changed, this implies that either 
			// case 1: we were at a command and we encounter a whitespace (-n \t)
			// case 2: we were at a whitespace and encounter a command (\t -n)

			if(space == false){
				// case 2, we found another command and create a new command
				*dindex = new char[end - start + 1];
				strncpy(*dindex, (source+start), end - start);
				(*dindex)[end - start] = '\0';

				// index des to the next position 
				dindex++;
			}			
			start = end;
			space = isspace(*sindex);
		}
		sindex++;
		end++;
	}	
	*dindex = new char[end - start + 1];
	strncpy(*dindex, (source+start), end - start);
	(*dindex)[end - start] = '\0';
}

// Default Constructor
Parameters::Parameters(){
	nD = 20;               // disparity levels (d = 0 .. nD-1)
	nP = 0;					// Number of pairwise.
	nF = 0;					// Number of scaling factor.
	nU = 0;
    ignoreVal = -1;        // Value of ground truth disp image to ignore
    birchfield = false;   // use Birchfield/Tomasi costs
    squaredDiffs = false; // use squared differences (absolute differences by default)
    truncDiffs = 255;      // truncated differences (before squaring), by default not
    rate = 0.001f;
    gaussSigma = 0;      // Scale of gaussian  regularizer (zero for none)
    factorU;              // Scaling factor for learning rate
    maxiter = 100;          // maximum number of iterations for gradient descent
    closeEnoughPercent = 2.0f; // when to stop GC or BP iterations
    outscale = -1;         // scale factor for disparities; -1 means full range 255.0/(nD-1)
    verbose = 1;               // print messages to stderr
}

// Parse and store parameters.
void Parameters::parseCommandLine(int argc, char **argv){
	// parse command-line arguments using the "getopt" utility
	int o;
	while ((o = getopt(argc, argv, "ed:n:bst:u:p:v:g:r:f:i:x:o:z:q")) != -1){
		switch (o) {
			case 'n': nD = atoi(optarg); break;
			case 'd': ignoreVal = atoi(optarg);
			case 'b': birchfield = 1; break;
			case 's': squaredDiffs = 1; break;
			case 't': truncDiffs = atoi(optarg); break;
			case 'u': thetaU.push_back((float)atof(optarg)); break;
			case 'v': thetaV.push_back((float)atof(optarg)); break;
			case 'p': thetaP.push_back((float)atof(optarg)); break;
			case 'g': gradThreshVec.push_back(atoi(optarg)); break;
			case 'r': rate = (float)atof(optarg); break;
			case 'f': factorU.push_back((float)atof(optarg)); break;
			case 'i': maxiter = atoi(optarg); break;
			case 'x': closeEnoughPercent = (float)atof(optarg); break;
			case 'o': outscale = atoi(optarg); break;
			case 'q': verbose = 0; break;
			case 'z': gaussSigma = (float)atof(optarg); break;
			default: 
				fprintf(stderr, "Ignoring unrecognized option %s\n", argv[optind-1]);
		}
	}
	if (optind != argc-4) {
		fprintf(stderr, usage, argv[0]);
		exit(1);
    }

	im1name = argv[optind++];
    im2name = argv[optind++];
    truedispname = argv[optind++];
    outstem = argv[optind++];

	if (outscale < 0)
		outscale = 255 / (nD-1);

	nU = (int)thetaU.size(); // number of data parameters (0 or 1 for now)
	nF = (int) factorU.size();
	nP = (int)thetaP.size();
	nT = (int)gradThreshVec.size(); // number of thresholds
	nV = (int)thetaV.size(); // number of smoothness parameters

	if(nF != nU){
		std::cout << "Scaling vector must have the same size of thetaU\n";
		factorU.resize(nU);
		for(fvec::iterator b = factorU.begin(); b != factorU.end(); ++b){
			*b = 1.0f;
		}	
	}
}

// Overload output operator.
ostream& operator<<(ostream& o, const Parameters& d){
	o << "disparity levels: " << d.nD << endl;
	o << "number of pairwise: " << d.nP << endl;
	o << "number of scaling factor: " << d.nF << endl;
	o << "value of ground truth disp image to ignore: " << d.ignoreVal << endl;
	o << "use Birchfield/Tomasi costs: " << d.birchfield << endl;
	o << "use squared differences: " << d.squaredDiffs << endl;
	o << "truncated differences: " << d.truncDiffs << endl;
	o << "rate: " << d.rate << endl;
	o << "Scale of gaussian  regularizer: " << d.gaussSigma << endl;
	o << "size of scaling factor for learning rate: " << d.factorU.size() << endl;
	o << "maximum number of iterations for gradient descent: " << d.maxiter << endl;
	o << "when to stop GC or BP iterations: " << d.closeEnoughPercent << endl;
	o << "scale factor for disparities; -1 means full range 255.0/(nD-1): " << d.outscale << endl;
	o << "print messages to stderr: " << d.verbose << endl;
	
	o << "image 1 name: " << d.im1name << endl;
	o << "image 2 name: " << d.im2name << endl;
	o << "true dispity image name: " << d.truedispname << endl;
	o << "print messages to stderr: " << d.outstem << endl;
	return o;
}

// Read parameter file.
void Parameters::readParameterFile(char *fileName){
	assert(fileName != '\0'); // null pointer

	FILE *fp;
	long len;
	char *buf;

	// Read the whole parameter file at once.
	fp=fopen(fileName,"rb");
	assert(fp != '\0'); // quit if cannot open.
	fseek(fp,0,SEEK_END); //go to end.
	len=ftell(fp); //get position at end (length)
	fseek(fp,0,SEEK_SET); //go to beg.
	buf = new char[len + 1]; //malloc buffer
	buf[len] = '\0';
	fread(buf,len,1,fp); //read into buffer
	fclose(fp);

	// Get the number of parameters.
	int pLen = getNumParameters(buf);

	// getopt skip the first element, we need to create 2 more than
	// the number of element.  1 for getopt and 1 to end the parameter.
	char** parameters = new char*[pLen+1];
	
	// convert single dimension to 2 dimensional array.
	parameters[0] = '\0';
	convertTo2DArray(parameters+1, buf);
	parseCommandLine(pLen+1,parameters);

	for(int i = 0; i < pLen + 1; ++i){
		if(parameters[i] != '\0'){
			delete[] parameters[i];
		}
	}

	//delete[] parameters;
	delete[] (buf);
}



int Parameters::getDisparityLevel() const{
	return nD;
}
	
int Parameters::getNumberPairwise() const{
	return nP;
}
	
int Parameters::getScaleFactor() const{
	return nF;
}

int Parameters::getOutScale() const{
	return outscale;
}
	
string Parameters::getImage1Name() const{
	return im1name;
}
	
string Parameters::getImage2Name() const{
	return im2name;
}
	
string Parameters::getTrueDispName() const{
	return truedispname;
}
	
string Parameters::getOutStem() const{
	return outstem;
}


unsigned int Parameters::getThetaUCount() const{
	return nU;
}

const fvec& Parameters::getThetaU() const{
	return thetaU;
}