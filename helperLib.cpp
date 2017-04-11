#include "helperLib.h"
#include <vector>
#include <algorithm>

using std::vector;
using std::sort;

// Read image from a directory.
int readImages(char *dirName, vector<CByteImage>& ims, bool testImage){
  /* check command line arguments */
  if (dirName == '\0') {
    fprintf (stderr, "invalid directory");
    return EXIT_FAILURE;
  }

  DIR *dir;
  struct dirent *ent;
  string name, fullPath;
  vector<string> files;

  /* open directory stream */
  dir = opendir (dirName);
  if (dir != NULL) {
	  /* print all the files and directories within directory */
      while ((ent = readdir (dir)) != NULL) {
		  name = ent->d_name ;
		  if(name == "..")
			  continue;
		  else if(name == ".")
			  continue;
		  else if( std::find(name.begin(), name.end(), '.') != name.end()){
			  fullPath = dirName;
			  if( fullPath[fullPath.length() -1] != dirs)
				  fullPath += dirs;
			  fullPath += name;
			  files.push_back(fullPath);
		  } else{ // a directory or unknow format. 
		  }  
      }
	sort(files.begin(), files.end());

	string testString = dirName;
	testString += dirs;
	testString += "test.png";
	if(testImage == 1){
		vector<string>::iterator it = std::find(files.begin(), files.end(), testString);
		if(it == files.end() ){
			throw CError("Directory does not contain test.png image \n");
		}
		CByteImage imTemp;
		ReadImage(imTemp, it->c_str());
		ims.push_back(imTemp);
		files.erase(it);
	}

	for(vector<string>::iterator b = files.begin(); b < files.end(); ++b){
		cout << *b << "\n";
		CByteImage imTemp;
		ReadImage(imTemp, b->c_str());
		ims.push_back(imTemp);
	}
	
 
	  closedir (dir);  
  } else {
      /* could not open directory */
      perror ("");
      return EXIT_FAILURE;  
  }
  return EXIT_SUCCESS;
}
