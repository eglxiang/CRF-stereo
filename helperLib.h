#ifndef HELPERLIB_H
#define HELPERLIB_H

#include <vector>
#include <iostream>
#include <algorithm>
#include "crfstereo.h"
#include "expectations.h"
#include "evaldisps.h"

#ifdef _MSC_VER        // Visual C++
#include "dirent.h"
const char dirs = '\\';
#else
const char dirs = '/';
# include <dirent.h>
#endif

using std::vector;

// read images from a directory.
int readImages(char *, vector<CByteImage>&, bool);

#endif
