/*
 * Config.h
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#ifndef CONFIG_H_
#define CONFIG_H_

// Common headers:

//multi-thread library OPENMP
#include <omp.h>

// C headers:
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>

// C++ headers:
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <string>
#include <vector>
#include <utility>
#include <limits>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <streambuf>
#include <map>


using namespace std;

typedef unsigned char UINT8;
typedef unsigned short UINT16;
typedef short INT16;
typedef unsigned long UINT32;
typedef long INT32;
typedef unsigned long long UINT64;
typedef long long INT64;

class Config {

private:
	static string operationCode;
	static string queryFilenameList;
	static vector<string> subjectFilenameList;

public:
	Config();
	~Config();
	static bool setConfig(int argc, char **argv);
	static string getOperation();
	static vector<string> getSubjectDatasetFilenames(); // subject is dataset B streamed from hard drive
	static string getQueryDatasetFilename();	// query is database A loaded to memory
	static int getminimumOverlapLength();
	static int getSingleKeyHashTableKeyLength();
	static int getDoubleKeyHashTableKeyLength();



};

#endif /* CONFIG_H_ */
