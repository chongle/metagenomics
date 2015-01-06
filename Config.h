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
	static string operationCode="";//"ConstructGraph","MergePairedEndReads","CorrectErrors"
	static string queryFilename = "";
	static vector<string> subjectFilenameList;

	static UINT16 minimumOverlapLength=40;
	static UINT16 hashKeyLength = 32;
	static UINT64 streamChunkSize = 100;
	static bool isSingleKey = true;

	static bool isSingleEnd = true; //subject is always treated as single end, while the query is only paired end when using "MergePairedEndReads"
	static bool isFilter = false;

	static bool perfectMatch = false;
	static UINT16 maxMismatch = 1; //only valid when perfectMatch is set to false
	static UINT16 maxIndel = 0; //only valid when perfectMatch is set to false

	static void printHelp();


public:
	Config();
	~Config();
	static bool setConfig(int argc, char **argv);
	static string getOperation();
	static vector<string> getSubjectDatasetFilenames(); // subject is dataset B streamed from hard drive
	static string getQueryDatasetFilename();	// query is database A loaded to memory
	static UINT16 getminimumOverlapLength();
	static UINT16 getHashKeyLength();
	static UINT64 getStreamChunkSize();

	static bool isSingleKeyHashTable();
	static bool allowMismatch();
	static bool isQuerySingleEndReads();

	static UINT16 getMaxMismatch();
	static UINT16 getMaxIndel();
	static bool isFilterOrAlign();



};

#endif /* CONFIG_H_ */
