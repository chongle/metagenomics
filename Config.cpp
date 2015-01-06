/*
 * Config.cpp
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#include "Config.h"

Config::Config() {
	// TODO Auto-generated constructor stub

}

Config::~Config() {
	// TODO Auto-generated destructor stub
}

static void Config::printHelp()
{
    std::cout << std::endl
              << "  Usage:" << std::endl
              << "    align [OPTION]...<PARAM>..." << std::endl
              << std::endl
              << "  <PARAM>" << std::endl
			  << "    -o/--operation\t ConstructGraph,MergePairedEndReads or CorrectErrors" << std::endl
              << "    -q/--query\t query file name" << std::endl  // file in fasta/fastq format.
              << "    -s/--subject\t subject file names (comma separated)" << std::endl  // Single-end file name list
              << "    --single/--double\t single hash table or double hash table method" << std::endl

              << "    -l\t minimum overlap length (default:40)" << std::endl
              << "    -k\t hash table key length (default:32)" << std::endl

              << std::endl
              << "  [OPTION]" << std::endl
              << "    -h/--help" << std::endl
			  << "    -f/--filter" << std::endl
			  << "    -a/--align" << std::endl
			  << std::endl;
}

static bool Config::setConfig(int argc, char **argv)
{
	vector<string> argumentsList;
	cout << "PRINTING ARGUMENTS" << endl;
	for(UINT64 i = 0; i < argc; i++)
	{
		cout << argv[i] << ' ';
	}
	cout << endl;
	while(argc--)
			argumentsList.push_back(*argv++);

	if(argumentsList.size() == 1)
	{
        Config::printHelp();
		return false;
	}

	for(UINT64 i = 1; i <= argumentsList.size()-1; i++)
	{
		if (argumentsList[i] == "-h" || argumentsList[i] == "--help")
	    {
	        Config::printHelp();
			return false;
	    }
		if (argumentsList[i] == "-f" || argumentsList[i] == "--filter")
	    {
	        Config::isFilter = true;

	    }
		if (argumentsList[i] == "-a" || argumentsList[i] == "--align")
	    {
	        Config::isFilter = false;

	    }
		if (argumentsList[i] == "-o" || argumentsList[i] == "--operation")
	    {
	        Config::operationCode = argumentsList[++i];

	    }
		else if(argumentsList[i] == "-s" || argumentsList[i] == "--subject")
		{
			string inputFilenames=argumentsList[++i];
            stringstream ss(inputFilenames);
            string item;

			while (getline(ss, item, ','))
			{
				Config::subjectFilenameList.push_back(item);
			}
		}
		else if(argumentsList[i] == "-q" || argumentsList[i] == "--query")
		{
			string inputFilename=argumentsList[++i];
			Config::queryFilename = inputFilename;

		}

		else if (argumentsList[i] == "-l")
			Config::minimumOverlapLength = atoi(argumentsList[++i].c_str());
		else if (argumentsList[i] == "-k")
			Config::hashKeyLength = atoi(argumentsList[++i].c_str());
		else if (argumentsList[i] == "--single")
			Config::isSingleKey = true;
		else if (argumentsList[i] == "--double")
			Config::isSingleKey = false;
		else
		{
          Config::printHelp();
          return false;
		}
	}



    if( (Config::subjectFilenameList.size() == 0) || Config::queryFilename == "")
    {
        cout << "missed -subject or -query input files!" << std::endl;
        printHelp();
		return false;
    }
    if(Config::operationCode == "")
    {
    	cout<< "operation code is missing" << std::endl;
    	printHelp();
    	return false;
    }
    return true;
}

static string Config::getOperation()
{
	return Config::operationCode;
}
static vector<string> Config::getSubjectDatasetFilenames()
{
	return Config::subjectFilenameList;
}
static string Config::getQueryDatasetFilename()
{
	return Config::queryFilename;
}
static UINT16 Config::getminimumOverlapLength()
{
	return Config::minimumOverlapLength;
}
static UINT16 Config::getHashKeyLength()
{
	return Config::hashKeyLength;
}
static UINT64 Config::getStreamChunkSize()
{
	return Config::streamChunkSize;
}

static bool Config::isSingleKeyHashTable()
{
	return Config::isSingleKey;
}
static bool Config::allowMismatch()
{
	return !Config::perfectMatch;
}
static bool Config::isQuerySingleEndReads()
{
	return Config::isSingleEnd;
}

static UINT16 Config::getMaxMismatch()
{
	return Config::maxMismatch;
}
static UINT16 Config::getMaxIndel()
{
	return Config::maxIndel;
}
static bool Config::isFilterOrAlign()
{
	return Config::isFilter;
}
