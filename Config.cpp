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

//Here is to create and initialize the static members
string Config::operationCode="";//"ConstructGraph","MergePairedEndReads","CorrectErrors"
string Config::queryFilename = "";
string Config::hashtabletype = "single";

vector<string> Config::subjectFilenameList;
string Config::outputfilename = "out.txt";

UINT16 Config::minimumOverlapLength=40;
UINT16 Config::hashKeyLength = 39;//key needs to be smaller than minimumoverlap, if double key, 2*key<minoverlap
UINT16 Config::hashKeyLength_left = 19;
UINT16 Config::hashKeyLength_right = 20;
UINT64 Config::streamChunkSize = 400;
bool Config::isSingleKey = true;
bool Config::useID = false;

bool Config::isSingleEnd = true; //subject is always treated as single end, while the query is only paired end when using "MergePairedEndReads"
bool Config::isFilter = false;

bool Config::perfectMatch = false;
UINT16 Config::maxMismatch = 1; //only valid when perfectMatch is set to false
UINT16 Config::maxIndel = 0; //only valid when perfectMatch is set to false

UINT16 Config::numberOfThreads = 1;




void Config::printHelp()
{
    std::cout << std::endl
              << "  Usage:" << std::endl
              << "    align_test [OPTION]...<PARAM>..." << std::endl
              << std::endl
              << "  <PARAM>" << std::endl
			  << "    -i/--instance\t RemoveContainedReads,ConstructGraph,MergePairedEndReads or CorrectErrors" << std::endl
              << "    -q/--query\t query file name" << std::endl  // file in fasta/fastq format.
              << "    -s/--subject\t subject file name(s) (comma separated)" << std::endl  // Single-end file name list
              << "    -ht/--hashtable single/double/omega\t single hash table or double hash table method or omega hash table method" << std::endl
              << "    single hash table method is default setting except that in RemoveContainedReads omega is default" << std::endl

              << "    single : please set up the key length -k" << std::endl
              << "    double : please set up the left key length -lk and right key length -rk" << std::endl
              << "    omega : key length will be minimum_overlap-1, will ignore -m setting because no mismatch is allowed." << std::endl

              << "    -l\t minimum overlap length (default:40)" << std::endl
              << "    -k\t single hash table key length (default:39)" << std::endl
              << "    -lk\t double hash table left key length (default:19)" << std::endl
              << "    -rk\t double hash table right key length (default:20)" << std::endl
			  << "    -m\t maximum allowed mismatch (default:1)" << std::endl
			  << "    -t\t number of threads (default:1 [single thread])" << std::endl
			  << "    -z\t stream chunk size of subject read file (default:400)" << std::endl
			  << "    -o/--out\t output file name (default:out.txt)" << std::endl

              << std::endl
              << "  [OPTION]" << std::endl
              << "    -h/--help\t only print out the help contents" << std::endl
			  << "    -id/--ID\t use and print IDs in the fasta file rather than the names" << std::endl
//			  << "    -f/--filter" << std::endl
//			  << "    -a/--align" << std::endl
			  << std::endl

	<< "Example: ./align_test -i RemoveContainedReads --ID --subject sreads1.fasta,sreads2.fasta --query qreads.fasta --out cleanreads.fasta -l 40 -t 4" << std::endl
	<< "Example: ./align_test -i ConstructGraph -ht omega --subject sreads1.fasta,sreads2.fasta --query qreads.fasta --out outgraph.txt -l 40 -t 4" << std::endl
    << "Example: ./align_test -i ConstructGraph -ht single --subject sreads1.fasta,sreads2.fasta --query qreads.fasta --out outgraph.txt -l 40 -m 1 -k 32 -t 4 -z 1000" << std::endl
	<< "Example: ./align_test -i ConstructGraph -ht double --subject sreads1.fasta,sreads2.fasta --query qreads.fasta --out outgraph.txt -l 40 -m 1 -lk 16 -rk 16 -t 4" << std::endl

    << "Example: ./align_test -i MergePairedEndReads -ht double --subject sreads1.fasta,sreads2.fasta --query qreads.fasta --out outgraph.txt -l 40 -m 1 -lk 16 -rk 16 -t 4" << std::endl;

}

bool Config::setConfig(int argc, char **argv)
{
	Config::subjectFilenameList.clear();

	vector<string> argumentsList;
	cout << "PRINTING ARGUMENTS" << endl;
	for(int i = 0; i < argc; i++)
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

		if (argumentsList[i] == "-id" || argumentsList[i] == "--ID")
	    {
	        Config::useID = true;

	    }
		else if (argumentsList[i] == "-i" || argumentsList[i] == "--instance")
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
		else if(argumentsList[i] == "--out")
		{
			string outputFilename=argumentsList[++i];
			Config::outputfilename = outputFilename;

		}
		else if(argumentsList[i] == "-ht" || argumentsList[i] == "--hashtable")
		{
			string type=argumentsList[++i];
			Config::hashtabletype = type;

		}

		else if (argumentsList[i] == "-l")
			Config::minimumOverlapLength = atoi(argumentsList[++i].c_str());
		else if (argumentsList[i] == "-k")
			Config::hashKeyLength = atoi(argumentsList[++i].c_str());
		else if (argumentsList[i] == "-lk")
			Config::hashKeyLength_left = atoi(argumentsList[++i].c_str());
		else if (argumentsList[i] == "-rk")
			Config::hashKeyLength_right = atoi(argumentsList[++i].c_str());
		else if (argumentsList[i] == "-m")
					Config::maxMismatch = atoi(argumentsList[++i].c_str());
		else if (argumentsList[i] == "-t")
					Config::numberOfThreads = atoi(argumentsList[++i].c_str());
		else if (argumentsList[i] == "-z")
					Config::streamChunkSize = atoi(argumentsList[++i].c_str());
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

string Config::getOperation()
{
	return Config::operationCode;
}
vector<string> Config::getSubjectDatasetFilenames()
{
	return Config::subjectFilenameList;
}
string Config::getQueryDatasetFilename()
{
	return Config::queryFilename;
}
UINT16 Config::getminimumOverlapLength()
{
	return Config::minimumOverlapLength;
}
UINT16 Config::getHashKeyLength()
{
	return Config::hashKeyLength;
}
UINT16 Config::getHashLeftKeyLength()
{
	return Config::hashKeyLength_left;
}
UINT16 Config::getHashRightKeyLength()
{
	return Config::hashKeyLength_right;
}
UINT64 Config::getStreamChunkSize()
{
	return Config::streamChunkSize;
}

bool Config::isSingleKeyHashTable()
{
	return Config::isSingleKey;
}
string Config::getHashTableType()
{
	return Config::hashtabletype;
}
bool Config::allowMismatch()
{
	return !Config::perfectMatch;
}
bool Config::isQuerySingleEndReads()
{
	return Config::isSingleEnd;
}

UINT16 Config::getMaxMismatch()
{
	return Config::maxMismatch;
}
UINT16 Config::getMaxIndel()
{
	return Config::maxIndel;
}
bool Config::isFilterOrAlign()
{
	return Config::isFilter;
}
