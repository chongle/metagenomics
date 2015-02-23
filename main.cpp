/*
 * main.cpp
 *
 *  Created on: Feb 11, 2015
 *      Author: qy2
 */
/*
 * main.cpp
 *
 * Created on: April 22, 2013
 * Author: Md. Bahlul Haider
 * Version: 1.0.1 (Beta) - Tae-Hyuk (Ted) Ahn: Validate input arguments
 */


#include "Common.h"
#include "Read.h"
#include "Dataset.h"
#include "HashTable.h"
#include "Edge.h"
#include "OverlapGraph.h"


//=============================================================================
// Help usage
//=============================================================================
void usage();


//=============================================================================
// Help usage
//=============================================================================
void usage()
{
    std::cout << std::endl
              << "  Usage:" << std::endl
              << "    omega [OPTION]...<PARAM>..." << std::endl
              << std::endl
              << "  <PARAM>" << std::endl
              << "    -pe\tpaired-end file names (comma separated)" << std::endl  // Paired-end file name in fasta/fastq format. mate pairs should be one after another in the file.
              << "    -se\tsingle-end file names (comma separated)" << std::endl  // Single-end file name in fasta/fastq format.
              << "    -l\tminimum overlap length" << std::endl  // Minimum overlap length for two reads to overlap in the overlap graph.
              << std::endl
              << "  [OPTION]" << std::endl
              << "    -h/--help" << std::endl
              << "    -f\tAll file name prefix (default: output)" << std::endl    // all output file with have this name with different extensions.
              << "    -s\tstart from unitig graph" << std::endl // -s means that the program will build the graph. Otherwise it will load the graph from the unitig graph file.
              << std::endl;
}


void parseArguments(int argc, char **argv, vector<string> & pairedEndFileNames, vector<string> & singleEndFileNames,string & allFileName, UINT64 & minimumOverlapLength, bool & startFromUnitigGraph);

int main(int argc, char **argv)
{
	CLOCKSTART;

	MEMORYSTART;
	UINT64 minimumOverlapLength;
	vector<string> pairedEndFileNames, singleEndFileNames;
	UINT64 counter;
	UINT64 iteration = 0;
	string allFileName;
	bool startFromUnitigGraph = false;
	parseArguments(argc, argv, pairedEndFileNames, singleEndFileNames, allFileName, minimumOverlapLength, startFromUnitigGraph);
	Dataset *dataSet = new Dataset(pairedEndFileNames, singleEndFileNames, minimumOverlapLength);
	MEMORYSTOP;

	//YAO 	OverlapGraph *overlapGraph;

	//YAO 	if(startFromUnitigGraph) // Read the graph from unitig file
	//YAO 	{
	//YAO 		overlapGraph = new OverlapGraph();
	//YAO 		overlapGraph->setDataset(dataSet);
	//YAO 		overlapGraph->readGraphFromFile(allFileName+"_unitig.unitig");
	//YAO 		overlapGraph->sortEdges();
	//YAO 	}
	//YAO 	else // Build the graph from the scratch
	{
				CLOCKSTART;
				MEMORYSTART;
				HashTable *hashTable=new HashTable();
				hashTable->insertDataset(dataSet, minimumOverlapLength);
				MEMORYSTOP;
				CLOCKSTOP;
				{
				CLOCKSTART;
				MEMORYSTART;
				OverlapGraph* overlapGraph=new OverlapGraph(hashTable); //hashTable deleted by this function after building the graph

		 		delete overlapGraph;
		 		MEMORYSTOP;
		 		CLOCKSTOP;
				}
				delete hashTable;
		//YAO 		dataSet->saveReads(allFileName+"_sortedReads.fasta");
		//YAO 		overlapGraph->sortEdges();
		//YAO 		overlapGraph->saveGraphToFile(allFileName+"_unitig.unitig");
	}

/*YAO
	overlapGraph->calculateFlow(allFileName+"_flow.input", allFileName+"_flow.output");
	cout << "nodes: " << overlapGraph->getNumberOfNodes() << " edges: " << overlapGraph->getNumberOfEdges() << endl;
	overlapGraph->printGraph(allFileName+"_graph1.gdl", allFileName+"_contigs1.fasta");

	overlapGraph->removeAllSimpleEdgesWithoutFlow();

	overlapGraph->calculateMeanAndSdOfInsertSize();



	do
	{
		// Mate pair paths are used to simplify the graph in this step
		cout << endl;
		cout << "===============================================================================================================================================" <<endl;
		cout << "FIRST LOOP ITERATION " << ++iteration << endl;
		cout << "===============================================================================================================================================" <<endl;
		overlapGraph->simplifyGraph();
		counter = overlapGraph->findSupportByMatepairsAndMerge();
	} while (counter > 0 && iteration < loopLimit); // To avoid infinite loops

	overlapGraph->printGraph(allFileName+"_graph2.gdl", allFileName+"_contigs2.fasta");

	iteration = 0;
	do
	{
		// Scaffolder
		cout << endl;
		cout << "===============================================================================================================================================" <<endl;
		cout << "SECOND LOOP ITERATION " << ++iteration << endl;
		cout << "===============================================================================================================================================" <<endl;
		overlapGraph->simplifyGraph();
		counter = overlapGraph->scaffolder();

	} while (counter > 0 && iteration < loopLimit);// To avoid infinite loops

	overlapGraph->printGraph(allFileName+"_graph3.gdl", allFileName+"_contigs3.fasta");

	iteration = 0;
	do
	{
		// Coverage depth information is used to resolve ambiguity
		cout << endl;
		cout << "===============================================================================================================================================" <<endl;
		cout << "THIRD LOOP ITERATION " << ++iteration << endl;
		cout << "===============================================================================================================================================" <<endl;
		overlapGraph->simplifyGraph();
		counter = overlapGraph->resolveNodes();

	} while (counter > 0 && iteration < loopLimit);// To avoid infinite loops

	overlapGraph->printGraph(allFileName+"_graph4.gdl", allFileName+"_contigs4.fasta");
*/
	delete dataSet;
//	delete overlapGraph;
	CLOCKSTOP;
}



/**********************************************************************************************************************
	Parse the input arguments
**********************************************************************************************************************/

void parseArguments(int argc, char **argv, vector<string> & pairedEndFileNames, vector<string> & singleEndFileNames, string & allFileName, UINT64 & minimumOverlapLength, bool & startFromUnitigGraph)
{
	allFileName = "output";
	minimumOverlapLength = 0;
	startFromUnitigGraph = false;
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
        usage();
		exit(0);
	}

	for(UINT64 i = 1; i <= argumentsList.size()-1; i++)
	{
		if (argumentsList[i] == "-h" || argumentsList[i] == "--help")
	    {
            usage();
		    exit(0);
	    }
		else if(argumentsList[i] == "-pe")
		{
			string inputFilenames=argumentsList[++i];
            stringstream ss(inputFilenames);
            string item;

			while (getline(ss, item, ','))
			{
				pairedEndFileNames.push_back(item);
			}
		}
		else if(argumentsList[i] == "-se")
		{
			string inputFilenames=argumentsList[++i];
            stringstream ss(inputFilenames);
            string item;

			while (getline(ss, item, ','))
			{
				singleEndFileNames.push_back(item);
			}
		}
		else if (argumentsList[i] == "-f")
			allFileName = argumentsList[++i];
		else if (argumentsList[i] == "-l")
			minimumOverlapLength = atoi(argumentsList[++i].c_str());
		else if (argumentsList[i] == "-s")
			startFromUnitigGraph = true;
		else
		{
            usage();
		    exit(0);
		}
	}

    if(minimumOverlapLength == 0)
    {
        cout << "missed -l option!" << std::endl;
        usage();
		exit(0);
    }

    if( (pairedEndFileNames.size() == 0) && (singleEndFileNames.size() == 0) )
    {
        cout << "missed -pe or -se input files!" << std::endl;
        usage();
		exit(0);
    }
}




