#include "Config.h"
#include "QueryDataset.h"
#include "HashTableMethod.h"
#include "SingleKeyHashTable.h"
#include "DoubleKeyHashTable.h"
#include "OmegaHashTable.h"
#include "OmegaGraphConstructor.h"
#include "OverlapGraphConstructor.h"
#include "QueryDatasetFilter.h"
//#include "OverlapGraphConstructor.h"
//#include "PairedReadMerger.h"
//#include "ErrorCorrector.h"


int main(int argc, char **argv)
{


	CLOCKSTART;


	// parse command line options:
	// align.exe -a -q query.fasta -s subject1.fasta ... -o CorrectErrors --single
	// align.exe -a --query query.fasta --subject subject1.fasta ... --operation CorrectErrors --double
	if(!Config::setConfig(argc, argv)){
		cout << "Error: wrong configurations" << endl;
		return false;
	}


//--------

	QueryDataset* queryDataset = new QueryDataset(Config::getQueryDatasetFilename());
	{
	CLOCKSTART;
	MEMORYSTART;

	if (!queryDataset->buildDataset()){
		cout << "Error: cannot build query dataset" << endl;
		return false;
	}
	MEMORYSTOP;
	CLOCKSTOP;

	}

//--------

	HashTableMethod* singleKeyHashTable = new SingleKeyHashTable(queryDataset);
	{
	CLOCKSTART;
	MEMORYSTART;
	cout << "-------singlekeyhashtable------" << endl;
	singleKeyHashTable->createHashTables();
	singleKeyHashTable->insertQueryDataset(queryDataset);
	MEMORYSTOP;
	CLOCKSTOP;

		{
		CLOCKSTART;
		MEMORYSTART;
		cout << "-------overlapgraph constructing------" << endl;
		OverlapGraphConstructor * overlapgraphConstructor = new OverlapGraphConstructor(singleKeyHashTable);
		overlapgraphConstructor->start();
		MEMORYSTOP;
		CLOCKSTOP;
		}
	}



//---------
	/*
	HashTableMethod* doubleKeyHashTable = new DoubleKeyHashTable(queryDataset);
	{
	CLOCKSTART;
	MEMORYSTART;
	cout << "-------doublekeyhashtable------" << endl;
	doubleKeyHashTable->createHashTables();
	doubleKeyHashTable->insertQueryDataset(queryDataset);
	MEMORYSTOP;
	CLOCKSTOP;
	}
	*/

//---------
/*
	OmegaHashTable *omegaHashTable=new OmegaHashTable();
	{
	CLOCKSTART;
	MEMORYSTART;
	cout << "-------omegahashtable------" << endl;
	omegaHashTable->insertDataset(queryDataset, Config::minimumOverlapLength);
	MEMORYSTOP;
	CLOCKSTOP;
	}
	OmegaGraphConstructor * omegaGraph = new OmegaGraphConstructor(omegaHashTable);
	omegaGraph->start();
	delete omegaGraph;
*/
//---------
/*
	OmegaHashTable *omegaHashTable=new OmegaHashTable();
	{
	CLOCKSTART;
	MEMORYSTART;
	cout << "-------remove contained reads omegahashtable------" << endl;
	omegaHashTable->insertDataset_half(queryDataset, Config::minimumOverlapLength);
	MEMORYSTOP;
	CLOCKSTOP;
	}
	QueryDatasetFilter * filter = new QueryDatasetFilter(omegaHashTable);
	filter->start();
	delete filter;
*/

//	delete omegaHashTable;
	delete queryDataset;
//	delete singleKeyHashTable;
//	delete doubleKeyHashTable;

	CLOCKSTOP;

/*
	if (Config::getOperation() == "ConstructGraph"){

		// start graph construction
		OverlapGraphConstructor* graphconstructor = new OverlapGraphConstructor();
		if(!graphconstructor->start())
		{
			cout<< "Failed to run OverlapGraphConstructor start()";
		}
		delete graphconstructor;

	}
	else if (Config::getOperation() == "MergePairedEndReads"){
		// start merge reads
		PairedReadMerger* merger = new PairedReadMerger();
		if(!merger->start())
		{
			cout<< "Failed to run PairedReadMerger start()";
		}
		delete merger;

	}
	else if (Config::getOperation() == "CorrectErrors"){
		// start correct errors
		ErrorCorrector* errorCorrector = new ErrorCorrector();
		if(!errorCorrector->start())
		{
			cout<< "Failed to run ErrorCorrector start()";
		}
		delete errorCorrector;
	}
	else{
		cout << "Error: invalid operation type."<< endl;
		return false;
	}
*/

//    int nthreads = omp_get_num_threads();
//    cout<<"Number of threads = "<<nthreads<<endl;
	return true;
}



