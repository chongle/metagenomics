#include "Config.h"
#include "QueryDataset.h"
//#include "SingleKeyHashTable.h"
#include "OmegaHashTable.h"
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
//---------
	/*
	SingleKeyHashTable* singleKeyHashTable = new SingleKeyHashTable();
	{
	CLOCKSTART;
	MEMORYSTART;

	if (!singleKeyHashTable->insertQueryDataset(queryDataset)){
		cout << "Error: cannot build hash table" << endl;
		return false;
	}
	MEMORYSTOP;
	CLOCKSTOP;
	}*/

//---------
	OmegaHashTable *omegaHashTable=new OmegaHashTable();
	{
	CLOCKSTART;
	MEMORYSTART;
	omegaHashTable->insertDataset(queryDataset, Config::minimumOverlapLength);
	MEMORYSTOP;
	CLOCKSTOP;
	}

	MEMORYSTART;
	delete omegaHashTable;
	delete queryDataset;
//	delete singleKeyHashTable;
	MEMORYSTOP;
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



