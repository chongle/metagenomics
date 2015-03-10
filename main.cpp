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
	if(!Config::setConfig(argc, argv)){
		cout << "Error: wrong configurations" << endl;
		return false;
	}


//--------generate data sets-----------

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

//--------single key hash table method------
if(Config::getOperation()=="ConstructGraph" && Config::getHashTableType() == "single")
{
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
		delete overlapgraphConstructor;
		}
	}
	delete singleKeyHashTable;
}


//---------double key hash table method-------
if(Config::getOperation()=="ConstructGraph" && Config::getHashTableType() == "double")
{
	HashTableMethod* doubleKeyHashTable = new DoubleKeyHashTable(queryDataset);
	{
	CLOCKSTART;
	MEMORYSTART;
	cout << "-------doublekeyhashtable------" << endl;
	doubleKeyHashTable->createHashTables();
	doubleKeyHashTable->insertQueryDataset(queryDataset);
	MEMORYSTOP;
	CLOCKSTOP;
	{
	CLOCKSTART;
	MEMORYSTART;
	cout << "-------overlapgraph constructing------" << endl;
	OverlapGraphConstructor * overlapgraphConstructor = new OverlapGraphConstructor(doubleKeyHashTable);
	overlapgraphConstructor->start();
	MEMORYSTOP;
	CLOCKSTOP;
	delete overlapgraphConstructor;
	}
	}
	delete doubleKeyHashTable;
}

//---------omega hash table method------
if(Config::getOperation()=="ConstructGraph" && Config::getHashTableType() == "omega")
{
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
	delete omegaHashTable;
}
//---------using omega hash table to filter the identical and contained reads---
if(Config::getOperation()=="RemoveContainedReads")
{
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
	delete omegaHashTable;
}


	delete queryDataset;


	CLOCKSTOP;


	return true;
}



