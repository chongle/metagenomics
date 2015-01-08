/*
 * SingleKeyHashTable.cpp
 *
 *  Created on: Jan 6, 2015
 *      Author: qy2
 */

#include "SingleKeyHashTable.h"

SingleKeyHashTable::SingleKeyHashTable() {
	// TODO Auto-generated constructor stub
	this->minimumOverlapLength = Config::minimumOverlapLength;
	this->hashKeyLength = Config::hashKeyLength;
}

SingleKeyHashTable::~SingleKeyHashTable() {
	// TODO Auto-generated destructor stub
}

string SingleKeyHashTable::getReadSubstring(string mode, UINT64 readID)
{
	string subString="";
	string readString="";
	QueryRead * read = queryDataSet->getReadFromID(readID);
	if(mode=="forwardprefix")
	{

		if(read!=NULL)
		{
			readString = read->getSequence();
			subString = readString.substr(0,hashKeyLength);
		}
	}
	else if(mode == "forwardsuffix")
	{

		if(read!=NULL)
		{
			readString = read->getSequence();
			subString = readString.substr(readString.length() - hashKeyLength, hashKeyLength);
		}

	}
	else if(mode == "reverseprefix")
	{
		if(read!=NULL)
		{
			readString = read->reverseComplement();
			subString = readString.substr(0,hashKeyLength);
		}
	}
	else if(mode == "reversesuffix")
	{
		if(read!=NULL)
		{
			readString = read->reverseComplement();
			subString = readString.substr(readString.length() - hashKeyLength, hashKeyLength);
		}
	}
	else cout<< "mode is illegal"<<endl;
	return subString;
}

bool SingleKeyHashTable::InitializeAllHashTables()
{
	if(hashTableNameList.empty() || queryDataSet == NULL) return false;

	for(int i = 0; i< hashTableNameList.size(); i++)
	{
		string stringmode = hashTableNameList.at(i);
		HashTable *currentHashtable = new HashTable(this->hashKeyLength);
		currentHashtable->InitializeWithDataSize(queryDataSet->getNumberOfUniqueReads());
		hashTableMap.insert(std::pair<string, HashTable*>(stringmode, currentHashtable));
	}
	return true;
}

bool SingleKeyHashTable::insertQueryDataset(QueryDataset* querydataset)
{

	queryDataSet = querydataset;
	UINT64 datasetsize = queryDataSet->getNumberOfUniqueReads();
	hashTableNameList.push_back("forwardprefix");
	hashTableNameList.push_back("forwardsuffix");
	hashTableNameList.push_back("reverseprefix");
	hashTableNameList.push_back("reversesuffix");
	InitializeAllHashTables();

#pragma omp parallel for
	{
		for(int i = 0; i< hashTableNameList.size(); i++)
		{
			string stringmode = hashTableNameList.at(i);
			UINT64 currentID = 1;
			while(currentID<=datasetsize)
			{
				if(currentID%1000000 == 0)
					cout << setw(10) << currentID << " reads inserted in the hash table. " << endl;
				QueryRead * queryread = queryDataSet->getReadFromID(currentID);
				insertQueryRead(queryread, stringmode);
				currentID++;
			}

		}
	}


	return true;
}

bool SingleKeyHashTable::insertQueryRead(QueryRead *read, string mode)
{
	UINT64 readID = read->getIdentifier();
	string keystring = getReadSubstring(mode,readID);
	HashTable * currentHashTable = hashTableMap[mode];
	return currentHashTable->insertIntoHashTable(keystring,readID);

}
