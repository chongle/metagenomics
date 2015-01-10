/*
 * DoubleKeyHashTable.cpp
 *
 *  Created on: Jan 9, 2015
 *      Author: qy2
 */

#include "DoubleKeyHashTable.h"

DoubleKeyHashTable::DoubleKeyHashTable() {
	// TODO Auto-generated constructor stub
	this->minimumOverlapLength = Config::minimumOverlapLength;
	this->hashKeyLength = Config::hashKeyLength;
}

DoubleKeyHashTable::~DoubleKeyHashTable() {
	// TODO Auto-generated destructor stub
}

string DoubleKeyHashTable::getReadSubstring(string mode, UINT64 readID)
{
	string subString="";
	string readString="";
	QueryRead * read = queryDataSet->getReadFromID(readID);
	if(mode=="forwardprefix1")
	{

		if(read!=NULL)
		{
			readString = read->getSequence();
			subString = readString.substr(0,hashKeyLength);
		}
	}
	else if(mode=="forwardprefix2")
	{

		if(read!=NULL)
		{
			readString = read->getSequence();
			subString = readString.substr(hashKeyLength,hashKeyLength);
		}
	}
	else if(mode == "forwardsuffix1")
	{

		if(read!=NULL)
		{
			readString = read->getSequence();
			subString = readString.substr(readString.length() - hashKeyLength - hashKeyLength, hashKeyLength);
		}

	}
	else if(mode == "forwardsuffix2")
	{

		if(read!=NULL)
		{
			readString = read->getSequence();
			subString = readString.substr(readString.length() - hashKeyLength, hashKeyLength);
		}

	}
	else if(mode == "reverseprefix1")
	{
		if(read!=NULL)
		{
			readString = read->reverseComplement();
			subString = readString.substr(0,hashKeyLength);
		}
	}
	else if(mode == "reverseprefix2")
	{
		if(read!=NULL)
		{
			readString = read->reverseComplement();
			subString = readString.substr(hashKeyLength,hashKeyLength);
		}
	}
	else if(mode == "reversesuffix1")
	{
		if(read!=NULL)
		{
			readString = read->reverseComplement();
			subString = readString.substr(readString.length() - hashKeyLength - hashKeyLength, hashKeyLength);
		}
	}
	else if(mode == "reversesuffix2")
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

bool DoubleKeyHashTable::InitializeAllHashTables()
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

bool DoubleKeyHashTable::insertQueryDataset(QueryDataset* querydataset)
{

	queryDataSet = querydataset;
	UINT64 datasetsize = queryDataSet->getNumberOfUniqueReads();
	hashTableNameList.push_back("forwardprefix1");
	hashTableNameList.push_back("forwardprefix2");
	hashTableNameList.push_back("forwardsuffix1");
	hashTableNameList.push_back("forwardsuffix2");
	hashTableNameList.push_back("reverseprefix1");
	hashTableNameList.push_back("reverseprefix2");
	hashTableNameList.push_back("reversesuffix1");
	hashTableNameList.push_back("reversesuffix2");
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

bool DoubleKeyHashTable::insertQueryRead(QueryRead *read, string mode)
{
	UINT64 readID = read->getIdentifier();
	string keystring = getReadSubstring(mode,readID);
	HashTable * currentHashTable = hashTableMap.at(mode);
	return currentHashTable->insertIntoHashTable(keystring,readID);

}
