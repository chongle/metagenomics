/*
 * DoubleKeyHashTable.cpp
 *
 *  Created on: Feb 15, 2015
 *      Author: qy2
 */

#include "DoubleKeyHashTable.h"

DoubleKeyHashTable::DoubleKeyHashTable() {
	// TODO Auto-generated constructor stub
	this->minimumOverlapLength = Config::minimumOverlapLength;
	this->hashKeyLength_left = Config::hashKeyLength_left;
	this->hashKeyLength_right = Config::hashKeyLength_right;
	this->maxMismatch = Config::maxMismatch;
	this->numberOfMode = 8;
	this->dataSet = NULL;
	this->hashTable = NULL;
	numberOfHashCollision = 0;
	maxSingleHashCollision = 0;
}

DoubleKeyHashTable::DoubleKeyHashTable(QueryDataset * qDataset) {
	// TODO Auto-generated constructor stub
	this->minimumOverlapLength = Config::minimumOverlapLength;
	this->hashKeyLength_left = Config::hashKeyLength_left;
	this->hashKeyLength_right = Config::hashKeyLength_right;
	this->maxMismatch = Config::maxMismatch;
	this->numberOfMode = 8;
	this->dataSet = qDataset;
	this->hashTable = NULL;
	numberOfHashCollision = 0;
	maxSingleHashCollision = 0;
}

DoubleKeyHashTable::~DoubleKeyHashTable() {
	// TODO Auto-generated destructor stub
	if(this->hashTable!=NULL)
	delete hashTable;
	this->dataSet = NULL;//we don't delete dataSet here
}

bool DoubleKeyHashTable::createHashTables()
{
	if(this->dataSet==NULL)
	{
		cout<<"no data set"<<endl;
		return false;
	}
	else if(this->hashTable!=NULL)
	{
		cout<<"Hash Table already exists."<<endl;
		return false;
	}
	else
	{
		this->hashTable = new HashTable(this->hashKeyLength_left+this->hashKeyLength_right,  this->dataSet, this->numberOfMode);

		return true;
	}
}

//left key
// 000 = 0 means prefix of the forward string.
// 001 = 1 means suffix of the forward string.
// 010 = 2 means prefix of the reverse string.
// 011 = 3 means suffix of the reverse string.
//right key
// 100 = 4 means prefix of the forward string.
// 101 = 5 means suffix of the forward string.
// 110 = 6 means prefix of the reverse string.
// 111 = 7 means suffix of the reverse string.
string DoubleKeyHashTable::getReadSubstring(UINT64 readID, int mode)
{
	QueryRead * read = this->dataSet->getReadFromID(readID);
	string str = (mode == 0 || mode == 1 || mode == 4 || mode == 5) ? read->getSequence() : read->reverseComplement();
	string subStr;
	switch(mode)
	{
	case 0:
	case 2:
		subStr = str.substr(0,this->hashKeyLength_left);
		break;
	case 1:
	case 3:
		subStr = str.substr(str.length() - this->hashKeyLength_left - this->hashKeyLength_right, this->hashKeyLength_left);
		break;
	case 4:
	case 6:
		subStr = str.substr(this->hashKeyLength_left,this->hashKeyLength_right);
		break;
	case 5:
	case 7:
		subStr = str.substr(str.length() - this->hashKeyLength_right, this->hashKeyLength_right);
		break;
	default:
		cout<<"wrong mode number"<<endl;
		subStr="";
		break;
	}
	return subStr;
}

bool DoubleKeyHashTable::insertQueryDataset(QueryDataset* d)
{
	if(this->hashTable==NULL)
	{
		cout<<"Hash Table hasn't been created yet."<<endl;
		return false;
	}
	else
	{
		UINT64 datasetsize = this->dataSet->getNumberOfUniqueReads();
//		omp_set_dynamic(0);
//		omp_set_num_threads(Config::numberOfThreads);
//		#pragma omp parallel
//			{
//		#pragma omp for schedule(dynamic)
		UINT64 currentID = 1;
		while(currentID<=datasetsize)
		{

			if(currentID%1000000 == 0)
				cout << setw(10) << currentID << " reads inserted in the hash table. " << endl;
			QueryRead * read = this->dataSet->getReadFromID(currentID);
			string forwardRead = read->getSequence();
			string reverseRead = read->reverseComplement();
			string prefixForward_left = forwardRead.substr(0,this->hashKeyLength_left);
			string prefixForward_right = forwardRead.substr(this->hashKeyLength_left,this->hashKeyLength_right);
			string suffixForward_left = forwardRead.substr(forwardRead.length() - this->hashKeyLength_left - this->hashKeyLength_right,this->hashKeyLength_left);
			string suffixForward_right = forwardRead.substr(forwardRead.length() - this->hashKeyLength_right, this->hashKeyLength_right);
			string prefixReverse_left = reverseRead.substr(0,this->hashKeyLength_left);
			string prefixReverse_right = reverseRead.substr(this->hashKeyLength_left,this->hashKeyLength_right);
			string suffixReverse_left = reverseRead.substr(reverseRead.length() - this->hashKeyLength_left - this->hashKeyLength_right,this->hashKeyLength_left);
			string suffixReverse_right = reverseRead.substr(reverseRead.length() - this->hashKeyLength_right, this->hashKeyLength_right);

			insertQueryRead(read, prefixForward_left, 0);
			insertQueryRead(read, suffixForward_left, 1);
			insertQueryRead(read, prefixReverse_left, 2);
			insertQueryRead(read, suffixReverse_left, 3);
			insertQueryRead(read, prefixForward_right, 4);
			insertQueryRead(read, suffixForward_right, 5);
			insertQueryRead(read, prefixReverse_right, 6);
			insertQueryRead(read, suffixReverse_right, 7);
			currentID++;
		}
		cout<<"Hash Table "<<" maximum collision number is: "<< this->numberOfHashCollision<<endl;
		cout<<"Hash Table "<<" maximum single read collision number is: "<< this->maxSingleHashCollision<<endl;


//			}//end of parallel

		return true;

	}
}
bool DoubleKeyHashTable::insertQueryRead(QueryRead *read, string subString, int mode)
{
	UINT64 currentCollision =0;

	UINT64 index = this->hashTable->hashFunction(subString);
	while(!hashTable->isEmptyAt(index))
	{
		map<int,vector<UINT64>*>::iterator p=hashTable->getDataVectorsAt(index)->begin();
		int keymode = p->first;
		UINT64 keyreadID = p->second->at(0);

		string keyStr = this->getReadSubstring(keyreadID,keymode);

		if(keyStr == subString)
				break;
		numberOfHashCollision++;
		currentCollision++;
		index = (index == hashTable->getHashTableSize() - 1) ? 0: index + 1; 	// Increment the index
	}
	hashTable->insertValueAt(index,mode,this->numberOfMode,read->getIdentifier());							// Add the string in the list.

	if(currentCollision> this->maxSingleHashCollision)
		this->maxSingleHashCollision = currentCollision;
	if(currentCollision > 1000)
	{
		cout << currentCollision << " collisions for read " << read->getIdentifier() << " " << subString << " " << mode << endl;
	}
	return true;
}
