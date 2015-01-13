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
	this->maxMismatch = Config::maxMismatch;
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
	modeList.push_back("forwardprefix");
	modeList.push_back("forwardsuffix");
	modeList.push_back("reverseprefix");
	modeList.push_back("reversesuffix");
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


bool DoubleKeyHashTable::doAlignment(Alignment* align, string mode, int subjectStart)
{
}

bool DoubleKeyHashTable::subjectWindowRange(int& startpoint, int& stoppoint, string mode, string& subjectRead)
{
	if(mode=="forwardprefix" || mode == "reverseprefix")
	{
		startpoint = 0;
		stoppoint = subjectRead.length()-this->minimumOverlapLength;
	}
	else if(mode == "forwardsuffix" || mode == "reversesuffix")
	{
		startpoint = this->minimumOverlapLength - hashKeyLength;
		stoppoint = subjectRead.length()-hashKeyLength;

	}
	else return false;

	return true;
}

bool DoubleKeyHashTable::doubleKeySearch(edge & Edge)
{
	string subjectRead = Edge.subjectReadSequence;


	for(int i =0; i<this->modeList.size();i++)
	{


	string mode = this->modeList.at(i);

	int startpoint,stoppoint;
	if(subjectWindowRange(startpoint, stoppoint, mode, subjectRead)) return false;//guarantee it meets the minimum overlaplength requirement for alignments.

	for(int j=startpoint;j<=stoppoint;j++)
	{
	string subStringLeft = subjectRead.substr(startpoint, hashKeyLength);
	string subStringRight = subjectRead.substr(startpoint+hashKeyLength, hashKeyLength);
	string modeLeft = mode+"1";
	string modeRight = mode+"2";
	vector<UINT64> LeftIDList = hashTableMap.at(modeLeft)->getReadIDListOfReads(subStringLeft);
	vector<UINT64> RightIDList = hashTableMap.at(modeRight)->getReadIDListOfReads(subStringRight);
	vector<UINT64> Key1OnlyList = new vector<UINT64>;
	vector<UINT64> Key2OnlyList = new vector<UINT64>;
	vector<UINT64> BothKeyList = new vector<UINT64>;
	wennDiagramTwoLists(LeftIDList,RightIDList,Key1OnlyList, Key2OnlyList, BothKeyList);

	for(int k=0;k<BothKeyList.size();k++)
	{
		UINT64 currentID = BothKeyList.at(k);
		QueryRead* queryRead = queryDataSet->getReadFromID(currentID);
		string querySequence = queryRead->getSequence();
		string queryName = queryRead->getName();
		bool alignFlag = true;
		if(queryName>Edge.subjectReadName)alignFlag = false; //only align when queryName is alphabetical smaller than subject, removing half of the pair-wise alignment redundancy
		if(alignFlag)
		{
			Alignment* align = new Alignment();
			align->subjectReadName = Edge.subjectReadName;
			align->subjectReadSequence = subjectRead;
			align->queryRead = queryRead;
			if(doAlignment(align, mode, startpoint)==false) delete align;
			else Edge.alignmentList.push_back(align);

		}
	}
	}

	}
}

bool DoubleKeyHashTable::wennDiagramTwoLists(vector<UINT64> list1, vector<UINT64> list2, vector<UINT64>& list1only, vector<UINT64>& list2only, vector<UINT64>& list12)//sorted from smaller to larger ID in the list
{
	int i,j;
	i=0;j=0;
	while(i<list1.size()&&j<list2.size())
	{
		if(list1.at(i)<list2.at(i))
		{
			UINT64 smallvalue = list1.at(i);
			list1.push_back(smallvalue);
			i++;
		}
		else if(list1.at(i)>list2.at(i))
		{
			UINT64 smallvalue = list2.at(i);
			list2.push_back(smallvalue);
			j++;
		}
		else //same value in both list
		{
			UINT64 value = list1.at(i);
			list12.push_back(value);
			i++;j++;
		}
	}

	while(i<list1.size())
	{
		UINT64 smallvalue = list1.at(i);
		list1.push_back(smallvalue);
		i++;
	}

	while(j<list2.size())
	{
		UINT64 smallvalue = list2.at(i);
		list2.push_back(smallvalue);
		j++;
	}
	return true;
}
