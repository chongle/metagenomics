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
    hashTableNameList.clear();
	modeList.clear();
    hashTableMap.clear();
	queryDataSet = NULL;
}

DoubleKeyHashTable::~DoubleKeyHashTable() {
	// TODO Auto-generated destructor stub
	for(int i = 0; i< hashTableNameList.size(); i++)
	{
		string stringmode = hashTableNameList.at(i);
		delete hashTableMap.at(stringmode);
	}
	hashTableMap.clear();
	hashTableNameList.clear();
	modeList.clear();
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

	for(unsigned int i = 0; i< hashTableNameList.size(); i++)
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

	omp_set_dynamic(0);
	omp_set_num_threads(Config::numberOfThreads);
#pragma omp parallel
	{
#pragma omp for
		for(unsigned int i = 0; i< hashTableNameList.size(); i++)
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
	}//end of parallel


	return true;
}

bool DoubleKeyHashTable::insertQueryRead(QueryRead *read, string mode)
{
	UINT64 readID = read->getIdentifier();
	string keystring = getReadSubstring(mode,readID);
	HashTable * currentHashTable = hashTableMap.at(mode);
	return currentHashTable->insertIntoHashTable(keystring,readID);

}




bool DoubleKeyHashTable::subjectWindowRange(int& startpoint, int& stoppoint, string mode, string& subjectRead)
{
	if(mode=="forwardprefix" || mode == "reverseprefix")
	{
		startpoint = 0;// can from -hashKeyLength
		stoppoint = subjectRead.length()-this->minimumOverlapLength;
	}
	else if(mode == "forwardsuffix" || mode == "reversesuffix")
	{
		startpoint = this->minimumOverlapLength - hashKeyLength*2;
		stoppoint = subjectRead.length()-hashKeyLength*2; // can be subjectRead.length()-hashKeyLength

	}
	else return false;

	return true;
}

bool DoubleKeyHashTable::doubleKeySearch(edge & Edge)
{
	string subjectRead = Edge.subjectReadSequence;


	for(unsigned int i =0; i<this->modeList.size();i++)
	{
	string mode = this->modeList.at(i);

	int startpoint,stoppoint;
	if(!subjectWindowRange(startpoint, stoppoint, mode, subjectRead)) return false;//guarantee it meets the minimum overlaplength requirement for alignments.

	for(int j=startpoint;j<=stoppoint;j++)
	{

	string subStringLeft = subjectRead.substr(j, hashKeyLength);
	string subStringRight = subjectRead.substr(j+hashKeyLength, hashKeyLength);
	string modeLeft = mode+"1";
	string modeRight = mode+"2";
	vector<UINT64>* LeftIDList = hashTableMap.at(modeLeft)->getReadIDListOfReads(subStringLeft);
	vector<UINT64>* RightIDList = hashTableMap.at(modeRight)->getReadIDListOfReads(subStringRight);
	vector<UINT64>* Key1OnlyList = new vector<UINT64>;
	vector<UINT64>* Key2OnlyList = new vector<UINT64>;
	vector<UINT64>* BothKeyList = new vector<UINT64>;
	Key1OnlyList->clear();
	Key2OnlyList->clear();
	BothKeyList->clear();

	wennDiagramTwoLists(LeftIDList,RightIDList,Key1OnlyList, Key2OnlyList, BothKeyList);


	for(unsigned int k=0;k<BothKeyList->size();k++)
	{
		UINT64 currentID = BothKeyList->at(k);
		QueryRead* queryRead = queryDataSet->getReadFromID(currentID);
		string querySequence = queryRead->getSequence();
		string queryName = queryRead->getName();
		bool alignFlag = true;
		if(queryName>=Edge.subjectReadName)alignFlag = false; //only align when queryName is alphabetical smaller than subject, removing half of the pair-wise alignment redundancy
		if(alignFlag)
		{
			Alignment* align = new Alignment();
			align->subjectReadName = Edge.subjectReadName;
			align->subjectReadSequence = subjectRead;
			align->queryRead = queryRead;
			if(createAlignment(align, j, j+this->hashKeyLength, j+this->hashKeyLength+this->hashKeyLength-1, mode, 3)==false) delete align;
			else Edge.alignmentList.push_back(align);

		}
	}

	for(unsigned int k=0;k<Key1OnlyList->size();k++)
	{
		UINT64 currentID = Key1OnlyList->at(k);
		QueryRead* queryRead = queryDataSet->getReadFromID(currentID);
		string querySequence = queryRead->getSequence();
		string queryName = queryRead->getName();
		bool alignFlag = true;
		if(queryName>=Edge.subjectReadName)alignFlag = false; //only align when queryName is alphabetical smaller than subject, removing half of the pair-wise alignment redundancy
		if(alignFlag)
		{
			Alignment* align = new Alignment();
			align->subjectReadName = Edge.subjectReadName;
			align->subjectReadSequence = subjectRead;
			align->queryRead = queryRead;
			if(createAlignment(align, j, j+this->hashKeyLength,j+this->hashKeyLength+this->hashKeyLength-1, mode, 1)==false) delete align;
			else Edge.alignmentList.push_back(align);

		}
	}
	for(unsigned int k=0;k<Key2OnlyList->size();k++)
	{
		UINT64 currentID = Key2OnlyList->at(k);
		QueryRead* queryRead = queryDataSet->getReadFromID(currentID);
		string querySequence = queryRead->getSequence();
		string queryName = queryRead->getName();
		bool alignFlag = true;
		if(queryName>=Edge.subjectReadName)alignFlag = false; //only align when queryName is alphabetical smaller than subject, removing half of the pair-wise alignment redundancy
		if(alignFlag)
		{
			Alignment* align = new Alignment();
			align->subjectReadName = Edge.subjectReadName;
			align->subjectReadSequence = subjectRead;
			align->queryRead = queryRead;
			if(createAlignment(align, j, j+this->hashKeyLength,j+this->hashKeyLength+this->hashKeyLength-1, mode, 2)==false) delete align;
			else Edge.alignmentList.push_back(align);

		}
	}

	}//loop for window j

	}//loop of i for each mode
	return true;
}

bool DoubleKeyHashTable::doubleKeySearch(SubjectAlignment & subjectAlign)
{
	string subjectRead = subjectAlign.subjectReadSequence;


	for(unsigned int i =0; i<this->modeList.size();i++)
	{
	string mode = this->modeList.at(i);

	int startpoint,stoppoint;
	if(!subjectWindowRange(startpoint, stoppoint, mode, subjectRead)) return false;//guarantee it meets the minimum overlaplength requirement for alignments.

	for(int j=startpoint;j<=stoppoint;j++)
	{
	string subStringLeft = subjectRead.substr(j, hashKeyLength);
	string subStringRight = subjectRead.substr(j+hashKeyLength, hashKeyLength);
	string modeLeft = mode+"1";
	string modeRight = mode+"2";
	vector<UINT64>* LeftIDList = hashTableMap.at(modeLeft)->getReadIDListOfReads(subStringLeft);
	vector<UINT64>* RightIDList = hashTableMap.at(modeRight)->getReadIDListOfReads(subStringRight);
	vector<UINT64>* Key1OnlyList = new vector<UINT64>;
	vector<UINT64>* Key2OnlyList = new vector<UINT64>;
	vector<UINT64>* BothKeyList = new vector<UINT64>;
	wennDiagramTwoLists(LeftIDList,RightIDList,Key1OnlyList, Key2OnlyList, BothKeyList);

	for(unsigned int k=0;k<BothKeyList->size();k++)
	{
		UINT64 currentID = BothKeyList->at(k);
		QueryRead* queryRead = queryDataSet->getReadFromID(currentID);
		string querySequence = queryRead->getSequence();
		string queryName = queryRead->getName();
		bool alignFlag = true;
		if(queryName==subjectAlign.subjectReadName)alignFlag = false; //only align when queryName is alphabetical smaller than subject, removing half of the pair-wise alignment redundancy
		if(alignFlag)
		{
			Alignment* align = new Alignment();
			align->subjectReadName = subjectAlign.subjectReadName;
			align->subjectReadSequence = subjectRead;
			align->queryRead = queryRead;
			if(createAlignment(align, j, j+this->hashKeyLength, j+this->hashKeyLength+this->hashKeyLength-1, mode, 3)==false) delete align;
			else subjectAlign.queryAlignmentList.push_back(align);

		}
	}

	for(unsigned int k=0;k<Key1OnlyList->size();k++)
	{
		UINT64 currentID = Key1OnlyList->at(k);
		QueryRead* queryRead = queryDataSet->getReadFromID(currentID);
		string querySequence = queryRead->getSequence();
		string queryName = queryRead->getName();
		bool alignFlag = true;
		if(queryName==subjectAlign.subjectReadName)alignFlag = false; //only align when queryName is alphabetical smaller than subject, removing half of the pair-wise alignment redundancy
		if(alignFlag)
		{
			Alignment* align = new Alignment();
			align->subjectReadName = subjectAlign.subjectReadName;
			align->subjectReadSequence = subjectRead;
			align->queryRead = queryRead;
			if(createAlignment(align, j, j+this->hashKeyLength,j+this->hashKeyLength+this->hashKeyLength-1, mode, 1)==false) delete align;
			else subjectAlign.queryAlignmentList.push_back(align);

		}
	}
	for(unsigned int k=0;k<Key2OnlyList->size();k++)
	{
		UINT64 currentID = Key2OnlyList->at(k);
		QueryRead* queryRead = queryDataSet->getReadFromID(currentID);
		string querySequence = queryRead->getSequence();
		string queryName = queryRead->getName();
		bool alignFlag = true;
		if(queryName==subjectAlign.subjectReadName)alignFlag = false; //only align when queryName is alphabetical smaller than subject, removing half of the pair-wise alignment redundancy
		if(alignFlag)
		{
			Alignment* align = new Alignment();
			align->subjectReadName = subjectAlign.subjectReadName;
			align->subjectReadSequence = subjectRead;
			align->queryRead = queryRead;
			if(createAlignment(align, j, j+this->hashKeyLength,j+this->hashKeyLength+this->hashKeyLength-1, mode, 2)==false) delete align;
			else subjectAlign.queryAlignmentList.push_back(align);

		}
	}

	Key1OnlyList->clear();
	Key2OnlyList->clear();
	BothKeyList->clear();
	delete Key1OnlyList;
	delete Key2OnlyList;
	delete BothKeyList;
	}//loop for window j

	}//loop of i for each mode
	return true;
}

//This function is designed for dynamic key range.
//keymatchmode = 1(only left key matched),2(only right key matched),3(both matched)
bool DoubleKeyHashTable::createAlignment(Alignment* subjectAlignment, int leftKeyPoint, int rightKeyPoint, int keyEndPoint, string mode, int keymatchmode)
{
	string subjectString=subjectAlignment->subjectReadSequence; // Get the forward of read1
	string queryString="";
	if(mode=="forwardprefix" || mode=="forwardsuffix")
	{
		queryString = subjectAlignment->queryRead->getSequence();
		subjectAlignment->queryOrientation = true;
	}
	else if(mode=="reverseprefix" || mode=="reversesuffix")
	{
		queryString = subjectAlignment->queryRead->reverseComplement();
		subjectAlignment->queryOrientation = false;
	}
	else return false;

	int leftKeyLength = rightKeyPoint-leftKeyPoint;
	int rightKeyLength = keyEndPoint-rightKeyPoint+1;

	bool isContainedAlign;
	int maxCompareLength;
	int offset;
	int seedStart, seedEnd;
	string subSubject, subQuery;
	if(mode == "forwardprefix" || mode=="reverseprefix")
	{
		subjectAlignment->subjectStart = 0-leftKeyPoint;
		subjectAlignment->subjectEnd = subjectAlignment->subjectStart + subjectAlignment->subjectReadSequence.length()-1;
		subjectAlignment->queryEnd = subjectAlignment->queryRead->getReadLength()-1;
		isContainedAlign = subjectAlignment->subjectReadSequence.length()- leftKeyPoint  >= subjectAlignment->queryRead->getReadLength();
		if(!isContainedAlign)
		{
			maxCompareLength = subjectAlignment->subjectReadSequence.length()- leftKeyPoint;
			offset = 0;
		}
		else
		{
			maxCompareLength = subjectAlignment->queryRead->getReadLength();
			offset = 0;
		}
		subSubject = subjectString.substr(0-subjectAlignment->subjectStart,maxCompareLength);
		subQuery = queryString.substr(0,maxCompareLength);
		//mark the perfect matched region
		if(keymatchmode ==1)
		{
			seedStart = 0;
			seedEnd = leftKeyLength-1;
		}
		else if(keymatchmode == 2)
		{
			seedStart = leftKeyLength;
			seedEnd = leftKeyLength+rightKeyLength-1;
		}
		else if(keymatchmode == 3)
		{
			seedStart = 0;
			seedEnd = leftKeyLength+rightKeyLength-1;
		}
		else return false;

	}
	else if (mode == "forwardsuffix" || mode=="reversesuffix")
	{
		subjectAlignment->subjectStart = subjectAlignment->queryRead->getReadLength()-1-keyEndPoint;
		subjectAlignment->subjectEnd = subjectAlignment->subjectStart + subjectAlignment->subjectReadSequence.length()-1;
		subjectAlignment->queryEnd = subjectAlignment->queryRead->getReadLength()-1;
		isContainedAlign = subjectAlignment->queryRead->getReadLength()-(leftKeyLength+rightKeyLength) <= leftKeyPoint;
		if(!isContainedAlign)
		{
			maxCompareLength = keyEndPoint + 1;
			offset = subjectAlignment->queryRead->getReadLength()-1-keyEndPoint;
		}
		else
		{
			maxCompareLength = subjectAlignment->queryRead->getReadLength();
			offset = 0;
		}
		subSubject = subjectString.substr(queryString.length()-maxCompareLength-subjectAlignment->subjectStart,maxCompareLength);
		subQuery = queryString.substr(queryString.length()-maxCompareLength,maxCompareLength);
		//int a = subSubject.length();
		int subQueryLen = subQuery.length();
		if(keymatchmode ==1)
		{
			seedStart = subQueryLen-rightKeyLength-leftKeyLength;
			seedEnd = subQueryLen-rightKeyLength - 1;
		}
		else if(keymatchmode == 2)
		{
			seedStart = subQueryLen-rightKeyLength;
			seedEnd = subQueryLen-1;
		}
		else if(keymatchmode == 3)
		{
			seedStart = subQueryLen-rightKeyLength-leftKeyLength;
			seedEnd = subQueryLen-1;

		}
		else return false;
	}
	else return false;

	return alignSubStrings(subjectAlignment, subSubject, subQuery, offset, seedStart, seedEnd);
}

bool DoubleKeyHashTable::alignSubStrings(Alignment* subjectAlignment, string& subSubject, string& subQuery, int offset, int seedStart, int seedEnd)
{
	int currentMismatchCount=0;
	int i = 0;
	while(i<subSubject.size())
	{
		if(i==seedStart)
		{
			i = seedEnd;
		}
		else
		{
			if(subSubject.at(i)!=subQuery.at(i))
			{
				currentMismatchCount++;
				if(currentMismatchCount>maxMismatch)return false;
				subjectAlignment->editInfor.insert(std::pair<int, char>(offset+i, subSubject.at(i)));
			}

		}
		i++;
	}
	return true;
}

bool DoubleKeyHashTable::wennDiagramTwoLists(vector<UINT64>* list1, vector<UINT64>* list2, vector<UINT64>* list1only, vector<UINT64>* list2only, vector<UINT64>* list12)//sorted from smaller to larger ID in the list
{
	unsigned int i,j;
	i=0;j=0;
	while(list1!=NULL&&list2!=NULL&& i<list1->size()&&j<list2->size())
	{
		if(list1->at(i)<list2->at(j))
		{
			UINT64 smallvalue = list1->at(i);
			list1only->push_back(smallvalue);
			i++;
		}
		else if(list1->at(i)>list2->at(j))
		{
			UINT64 smallvalue = list2->at(j);
			list2only->push_back(smallvalue);
			j++;
		}
		else //same value in both list
		{
			UINT64 value = list1->at(i);
			list12->push_back(value);
			i++;j++;
		}
	}

	while(list1!=NULL && i<list1->size())
	{
		UINT64 smallvalue = list1->at(i);
		list1only->push_back(smallvalue);
		i++;
	}

	while(list2!=NULL && j<list2->size())
	{
		UINT64 smallvalue = list2->at(j);
		list2only->push_back(smallvalue);
		j++;
	}
	return true;
}
