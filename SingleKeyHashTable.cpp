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
	this->maxMismatch = Config::maxMismatch;
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
	HashTable * currentHashTable = hashTableMap.at(mode);
	return currentHashTable->insertIntoHashTable(keystring,readID);

}

bool SingleKeyHashTable::doAlignment(Alignment* align, string mode, int subjectStart)
{

	if(mode=="forwardprefix")
	{
		align->queryOrientation = true;
		if(align->subjectReadSequence.length()- subjectStart - hashKeyLength >= align->queryRead->getSequence().length() - hashKeyLength) // The overlap must continue till the end.
			checkForContainedAlignment(align, mode, subjectStart);
		string restSubject = align->subjectReadSequence.substr(subjectStart + hashKeyLength, align->subjectReadSequence.length()-(subjectStart + hashKeyLength));
		string restQuery = align->queryRead->getSequence().substr(hashKeyLength,  restSubject.length());
		int currentMismatchCount=0;
		for(int i=0; i<restQuery.length();i++)
		{
			if(restQuery.at(i)!=restSubject.at(i))
			{
				currentMismatchCount++;
				if(currentMismatchCount>maxMismatch)return false;
				align->editInfor.insert(std::pair<int, char>(hashKeyLength+i, restSubject.at(i)));
			}
		}
		align->subjectStart = -subjectStart;
		align->subjectEnd = align->subjectStart + align->subjectReadSequence.length()-1;
		align->queryEnd = align->queryRead->getReadLength()-1;
	}
	else if(mode == "forwardsuffix")
	{
		align->queryOrientation = true;
		if(align->queryRead->getSequence().length()-hashKeyLength <= subjectStart)
			checkForContainedAlignment(align, mode, subjectStart);
		string restSubject = align->subjectReadSequence.substr(0, subjectStart);
		string restQuery = align->queryRead->getSequence().substr(align->queryRead->getReadLength()-hashKeyLength-subjectStart, restSubject.length());
		int currentMismatchCount=0;
		for(int i=0; i<restQuery.length();i++)
		{
			if(restQuery.at(i)!=restSubject.at(i))
			{
				currentMismatchCount++;
				if(currentMismatchCount>maxMismatch)return false;
				align->editInfor.insert(std::pair<int, char>(align->queryRead->getReadLength()-hashKeyLength-subjectStart+i, restSubject.at(i)));
			}
		}
		align->subjectStart = align->queryRead->getReadLength()-hashKeyLength-subjectStart;
		align->subjectEnd = align->subjectStart + align->subjectReadSequence.length()-1;
		align->queryEnd = align->queryRead->getReadLength()-1;
	}
	else if(mode == "reverseprefix")
	{
		align->queryOrientation = false;
		if(align->subjectReadSequence.length()- subjectStart - hashKeyLength >= align->queryRead->getSequence().length() - hashKeyLength) // The overlap must continue till the end.
			checkForContainedAlignment(align, mode, subjectStart);
		string restSubject = align->subjectReadSequence.substr(subjectStart + hashKeyLength, align->subjectReadSequence.length()-(subjectStart + hashKeyLength));
		string restQuery = align->queryRead->reverseComplement().substr(hashKeyLength,  restSubject.length());
		int currentMismatchCount=0;
		for(int i=0; i<restQuery.length();i++)
		{
			if(restQuery.at(i)!=restSubject.at(i))
			{
				currentMismatchCount++;
				if(currentMismatchCount>maxMismatch)return false;
				align->editInfor.insert(std::pair<int, char>(hashKeyLength+i, restSubject.at(i)));
			}
		}
		align->subjectStart = -subjectStart;
		align->subjectEnd = align->subjectStart + align->subjectReadSequence.length()-1;
		align->queryEnd = align->queryRead->getReadLength()-1;
	}
	else if(mode == "reversesuffix")
	{
		align->queryOrientation = false;
		if(align->queryRead->getSequence().length()-hashKeyLength <= subjectStart)
			checkForContainedAlignment(align, mode, subjectStart);
		string restSubject = align->subjectReadSequence.substr(0, subjectStart);
		string restQuery = align->queryRead->reverseComplement().substr(align->queryRead->getReadLength()-hashKeyLength-subjectStart, restSubject.length());
		int currentMismatchCount=0;
		for(int i=0; i<restQuery.length();i++)
		{
			if(restQuery.at(i)!=restSubject.at(i))
			{
				currentMismatchCount++;
				if(currentMismatchCount>maxMismatch)return false;
				align->editInfor.insert(std::pair<int, char>(align->queryRead->getReadLength()-hashKeyLength-subjectStart+i, restSubject.at(i)));
			}
		}
		align->subjectStart = align->queryRead->getReadLength()-hashKeyLength-subjectStart;
		align->subjectEnd = align->subjectStart + align->subjectReadSequence.length()-1;
		align->queryEnd = align->queryRead->getReadLength()-1;
	}
	else return false;

	return true;
}

//the choice of the start and stop position should meet the minimum overlap requirement.
bool SingleKeyHashTable::subjectWindowRange(int& startpoint, int& stoppoint, string mode, string& subjectRead)
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

bool SingleKeyHashTable::checkForContainedAlignment(Alignment* align, string mode, int subjectStart)
{
	string subjectString=align->subjectReadSequence; // Get the forward of read1
	string queryString="";
//	string queryString = (mode=="forwardprefix" || mode=="forwardsuffix") ? align->queryRead->getSequence() : align->queryRead->reverseComplement(); // Get the string in read2 based on the orientation.
	if(mode=="forwardprefix" || mode=="forwardsuffix")
	{
		queryString = align->queryRead->getSequence();
		align->queryOrientation = true;
	}
	else if(mode=="reverseprefix" || mode=="reversesuffix")
	{
		queryString = align->queryRead->reverseComplement();
		align->queryOrientation = false;
	}
	else return false;

	if(mode=="forwardprefix" || mode=="reverseprefix")
									// mode = forwardprefix
									//   >--------MMMMMMMMMMMMMMM*******------> subject read1      M means match found by hash table
									//            MMMMMMMMMMMMMMM*******>       query read2      * means we need to check these characters for match
									//				OR
									// mode = reverseprefix
									//	 >---*****MMMMMMMMMMMMMMM*******------> subject read1
									//		      MMMMMMMMMMMMMMM*******<	    Reverese complement of query read2
	{
		int restSubjectLength = subjectString.length() - subjectStart - this->hashKeyLength; 	// This is the remaining of read1
		int restQueryLength = queryString.length() - this->hashKeyLength; 	// This is the remaining of read2
		if(restSubjectLength >= restQueryLength)
		{
			string restSubject = subjectString.substr(subjectStart + hashKeyLength, restQueryLength);
			string restQuery = queryString.substr(hashKeyLength,  restQueryLength);
			int currentMismatchCount=0;
			for(int i=0; i<restQuery.length();i++)
			{
				if(restQuery.at(i)!=restSubject.at(i))
				{
					currentMismatchCount++;
					if(currentMismatchCount>maxMismatch)return false;
					align->editInfor.insert(std::pair<int, char>(hashKeyLength+i, restSubject.at(i)));
				}
			}
			align->subjectStart = -subjectStart;
			align->subjectEnd = align->subjectStart + align->subjectReadSequence.length()-1;
			align->queryEnd = align->queryRead->getReadLength()-1;

		}
	}
	else if(mode=="forwardsuffix" || mode=="reversesuffix")
									// mode = forwardsuffix
									//   >---*****MMMMMMMMMMMMMMM-------------> subject read1      M means match found by hash table
									//      >*****MMMMMMMMMMMMMMM       		query read2      * means we need to check these characters for match
									//				OR
									// mode = reversesuffix
									//	 >---*****MMMMMMMMMMMMMMM-------------> subject read1
									//		<*****MMMMMMMMMMMMMMM				Reverse Complement of query Read2
	{
		int restSubjectLength = subjectStart;
		int restQueryLength = queryString.length() - this->hashKeyLength;
		if(restSubjectLength >= restQueryLength)
		{
			string restSubject = subjectString.substr(subjectStart-restQueryLength, restQueryLength);
			string restQuery = queryString.substr(0, restQueryLength);
			int currentMismatchCount=0;
			for(int i=0; i<restQuery.length();i++)
			{
				if(restQuery.at(i)!=restSubject.at(i))
				{
					currentMismatchCount++;
					if(currentMismatchCount>maxMismatch)return false;
					align->editInfor.insert(std::pair<int, char>(i, restSubject.at(i)));
				}
			}
			align->subjectStart = align->queryRead->getReadLength()-hashKeyLength-subjectStart;
			align->subjectEnd = align->subjectStart + align->subjectReadSequence.length()-1;
			align->queryEnd = align->queryRead->getReadLength()-1;
		}
	}
	else return false;

	return true;

}

bool SingleKeyHashTable::singleKeySearch(edge & Edge)
{

	string subjectRead = Edge.subjectReadSequence;
	int startpoint = this->minimumOverlapLength-this->hashKeyLength;
	int stoppoint = subjectRead.length()-this->minimumOverlapLength;

	for(int i =0; i<this->hashTableNameList.size();i++)
	{


	string modestring = this->hashTableNameList.at(i);

	int startpoint,stoppoint;
	if(subjectWindowRange(startpoint, stoppoint, modestring, subjectRead)) return false;//guarantee it meets the minimum overlaplength requirement for alignments.

	for(int j=startpoint;j<=stoppoint;j++)
	{
	string subString = subjectRead.substr(startpoint, hashKeyLength);
	vector<UINT64> currentIDList = hashTableMap.at(modestring)->getReadIDListOfReads(subString);
	for(int k=0;k<currentIDList.size();k++)
	{
		UINT64 currentID = currentIDList.at(k);
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
			if(doAlignment(align, modestring, startpoint)==false) delete align;
			else Edge.alignmentList.push_back(align);

		}
	}
	}

	}

}
