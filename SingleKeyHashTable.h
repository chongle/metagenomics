/*
 * SingleKeyHashTable.h
 *
 *  Created on: Jan 6, 2015
 *      Author: qy2
 */

#ifndef SINGLEKEYHASHTABLE_H_
#define SINGLEKEYHASHTABLE_H_

#include "Config.h"
#include "HashTable.h"
#include "QueryDataset.h"
#include "Alignment.h"
//#include "AlignmentPairedEnd.h"

class SingleKeyHashTable {

	// 00 = 0 means prefix of the forward string.
	// 01 = 1 means suffix of the forward string.
	// 10 = 2 means prefix of the reverse string.
	// 11 = 3 means suffix of the reverse string.
	vector<string> hashTableNameList;	//forwardprefix, forwardsuffix, reverseprefix, reversesuffix
	map<string, HashTable*> hashTableMap;
	UINT16 minimumOverlapLength;
	UINT16 hashKeyLength;
	int maxMismatch;
	QueryDataset * queryDataSet;



	string getReadSubstring(string mode, UINT64 readID);// mode ={forwardprefix, forwardsuffix, reverseprefix, reversesuffix}
//	bool doAlignment(Alignment* align, string mode, int subjectStart);
//	bool checkForContainedAlignment(Alignment* align, string mode, int subjectStart);
//	bool subjectWindowRange(int& startpoint, int& stoppoint, string mode, string& subjectRead);

public:
	SingleKeyHashTable();
	virtual ~SingleKeyHashTable();
	bool InitializeAllHashTables();
	bool insertQueryDataset(QueryDataset* d);
	bool insertQueryRead(QueryRead *read, string mode);// mode ={forwardprefix, forwardsuffix, reverseprefix, reversesuffix}
//	bool singleKeySearch(edge & Edge);
//	bool singleKeySearch(SubjectAlignment & subjectAlign);
//	bool singleKeySearch(SubjectAlignmentPairedEnd & subjectAlignment);
};

#endif /* SINGLEKEYHASHTABLE_H_ */
