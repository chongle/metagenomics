/*
 * DoubleKeyHashTable.h
 *
 *  Created on: Jan 9, 2015
 *      Author: qy2
 */

#ifndef DOUBLEKEYHASHTABLE_H_
#define DOUBLEKEYHASHTABLE_H_

#include "Config.h"

class DoubleKeyHashTable {

	//1 is the key in the front, 2 is the following; or 1 is on the left side, 2 is on the right side
    vector<string> hashTableNameList;	//forwardprefix1,forwardprefix2, forwardsuffix1,forwardsuffix2, reverseprefix1, reverseprefix2, reversesuffix1, reversesuffix2
	map<string, HashTable*> hashTableMap;
	UINT16 minimumOverlapLength;
	UINT16 hashKeyLength;
	QueryDataset * queryDataSet = NULL;



	string getReadSubstring(string mode, UINT64 readID);// mode ={forwardprefix1,forwardprefix2, forwardsuffix1,forwardsuffix2, reverseprefix1, reverseprefix2, reversesuffix1, reversesuffix2}
public:
	DoubleKeyHashTable();
	virtual ~DoubleKeyHashTable();
	bool InitializeAllHashTables();
	bool insertQueryDataset(QueryDataset* d);
	bool insertQueryRead(QueryRead *read, string mode);// mode ={forwardprefix1,forwardprefix2, forwardsuffix1,forwardsuffix2, reverseprefix1, reverseprefix2, reversesuffix1, reversesuffix2}
	bool doubleKeySearch(edge & Edge);
	bool doubleKeySearch(SubjectAlignment & subject);
	bool doubleKeySearch(SubjectAlignmentPairedEnd & subjectAlignment);

};

#endif /* DOUBLEKEYHASHTABLE_H_ */
