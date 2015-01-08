/*
 * SingleKeyHashTable.h
 *
 *  Created on: Jan 6, 2015
 *      Author: qy2
 */

#ifndef SINGLEKEYHASHTABLE_H_
#define SINGLEKEYHASHTABLE_H_

#include "Config.h"

class SingleKeyHashTable {

	// 00 = 0 means prefix of the forward string.
	// 01 = 1 means suffix of the forward string.
	// 10 = 2 means prefix of the reverse string.
	// 11 = 3 means suffix of the reverse string.
	vector<string> hashTableNameList;	//forwardprefix, forwardsuffix, reverseprefix, reversesuffix
	map<string, HashTable*> hashTableMap;
	UINT16 minimumOverlapLength;
	UINT16 hashKeyLength;
	QueryDataset * queryDataSet = NULL;



	string getReadSubstring(string mode, UINT64 readID);// mode ={forwardprefix, forwardsuffix, reverseprefix, reversesuffix}
public:
	SingleKeyHashTable();
	virtual ~SingleKeyHashTable();
	bool InitializeAllHashTables();
	bool insertQueryDataset(QueryDataset* d);
	bool insertQueryRead(QueryRead *read, string mode);// mode ={forwardprefix, forwardsuffix, reverseprefix, reversesuffix}
};

#endif /* SINGLEKEYHASHTABLE_H_ */
