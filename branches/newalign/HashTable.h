/*
 * HashTable.h
 *
 *  Created on: Jan 6, 2015
 *      Author: qy2
 */

#ifndef HASHTABLE_H_
#define HASHTABLE_H_

#include "Config.h"

class DataVector{

public:
	vector<UINT64> * dataList;
	string keystring;
DataVector(){dataList = NULL;keystring = "";};
~DataVector();
};

class HashTable {
	string hashTableName;

	UINT64 hashTableSize; 						// Size of hash table, primer number usually.
	UINT16 hashKeyLength;					// Length of the substring of the reads to hash. It's also the key length.
												// It has the 64 bps limit. See the comments in the hashFunction().
	UINT64 numberOfHashCollision;				// Counted total number of hash collisions. For debugging only.
	UINT64 maxSingleHashCollision;				// Counted maximal number of hash collision for a single case.
	vector <DataVector *> *hashTable; 		// Main structure of hashtable, storing read identifiers.

	UINT64 hashFunction(const string & subString);
	UINT64 getPrimeLargerThanNumber(UINT64 number);
	void setHashTableSizeAndInitialize(UINT64 size);
	bool findInsertIndex(string keyString, UINT64& hashTableIndex);
	bool getIndexFromHashTable(string subString, UINT64& hashTableIndex);
	string getHashKeyFromIndex(UINT64 hashTableIndex);

public:
	HashTable();
	HashTable(UINT16 keylength);
	HashTable(UINT16 keylength, string name);
	HashTable(UINT16 keylength,  UINT64 size);
	virtual ~HashTable();
	void InitializeWithDataSize(UINT64 dataSetSize);//hash table will automatically determine a good hash table size
	bool insertIntoHashTable(string keyString, UINT64 readID);//handle collision internally
	vector<UINT64> * getReadIDListOfReads(string subString);



};

#endif /* HASHTABLE_H_ */
