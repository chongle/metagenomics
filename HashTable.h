/*
 * HashTable.h
 *
 *  Created on: Jan 6, 2015
 *      Author: qy2
 */

#ifndef HASHTABLE_H_
#define HASHTABLE_H_

#include "Config.h"

class HashTable {

	UINT64 hashTableSize; 						// Size of hash table, primer number usually.
	UINT16 hashKeyLength;					// Length of the substring of the reads to hash. It's also the key length.
												// It has the 64 bps limit. See the comments in the hashFunction().
	UINT64 numberOfHashCollision;				// Counted total number of hash collisions. For debugging only.
	UINT64 maxSingleHashCollision;				// Counted maximal number of hash collision for a single case.
	vector <vector<UINT64> *> *hashTable; 		// Main structure of hashtable, storing read identifiers.

	UINT64 hashFunction(const string & subString);
	UINT64 getPrimeLargerThanNumber(UINT64 number);
	void setHashTableSizeAndInitialize(UINT64 size);
	//UINT64 getIndexFromHashTable(string String);
	//bool getLinkedListFromHashTable();
	//bool getHashKeyFromIndex();

public:
	HashTable();
	HashTable(UINT16 keylength,  UINT64 size);
	virtual ~HashTable();

	//bool insertIntoHashTable();



};

#endif /* HASHTABLE_H_ */
