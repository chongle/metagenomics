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
	map<string, HashTable*> hashTableMap;
	UINT16 minimumOverlapLength;
	UINT16 hashKeyLength;
public:
	SingleKeyHashTable();
	virtual ~SingleKeyHashTable();
	bool InitializeAllHashTables();
};

#endif /* SINGLEKEYHASHTABLE_H_ */
