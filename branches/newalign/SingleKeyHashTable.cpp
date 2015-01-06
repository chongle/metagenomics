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
}

SingleKeyHashTable::~SingleKeyHashTable() {
	// TODO Auto-generated destructor stub
}

bool InitializeAllHashTables()
