/*
 * ErrorCorrector.h
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#ifndef ERRORCORRECTOR_H_
#define ERRORCORRECTOR_H_

#include "Config.h"
#include "QueryDataset.h"
#include "Alignment.h"



class ErrorCorrector {

	QueryDataset * queryDataset;
//	HashTable * hashTable;

public:
	ErrorCorrector();
	~ErrorCorrector();
	bool start();
	bool searchHashTable(SubjectAlignment & subject);

};

#endif /* ERRORCORRECTOR_H_ */
