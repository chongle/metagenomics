/*
 * ErrorCorrector.h
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#ifndef ERRORCORRECTOR_H_
#define ERRORCORRECTOR_H_

#include "Config.h"


class SubjectAlignment {
	// alignment input
	string subjectReadSequence;
	string subjectReadName;

	// alignment results
	vector<Alignment> queryAlignment;

};

class ErrorCorrector {

	QueryDataset * queryDataset;
	HashTable * hashTable;

public:
	ErrorCorrector();
	~ErrorCorrector();
	bool start();
	bool searchHashTable(SubjectAlignment & subject);

};

#endif /* ERRORCORRECTOR_H_ */
