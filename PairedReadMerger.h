/*
 * PairedReadMerger.h
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#ifndef PAIREDREADMERGER_H_
#define PAIREDREADMERGER_H_

#include "Config.h"
#include "QueryDataset.h"
#include "SubjectDataset.h"
#include "AlignmentPairedEnd.h"




class PairedReadMerger {

	QueryDataset * queryDataset;
//	HashTable * hashTable;
	SubjectDataset * subject;

//	search(SubjectAlignmentPairedEnd);

public:
	PairedReadMerger();
	~PairedReadMerger();

	bool start();
	bool searchHashTable(SubjectAlignmentPairedEnd & subjectAlignment);

};

#endif /* PAIREDREADMERGER_H_ */
