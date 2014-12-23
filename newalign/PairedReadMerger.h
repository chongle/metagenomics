/*
 * PairedReadMerger.h
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#ifndef PAIREDREADMERGER_H_
#define PAIREDREADMERGER_H_

#include "Config.h"

class AlignmentPairedEnd {

public:

	QueryReadPairedEnd * queryPairedEndRead;

	// this is in the same orientation as the query forward
	string gapSequence;

	// orientation of the query read
	bool queryOrientation;		// true for forward, false for reverse

	//  The following case will have a positive subject start position
	//  query:        XXXXXXMMMMMMMM
	//	subject:            MMMMMMMMXXXXXXX
	//  The following case will have a negative subject start position
	//  query:        MMMMMMXXXXXXXX
	//	subject: XXXXXMMMMMM
	int subjectStart;
	int leftQueryEnd;
	int rightQueryStart;
	int subjectEnd;
	int rightQueryEnd;

	// coordinates are defined by the query reads
	// all insertion and deletion is on the subject read
	// lower case char means substitution
	// upper case char means insertion
	// 'D' means deletion
	//
	map<int, char> editInfor;
};

class SubjectAlignmentPairedEnd {
	// alignment input
	string subjectReadSequence;
	string subjectReadName;

	// alignment results
	vector<AlignmentPairedEnd> queryAlignmentList;

};

class PairedReadMerger {

	QueryDataset * queryDataset;
	HashTable * hashTable;
	SubjectDataset * subject;

//	search(SubjectAlignmentPairedEnd);

public:
	PairedReadMerger();
	~PairedReadMerger();

	bool start();
	bool searchHashTable(SubjectAlignmentPairedEnd & subjectAlignment);

};

#endif /* PAIREDREADMERGER_H_ */
