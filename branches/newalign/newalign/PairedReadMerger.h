/*
 * PairedReadMerger.h
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#ifndef PAIREDREADMERGER_H_
#define PAIREDREADMERGER_H_


class AlignmentPairedEnd {

public:

	QueryReadPairedEnd * queryPairedEndRead;
	string subjectReadSequence;

	// results
	int orientation;
	int start;
	int stop;
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
};

#endif /* PAIREDREADMERGER_H_ */
