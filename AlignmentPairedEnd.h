/*
 * AlignmentPairedEnd.h
 *
 *  Created on: Jan 6, 2015
 *      Author: qy2
 */

#ifndef ALIGNMENTPAIREDEND_H_
#define ALIGNMENTPAIREDEND_H_

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
	// upper case char means substitution
	// lower case char means insertion
	// 'D' means deletion
	//
	map<int, char> editInfor;

	AlignmentPairedEnd();
    virtual ~AlignmentPairedEnd();
};
#endif /* ALIGNMENTPAIREDEND_H_ */
