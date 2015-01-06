/*
 * Alignment.h
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

#include "Config.h"


class Alignment {
public:
	Alignment();
	virtual ~Alignment();

	QueryRead * queryRead;

	// Don't need these two variables for error correction to save memory.
	string subjectReadSequence = "";
	string subjectReadName =  "";

	// orientation of the query read
	bool queryOrientation;		// true for forward, false for reverse

	// coordinates of the overlap alignment, which is defined by the query reads
	//  The following case will have a positive subject start position
	//  query:        XXXXXXMMMMMMMM
	//	subject:            MMMMMMMMXXXXXXX
	//  The following case will have a negative subject start position
	//  query:        MMMMMMXXXXXXXX
	//	subject: XXXXXMMMMMM
	int subjectStart;
	int queryEnd;
	int subjectEnd;

	// coordinates are defined by the query reads
	// all insertion and deletion is on the subject read
	// lower case char means substitution
	// upper case char means insertion
	// 'D' means deletion
	//
	map<int, char> editInfor;
};

#endif /* ALIGNMENT_H_ */
