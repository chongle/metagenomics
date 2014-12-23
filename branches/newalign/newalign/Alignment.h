/*
 * Alignment.h
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

include "Contig"


class Alignment {
public:
	Alignment();
	virtual ~Alignment();

	QueryRead * queryRead;
	string subjectReadSequence;
	string subjectReadName;
	int orientation;
	int start;
	int stop;
};

#endif /* ALIGNMENT_H_ */
