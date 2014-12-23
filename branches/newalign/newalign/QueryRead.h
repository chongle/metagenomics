/*
 * QueryRead.h
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#ifndef QUERYREAD_H_
#define QUERYREAD_H_

class QueryRead {
	vector<Alignment> queryAlignmentList;
	string readSequence;
	string readName;
public:
	QueryRead();
	virtual ~QueryRead();
	bool addAlignment(Alignment subjectAlignment);
	bool correctErrors();
};

#endif /* QUERYREAD_H_ */
