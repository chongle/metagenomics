/*
 * QueryReadPairedEnd.h
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#ifndef QUERYREADPAIREDEND_H_
#define QUERYREADPAIREDEND_H_

#include "Config.h"

class QueryReadPairedEnd {
public:
	QueryReadPairedEnd();
	virtual ~QueryReadPairedEnd();
	string readLeft;
	string readRight;
	string mergedRead;
	bool merge();
};

#endif /* QUERYREADPAIREDEND_H_ */
