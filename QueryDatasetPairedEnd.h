/*
 * QueryDatasetPairedEnd.h
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#ifndef QUERYDATASETPAIREDEND_H_
#define QUERYDATASETPAIREDEND_H_

#include "Config.h"

class QueryDatasetPairedEnd {
public:
	QueryDatasetPairedEnd();
	virtual ~QueryDatasetPairedEnd();
	vector<QueryReadPairedEnd *> queryReadList;

};

#endif /* QUERYDATASETPAIREDEND_H_ */
