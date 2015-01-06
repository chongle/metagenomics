/*
 * QueryDataset.h
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#ifndef QUERYDATASET_H_
#define QUERYDATASET_H_

#include "Config.h"

class QueryDataset {
public:
	QueryDataset();
	virtual ~QueryDataset();
	vector<QueryRead *> queryReadList;
	bool buildDataset(const string & sQueryFilename);

};

#endif /* QUERYDATASET_H_ */
