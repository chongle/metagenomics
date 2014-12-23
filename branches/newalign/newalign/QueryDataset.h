/*
 * QueryDataset.h
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#ifndef QUERYDATASET_H_
#define QUERYDATASET_H_

class QueryDataset {
public:
	QueryDataset();
	virtual ~QueryDataset();
	vector<QueryRead *> queryReadList;

};

#endif /* QUERYDATASET_H_ */
