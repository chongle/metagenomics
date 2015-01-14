/*
 * OverlapGraphConstructor.h
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#ifndef OVERLAPGRAPHCONSTRUCTOR_H_
#define OVERLAPGRAPHCONSTRUCTOR_H_

#include "Config.h"

class edge{
public:
	// alignment input
	string subjectReadSequence;
	string subjectReadName;



	// alignment results
	vector<Alignment*> alignmentList;

	// a list of query reads that are contained by the subject read or duplicate with the subject read
	// when a query is duplicate with a subject read, we compare their names and add the query read to the duplicate read list if its name is greater than the subject read name.
	vector<QueryRead *> containedOrDeplicateReadList;//It's all for contained read if the duplicated reads are removed from the beginning.


};

class OverlapGraphConstructor {
	QueryDataset * queryDataset;
	SingleKeyHashTable * singleKeyHashTable;
	DoubleKeyHashTable * doubleKeyHashTable;
	SubjectDataset * subject;
	bool isContainedAlignment(Alignment * subjectAlignment);
public:
	OverlapGraphConstructor();
	~OverlapGraphConstructor();
	bool start();
	bool searchHashTable(edge & currentEdge);

};

#endif /* OVERLAPGRAPHCONSTRUCTOR_H_ */
