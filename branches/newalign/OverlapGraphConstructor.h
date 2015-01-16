/*
 * OverlapGraphConstructor.h
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#ifndef OVERLAPGRAPHCONSTRUCTOR_H_
#define OVERLAPGRAPHCONSTRUCTOR_H_

#include "Config.h"
#include "SingleKeyHashTable.h"
#include "DoubleKeyHashTable.h"
#include "QueryDataset.h"
#include "SubjectDataset.h"
#include "Alignment.h"


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
	void printEdgesToFile(bool nonRemovedReads, string outFileName);

};

#endif /* OVERLAPGRAPHCONSTRUCTOR_H_ */
