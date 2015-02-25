/*
 * OmegaGraphConstructor.h
 *
 *  Created on: Feb 16, 2015
 *      Author: qy2
 */

#ifndef OMEGAGRAPHCONSTRUCTOR_H_
#define OMEGAGRAPHCONSTRUCTOR_H_

#include "Config.h"
#include "OmegaHashTable.h"
#include "SubjectDataset.h"
class OmegaGraphConstructor {
	OmegaHashTable * omegaHashTable;
	UINT64 totaledgenumber;
	ofstream filePointer;
public:
	OmegaGraphConstructor(OmegaHashTable * omegaHashTable);
	~OmegaGraphConstructor();
	bool start();
	bool searchHashTable(SubjectEdge * subjectEdge);
	bool checkOverlap(SubjectRead *read1, QueryRead *read2, UINT64 orient, UINT64 start);
	bool checkOverlap(QueryRead *read1, QueryRead *read2, UINT64 orient, UINT64 start);
	bool insertAllEdgesOfRead(UINT64 readNumber, vector<int> * exploredReads);
	bool printEdgesToFile(string outFileName);
};

#endif /* OMEGAGRAPHCONSTRUCTOR_H_ */
