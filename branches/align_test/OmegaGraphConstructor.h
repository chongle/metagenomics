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
public:
	OmegaGraphConstructor(OmegaHashTable * omegaHashTable);
	~OmegaGraphConstructor();
	bool start();
	bool searchHashTable(SubjectEdge * subjectEdge);
};

#endif /* OMEGAGRAPHCONSTRUCTOR_H_ */
