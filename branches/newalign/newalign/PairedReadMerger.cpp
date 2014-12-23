/*
 * PairedReadMerger.cpp
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#include "PairedReadMerger.h"

PairedReadMerger::PairedReadMerger() {
	// TODO Auto-generated constructor stub

}

PairedReadMerger::~PairedReadMerger() {
	// TODO Auto-generated destructor stub
}

bool PairedReadMerger::start() {
	if (!queryDatasetPairedEnd->buildDataset(Config::getQueryDatasetFilename())){
		cout << "Error: cannot build query dataset" << endl;
		return false;
	}
	if (!hashTable->insertDatasetPairedEnd(queryDataset)){
		cout << "Error: cannot build hash table" << endl;
		return false;
	}

	vector<SubjectAlignmentPairedEnd> subjectList;

	while(subject->loadNextChunk(subjectList)){
		// #pragma omp parallel for
		for(int i = 0; i < subjectList.size(); i++){
				searchHashTable(subjectList[i]);

		}
		// end of omp parallel

		for(int i = 0; i < subjectList.size(); i++){
			for(int j = 0; j < subjectList[i].queryAlignmentList.size(); j++){
				subjectList[i].queryAlignmentList[j].queryRead->addAlignment( subjectList[i].queryAlignmentList[j] );
			}
		}
		subjectList.clear();
	}

	// #pragma omp parallel for
	for(int i = 0; i < queryDataset->queryReadList.size(); i++){
		queryDataset->queryReadList[i]->merge();

	}
	// end of omp parallel

	queryDataset->write2File();

}
