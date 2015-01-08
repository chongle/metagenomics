/*
 * ErrorCorrector.cpp
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#include "ErrorCorrector.h"

ErrorCorrector::ErrorCorrector() {
	// TODO Auto-generated constructor stub

}

ErrorCorrector::~ErrorCorrector() {


}

bool ErrorCorrector::start() {

	if (!queryDataset->buildDataset(Config::getQueryDatasetFilename())){
		cout << "Error: cannot build query dataset: " << Config::getQueryDatasetFilename() << endl;
		return false;
	}

	if (!hashTable->insertDataset(queryDataset)){
		cout << "Error: cannot build hash table" << endl;
		return false;
	}

	SubjectDataset * subject = new SubjectDataset;
	if(!subject->setFilename(Config::getSubjectDatasetFilenames())){
		cout << "Error: cannot load subject dataset" << endl;
		return false;
	}

	vector<SubjectAlignment*> subjectList;

	while(subject->loadNextChunk(subjectList)){
		// #pragma omp parallel for
		for(int i = 0; i < subjectList.size(); i++){
				searchHashTable(subjectList[i])

		}
		// end of omp parallel

		for(int i = 0; i < subjectList.size(); i++){
			for(int j = 0; j < subjectList[i].queryAlignment.size(); j++){
				subjectList[i].queryAlignment[j].queryRead->addAlignment( subjectList[i].queryAlignment[j] );
			}
		}
		subjectList.clear();
	}

	// #pragma omp parallel for
	for(int i = 0; i < queryDataset->queryReadList.size(); i++){
		queryDataset->queryReadList[i]->correctErrors();

	}
	// end of omp parallel

	queryDataset->write2File();

	delete subject;

}


