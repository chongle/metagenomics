/*
 * OverlapGraphConstructor.cpp
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#include "OverlapGraphConstructor.h"

OverlapGraphConstructor::OverlapGraphConstructor() {
	// TODO Auto-generated constructor stub

}

OverlapGraphConstructor::~OverlapGraphConstructor() {
	// TODO Auto-generated destructor stub
}


bool OverlapGraphConstructor::start() {
	if (!queryDataset->buildDataset(Config::getQueryDatasetFilename())){
		cout << "Error: cannot build query dataset" << endl;
		return false;
	}
	if (!hashTable->insertDataset(queryDataset)){
		cout << "Error: cannot build hash table" << endl;
		return false;
	}
	vector<edge> edgeList;

	while(subject->loadNextChunk(edgeList)){
		// #pragma omp parallel for
		for(int i = 0; i < edgeList.size(); i++){
				hashTable->search(edgeList[i]);

		}
		// end of omp parallel
		for(int i = 0; i < edgeList.size(); i++){
				edgeList[i].print2File();

		}

		edgeList.clear();

	}
}
