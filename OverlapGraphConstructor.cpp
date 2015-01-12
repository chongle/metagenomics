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

bool OverlapGraphConstructor::searchHashTable(edge & currentEdge)
{
	if(Config::isSingleKey)
	{

		return singleKeyHashTable->singleKeySearch(currentEdge);
	}
	else
	{
		return doubleKeyHashTable->doubleKeySearch(currentEdge);
	}
}
bool OverlapGraphConstructor::start() {
	if (!queryDataset->buildDataset(Config::getQueryDatasetFilename())){
		cout << "Error: cannot build query dataset" << endl;
		return false;
	}


	if (!singleKeyHashTable->insertQueryDataset(queryDataset)){
		cout << "Error: cannot build hash table" << endl;
		return false;
	}
	vector<edge*> edgeList;

	while(subject->loadNextChunk(edgeList)){
	#pragma omp parallel for
	{
		for(int i = 0; i < edgeList.size(); i++){
				searchHashTable(*edgeList[i]);

		}
	}
		// end of omp parallel


	//consider contained reads or not???
		for(int i = 0; i < edgeList.size(); i++){
			for(int j = 0; j < edgeList[i]->alignmentList.size(); j++){
				edgeList[i].alignmentList[j].queryRead->addAlignment( edgeList[i].alignmentList[j] );
			}
			for(int j = 0; j < edgeList[i].containedOrDeplicateReadList.size(); j++){
				edgeList[i].containedOrDeplicateReadList[j]->flag4Removal = true;
			}
		}

		edgeList.clear();

	}

	for(int i = 0; i < queryDataset->queryReadList.size(); i++){
		if(!queryReadList[i]->flag4Removal){
			queryDataset->queryReadList[i]->printEdges2File();
		}
	}
	return true;
}



