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


	//filter out contained alignment and mark the contained reads
		for(int i = 0; i < edgeList.size(); i++){
			for(int j = 0; j < edgeList[i]->alignmentList.size(); j++){
				if(!isContainedAlignment(edgeList[i]->alignmentList[j])){
				edgeList[i]->alignmentList[j]->queryRead->addAlignment( edgeList[i]->alignmentList[j] );//only add non-contained alignment
				}
				else if(edgeList[i]->subjectReadSequence.length()>=edgeList[i]->alignmentList[j]->queryRead->getReadLength())//query read is contained in subject read
				{
					edgeList[i]->containedOrDeplicateReadList.push_back(edgeList[i]->alignmentList[j]->queryRead);
					edgeList[i]->alignmentList[j]->queryRead->flag4Removal = true;
				}
				else //subject is contained in query read
				{
					//it's almost impossible to know the subject is a contained read, unless it's the terminal case. Accept this reality.
				}
			}

		}

		edgeList.clear();

	}

	for(int i = 0; i < queryDataset->queryReadList.size(); i++){
		if(!queryDataset->queryReadList[i]->flag4Removal){
			queryDataset->queryReadList[i]->printEdges2File();
		}
	}
	return true;
}

//either subject is contained in query or query is contained in subject
bool OverlapGraphConstructor::isContainedAlignment(Alignment * subjectAlignment)
{
	if(subjectAlignment->subjectStart<=0&&subjectAlignment->subjectEnd>=subjectAlignment->queryEnd&&subjectAlignment->editInfor.empty())
		return true;
	else if(subjectAlignment->subjectStart>=0&&subjectAlignment->subjectEnd<=subjectAlignment->queryEnd&&subjectAlignment->editInfor.empty())
		return true;
	else return false;
}

