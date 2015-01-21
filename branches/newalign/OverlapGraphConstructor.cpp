/*
 * OverlapGraphConstructor.cpp
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#include "OverlapGraphConstructor.h"

OverlapGraphConstructor::OverlapGraphConstructor() {
	// TODO Auto-generated constructor stub
	queryDataset = new QueryDataset();
	subject = new SubjectDataset();

	if(Config::isSingleKey)
	{
		singleKeyHashTable = new SingleKeyHashTable();
	}
	else
	{
		doubleKeyHashTable = new DoubleKeyHashTable();
	}

}

OverlapGraphConstructor::~OverlapGraphConstructor() {
	// TODO Auto-generated destructor stub
	delete queryDataset;
	delete subject;
	if(Config::isSingleKey)
	{
		delete singleKeyHashTable;
	}
	else
	{
		delete doubleKeyHashTable;
	}
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

	if(Config::isSingleKey)
	{
		if (!singleKeyHashTable->insertQueryDataset(queryDataset)){
			cout << "Error: cannot build hash table" << endl;
			return false;
		}
	}
	else
	{
		if (!doubleKeyHashTable->insertQueryDataset(queryDataset)){
			cout << "Error: cannot build hash table" << endl;
			return false;
		}
	}

	vector<edge*> edgeList;
	edgeList.clear();

	subject->setFilenameList(Config::subjectFilenameList);
	while(subject->loadNextChunk(edgeList))
	{
//	pragma omp parallel
	{
//		pragma omp for
		for(UINT64 i = 0; i < edgeList.size(); i++)
				searchHashTable(*edgeList[i]);


	}
		// end of omp parallel


	//filter out contained alignment and mark the contained reads
		for(UINT64 i = 0; i < edgeList.size(); i++){cout<<i<<endl;
			for(UINT16 j = 0; j < edgeList[i]->alignmentList.size(); j++){
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

		for(UINT64 i=0;i<edgeList.size();i++)
			delete edgeList.at(i);
		edgeList.clear();

	}// end of while chunk loop
	printEdgesToFile(true, "ourgraph.txt");

	return true;
}

void OverlapGraphConstructor::printEdgesToFile(bool nonRemovedReads, string outFileName)
{
	for(UINT64 i = 0; i < queryDataset->queryReadList.size(); i++){
		if(nonRemovedReads)
		{
		if(!queryDataset->queryReadList[i]->flag4Removal)
		{
			queryDataset->queryReadList[i]->printAlignmentToFile(outFileName);
		}
		}
		else
		{
			queryDataset->queryReadList[i]->printAlignmentToFile(outFileName);

		}

	}
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

