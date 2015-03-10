/*
 * OverlapGraphConstructor.cpp
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#include "OverlapGraphConstructor.h"

OverlapGraphConstructor::OverlapGraphConstructor(HashTableMethod * hashTableMethod) {
	// TODO Auto-generated constructor stub
	this->totaledgenumber = 0;
	this->minimumOverlapLength = Config::minimumOverlapLength;
	this->hashTableMethod = hashTableMethod;

}

OverlapGraphConstructor::~OverlapGraphConstructor() {
	// TODO Auto-generated destructor stub

}


bool OverlapGraphConstructor::start()
{

	CLOCKSTART;
	MEMORYSTART;
	vector<SubjectRead *>* subjectReadList= new vector<SubjectRead *>();
	subjectReadList->clear();

	SubjectDataset *subjectDataset = new SubjectDataset();
	subjectDataset->setFilenameList(Config::subjectFilenameList);
	while(subjectDataset->loadNextChunk(subjectReadList))
	{
		vector<SubjectEdge*>* subjectEdgeList = new  vector<SubjectEdge*>();
		for(UINT64 i = 0; i < subjectReadList->size(); i++)
		{
			SubjectRead * sRead = subjectReadList->at(i);
			SubjectEdge * sEdge = new SubjectEdge(sRead);
			subjectEdgeList->push_back(sEdge);
		}
omp_set_dynamic(0);
omp_set_num_threads(Config::numberOfThreads);
#pragma omp parallel
	{
#pragma omp for schedule(dynamic)
		for(UINT64 i = 0; i < subjectEdgeList->size(); i++)
		{
			SubjectEdge * subjectEdge = subjectEdgeList->at(i);
			this->hashTableMethod->searchHashTable(subjectEdge);
		}

	}
		// end of omp parallel


	//clear the unused subject reads
		for(UINT64 i = 0; i < subjectEdgeList->size(); i++)
		{
			SubjectEdge * subjectEdge = subjectEdgeList->at(i);
			if(subjectEdge->alignmentList==NULL)
			{
				delete subjectEdge->subjectRead;
				subjectEdge->subjectRead = NULL;
			}
			else
			{
				//add alignments to the query reads before the edges are destroyed.
				for(UINT16 j = 0; j < subjectEdge->alignmentList->size(); j++)
				{
					 Alignment * alignment= subjectEdge->alignmentList->at(j);
					 alignment->queryRead->addAlignment(alignment);
					 this->totaledgenumber++;

				}
				 subjectDataset->addSubjectRead(subjectEdge->subjectRead);
//				 cout<<"add subject:"<<subjectEdge->subjectRead->getName()<<endl;
			}
			delete subjectEdge;

		}

		//
		subjectEdgeList->clear();
		delete subjectEdgeList;
		//clear the subject read vector for the next round
		subjectReadList->clear();
		subjectReadList->resize(0);

	}// end of while chunk loop

	printEdgesToFile(Config::outputfilename);
	delete subjectReadList;
	delete subjectDataset;
	MEMORYSTOP;
	CLOCKSTOP;
cout<<"edge: "<<this->totaledgenumber<<endl;

	return true;

}

bool OverlapGraphConstructor::printEdgesToFile(string outFileName)
{
	ofstream filePointer;
	filePointer.open(outFileName.c_str());
	if(filePointer == NULL)
	{
		cout<<"Unable to open file: "<<endl;
		return false;
	}
	QueryDataset * dataset = this->hashTableMethod->getDataset();

	for(UINT64 i = 1; i <= dataset->getNumberOfUniqueReads(); i++)
		dataset->getReadFromID(i)->printAlignmentToFile(filePointer);





	filePointer.close();
	return true;
}


//either subject is contained in query or query is contained in subject
bool OverlapGraphConstructor::isContainedAlignment(Alignment * subjectAlignment)
{
	if(subjectAlignment->subjectStart<=0&&subjectAlignment->subjectEnd>=subjectAlignment->queryEnd)
		return true;
	else if(subjectAlignment->subjectStart>=0&&subjectAlignment->subjectEnd<=subjectAlignment->queryEnd)
		return true;
	else return false;
}


