/*
 * OmegaGraphConstructor.cpp
 *
 *  Created on: Feb 16, 2015
 *      Author: qy2
 */

#include "OmegaGraphConstructor.h"

OmegaGraphConstructor::OmegaGraphConstructor(OmegaHashTable * omegaHashTable)
{
	// TODO Auto-generated constructor stub
	this->omegaHashTable = omegaHashTable;
}

OmegaGraphConstructor::~OmegaGraphConstructor() {
	// TODO Auto-generated destructor stub
}

bool OmegaGraphConstructor::start()
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
			this->searchHashTable(subjectEdge);
		}

//		  int tid = omp_get_thread_num();
//		  if (tid == 0)
//		    {
//		    int nthreads = omp_get_num_threads();
//		    printf("Number of threads = %d\n", nthreads);
//		    }

	}
		// end of omp parallel


	//filter out contained alignment and mark the contained reads
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
				}
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
	MEMORYSTOP;
	CLOCKSTOP;

	return true;
}

bool OmegaGraphConstructor::searchHashTable(SubjectEdge * subjectEdge)
{

}
