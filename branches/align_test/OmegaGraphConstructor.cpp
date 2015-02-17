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
	this->totaledgenumber = 0;
}

OmegaGraphConstructor::~OmegaGraphConstructor() {
	// TODO Auto-generated destructor stub
}

bool OmegaGraphConstructor::start()
{

	filePointer.open("tempoutput.txt");
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
cout<<"edge: "<<this->totaledgenumber<<endl;
filePointer.close();
	return true;
}

bool OmegaGraphConstructor::searchHashTable(SubjectEdge * subjectEdge)
{
	SubjectRead *sRead = subjectEdge->subjectRead; 	// Get the current read subject read.
	string subjectReadString = sRead->getSequence(); 		// Get the forward string of subject read.
	string subString;
	for(UINT64 j = 0; j <= sRead->getReadLength()-this->omegaHashTable->getHashStringLength(); j++)
	{
		subString = subjectReadString.substr(j,this->omegaHashTable->getHashStringLength());
		vector<UINT64> * listOfReads=this->omegaHashTable->getListOfReads(subString); // Search the string in the hash table.
		if(!listOfReads->empty()) // If there are some reads that contain subString as prefix or suffix of the read or their reverse complement
		{
			for(UINT64 k = 0; k < listOfReads->size(); k++) // For each such reads.
			{
				UINT64 data = listOfReads->at(k);			// We used bit operations in the hash table. Most significant 2 bits store orientation and least significant 62 bits store read ID.
				UINT16 overlapOffset;
				UINT8 orientation;
				QueryRead *qRead = this->omegaHashTable->getDataset()->getReadFromID(data & 0X3FFFFFFFFFFFFFFF); 	// Least significant 62 bits store the read number.
				if(sRead->getName()<qRead->getName()) 			// No need to discover the same edge again. Only need to explore half of the combinations
				{
//					cout << sRead->getName() << " <<<<<< " << qRead->getName() << " ????? " <<j<<" and " << k << endl;
					continue;
				}
				else
				{
//					cout << sRead->getName() << " >>> " << qRead->getName() << " ????? " <<j<<" and " << k << endl;

				}
				if(sRead->flag4Removal==false && qRead->flag4Removal==false && checkOverlap(sRead,qRead,(data >> 62),j)) // Both read need to be non contained.
				{
					int orientation = data >> 62;
					cout<<sRead->getName() <<" "<< " >>> " << qRead->getName() <<" orientation: "<<orientation<<" position: "<<j<<endl;

//					cout << "S  " << sRead->getSequence() << endl;
//					cout << "Q  " << qRead->getSequence() << endl;
					std::stringstream sstm;
					if(sRead->getSequence()<qRead->getSequence())
					sstm<<" orientation: "<<orientation<<" position: "<<j<<"||||"<<sRead->getSequence()<<"||||"<<qRead->getSequence()<<endl;
					else
					sstm<<" orientation: "<<orientation<<" position: "<<j<<"||||"<<qRead->getSequence()<<"||||"<<sRead->getSequence()<<endl;
					filePointer<<sstm.str();
					this->totaledgenumber++;
/*					switch (data >> 62) // Most significant 2 bit represents  00 - prefix forward, 01 - suffix forward, 10 -  prefix reverse, 11 -  suffix reverse.
					{
						case 0: orientation = 3; overlapOffset = read1->getReadLength() - j; break; 				// 3 = r1>------->r2
						case 1: orientation = 0; overlapOffset = hashTable->getHashStringLength() + j; break; 		// 0 = r1<-------<r2
						case 2: orientation = 2; overlapOffset = read1->getReadLength() - j; break; 				// 2 = r1>-------<r2
						case 3: orientation = 1; overlapOffset = hashTable->getHashStringLength() + j; break; 		// 1 = r2<------->r2
					}
					insertEdge(read1,read2,orientation,read1->getStringForward().length()-overlapOffset); 			// Insert the edge in the graph.
				*/
				}
			}
		}
	}
//YAO	if(graph->at(readNumber)->size() != 0)
//YAO		sort(graph->at(readNumber)->begin(),graph->at(readNumber)->end(), compareEdges); // Sort the list of edges of the current node according to the overlap offset (ascending).
	return true;
}

/**********************************************************************************************************************
	Checks if two read overlaps.
	Hash table search found that a proper substring of read1 is a prefix or suffix of read2 or reverse complement of
	read2 (depents on the orient).
	Orientation 0 means prefix of forward of the read2
	Orientation 1 means suffix of forward of the read2
	Orientation 2 means prefix of reverse of the read2
	Orientation 3 means prefix of reverse of the read2
	We need to check if the remaining of the stings match to see if read1 and read2 overlap.
**********************************************************************************************************************/
bool OmegaGraphConstructor::checkOverlap(SubjectRead *read1, QueryRead *read2, UINT64 orient, UINT64 start)
{
	string string1=read1->getSequence(); // Get the forward string of read1
	UINT64 hashStringLength = this->omegaHashTable->getHashStringLength();
	string string2 = (orient == 0 || orient== 1) ? read2->getSequence() : read2->reverseComplement(); // Get the string from read2 according to orient.
	if(orient == 0 || orient == 2)		// orient 0
										//   >--------MMMMMMMMMMMMMMM*************> 			read1      M means match found by hash table
										//            MMMMMMMMMMMMMMM*************------->      read2      * means we need to check these characters for match
										//				OR
										// orient 2
										//	 >---*****MMMMMMMMMMMMMMM*************> 			read1
										//		      MMMMMMMMMMMMMMM*************-------<	    Reverese complement of read2
	{
		if(string1.length()- start - hashStringLength >= string2.length() - hashStringLength) // The overlap must continue till the end.
			return false;
		return string1.substr(start + hashStringLength, string1.length()-(start + hashStringLength)) == string2.substr(hashStringLength,  string1.length()-(start + hashStringLength)); // If the remaining strings match.
	}
	else								// orient 1
										//   	>********MMMMMMMMMMMMMMM-------------> 			read1      M means match found by hash table
										//  >----********MMMMMMMMMMMMMMM       		    		read2      * means we need to check these characters for match
										//				OR
										// orient 3
										//	 	>********MMMMMMMMMMMMMMM-------------> 			read1
										//	<----********MMMMMMMMMMMMMMM						Reverse Complement of Read2
	{
		if(string2.length()-hashStringLength < start)
			return false;
		return string1.substr(0, start) == string2.substr(string2.length()-hashStringLength-start, start); // If the remaining strings match.
	}
}

