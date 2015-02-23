/*
 * QueryDatasetFilter.cpp
 *
 *  Created on: Feb 23, 2015
 *      Author: qy2
 */

#include "QueryDatasetFilter.h"

QueryDatasetFilter::QueryDatasetFilter(OmegaHashTable * omegaHashTable) {
	// TODO Auto-generated constructor stub
	this->omegaHashTable = omegaHashTable;
	this->dataset = this->omegaHashTable->getDataset();

	this->tagList =  NULL;

	this->superNameList = NULL;

	this->superReadLength = NULL;
}
/*
QueryDatasetFilter::~QueryDatasetFilter() {
	// TODO Auto-generated destructor stub
	if(this->tagList!=NULL)
	{
		this->tagList->clear();
		delete this->tagList;
	}
	if(this->superNameList!=NULL)
	{
		this->superNameList->clear();
		delete this->superNameList;
	}
	if(this->superReadLength!=NULL)
	{
		this->superReadLength->clear();
		delete this->superReadLength;
	}
}

bool QueryDatasetFilter::start()
{
	filePointer.open("tempoutput.txt");
	CLOCKSTART;
	MEMORYSTART;

	if(this->tagList==NULL)
	{
		this->tagList = new vector<UINT8>;
		this->tagList->reserve(dataset->getNumberOfUniqueReads()+1);
	}
	if(this->superNameList==NULL)
	{
		this->superNameList = new vector<string&>;
		this->superNameList->reserve(dataset->getNumberOfUniqueReads()+1);
	}
	if(this->superReadLength==NULL)
	{
		this->superReadLength = new vector<int>;
		this->superReadLength->reserve(dataset->getNumberOfUniqueReads()+1);
	}


	string emptystring = 0;
	for(UINT64 i = 0; i <= dataset->getNumberOfUniqueReads(); i++) // Initialization
	{
		this->tagList->push_back(0);
		this->superNameList->push_back(emptystring);
		this->superReadLength->push_back(0);
	}

	vector<SubjectRead *>* subjectReadList= new vector<SubjectRead *>();
	subjectReadList->clear();

	SubjectDataset *subjectDataset = new SubjectDataset();
	subjectDataset->setFilenameList(Config::subjectFilenameList);
	while(subjectDataset->loadNextChunk(subjectReadList,this->omegaHashTable->getDataset()))
	{
		vector<SubjectEdge*>* subjectEdgeList = new  vector<SubjectEdge*>();
		for(UINT64 i = 0; i < subjectReadList->size(); i++)
		{
			SubjectRead * subjectRead = subjectReadList->at(i);
			SubjectEdge * subjectEdge = new SubjectEdge(subjectRead);
			subjectEdgeList->push_back(subjectEdge);
		}
//omp_set_dynamic(0);
//omp_set_num_threads(Config::numberOfThreads);
//#pragma omp parallel
	{
//#pragma omp for schedule(dynamic)
		for(UINT64 i = 0; i < subjectEdgeList->size(); i++)
		{
			SubjectEdge * subjectEdge = subjectEdgeList->at(i);
			this->searchHashTable(subjectEdge);
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
//				UINT16 overlapOffset;
//				UINT8 orientation;
				QueryRead *qRead = this->omegaHashTable->getDataset()->getReadFromID(data & 0X3FFFFFFFFFFFFFFF); 	// Least significant 62 bits store the read number.

				if(checkIdenticalRead(sRead,qRead,(data >> 62),j))
				{
					string subjectName = sRead->getName();
					if(subjectName>=qRead->getName())
				}

				if(checkOverlapForContainedRead(sRead,qRead,(data >> 62),j)) // Both read need to be non contained.
				{
//YAO				int orientation = data >> 62;
//YAO					cout<<sRead->getName() <<" "<< " >>> " << qRead->getName() <<" orientation: "<<orientation<<" position: "<<j<<endl;

//YAO					cout << "S  " << sRead->getSequence() << endl;
//YAO				cout << "Q  " << qRead->getSequence() << endl;
//YAO					std::stringstream sstm;
//YAO					if(sRead->getSequence()<qRead->getSequence())
//YAO					sstm<<" orientation: "<<orientation<<" position: "<<j<<"||||"<<sRead->getSequence()<<"||||"<<qRead->getSequence()<<endl;
//YAO					else
//YAO					sstm<<" orientation: "<<orientation<<" position: "<<j<<"||||"<<qRead->getSequence()<<"||||"<<sRead->getSequence()<<endl;
//YAO				filePointer<<sstm.str();
					this->totaledgenumber++;
					switch (data >> 62) // Most significant 2 bit represents  00 - prefix forward, 01 - suffix forward, 10 -  prefix reverse, 11 -  suffix reverse.
					{
						case 0: orientation = 3; overlapOffset = read1->getReadLength() - j; break; 				// 3 = r1>------->r2
						case 1: orientation = 0; overlapOffset = hashTable->getHashStringLength() + j; break; 		// 0 = r1<-------<r2
						case 2: orientation = 2; overlapOffset = read1->getReadLength() - j; break; 				// 2 = r1>-------<r2
						case 3: orientation = 1; overlapOffset = hashTable->getHashStringLength() + j; break; 		// 1 = r2<------->r2
					}
					insertEdge(read1,read2,orientation,read1->getStringForward().length()-overlapOffset); 			// Insert the edge in the graph.

				}
			}
		}
	}
//YAO	if(graph->at(readNumber)->size() != 0)
//YAO		sort(graph->at(readNumber)->begin(),graph->at(readNumber)->end(), compareEdges); // Sort the list of edges of the current node according to the overlap offset (ascending).
	return true;
}
bool QueryDatasetFilter::checkIdenticalRead(SubjectRead *read1, QueryRead *read2, UINT64 orient, UINT64 start)
{
	if(start!=0)return false;
	string string1=read1->getSequence();
	string string2 = (orient == 0 || orient== 1) ? read2->getSequence() : read2->reverseComplement();
	if(string1 == string2) return true;
}

*/

/**********************************************************************************************************************
	Hash table search found that a proper substring of read1 is a prefix or suffix of read2 or reverse complement of
	read2 (depents on the orient).
	orient 0 means prefix of forward of the read2
	orient 1 means suffix of forward of the read2
	orient 2 means prefix of reverse of the read2
	orient 3 means prefix of reverse of the read2
	We need to check if the remaining of the stings match to see if read2 is contained in read1.
**********************************************************************************************************************/
bool QueryDatasetFilter::checkOverlapForContainedRead(SubjectRead *read1, QueryRead *read2, UINT64 orient, UINT64 start)
{
	string string1=read1->getSequence(); // Get the forward of read1
	UINT64 hashStringLength = this->omegaHashTable->getHashStringLength(), lengthRemaining1, lengthRemaining2;
	string string2 = (orient == 0 || orient== 1) ? read2->getSequence() : read2->reverseComplement(); // Get the string in read2 based on the orientation.
	if(orient == 0 || orient == 2)
									// orient 0
									//   >--------MMMMMMMMMMMMMMM*******------> read1      M means match found by hash table
									//            MMMMMMMMMMMMMMM*******>       read2      * means we need to check these characters for match
									//				OR
									// orient 2
									//	 >---*****MMMMMMMMMMMMMMM*******------> read1
									//		      MMMMMMMMMMMMMMM*******<	    Reverese complement of read2
	{
		lengthRemaining1 = string1.length() - start - hashStringLength; 	// This is the remaining of read1
		lengthRemaining2 = string2.length() - hashStringLength; 	// This is the remaining of read2
		if(lengthRemaining1 >= lengthRemaining2)
		{
			return string1.substr(start + hashStringLength, lengthRemaining2) == string2.substr(hashStringLength, lengthRemaining2); // If the remaining of the string match, then read2 is contained in read1
		}
	}
	else							// orient 1
									//   >---*****MMMMMMMMMMMMMMM-------------> read1      M means match found by hash table
									//      >*****MMMMMMMMMMMMMMM       		read2      * means we need to check these characters for match
									//				OR
									// orient 3
									//	 >---*****MMMMMMMMMMMMMMM-------------> read1
									//		<*****MMMMMMMMMMMMMMM				Reverse Complement of Read2
	{
		lengthRemaining1 = start;
		lengthRemaining2 = string2.length() - hashStringLength;
		if(lengthRemaining1 >= lengthRemaining2)
		{
			return string1.substr(start - lengthRemaining2, lengthRemaining2) == string2.substr(0, lengthRemaining2); // If the remaining of the string match, then read2 is contained in read1
		}
	}
	return false;

}
