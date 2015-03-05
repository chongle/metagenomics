/*
 * Alignment.h
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

#include "Config.h"
#include "QueryRead.h"
#include "SubjectRead.h"

class QueryRead;




class Alignment {



public:
	Alignment();
	Alignment(SubjectRead* sRead);
	Alignment(QueryRead* qRead);
	Alignment(SubjectRead* sRead, QueryRead* qRead);
	~Alignment();

	/*
	// Don't need these two variables for error correction to save memory.
	string subjectReadSequence;
	string subjectReadName;
	*/
	SubjectRead * subjectRead;
	//alignment is subject to the subject reads, which is on top of the orientation graph

	QueryRead * queryRead;

	// orientation of the query read
	bool queryOrientation;		// true for forward, false for reverse

//******OMEGA original definition for the alignment orientation******
	// orient 0
	//   >--------MMMMMMMMMMMMMMM*************> 			read1      M means match found by hash table
	//            MMMMMMMMMMMMMMM*************------->      read2      * means we need to check these characters for match

	// orient 2
	//	 >---*****MMMMMMMMMMMMMMM*************> 			read1
	//		      MMMMMMMMMMMMMMM*************-------<	    Reverese complement of read2

	// orient 1
	//   	>********MMMMMMMMMMMMMMM-------------> 			read1      M means match found by hash table
	//  >----********MMMMMMMMMMMMMMM       		    		read2      * means we need to check these characters for match

	// orient 3
	//	 	>********MMMMMMMMMMMMMMM-------------> 			read1
	//	<----********MMMMMMMMMMMMMMM						Reverse Complement of Read2

//*******************************************************************


	// coordinates of the overlap alignment, which is defined by the query reads
	//  The following case will have a positive subject start position
	//  query:        XXXXXXMMMMMMMM
	//	subject:            MMMMMMMMXXXXXXX
	//  The following case will have a negative subject start position
	//  query:        MMMMMMXXXXXXXX
	//	subject: XXXXXMMMMMM
	int subjectStart;
	int queryEnd;
	int subjectEnd;

	UINT8 containedTag; //0-not contained alignment, 1-subject is contained in query 2-query is contained in subject



	int orientationTranslate();
	bool isContainedAlignment();
	int getEditDistance(); //return the size of the editInfor


	// coordinates are defined by the query reads
	// all insertion and deletion is on the subject read
	// upper case char means substitution
	// lower case char means insertion
	// 'D' means deletion
	//
	map<int, char>* editInfor;
};

class Edge{
public:
	QueryRead* source;
	QueryRead* destination;
	Edge(){source = NULL; destination = NULL;};
	~Edge(){};
};

class ContainedAlignment {



public:

	ContainedAlignment(SubjectRead* sRead, QueryRead* qRead);
	~ContainedAlignment();

	/*
	// Don't need these two variables for error correction to save memory.
	string subjectReadSequence;
	string subjectReadName;
	*/
	SubjectRead * subjectRead;
	//alignment is subject to the subject reads, which is on top of the orientation graph

	QueryRead * queryRead;
};


//class edge is used by OverlapGraphConstructor
class SubjectEdge{
public:
	// alignment input
	SubjectRead* subjectRead;

	// alignment results
	vector<Alignment*>* alignmentList;
	vector<ContainedAlignment*>* contained_alignmentList;

	// a list of query reads that are duplicated with the subject read
	// when a query is duplicate with a subject read, we compare their names and add the query read to the duplicate read list if its name is greater than the subject read name.
	vector<QueryRead *>* DuplicateReadList;
	SubjectEdge(SubjectRead* sRead);
	~SubjectEdge();
	bool addAlignment(Alignment* subjectAlignment);
	bool addContainedAlignment(ContainedAlignment* subjectAlignment);
	bool addDuplicateList(QueryRead* queryRead);

};

#endif /* ALIGNMENT_H_ */
