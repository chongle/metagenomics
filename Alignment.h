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

class QueryRead;




class Alignment {
public:
	Alignment();
	virtual ~Alignment();

	QueryRead * queryRead;

	// Don't need these two variables for error correction to save memory.
	string subjectReadSequence;
	string subjectReadName;

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

	int orientationTranslate();

	// coordinates are defined by the query reads
	// all insertion and deletion is on the subject read
	// upper case char means substitution
	// lower case char means insertion
	// 'D' means deletion
	//
	map<int, char> editInfor;
};


//class edge is used by OverlapGraphConstructor
class edge{
public:
	// alignment input
	string subjectReadSequence;
	string subjectReadName;



	// alignment results
	vector<Alignment*> alignmentList;

	// a list of query reads that are contained by the subject read or duplicate with the subject read
	// when a query is duplicate with a subject read, we compare their names and add the query read to the duplicate read list if its name is greater than the subject read name.
	vector<QueryRead *> containedOrDeplicateReadList;//It's all for contained read if the duplicated reads are removed from the beginning.
	edge();
	~edge();

};


//class SubjectAlignment is used by ErrorCorrector
class SubjectAlignment {
public:
	// alignment input
	string subjectReadSequence;
	string subjectReadName;

	// alignment results
	vector<Alignment*> queryAlignmentList;
	SubjectAlignment();
	~SubjectAlignment();

};

#endif /* ALIGNMENT_H_ */
