/*
 * Alignment.cpp
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#include "Alignment.h"

edge::edge()
{
	subjectReadSequence="";
	subjectReadName="";
	containedOrDeplicateReadList.clear();
	alignmentList.clear();
}
edge::~edge()
{
	containedOrDeplicateReadList.clear();
	alignmentList.clear();
}

SubjectAlignment::SubjectAlignment()
{
	subjectReadSequence="";
	subjectReadName="";
	queryAlignmentList.clear();
}
SubjectAlignment::~SubjectAlignment()
{
	queryAlignmentList.clear();
}

Alignment::Alignment()
{
	// TODO Auto-generated constructor stub
	this->subjectReadName = "";
	this->subjectReadSequence = "";
	queryRead = NULL;
	queryOrientation = true;
	subjectStart = 0;
	queryEnd = 0;
	subjectEnd = 0;

}

Alignment::~Alignment() {
	// TODO Auto-generated destructor stub
	editInfor.clear();
	queryRead = NULL;
}

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


//in our coordinate system, subject is always on the top since it's using sliding window and always going forward,
//the query strand can be forward or reversed, since both cases are hashed in the hash table.
int Alignment::orientationTranslate()
{
	if(this->subjectStart>0)
	{
		if(this->queryOrientation==true) return 1; //read1 = subject, read2 = query
		else return 3; //read1 = subject, read2 = query
	}
	else if(this->subjectStart<0)
	{
		if(this->queryOrientation==true) return 0; //read1 = subject, read2 = query
		else return 2; //read1=subject, read2 = query
	}
	else return -1;
}

