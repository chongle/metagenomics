/*
 * Alignment.cpp
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#include "Alignment.h"



Alignment::Alignment()
{
	// TODO Auto-generated constructor stub
	subjectRead = NULL;
	queryRead = NULL;
	queryOrientation = true;
	subjectStart = 0;
	queryEnd = 0;
	subjectEnd = 0;
	editInfor = NULL;
	containedTag =0;

}
Alignment::Alignment(SubjectRead* sRead)
{
	// TODO Auto-generated constructor stub
	subjectRead = sRead;
	queryRead = NULL;
	queryOrientation = true;
	subjectStart = 0;
	queryEnd = 0;
	subjectEnd = 0;
	editInfor = NULL;
	containedTag =0;

}
Alignment::Alignment(QueryRead* qRead)
{
	// TODO Auto-generated constructor stub
	subjectRead = NULL;
	queryRead = qRead;
	queryOrientation = true;
	subjectStart = 0;
	queryEnd = 0;
	subjectEnd = 0;
	editInfor = NULL;
	containedTag =0;

}

Alignment::Alignment(SubjectRead* sRead, QueryRead* qRead)
{
	// TODO Auto-generated constructor stub
	subjectRead = sRead;
	queryRead = qRead;
	queryOrientation = true;
	subjectStart = 0;
	queryEnd = 0;
	subjectEnd = 0;
	editInfor = NULL;
	containedTag =0;

}

Alignment::~Alignment() {
	// TODO Auto-generated destructor stub
	if(editInfor!=NULL)
	{
		(*editInfor).clear();
		delete editInfor;
	}
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

//either subject is contained in query or query is contained in subject
bool Alignment::isContainedAlignment()
{
	if(this->subjectStart<=0&&this->subjectEnd>=this->queryEnd)
	{
		containedTag =2;
		return true;
	}
	else if(this->subjectStart>=0&&this->subjectEnd<=this->queryEnd)
	{
		containedTag =1;
		return true;
	}
	else return false;
}

int Alignment::getEditDistance()
{
	int size;
	if(editInfor!=NULL)
	 size = editInfor->size();
	else size = 0;
	return size;
}
