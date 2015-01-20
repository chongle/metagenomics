/*
 * AlignmentPairedEnd.cpp
 *
 *  Created on: Jan 6, 2015
 *      Author: qy2
 */

#include "AlignmentPairedEnd.h"

SubjectAlignmentPairedEnd::SubjectAlignmentPairedEnd()
{
	subjectReadSequence="";
	subjectReadName="";
	queryAlignmentList.clear();
}
SubjectAlignmentPairedEnd::~SubjectAlignmentPairedEnd()
{
	queryAlignmentList.clear();
}

AlignmentPairedEnd::AlignmentPairedEnd() {
	// TODO Auto-generated constructor stub

}

AlignmentPairedEnd::~AlignmentPairedEnd() {
	// TODO Auto-generated destructor stub
}

