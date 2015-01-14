/*
 * QueryRead.cpp
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#include "QueryRead.h"

QueryRead::QueryRead() {

	readSequence = "";
	readName = "";
	readID = 0;
	frequency = 0;
	flag4Removal = false;
}

QueryRead::QueryRead(string& sequence, string& name)
{
	readSequence = sequence;
	readName = name;
	readID = 0;
	frequency = 1;
	flag4Removal = false;
}

QueryRead::~QueryRead() {
	// TODO Auto-generated destructor stub
}

bool QueryRead::needRemoval()
{
	return flag4Removal;
}
void QueryRead::setName(string & name)
{
	readName = name;
}
void QueryRead::setSequence(string & sequence)
{
	readSequence = sequence;
}
void QueryRead::setFrequency(UINT32 number)
{
	frequency = number;
}
void QueryRead::setIdentifier(UINT64 id)
{
	readID = id;
}
string QueryRead::getName()
{
	return readName;
}
string QueryRead::getSequence()
{
	return readSequence;
}

UINT64 QueryRead::getIdentifier()
{
	return readID;
}
UINT32 QueryRead::getFrequency()
{
	return frequency;
}
UINT32 QueryRead::getReadLength()
{
	return this->readSequence.length();
}

bool QueryRead::addAlignment(Alignment* subjectAlignment)
{
	this->queryAlignmentList.push_back(subjectAlignment);
	return true;
}
static string QueryRead::reverseComplement(const string & read)
{
	UINT64 stringLength = read.length();
	string reverse(stringLength,'0');
	for(UINT64 i = 0;i < stringLength; i++)						// Then complement the string. Change A to T, C to G, G to C and T to A.
	{
		if( read[i] & 0X02 ) // C <==> G
			reverse.at(stringLength -  i - 1 ) = read[i] ^ 0X04;
		else // A <==> T
			reverse.at(stringLength -  i - 1 ) = read[i] ^ 0X15;
	}
	return reverse; // return the reverse complement as a string
}
string QueryRead::reverseComplement()
{
	return QueryRead::reverseComplement(this->readSequence);
}

bool QueryRead::correctErrors(){

	for(int i = 0; i < queryAlignmentList.size(); i++){

	}
}
