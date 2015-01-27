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
	queryAlignmentList.clear();
}

QueryRead::QueryRead(string& sequence, string& name)
{
	readSequence = sequence;
	readName = name;
	readID = 0;
	frequency = 1;
	flag4Removal = false;
	queryAlignmentList.clear();
}

QueryRead::~QueryRead() {
	// TODO Auto-generated destructor stub
	for(int i = 0; i< queryAlignmentList.size();i++)
		delete queryAlignmentList.at(i);
	queryAlignmentList.clear();
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
string QueryRead::reverseComplement(const string & read)
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

bool wayToSort(int i, int j) { return i < j; }

bool QueryRead::printAlignmentToFile(ofstream & filePointer)
{


	for(int i =0; i< this->queryAlignmentList.size();i++)
	{
	string outputString;
	Alignment *align = this->queryAlignmentList[i];
	/*SourceVertexId       DestinationVertexId     Properties

	where properties should have

	orientation, overlap length, substitutions, edits, length1, start1,
	stop1, length2, start2, stop2, error info
	*/
	vector<int> coordinates;
	coordinates.push_back(align->queryEnd);
	coordinates.push_back(align->subjectStart);
	coordinates.push_back(align->subjectEnd);
	coordinates.push_back(0);
	sort(coordinates.begin(), coordinates.end(), wayToSort);
	int overlapLength = coordinates.at(2)-coordinates.at(1)+1;
	int subStart = coordinates.at(1)-align->subjectStart;
	int subEnd = coordinates.at(2)-align->subjectStart;
	std::stringstream sstm;
	// subject print first for good orientation
	sstm <<  align->subjectReadName + "\t" + this->readName + "\t" <<align->orientationTranslate() <<
			"," << overlapLength << "," << align->editInfor.size() << "," << align->editInfor.size()
			<< "," << align->subjectReadSequence.length() << "," << subStart << "," << subEnd <<","
			<< this->getReadLength() << "," << coordinates.at(1)<< "," << coordinates.at(2)<< ", no error info \r\n"<<endl;
	//query print first, but orientation is wrong then
	/*
	sstm << this->readName + "\t" + align->subjectReadName + "\t" <<align->orientationTranslate() <<
			"," << overlapLength << "," << align->editInfor.size() << "," << align->editInfor.size()
			<< "," << this->getReadLength() <<  "," << coordinates.at(1)<<  "," << coordinates.at(2) << ","
			<<align->subjectReadSequence.length() << "," << subStart << "," << subEnd << ", no error info \r\n"<<endl;
			*/
	outputString = sstm.str();
	filePointer<<outputString;
	}

	return true;
}
/*
bool QueryRead::correctErrors(int minimumDepth, double cutoff){

	for(int i = 0; i < queryAlignmentList.size(); i++){

	}
}
*/
