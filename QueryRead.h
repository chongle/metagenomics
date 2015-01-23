/*
 * QueryRead.h
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#ifndef QUERYREAD_H_
#define QUERYREAD_H_

#include "Config.h"
#include "Alignment.h"
class Alignment;
class QueryRead {
	vector<Alignment*> queryAlignmentList;
	string readSequence;
	string readName;
	UINT64 readID; 						// Unique Identification of the read. start from one. zero means a new read.
	UINT32 frequency; 						// Frequency of the read. Number of times this read is present in the dataset. Used in some statistical analysis.
	string correctedRead;


public:
	bool flag4Removal; // this is a contained or duplicate read
	QueryRead();
	QueryRead(string & sequence, string & name);
	virtual ~QueryRead();
	bool addAlignment(Alignment* subjectAlignment);
	bool correctErrors();
	bool needRemoval();
	void setName(string & name);
	void setSequence(string & sequence);
	void setFrequency(UINT32 number);
	void setIdentifier(UINT64 id);
	string getName();
	string getSequence();
	UINT64 getIdentifier();
	UINT32 getFrequency();
	UINT32 getReadLength();
	static string reverseComplement(const string & read);
	string reverseComplement();
	bool printAlignmentToFile(ofstream & filePointer);
};

#endif /* QUERYREAD_H_ */
