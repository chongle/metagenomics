/*
 * OverlapGraphConstructor.h
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#ifndef OVERLAPGRAPHCONSTRUCTOR_H_
#define OVERLAPGRAPHCONSTRUCTOR_H_

class edge{
public:
	// alignment input
	string subjectReadSequence;
	string subjectReadName;



	// alignment results
	vector<string> queryReadName;
	vector<int> orientation;
	vector<int> overlapLength;


};

class OverlapGraphConstructor {
	QueryDataset * queryDataset;
	HashTable * hashTable;
	SubjectDataset * subject;
public:
	OverlapGraphConstructor();
	~OverlapGraphConstructor();
	bool start();

};

#endif /* OVERLAPGRAPHCONSTRUCTOR_H_ */
