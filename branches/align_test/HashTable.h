/*
 * HashTable.h
 *
 *  Created on: Jan 6, 2015
 *      Author: qy2
 */

#ifndef HASHTABLE_H_
#define HASHTABLE_H_

#include "Config.h"
#include "QueryDataset.h"
/*
class DataVector
{
public:
	vector<UINT64>** dataList;
	//	string keystring;
	DataVector(int n)
	{dataList = new vector<UINT64>*[n]();
	for(int i=0;i<n;i++)
		dataList[i] = NULL;
	};
	~DataVector(){delete[] dataList;};
};
*/

//Since there is no internal keystring stored, we cannot handle the collision internally
class HashTable {

//there are three variables which are potentially useless which can be taken off.
//hashKeyLength
//numberOfHashCollision
//maxSingleHashCollision
//because the collsion is not handled here.
//and the key squence and extraction step is not handled here.
private:
	QueryDataset * dataSet;							// Pointer of the dataset.
	string hashTableName;

	UINT64 hashTableSize; 					// Size of hash table, primer number usually.
	UINT16 hashKeyLength;					// Length of the substring of the reads to hash. It's also the key length.
												// It has the 64 bps limit. See the comments in the hashFunction().


	int dataVectorSize;
//	vector <DataVector *> *hashTable; 		// Main structure of hashtable, storing read identifiers.
	vector <map<int,vector<UINT64>*> *> *hashTable; 		// Main structure of hashtable, storing read identifiers.


	UINT64 getPrimeLargerThanNumber(UINT64 number);
	void setHashTableSizeAndInitialize(UINT64 size);


public:
	UINT64 numberOfHashCollision;				// Counted total number of hash collisions. For debugging only.
	UINT64 maxSingleHashCollision;				// Counted maximal number of hash collision for a single case.

	HashTable(int insertSize);
	HashTable(UINT16 keylength,int insertSize);
	HashTable(UINT16 keylength, string name, int insertSize);
	HashTable(UINT16 keylength,  UINT64 size, int insertSize);
	HashTable(UINT16 keylength,  QueryDataset * dataset, int insertSize);
	~HashTable();
	UINT64 hashFunction(const string & subString);
	void InitializeWithDataSize(UINT64 dataSetSize);//hash table will automatically determine a good hash table size

	UINT64 getMaximumHashTableCollision();
	UINT64 getTotalHashTableCollision();
	UINT64 getHashTableSize();


	bool isEmptyAt(UINT64 hashTableIndex);
	map<int,vector<UINT64>*>* getDataVectorsAt(UINT64 hashTableIndex);
	vector<UINT64> * getReadIDListAt(UINT64 hashTableIndex, int dataVectorNum);
	bool insertValueAt(UINT64 hashTableIndex, int dataVectorNum, int dataVecorSize, UINT64 value);



};

#endif /* HASHTABLE_H_ */
