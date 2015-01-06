/*
 * QueryDataset.cpp
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#include "QueryDataset.h"

QueryDataset::QueryDataset() {
	// TODO Auto-generated constructor stub

}

QueryDataset::~QueryDataset() {
	for(int i = 0; i < queryReadList.size(); i++){
		delete queryReadList[i];
	}
}

bool QueryDataset::buildDataset(const string & sQueryFilename)
{
	// get from Dataset.cpp in Omega

	// queryReadList
}
