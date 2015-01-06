/*
 * SubjectDataset.h
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#ifndef SUBJECTDATASET_H_
#define SUBJECTDATASET_H_

#include "Config.h"

class SubjectDataset {

	vector<string> sFilenameList;

	int currentFileIndex; // the index of the file in sFilenameList, which is currently be processinged loadNextChunk

	ifstream currentFileStreamer;


public:
	SubjectDataset();
	~SubjectDataset();


	bool setFilenameList(const vector<string> & sFilenames);
	// load the next chunk and populate currentChunk
	// return false if there is no reads in the file.
	bool loadNextChunk(vector<SubjectAlignment> & chunkAlignment);
	bool clearCurrentChunk();
};

#endif /* SUBJECTDATASET_H_ */
