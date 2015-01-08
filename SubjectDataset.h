/*
 * SubjectDataset.h
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#ifndef SUBJECTDATASET_H_
#define SUBJECTDATASET_H_

#include "Config.h"

enum FileType {FASTA, FASTQ, UNDEFINED};


class SubjectDataset {

	vector<string> subjectFilenameList;

	UINT64 chunkSize;

	int currentFileIndex; // the index of the file in sFilenameList, which is currently be processinged loadNextChunk
	FileType currentFileType = UNDEFINED;

	ifstream currentFileStreamer;
	UINT64 currentChunkSize;

	bool outOfDataFlag;
	bool newFileFlag;


public:
	SubjectDataset();
	~SubjectDataset();


	bool setFilenameList(const vector<string> & sFilenames);
	// load the next chunk and populate currentChunk
	// return false if there is no reads in the file.
	bool loadNextChunk(vector<SubjectAlignment*> & chunkData);
	bool loadNextChunk(vector<edge*> & chunkData);
	bool loadNextChunk(vector<SubjectAlignmentPairedEnd*> & chunkData);
//	bool clearCurrentChunk();
};

#endif /* SUBJECTDATASET_H_ */