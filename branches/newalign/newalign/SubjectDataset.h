/*
 * SubjectDataset.h
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#ifndef SUBJECTDATASET_H_
#define SUBJECTDATASET_H_

class SubjectDataset {
public:
	SubjectDataset();
	~SubjectDataset();


	bool setFilename(string sFilename);
	// load the next chunk and populate currentChunk
	// return false if there is no reads in the file.
	bool loadNextChunk();
	bool clearCurrentChunk();
};

#endif /* SUBJECTDATASET_H_ */
