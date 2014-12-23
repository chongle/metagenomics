/*
 * SubjectDataset.cpp
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#include "SubjectDataset.h"

SubjectDataset::SubjectDataset() {
	// TODO Auto-generated constructor stub
	currentFileIndex = 0;

}

SubjectDataset::~SubjectDataset() {
	// TODO Auto-generated destructor stub
}

bool SubjectDataset::setFilenameList(const vector<string> & sFilenames){

	sFilenameList = sFilenames;
	currentFileIndex = 0;

}


bool SubjectDataset::loadNextChunk(vector<SubjectAlignment> & chunkAlignment){


	int k = 0;
	while (k < Config::getChunkSize()){
		string currentLine = currentFileStreamer.getline();

		// if this is end of the file, break out this while loop

		k++;
	}


}
