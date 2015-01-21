/*
 * SubjectDataset.cpp
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#include "SubjectDataset.h"
#include "QueryDataset.h"

SubjectDataset::SubjectDataset() {
	// TODO Auto-generated constructor stub
	currentFileIndex = 0;
	currentFileType = UNDEFINED;
	chunkSize = Config::getStreamChunkSize();
	outOfDataFlag = false;
	newFileFlag = true;
	currentChunkSize = 0;
	subjectFilenameList.clear();


}

SubjectDataset::~SubjectDataset() {
	// TODO Auto-generated destructor stub
	subjectFilenameList.clear();
}

bool SubjectDataset::setFilenameList(const vector<string> & sFilenames){

	subjectFilenameList = sFilenames;
	currentFileIndex = 0;
	return true;

}


bool SubjectDataset::loadNextChunk(vector<SubjectAlignment*> & chunkData)
{
	if(outOfDataFlag) return false;

	if(newFileFlag)
	{
		string currentFileName = subjectFilenameList.at(currentFileIndex);

		currentFileStreamer.open (currentFileName.c_str());
		if(currentFileStreamer == NULL)
			cout<<"Unable to open file: " << currentFileStreamer <<endl;

		if(currentFileType == UNDEFINED)
		{
			string text;
			getline (currentFileStreamer,text);
			if(text[0] == '>')
				currentFileType = FASTA;
			else if(text[0] == '@')
				currentFileType = FASTQ;
			else
				cout<< "Unknown input file format."<<endl;
			currentFileStreamer.seekg(0, ios::beg);
		}
		newFileFlag = false;
	}

	currentChunkSize = 0;
	for(UINT64 k=0;k < chunkSize;k++)
	{
		// if this is end of the file, break out this while loop
		if(currentFileStreamer.eof())
		{
			currentFileIndex++;
			if(currentFileIndex>=subjectFilenameList.size())
				outOfDataFlag = true;
			else newFileFlag = true;
			currentFileStreamer.close();
			break;
		}

		vector<string> line;
		string line1,line0, text, line1ReverseComplement,line2ReverseComplement;
		if(currentFileType == FASTA) 			// Fasta file
		{
			getline(currentFileStreamer,text);
			line.push_back(text);
			getline(currentFileStreamer,text,'>');
			line.push_back(text);

			line.at(1).erase(std::remove(line.at(1).begin(), line.at(1).end(), '\n'), line.at(1).end());
			line1 = line.at(1);								// The first string is in the 2nd line.

			line0 = line.at(0);                              //title line is the very first line
		}
		else if(currentFileType == FASTQ) 					// Fastq file.
		{
			for(UINT64 i = 0; i < 4; i++) 	// Read the remaining 3 lines. Total of 4 lines represent one sequence in a fastq file.
			{
				getline (currentFileStreamer,text);
				line.push_back(text);
			}
			line1 = line.at(1); 			// The first string is in the 2nd line.
			line0 = line.at(0);             //title line is the very first line
		}

		//remove the > if it appears in the title or name
		string readname="";
		if(line0[0] == '>')
			readname = line0.substr(1);
		else readname = line0;

		for (std::string::iterator p = line1.begin(); line1.end() != p; ++p) // Change the case
		    *p = toupper(*p);
		if(line1.length() > Config::getminimumOverlapLength() && QueryDataset::qualityFilter(line1) ) // Test the read is of good quality.
		{


			SubjectAlignment* subAlign = new SubjectAlignment();
			subAlign->subjectReadSequence = line1;
			subAlign->subjectReadName = readname;
			chunkData.push_back(subAlign);
			currentChunkSize++;
		}



		line.clear();

	}

	return true;
}

bool SubjectDataset::loadNextChunk(vector<edge*> & chunkData)
{
	if(outOfDataFlag) return false;

	if(newFileFlag)
	{
		string currentFileName = subjectFilenameList.at(currentFileIndex);

		currentFileStreamer.open (currentFileName.c_str());
		if(currentFileStreamer == NULL)
			cout<<"Unable to open file: " << currentFileStreamer <<endl;

		if(currentFileType == UNDEFINED)
		{
			string text;
			getline (currentFileStreamer,text);
			if(text[0] == '>')
				currentFileType = FASTA;
			else if(text[0] == '@')
				currentFileType = FASTQ;
			else
				cout<< "Unknown input file format."<<endl;
			currentFileStreamer.seekg(0, ios::beg);
		}
		newFileFlag = false;
	}

	currentChunkSize = 0;
	for(UINT64 k=0;k < chunkSize;k++)
	{
		// if this is end of the file, break out this while loop
		if(currentFileStreamer.eof())
		{
			currentFileIndex++;
			if(currentFileIndex>=subjectFilenameList.size())
				outOfDataFlag = true;
			else newFileFlag = true;
			currentFileStreamer.close();
			break;
		}

		vector<string> line;
		string line1,line0, text, line1ReverseComplement,line2ReverseComplement;
		if(currentFileType == FASTA) 			// Fasta file
		{
			getline(currentFileStreamer,text);
			line.push_back(text);
			getline(currentFileStreamer,text,'>');
			line.push_back(text);

			line.at(1).erase(std::remove(line.at(1).begin(), line.at(1).end(), '\n'), line.at(1).end());
			line1 = line.at(1);								// The first string is in the 2nd line.

			line0 = line.at(0);                              //title line is the very first line
		}
		else if(currentFileType == FASTQ) 					// Fastq file.
		{
			for(UINT64 i = 0; i < 4; i++) 	// Read the remaining 3 lines. Total of 4 lines represent one sequence in a fastq file.
			{
				getline (currentFileStreamer,text);
				line.push_back(text);
			}
			line1 = line.at(1); 			// The first string is in the 2nd line.
			line0 = line.at(0);             //title line is the very first line
		}

		//remove the > if it appears in the title or name
		string readname="";
		if(line0[0] == '>')
			readname = line0.substr(1);
		else readname = line0;

		for (std::string::iterator p = line1.begin(); line1.end() != p; ++p) // Change the case
		    *p = toupper(*p);
		if(line1.length() > Config::getminimumOverlapLength() && QueryDataset::qualityFilter(line1) ) // Test the read is of good quality.
		{


			edge* Edge = new edge();
			Edge->subjectReadSequence = line1;
			Edge->subjectReadName = readname;
			chunkData.push_back(Edge);
			currentChunkSize++;
		}



		line.clear();

	}

	return true;
}

bool SubjectDataset::loadNextChunk(vector<SubjectAlignmentPairedEnd*> & chunkData)
{
	if(outOfDataFlag) return false;

	if(newFileFlag)
	{
		string currentFileName = subjectFilenameList.at(currentFileIndex);

		currentFileStreamer.open (currentFileName.c_str());
		if(currentFileStreamer == NULL)
			cout<<"Unable to open file: " << currentFileStreamer <<endl;

		if(currentFileType == UNDEFINED)
		{
			string text;
			getline (currentFileStreamer,text);
			if(text[0] == '>')
				currentFileType = FASTA;
			else if(text[0] == '@')
				currentFileType = FASTQ;
			else
				cout<< "Unknown input file format."<<endl;
			currentFileStreamer.seekg(0, ios::beg);
		}
		newFileFlag = false;
	}

	currentChunkSize = 0;
	for(UINT64 k=0;k < chunkSize;k++)
	{
		// if this is end of the file, break out this while loop
		if(currentFileStreamer.eof())
		{
			currentFileIndex++;
			if(currentFileIndex>=subjectFilenameList.size())
				outOfDataFlag = true;
			else newFileFlag = true;
			currentFileStreamer.close();
			break;
		}

		vector<string> line;
		string line1,line0, text, line1ReverseComplement,line2ReverseComplement;
		if(currentFileType == FASTA) 			// Fasta file
		{
			getline(currentFileStreamer,text);
			line.push_back(text);
			getline(currentFileStreamer,text,'>');
			line.push_back(text);

			line.at(1).erase(std::remove(line.at(1).begin(), line.at(1).end(), '\n'), line.at(1).end());
			line1 = line.at(1);								// The first string is in the 2nd line.

			line0 = line.at(0);                              //title line is the very first line
		}
		else if(currentFileType == FASTQ) 					// Fastq file.
		{
			for(UINT64 i = 0; i < 4; i++) 	// Read the remaining 3 lines. Total of 4 lines represent one sequence in a fastq file.
			{
				getline (currentFileStreamer,text);
				line.push_back(text);
			}
			line1 = line.at(1); 			// The first string is in the 2nd line.
			line0 = line.at(0);             //title line is the very first line
		}

		//remove the > if it appears in the title or name
		string readname="";
		if(line0[0] == '>')
			readname = line0.substr(1);
		else readname = line0;

		for (std::string::iterator p = line1.begin(); line1.end() != p; ++p) // Change the case
		    *p = toupper(*p);
		if(line1.length() > Config::getminimumOverlapLength() && QueryDataset::qualityFilter(line1) ) // Test the read is of good quality.
		{


			SubjectAlignmentPairedEnd* subAlign = new SubjectAlignmentPairedEnd();
			subAlign->subjectReadSequence = line1;
			subAlign->subjectReadName = readname;
			chunkData.push_back(subAlign);
			currentChunkSize++;
		}



		line.clear();

	}

	return true;
}
