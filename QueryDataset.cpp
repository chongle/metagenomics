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

UINT64 QueryDataset::getNumberOfReads()
{
	return this->numberOfReads;
}
UINT64 QueryDataset::getNumberOfUniqueReads()
{
	return this->numberOfUniqueReads;
}

QueryRead * QueryDataset::getReadFromID(UINT64 ID)
{
	if(ID < 1 || ID > numberOfUniqueReads)	// ID outside the range.
	{

		cout << "ID " << ID << " out of bound."<<endl;
		return NULL;

	}
	else return queryReadList.at(ID - 1);
}

QueryRead * QueryDataset::getReadFromString(const string & read)
{
	UINT64 min = 0, max = getNumberOfUniqueReads()-1;
	string readReverse = QueryRead::reverseComplement(read);
	int comparator;
	if(read.compare(readReverse) < 0)
	{
		while (max >= min) 															// At first search for the forward string.
		{
			UINT64 mid = (min + max) / 2; 	// Determine which subarray to search.
			comparator = queryReadList.at(mid)->getSequence().compare(read.c_str());
			if(comparator == 0)
				return queryReadList.at(mid);
			else if (comparator < 0) 	// Change min index to search upper subarray.
				min = mid + 1;
			else if (comparator > 0) 	// Change max index to search lower subarray.
				max = mid - 1;
		}
	}
	else
	{
		while (max >= min) 																	// If forward string is not found then search for the reverse string
		{
			UINT64 mid = (min+max) / 2; 													// Determine which subarray to search
			comparator = queryReadList.at(mid)->getSequence().compare(readReverse.c_str());
			if( comparator == 0)
				return queryReadList.at(mid);
			else if (comparator < 0) 	// Change min index to search upper subarray.
				min = mid + 1;
			else if (comparator > 0) 	// Change max index to search lower subarray.
				max = mid - 1;
		}
	}
	cout<<"String not found "<<endl;
	return NULL;
}


static bool QueryDataset::qualityFilter(string & sequence)
{
	//Returns true if the read contains only {A,C,G,T} and does not contain more than 80% of the same nucleotide
	UINT64 cnt[4] = {0,0,0,0};
	UINT64 readLength = sequence.length();
	for(UINT64 i = 0; i < readLength; i++) // Count the number of A's, C's , G's and T's in the string.
	{
		if(sequence[i]!= 'A' && sequence[i] != 'C' && sequence[i] != 'G' && sequence[i] != 'T')
			return false;
		cnt[(sequence[i] >> 1) & 0X03]++; // Least significant 2nd and 3rd bits of ASCII value used here
	}
	UINT64 threshold = sequence.length()*.8;	// 80% of the length.
	if(cnt[0] >= threshold || cnt[1] >= threshold || cnt[2] >= threshold || cnt[3] >= threshold)
		return false;	// If 80% bases are the same base.
	return true;
}

bool QueryDataset::duplicateFilter()
{
	UINT64 j = 0;
	QueryRead *temp;
	for(UINT64 i = 0; i < queryReadList.size(); i++)		// Move the unique reads in the top of the sorted list. Store the frequencey of the duplicated reads.
	{
		if(queryReadList.at(j)->getSequence()!= queryReadList.at(i)->getSequence())
		{
			j++;
			temp = queryReadList.at(j);
			queryReadList.at(j) = queryReadList.at(i);
			queryReadList.at(i) = temp;
		}
		else if(i!=j)
		{
			queryReadList.at(j)->setFrequency(queryReadList.at(j)->getFrequency() + 1);
			//set the name as the alphabetical larger name
			string newname = queryReadList.at(i)->getName();
			if(queryReadList.at(j)->getName()< newname)
				queryReadList.at(j)->setName(newname);
		}
	}
	numberOfUniqueReads = j+1;
	cout <<"Number of unique reads: " << numberOfUniqueReads << endl;
	for(UINT64 i = 0 ; i < queryReadList.size(); i++) 		// Assign ID's to the reads.
	{
		if(i < getNumberOfUniqueReads())
			queryReadList.at(i)->setIdentifier(i + 1);
		else
			delete queryReadList.at(i); 					// Free the unused reads.
	}
	queryReadList.resize(numberOfUniqueReads);				//Resize the list.

	return true;
}

void QueryDataset::sortReads()
{
	sort(queryReadList.begin(),queryReadList.end(), compareReads);	// Sort the reads lexicographically.

}

bool QueryDataset::compareReads (QueryRead *read1, QueryRead *read2)
{
	return read1->getSequence() < read2->getSequence();
}



bool QueryDataset::buildDataset(const string & QueryFilename)
{
	string fileName = QueryFilename;
	cout << "Reading dataset: "  << " from file: " << fileName << endl;
	ifstream myFile;
	myFile.open (fileName.c_str());
	if(myFile == NULL)
		cout<<"Unable to open file: " << fileName <<endl;
	UINT64 goodReads = 0, badReads = 0;
	vector<string> line;
	string line1,line0, text, line1ReverseComplement,line2ReverseComplement;
	enum FileType {FASTA, FASTQ, UNDEFINED};
	FileType fileType = UNDEFINED;
	while(!myFile.eof())
	{
		if( (goodReads + badReads ) != 0 && (goodReads + badReads)%1000000 == 0)
			cout<< setw(10) << goodReads + badReads << " reads processed in dataset " << setw(2) << setw(10) << goodReads << " good reads." << setw(10) << badReads << " bad reads." << endl;
		if(fileType == UNDEFINED)
		{
			getline (myFile,text);
			if(text[0] == '>')
				fileType = FASTA;
			else if(text[0] == '@')
				fileType = FASTQ;
			else
				cout<< "Unknown input file format."<<endl;
			myFile.seekg(0, ios::beg);
		}
		line.clear();
		if(fileType == FASTA) 			// Fasta file
		{
			getline (myFile,text);
			line.push_back(text);
			getline (myFile,text,'>');
			line.push_back(text);

			line.at(1).erase(std::remove(line.at(1).begin(), line.at(1).end(), '\n'), line.at(1).end());
			line1 = line.at(1);								// The first string is in the 2nd line.

			line0 = line.at(0);                              //title line is the very first line
		}
		else if(fileType == FASTQ) 					// Fastq file.
		{
			for(UINT64 i = 0; i < 4; i++) 	// Read the remaining 3 lines. Total of 4 lines represent one sequence in a fastq file.
			{
				getline (myFile,text);
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
		if(line1.length() > Config::getminimumOverlapLength() && qualityFilter(line1) ) // Test the read is of good quality.
		{
			QueryRead *r1=new QueryRead();
			r1->setName(readname);
			line1ReverseComplement = r1->reverseComplement(line1);
			if(line1.compare(line1ReverseComplement) < 0) // Store lexicographically smaller between the read and its reverse complement.
				r1->setSequence(line1);
			else
				r1->setSequence(line1ReverseComplement);
			r1->setFrequency(1);

			UINT64 len = r1->getReadLength();

			if(len > longestReadLength)
				longestReadLength = len;
			if(len < shortestReadLength)
				shortestReadLength = len;

			queryReadList.push_back(r1);						// Store the first string in the dataset.
			numberOfReads++;							// Counter of the total number of reads.
			goodReads++;
		}
		else
			badReads++;
	}

	myFile.close();

    cout << "File name: " << fileName << endl;
	cout << setw(10) << goodReads << " good reads in current dataset."  << endl;
	cout << setw(10) << badReads << " bad reads in current dataset." << endl;
	cout << setw(10) << goodReads + badReads << " total reads in current dataset." << endl;
	cout << setw(10) << numberOfReads << " good reads in all datasets." << endl << endl;


	//sort the reads by its alphabetical order; and remove the duplicate reads so then the readID can be assigned uniquely
	sortReads();
	duplicateFilter();

	return true;
}
