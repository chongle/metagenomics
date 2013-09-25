/*
 * Reads.h
 *
 * Created on: April 22, 2013
 * Author: Md. Bahlul Haider
 */


#ifndef READ_H_
#define READ_H_

#include "Common.h"

class Edge;

// Ted: SVN test
struct MPlist
{
	UINT64 matePairID;				// This list stores the ID of the matepair of the current read.
									// CP: this is NOT a list, but a read ID, right?
									// BH: yes this in not a list. It contains only one read ID. We use a vector of the structure to store all the matepairs.
	UINT8 matePairOrientation; 		// List of mate pair orientation.
									// 0 = 00 means the reverse of this read and the reverse of the second read are matepairs.
									// 1 = 01 means the reverse of this read and the forward of the second read are matepairs.
									// 2 = 10 means the forward of this read and the reverse of the second read are matepairs.
									// 3 = 11 means the forward of this read and the forward of the second read are matepairs.
	UINT8 datasetNumber;			//Dataset number.
									// CP: different datasets have different insert sizes
};

/**********************************************************************************************************************
	Class to store a read.
**********************************************************************************************************************/

class Read
{
	private:
		UINT64 readNumber; 						// Unique Identification of the read.
		string read; 							// String representation of the read.
		// CP: Is readReverse stored in memory for every read? Why can't this be calculated on the flight to save memory?
		// BH: We need to find reverse complement of the reads in many places. That's whey I computed them once and saved them.
		string readReverse;
		UINT32 frequency; 						// Frequency of the read. Number of times this read is present in the dataset. Used in some statistical analysis.
		vector<MPlist> *matePairList;
		vector<Edge *> *listOfEdgesForward;   	// List of edges that contain the forward string of this read.
		// CP: How is the location defined? from where to the beginning of this read?
		// BH: Location is the distance from the first read on the edge.
		vector<UINT64> *locationOnEdgeForward;	// List of locations on the edges that contain the forward string of the current read.
		vector<Edge *> *listOfEdgesReverse;		// List of edges that contain the reverse string of this read.
		vector<UINT64> *locationOnEdgeReverse;	// List of locations on the edges that contain the reverse string of the current read.

		string reverseComplement(const string & read);

	public:
		// CP: can the two variables be removed?
		// BH: We used this variables in resolved Node function. Later we might want to use it in other functions. That's why I did not remove them.
		UINT64 coverageDepth;					// Estimated depth of coverage.  --> Not use anymore
		UINT64 locationInDataset;				// Not for debugging only. Will be removed.  --> Not use anymore
		bool isContainedRead;
		UINT64 superReadID;						// 0 = not a contained read
												// otherwise superReadID contains the ID of the uniqe super read.
		Read(void);								// Default constructor.
		Read(const string & s);					// Another constructor.
		~Read(void);							// Destructor.

		bool setRead(const string & s); 		// Set the read.
		bool setReadNumber(UINT64 id); 			// Set the read number.
		bool setFrequency(UINT32 freq);			// Set the ferquency of the read.

		string getStringForward(void){return read;} 									// Get the forward string of the current read.
		string getStringReverse(void){return readReverse;} 								// Get the reverse string of the current read.
		UINT16 getReadLength(void){return read.length();} 								// Get the length of the string in the current read.
		UINT64 getReadNumber(void) {return readNumber;} 								// Get the read number of the current read.
		UINT32 getFrequency(void) {return frequency;}									// Get the frequency of the current read.
		vector<MPlist> * getMatePairList(void) {return matePairList;} 					// Get the list of matepairs.
		vector<Edge *> * getListOfEdgesForward(void){return listOfEdgesForward;}		// Get the list of edges that contain the forward string of the current read.
		vector<UINT64> * getLocationOnEdgeForward(void){return locationOnEdgeForward;}	// Get the list of locations on the edges that contain the forward string of the current read.
		vector<Edge *> * getListOfEdgesReverse(void){return listOfEdgesReverse;}		// Get the list of edges that contain the reverse string of the current read.
		vector<UINT64> * getLocationOnEdgeReverse(void){return locationOnEdgeReverse;}	// Get the list of locations on the edges that contain the reverse string of the current read.
		bool addMatePair(UINT64 read, UINT8 orientation, UINT64 datasetNumber);				// Add a matepair in the list.

};

#endif /* READS_H_ */