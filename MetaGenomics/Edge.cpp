/*
 * Edge.cpp
 *
 * Created on: April 22, 2013
 * Author: Md. Bahlul Haider
 */

#include "Common.h"
#include "Edge.h"


/**********************************************************************************************************************
	Default Constructor
**********************************************************************************************************************/
Edge::Edge(void)
{

	// Initialize the variables.
	overlapOffset = 0;
	overlapOrientation = 10;
	transitiveRemovalFlag = false;
	flow = 0;
	coverageDepth = 0;
	edgeID = 0xFFFFFFFFFFFFFFFF;
/*
	listOfReads = new vector<UINT64>;
	listOfReads->resize(listOfReads->size());						// Resize to reduce space.

	listOfOverlapLengths = new vector<UINT16>;
	listOfOverlapLengths->resize(listOfOverlapLengths->size());		// Resize to reduce space.

	listOfOrientations = new vector<UINT8>;
	listOfOrientations->resize(listOfOrientations->size());			// Resize to reduce space.
*/
}



/**********************************************************************************************************************
	Another Constructor
**********************************************************************************************************************/
Edge::Edge(Read *from, Read *to, UINT64 orient, UINT64 overlapOffsetInput)
{
	makeEdge(from, to, orient, overlapOffsetInput);
}



/**********************************************************************************************************************
 	 Another Constructor
**********************************************************************************************************************/
Edge::Edge(Read *from, Read *to, UINT64 orient, UINT64 overlapOffsetInput, vector<UINT64> *listReads, vector<UINT16> *listOverlapOffsets, vector<UINT8> * listOrientations)
{
	makeEdge(from, to, orient, overlapOffsetInput, listReads, listOverlapOffsets, listOrientations);
}



/**********************************************************************************************************************
	Default destructor
**********************************************************************************************************************/
Edge::~Edge()
{
	// Free the memory used by the current edge.
	delete listOfReads;
	delete listOfOverlapOffsets;
	delete listOfOrientations;
}



/**********************************************************************************************************************
	Function to insert a simple edge
**********************************************************************************************************************/
bool Edge::makeEdge(Read *from, Read *to, UINT64 orient, UINT64 overlapOffsetInput)
{
	source = from;
	destination = to;
	overlapOrientation = orient;
	overlapOffset = overlapOffsetInput;
	hignCoverageAndMatepairFlag = false;
	edgeID = 0xFFFFFFFFFFFFFFFF;

	// Initialize variables.
	transitiveRemovalFlag = false;
	flow = 0;
	coverageDepth = 0;

	listOfReads = new vector<UINT64>;
	listOfReads->resize(listOfReads->size());						// Resize to reduce space.

	listOfOverlapOffsets = new vector<UINT16>;
	listOfOverlapOffsets->resize(listOfOverlapOffsets->size());		// Resize to reduce space.

	listOfOrientations = new vector<UINT8>;
	listOfOrientations->resize(listOfOrientations->size());			// Resize to reduce space.

	return true;
}



/**********************************************************************************************************************
	Function to add a composite edge
**********************************************************************************************************************/
bool Edge::makeEdge(Read *from, Read *to, UINT64 orient, UINT64 overlapOffsetInput,  vector<UINT64> *listReads, vector<UINT16> *listOverlapOffsets, vector<UINT8> * listOrientations)
{
	source = from;
	destination = to;
	overlapOrientation = orient;
	overlapOffset = overlapOffsetInput;
	transitiveRemovalFlag = false;
	hignCoverageAndMatepairFlag = false;
	edgeID = 0xFFFFFFFFFFFFFFFF;
	flow = 0;
	coverageDepth = 0;
	listOfReads = listReads;
	listOfOverlapOffsets = listOverlapOffsets;
	listOfOrientations = listOrientations;
	listOfReads->resize(listOfReads->size());					// Resize to reduce space.
	listOfOverlapOffsets->resize(listOfOverlapOffsets->size());	// Resize to reduce space.
	listOfOrientations->resize(listOfOrientations->size());		// Resize to reduce space.

	return true;
}



/**********************************************************************************************************************
	Function to set the pointer to the reverse edge;
**********************************************************************************************************************/
bool Edge::setReverseEdge(Edge * edge)
{
	reverseEdge = edge;
	return true;
}
