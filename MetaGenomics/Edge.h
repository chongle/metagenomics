/*
 * Edge.h
 *
 * Created on: April 22, 2013
 * Author: Md. Bahlul Haider
 */

#ifndef EDGE_H_
#define EDGE_H_
#include "Common.h"
#include "Read.h"

/**********************************************************************************************************************
	Class to store an edge.
**********************************************************************************************************************/
class Edge{
	private:
		// CP: How is it determined which read is the source and which read is the destination?
		// BH: An edge (u,v) is stored twice in the graph. In the list of edges of the node u, the node u is the source  and the edge is (u,v)
		// Similarly in the list of edges of the node v, the node v is the source and the edge is (v,u).
		// read u and v may be in multiple edges as  source/destination or inside read
		Read *source; 							// Source read u
		Read *destination; 						// Destination read v
		UINT8 overlapOrientation;				// Orientation of overlap:	go from u to v
												// 0 = u<-----------<v		reverse of u to reverse of v
												// 1 = u<----------->v		reverse of u to forward of v
												// 2 = u>-----------<v		forward of u to reverse of v
												// 3 = u>----------->v		forward of u to forware of v
		UINT64 overlapOffset;					// Length of the overlap.
												// overlap offset in the following example is 6
												// 	  012345678901234567890123456789
												//	u ACTTACGGGATTATACCATCGAGA
												//	v       GGGATTATACCATCGAGATTCAAT
												// CP: what about a composite edge? from the beginning of the first read to the beginning of the last read?
												// BH: Yes. Froom the beginning of the first read, u, to the beginning of the last read, v, for composite edge.

		vector<UINT64> * listOfReads; 			// List of ordered reads in the current edge. NOT including u and v.
		vector<UINT16> * listOfOverlapOffsets; 	// List of overlap offsets of the ordered reads in the current edge.
												// CP: this is relative to the intermediate previous read?
												// BH: Yes this is the offset from the previous read in the list.
												// The offset of the destination read is calculated as overlapOffset - sum(listOfOverlapOffsets[i]) or it can retrieved from the reverse edge
		vector<UINT8> * listOfOrientations;		// List of orientations of the ordered reads in the current edge.
												// orientation is 0 or 1. 0 means read's reverse complement read.getStringReverse(). 1 means read's forward string.
		Edge *reverseEdge;						// This points to the reverse edge in the overlap graph.
												// Each edge (u,v) is present twice in the graph.
												// Edge (u,v) is present in the list of edges of node u
												// Edge (v,u) is present in the list of edges of node v.
												// Edge (u,v) and (v,u) are called twin edge and are reverse of one another.'
												// the reverseEdge always have the same flow as the forward Edge
												// the listOfReads of the reverse edge is the same reads in reverse order as the forward edge
												// the listOfOrientations of the reverse edge is reverse complement


	public:
		bool transitiveRemovalFlag;							// Used to mark transitive edges.
		UINT16 flow;							// Store the flow in the current edge.
		UINT64 coverageDepth;					// Estimated depth of coverage.
		UINT64 SD;
		Edge(void);								// Default constructor.
		Edge(Read *from, Read *to, UINT64 orient, UINT64 overlapOffsetInput); 	// constructor for simple edge, length is the overlap length (i.e. length = 18 in the example above)
		Edge(Read *from, Read *to, UINT64 orient, UINT64 overlapOffsetInput, vector<UINT64> *listReads, vector<UINT16> *listOverlapOffsets, vector<UINT8> * listOrientations);  // constructor for composite edge
		~Edge();								// Destructor.
		bool makeEdge(Read *from, Read *to, UINT64 orient, UINT64 overlapOffsetInput);		// make a simple edge
		bool makeEdge(Read *from, Read *to, UINT64 orient, UINT64 overlapOffsetInput,  vector<UINT64> *listReads, vector<UINT16> *listOverlapOffsets, vector<UINT8> * listOrientations);
		string getStringInEdge(void); 			// return the string in the edge
		bool setReverseEdge(Edge * edge);		// Set the pointer to the reverse edge.
		Read * getSourceRead() {return source;}	// Get the read object of the source node.
		Read * getDestinationRead() {return destination; }	// Get the read object of the destination node.
		UINT8 getOrientation() {return overlapOrientation;}	// Return the orientation of the edge.
		UINT64 getOverlapOffset() {return overlapOffset;}	// Return the overlap offset.
		vector<UINT64> * getListOfReads() {return listOfReads;}		// Get the ordered list of reads in the current edge.
		vector<UINT16> * getListOfOverlapOffsets() {return listOfOverlapOffsets;} // Get the list of ordered offset.
		vector<UINT8> * getListOfOrientations() {return listOfOrientations;}	// Get the ordered orientation of the reads. 1 means forward. 0 means reverse.
		Edge * getReverseEdge() {return reverseEdge;}	// Return the pointer to the reverse edge.
};

#endif /* EDGE_H_ */
