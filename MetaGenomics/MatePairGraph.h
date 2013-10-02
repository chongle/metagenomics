/*
 * MatePairGraph.h
 *
 *  Created on: Sep 28, 2013
 *      Author: test
 */

#ifndef MATEPAIRGRAPH_H_
#define MATEPAIRGRAPH_H_

#include "Common.h"
#include "Read.h"
#include "Edge.h"
#include "OverlapGraph.h"

//forward declaration
class OverlapGraph;

enum MatePairOrientationType{
	RevRev = 0,		// source<-----------<destination		reverse of source to reverse of destination
	RevFwd = 1,		// source<----------->destination		reverse of source to forward of destination
	FwdRev = 2,		// source>-----------<destination		forward of source to reverse of destination
	FwdFwd = 3		// source>----------->destination		forward of source to forward of destination
};

class MatePairLinks{
	private:
		/*
		 *   Mate pairs are in forward - reverse orientation
		 *   fwdEdgeA and revEdgeA are twin edges from read u to read v
		 *   fwdEdgeB and revEdgeB are twin edges from read x to read y
		 *   fwdEdgeA is in front of fwdEdgeB in sequence
		 *   revEdgeB is in front of revEdgeA in sequence
		 *   u-----fwdEdgeA----v	+	x-----fwdEdgeB----y
		 *   u-----revEdgeA----v	+	x-----revEdgeB----y
		 *
		 *   gapDistance are from the end of read v to the beginning of read x
		 *   it's possible that gapDistance is negative, meaning the two edges overlap
		 */

		// source and destination are the linked edges.
		Edge * source;
		Edge * destination;
		MatePairOrientationType orientation;
		INT64 averageGapDistance;
		// source and revSource are twin edges. The edge with source read ID < destination read ID is designated as the source edge.
		Edge * revSource;
		// destination and revDestination are twin edges. The edge with source read ID < destination read ID is designated as the destination edge
		Edge * revDestination;

		/*
		 *   pairedReadsInSource[i] and pairedReadsInDestination[i] are a pair of reads in source edge and destination edge, respectively
		 *   matePairCount are the number of matepairs linking the two edges
		 *   matePairCount == pairedReadsInSource.size() == pairedReadsInDestination.size()
		 */
		UINT64 matePairCount;
		vector<Read *> *pairedReadsInSource;
		vector<Read *> *pairedReadsInDestination;
		vector<INT64> *gapDistance;		// gap distance according to this pair of reads
		vector<bool> areUniqueReads;	// are both reads in this pair unique to their respective edges

		// the two linked edges have a feasible path between them within the range of gap distance
		// Not sure if we need to check this or not.
		// can use the exploreGraph function in the Overlap graph. Simply give the last read
		// UINT64 exploreGraph(Edge* firstEdge, Edge * lastEdge, UINT64 distanceOnFirstEdge, UINT64 distanceOnLastEdge,
		// UINT64 datasetNumber, UINT64 level, vector <Edge *> &firstPath, vector <UINT64> &flags);
		bool hasFeasiblePath;

		// If the two edges have a small overlap and a small gap distance, then true
		// this is similar to a step in scaffolding
		// not sure if this is needed
		bool canBeMerged;


	public:
		bool isTransitive;
		MatePairLinks(){isTransitive = false;};
		~MatePairLinks(){}
		void setVariables(Edge *edge1, Edge *edge2, MatePairOrientationType orient, UINT64 support, INT64 averageGapDistance, vector<Read *> *pairedReadsInSource, vector<Read *> *pairedReadsInDestination, vector<INT64> *gapDistance)
		{source = edge1; revSource = edge1->getReverseEdge() ; destination = edge2; revDestination = edge2->getReverseEdge(); orientation = orient; matePairCount=support;this->averageGapDistance=averageGapDistance; this->pairedReadsInSource = pairedReadsInSource; this->pairedReadsInDestination = pairedReadsInDestination; this->gapDistance = gapDistance;}
		Edge * getSourceEdge(){return source;}
		Edge * getDestinationEdge(){return destination;}
		MatePairOrientationType getOrient() {return orientation;}
		UINT64 getSupport(){return matePairCount;}
		INT64 getAverageGapDistance(){return averageGapDistance;}
		vector<Read *> * getPairedReadsInSource() {return pairedReadsInSource;}
		vector<Read *> * getPairedReadsInDestination() {return pairedReadsInDestination;}
		vector<INT64> * getGapDistance() {return gapDistance;}

};

class MatePairGraph{
	private:
		vector< vector<MatePairLinks> > *linkList;
		vector< Edge *> listOfCompositeEdges;

		void addLink(MatePairLinks mpLink);
	//	vector<Edge *> * getListOfFeasibleEdges(const Edge *edge);

		void markTransitiveEdge();

	public:
		MatePairGraph();
		~MatePairGraph();

		// Populate linkList
		// Ideally this should work for both unitig graph and contig graph
		void buildMatePairGraph(Dataset * dataSet, OverlapGraph * overlapGraph);


		// Mark additional edges linked to already marked edges
		void markEdgesByMatePairs();
		void printMatePairLinkageGraph();
};

#endif /* MATEPAIRGRAPH_H_ */
