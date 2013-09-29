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

enum MatePairOrientationType{
	RevRev = 0,		// source<-----------<destination		reverse of source to reverse of destination
	RevFwd = 1,		// source<----------->destination		reverse of source to forward of destination
	FwdRev = 2,		// source>-----------<destination		forward of source to reverse of destination
	FwdFwd = 3		// source>----------->destination		forward of source to forward of destination
};

class LinkedEdgePair{
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

		// Not sure if source read ID are unique, we need an edge ID.
		// source and destination are the mate pair. The edge with smaller source ID is designated as the source
		Edge * source;
		Edge * destination;
		MatePairOrientationType orientation;
		INT32 gapDistance;
		// source and revSource are twin edges. The edge with source read ID < destination read ID is designated as the source edge.
		Edge * revSource;
		// destination and revDestination are twin edges. The edge with source read ID < destination read ID is designated as the destination edge
		Edge * revDestination;
		/*
		 *   pairedReadsInSource[i] and pairedReadsInDestination[i] are a pair of reads in source edge and destination edge, respectively
		 *   matePairCount are the number of matepairs linking the two edges
		 *   matePairCount == pairedReadsInSource.size() == pairedReadsInDestination.size()
		 */
		vector<Read *> pairedReadsInSource;
		vector<Read *> pairedReadsInDestination;
		UINT64 matePairCount;

	public:
		LinkedEdgePair();
		~LinkedEdgePair();
};

class MatePairGraph{
	private:
		vector<LinkedEdgePair> edgePairList;

		// expand the coverage
		// change EdgeMarkType of edges from Unvisited to UpstreamStop or DownstreamStop
		// return the number of edges with changed EdgeMarkType
		int expand();
		// confirm the coverage
		// change EdgeMarkType of edges from UpstreamStop to UpstreamGo or from DownstreamStop to DownstreamGo
		// return the number of edges with changed EdgeMarkType
		int lookBack();

	public:
		MatePairGraph();
		~MatePairGraph();

		// Populate edgePairList
		void buildMatePairGraph();

		/*
		 *   An edge has three mark types: Covered, Terminal, and Uncovered
		 *   Start from a set of edges marked as Covered and all the other edges marked as Uncovered
		 *   Mark all the Uncovered edges linked to Covered edges as Terminal
		 *   Mark all the Terminal edge
		 */
		void iterativelyExpandCoverage()
			{
			// expand until no edge's EdgeMarkType is changed by lookBack()
				do
				{
					expand();
				}
				while( lookBack() > 0 );
			}

};

#endif /* MATEPAIRGRAPH_H_ */
