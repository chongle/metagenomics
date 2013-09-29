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
		 *   gapDistance are from the end of read v to the beginning of read x
		 *   it's possible that gapDistance is negative, meaning the two edges overlap
		 */
		Edge * fwdEdgeA;
		Edge * revEdgeA;
		Edge * fwdEdgeB;
		Edge * revEdgeB;
		INT32 gapDistance;
		/*
		 *   readsEdgeA[i] and readsEdgeB[i] are a pair of reads in edge A and edge B, respectively
		 *   matePairCount are the number of matepairs linking the two edges
		 *   matePairCount == readsEdgeA.size() == readsEdgeB.size()
		 */
		vector<Read *> readsEdgeA;
		vector<Read *> readsEdgeB;
		UINT64 matePairCount;

	public:
		LinkedEdgePair();
		~LinkedEdgePair();
};

class MatePairGraph{
	private:
		vector<LinkedEdgePair> edgePairList;

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
		void expandCoveredEdges();

};

#endif /* MATEPAIRGRAPH_H_ */
