/*
 * OverlapGraph.h
 *
 * Created on: April 22, 2013
 * Author: Md. Bahlul Haider
 */


#ifndef OVERLAPGRAPH_H_
#define OVERLAPGRAPH_H_
#include "Common.h"
#include "Dataset.h"
#include "HashTable.h"
#include "MatePairGraph.h"
#include "Edge.h"

/**********************************************************************************************************************
	Class to store the overlap graph.
**********************************************************************************************************************/

enum nodeType {
	// Current node u is not explored yet. Meaning that there is no edge (u,v) in the graph.
	UNEXPLORED = 0,
	//  Current node u is explored. Meaning that all edges (u,v) are inserted in the dataset.
	EXPLORED = 1,
	// Meaning that all transitive edges (u,v) of current node u is marked and its neighbors transitive edges are also marked.
	// Now it is safe to remove any transitive edge (u,v) from node u.
	EXPLORED_AND_TRANSITIVE_EDGES_MARKED = 2
};

enum markType{
	// This type is used when we discover transitive edges.
	VACANT = 0,	// Initially all nodes are marked as VACANT.
	INPLAY = 1,	// then all the neighbors of the current node is marked as INPLAY
	ELIMINATED = 2	// Then neighbors' (marked as INPLAY) neighbors (marked as INPLAY) are marked as ELIMINATED
					// At the end all edges to the neighbor marked as ELIMINATED are marked as transitive edges and will be removed.
};

// This structure is used to store list of pair of edges and their support. Used in two function: 1. when we find path by mate-pairs 2. scaffolder.
struct pairedEdges
{

		/*
		 *  edge 1 is in front of edge 2 in the sequence
		 *  CP2: is this right? the reverse of edge2 are in the same orientation as edge1??
		 *  -----edge1----		-------------
		 *  --------------		-----edge2---
		 */
		Edge * edge1;
		Edge * edge2;
		UINT64 support;			// number of matepairs supporting these two edges
		UINT64 distance;		// sum of the end of read1 to the end of edge 1 and the beginning of read2 to the beginning of edge2
		bool isFreed;
		bool operator < (const pairedEdges& rhs) const
		{
		       return support > rhs.support;
		}

};

class OverlapGraph
{
	private:
		Dataset * dataSet; 											// Pointer to the dataset containing all the reads. this is NOT modified here

		UINT16 hashStringLength;									// length of hash string, got it from hashTable.getHashStringLength()

		// CP: the overlap graph is a vector of nodes (i.e. reads) with their vectors of incident edges
		// CP: the readNumber of a read is equal to the index of the node in this vector. Because readNumber starts from 1, the
		// CP2: how to iterate through all edges??
		vector< vector<Edge *> * > *graph;							// Adjacency list of the graph.
		UINT64 numberOfNodes;										// Number of nodes in the overlap graph.
		UINT64 numberOfEdges;										// Number of edges in the overlap graph.
		bool flowComputed;											// Flag to check wheather the flow is computed or not.

		// the vector corresponds to a list of datasets. These manages the insert sizes of different datasets
		vector<INT64> meanOfInsertSizes; 							// Mean of insert sizes. index is the datasetNumber
		vector<INT64> sdOfInsertSizes; 							// Standard deviation of insert sizes.
		INT64 longestMeanOfInsertSize;								// CP: the longest mean insert size out of all datasets: max(meanOfInsertSizes[i])
		UINT64 getMean(UINT64 datasetNumber) {return meanOfInsertSizes.at(datasetNumber);} 			// Get the mean of insert sizes.
		UINT64 getSD(UINT64 datasetNumber) {return sdOfInsertSizes.at(datasetNumber);} 				// Get standard deviation of insert sizes

		// Genome size estimation. Only for testing
		UINT64 estimatedGenomeSize;									// Estimated genome size. Works for isolated genome. Will not work for Metagenomics.
		UINT64 getEstimatedGenomeSize() {return estimatedGenomeSize;} 	// Get the estimated genome size.
		bool estimateGenomeSize(void);								// Estimate the genome size. Works only for isolated genome. Does not work for Metagenomics.

		// functions that used hash table to find overlaps and contained reads. they are called in buildOverlapGraphFromHashTable
		bool insertAllEdgesOfRead(const HashTable & hashTable, UINT64 readNumber, vector<nodeType> * exploredReads);	// Insert into the overlap graph all edges of a read.
		void markContainedReads(const HashTable & hashTable);								// Find superReads for each read and mark them as contained read.
		bool checkOverlap(const Read *read1, const Read *read2, UINT64 orient, UINT64 start); // Check overlap between two reads after a match is found using the hash table.
		bool checkOverlapForContainedRead(const Read *read1, const Read *read2, UINT64 orient, UINT64 start);
		UINT64 findOverlap(const string & string1, const string & string2);			// Find overlap length between two strings.

		// functions used in building the overlap graph
		bool markTransitiveEdges(UINT64 readNumber, vector<markType> * markedNodes); // Mark transitive edges of a read.
		bool removeTransitiveEdges(UINT64 readNumber);				// Remove all transitive edges from the overlap graph incident to a given read.
		bool insertEdge(Edge * edge); 								// Insert an edge in the overlap graph.
		bool insertEdge(Read *read1, Read *read2,  UINT8 orient, UINT16 overlapOffset); // Insert an edge in the overlap graph.

		// functions for merging two edges
		bool mergeEdges(Edge *edge1, Edge *edge2);					// Merge two edges in the  overlap graph.
		UINT8 mergedEdgeOrientation(const Edge *edge1, const Edge *edge2);		// Orientation of the edge when two edges are merged.
		UINT8 twinEdgeOrientation(UINT8 orientation);				// Orientation of the reverse edge.
		// When we merge two edges, we need to merge the list of ordered reads, their overlap offsets and orientations.
		bool mergeList(const Edge *edge1, const Edge *edge2, vector<UINT64> *listReads, vector<UINT16> *listOverlapOffsets, vector<UINT8> * ListOrientations);
		bool removeEdge(Edge *edge); 						// Remove an edge from the overlap graph.
		bool updateReadLocations(Edge *edge);				// Update the location of all the reads in the current edge. This function is called when a new edge is inserted.
		bool removeReadLocations(Edge *edge);				// Remove the location of all the reads from the current edge. This function is called when an edge is removed.

		Edge *findEdge(UINT64 source, UINT64 destination);			// Find an edge from source to destination in the overlap graph.
		bool isEdgePresent(UINT64 source, UINT64 destination);		// Check if an edge is present in the overlap graph between source and destination.
		bool mergeEdgesDisconnected(Edge *edge1, Edge *edge2, UINT64 gapLength); // Merge to edges that do not share any node.
		UINT8 mergedEdgeOrientationDisconnected(const Edge *edge1, const Edge *edge2);		// Orientation of the edge when two disconnected edge are merged.
		bool mergeListDisconnected(Edge *edge1, Edge *edge2, UINT64 overlapOffset,
				UINT64 gapLength, vector<UINT64> *listReads, vector<UINT16> *listOverlaps, vector<UINT8> * listOrientations);

		// Use paired end information to merge and scaffold
		// Find support between to reads in the the overlap graph. Find all pathas between two reads in the graph.
		// Only the subpath that are common in all such paths are considered to be supported.
		bool findPathBetweenMatepairs(const Read * read1, const Read * read2, UINT8 orient, UINT8 datasetNumbe, vector <Edge *> &copyOfPath, vector <UINT64> &copyOfFlags);
		// Starting from the firstEdge we try to reach lastEdge by distance mean + 3 *sd of the datasetNumber.
		UINT64 exploreGraph(Edge* firstEdge, Edge * lastEdge, UINT64 distanceOnFirstEdge, UINT64 distanceOnLastEdge,
				UINT64 datasetNumber, UINT64 level, vector <Edge *> &firstPath, vector <UINT64> &flags);



		// Functions related to flow analysis
		bool calculateBoundAndCost(const Edge *edge, INT64* FLOWLB, INT64* FLOWUB, INT64* COST); // Calculate bounds and costs of flow for minimum cost flow in the overlap graph.
		bool calculateBoundAndCost2(Edge *edge, INT64* FLOWLB, INT64* FLOWUB, INT64* COST); // Calculate bounds and costs of flow for minimum cost flow in the overlap graph.
		void markEdgeThatWillHaveOneFlow(void);

		// functions for simplifying the graph. Used in simplifyGraph()
		UINT64 removeDeadEndNodes(void); 					        // Remove dead-ends from the overlap graph.
		UINT64 contractCompositePaths(void); 						// Contract composite paths in the overlap graph.
		UINT64 removeSimilarEdges(void);							// remove multi-edges with similar strings
		UINT64 reduceTrees(void);									// Remove trees in the overlap graph.
		UINT64 reduceLoops(void);									// loops that can be traversed only one way
		UINT64 calculateEditDistance(const string & s1, const string & s2);	// Find the edit distance between two strings. Used by removeSimilarEdges()

		void getBaseByBaseCoverage(Edge *edge);						// Get the coverage Mean and SD of an edge. Only considering the unique reads.

		string getStringInEdge(const Edge *edge);							// Get the string in an edge by overlapping the ordered reads in the edge.

	public:

		OverlapGraph(void);											// Default constructor.
		~OverlapGraph();											// Destructor.

		// Initialize the overlap using the dataset and the hash table
		bool setDataset(Dataset *dataset){dataSet=dataset;  return true;}	// Set the dataset pointer.
		bool buildOverlapGraphFromHashTable(const HashTable & hashTable);			// Build the overlap graph using hashtable.

		// Save and reload the overlap graph
		void sortEdges();						 // Sort edges of each read based on ID of the destination read.
		bool readGraphFromFile(string fileName); // Read the graph from a file stored before.
		bool saveGraphToFile(string fileName);	 // Save the unitig graph to a file. This file can be used to reproduced the graph. Useful to avoid recomputing the graph.
		bool printGraph(string graphFileName, string contigFileName);	// Store the overlap graph for visual display and also store the contigs/scaffods in a file.

		bool calculateMeanAndSdOfInsertSize(void); 								// Estimate the mean and standard deviation of insert sizes.
		UINT64 getNumberOfEdges(void){return numberOfEdges;}		// Get the number of edges in the overlap graph.
		UINT64 getNumberOfNodes(void){return numberOfNodes;}		// Get the number of nodes in the overlap graph.

		bool simplifyGraph(void);									// Some simple simplification.
		UINT64 findSupportByMatepairsAndMerge(void);				// Using matepair information find which pair of edges are supported.
		UINT64 scaffolder(void);									// Construct scaffolds using matepair information.
		UINT64 scaffolderTemp(void);									// Construct scaffolds using matepair information.
		UINT64 resolveNodes(void);									// Merge edges based on coverage depth

		bool calculateFlow(string inputFileName, string outputFileName);									// Calculate the minimum cost flow of the overlap graph.
		bool calculateFlow2(string inputFileName, string outputFileName);


		vector<Edge *> * getListOfFeasibleEdges(const Edge *edge);
		UINT64 checkForScaffold(const Edge *edge1, const Edge *edge2, UINT64 *distance);

		UINT64 removeAllSimpleEdgesWithoutFlow();						// Not used. Only for testing.
		UINT64 removeAllEdgesWithoutFlow();						// Not used. Only for testing.
		vector< vector<Edge *> * > * getGraph() {return graph;};

};



#endif /* OVERLAPGRAPH_H_ */