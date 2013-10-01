/*
 * MatePairGraph.cpp
 *
 *  Created on: Sep 30, 2013
 *      Author: b72
 */

#include "Common.h"
#include "MatePairGraph.h"


MatePairGraph::MatePairGraph()
{


}
MatePairGraph::~MatePairGraph()
{
	delete linkList;
}

void MatePairGraph::buildMatePairGraph(Dataset * dataSet, OverlapGraph * overlapGraph)
{
	CLOCKSTART;
	vector< vector<Edge *> * > *graph = overlapGraph->getGraph();
	UINT64 ID = 1;
	Edge * edge = new Edge(); // empty edge at 0th location.
	listOfCompositeEdges.push_back(edge);
	// Get the list of composite edges.
	for(UINT64 i = 1; i < graph->size(); i++) // For each node
	{
		for(UINT64 j = 0 ; j < graph->at(i)->size(); j++) // for each edge
		{
			Edge *edge = graph->at(i)->at(j);
			UINT64 u = edge->getSourceRead()->getReadNumber(), v = edge->getDestinationRead()->getReadNumber();
			if(!edge->getListOfReads()->empty() && ( u < v || (u==v && edge < edge->getReverseEdge())) ) // composite edge
			{
				edge->setEdgeID(ID); // Forward edge to positive ID
				edge->getReverseEdge()->setEdgeID(-ID); // Reverse edge to negative ID
				listOfCompositeEdges.push_back(edge);	// List of all the composite edges.
				ID++;
			}
		}
	}
	cout << "Total Composite Edges: " << ID-1 << endl;
	linkList = new vector< vector<MatePairLinks> >;
	for(UINT64 i = 0; i< ID; i++)
	{
		vector<MatePairLinks> a;
		linkList->push_back(a);
	}

	for(UINT64 i = 1 ; i < listOfCompositeEdges.size(); i++) // For each composite edge in the graph
	{
		Edge *edge1 = listOfCompositeEdges.at(i);
		// Find the list of other edges that share unique matepairs with the current edge. Only check one of the endpoints of the current edge.
		vector<Edge *> *listOfFeasibleEdges = overlapGraph->getListOfFeasibleEdges(edge1);
		for(UINT64 j = 0; j < listOfFeasibleEdges->size(); j++ ) // Check the current edge vs the list of edges for suppor for scaffolder
		{


			Edge *edge2 =listOfFeasibleEdges->at(j);
			INT64 gapDistance;
			UINT64 support = 0;
			support = overlapGraph->checkForScaffold(edge1,edge2,&gapDistance); // check the support and distance of the forward edge
			if(support > 0)
			{
				INT64 ID2 = edge2->getEdgeID();
				MatePairLinks mpLinke1e2, mpLinke2e1;
				MatePairOrientationType oriente1e2;
				if(ID2 > 0)
				{
					oriente1e2 = FwdFwd;
				}
				else
				{
					oriente1e2 = FwdRev;
				}
				if(edge2->getEdgeID() < 0)
					edge2=edge2->getReverseEdge();
				mpLinke1e2.setVariables(edge1, edge2, oriente1e2, support, gapDistance);
				addLink(mpLinke1e2);
			}
		}
		delete listOfFeasibleEdges;	// Free the memory

		Edge *edge1Reverse = listOfCompositeEdges.at(i)->getReverseEdge();
		// Find the list of other edges that share unique matepairs with the current edge. Only check one of the endpoints of the current edge.
		vector<Edge *> *listOfFeasibleEdgesReverse = overlapGraph->getListOfFeasibleEdges(edge1Reverse);
		for(UINT64 j = 0; j < listOfFeasibleEdgesReverse->size(); j++ ) // Check the current edge vs the list of edges for suppor for scaffolder
		{
			Edge *edge2 =listOfFeasibleEdgesReverse->at(j);
			INT64 gapDistance;
			UINT64 support = 0;
			support = overlapGraph->checkForScaffold(edge1Reverse,edge2,&gapDistance); // check the support and distance of the reverse edge
			if(support > 0)
			{
				INT64 ID2 = edge2->getEdgeID();
				MatePairLinks mpLinke1e2, mpLinke2e1;
				MatePairOrientationType oriente1e2;
				if(ID2 > 0)
				{
					oriente1e2 = RevFwd;
				}
				else
				{
					oriente1e2 = RevRev;
				}
				if(edge2->getEdgeID() < 0)
					edge2=edge2->getReverseEdge();
				mpLinke1e2.setVariables(edge1Reverse->getReverseEdge(), edge2, oriente1e2, support, gapDistance);
				addLink(mpLinke1e2);
			}
		}
		delete listOfFeasibleEdgesReverse;	// Free the memory
	}
	CLOCKSTOP;
}

void MatePairGraph::addLink(MatePairLinks mpLink)
{
	INT64 ID =	mpLink.getSourceEdge()->getEdgeID();
	if(ID < 0) ID = -ID;
	linkList->at(ID).push_back(mpLink);
}



// Mark the transitive links.
void MatePairGraph::markTransitiveEdge()
{
	for(UINT64 i = 1; i < linkList->size(); i++) // For each composite edge e
	{
		if(linkList->at(i).size() > 0) // e has some linked edges.
		{
			for(UINT64 j = 0; j < linkList->at(i).size(); j++) // for the first linked edge e1
			{
				Edge *destinationEdge1 = linkList->at(i).at(j).getDestinationEdge();
				INT64 destinationEdge1ID = destinationEdge1->getEdgeID();
				MatePairOrientationType orient1 = linkList->at(i).at(j).getOrient();
				for(UINT64 k = 0; k < linkList->at(i).size(); k++)		// for the second linked edge e2
				{
					if(j==k) // if e1== e2
						continue;
					Edge *destinationEdge2 = linkList->at(i).at(k).getDestinationEdge();
					INT64 destinationEdge2ID = destinationEdge2->getEdgeID();
					MatePairOrientationType orient2 = linkList->at(i).at(k).getOrient();

					for(UINT64 l = 0; l <linkList->at(destinationEdge1ID).size() ; l++) // check for a link between e1 and e2
					{
						Edge *destinationEdge3 = linkList->at(destinationEdge1ID).at(l).getDestinationEdge();
						INT64 destinationEdge3ID = destinationEdge3->getEdgeID();
						MatePairOrientationType orient3 = linkList->at(destinationEdge1ID).at(l).getOrient();
						if(destinationEdge2ID == destinationEdge3ID) // we can go to same edge from another path;
						// TODO: We need to also check if we can to to that edge by following a node in the overlap graph.
						{
							// (orient1 & 1) == ((orient2&2)>>1) checks if
							// e->e1    e->e2
							// *Fwd and Fwd*
							// *Rev and Rev*

							// (orient1 & 2) | (orient2 & 1)) == orient3) checks
							// e->e1  e->e2  e1->e2
							// Fwd* + *Fwd = FwdFwd
							// Fwd* + *Rev = FwdRev
							// Rev* + *Fwd = RevFwd
							// Rev* + *Rev = RevRev
							if( ((orient1 & 1) == ((orient2&2)>>1) ) && (((orient1 & 2) | (orient2 & 1)) == orient3))
							{
								linkList->at(destinationEdge1ID).at(l).isTransitive = true;
								//cout << "Transitive Edge Found" << endl;
							}

						}
					}
				}
			}
		}
	}

}



// Mark the transitive links.
void MatePairGraph::markEdgesByMatePairs()
{
	CLOCKSTART;
	markTransitiveEdge();

	for(UINT64 i = 1; i < linkList->size(); i++) // For each composite edge e
	{
		if(linkList->at(i).size() > 0) // if the edge e has some matepair linkage
		{
			MatePairLinks mpLink = linkList->at(i).at(0);
			if(mpLink.getSourceEdge()->coverageDepth >=coverageDepthLB && mpLink.getSourceEdge()->coverageDepth <=coverageDepthUB) // If the coverage depth of e is within desired range. This edges are already marked for flow lb 1. We want to mark it unambiguous linked edges.
			{
				UINT64 fwdEdges = 0, revEdges=0;
				Edge * forwardLinkedEdge, *reverseLinkedEdge;
				for(UINT64 j = 0; j < linkList->at(i).size(); j++)
				{
					MatePairOrientationType orient = linkList->at(i).at(j).getOrient();
					bool isTransitive = linkList->at(i).at(j).isTransitive;
					if(!isTransitive)
					{
						if((orient & 2) == 2) // Count the forward edge links
						{
							forwardLinkedEdge = linkList->at(i).at(j).getDestinationEdge();
							fwdEdges++;
						}
						else	// count the reverse edge links
						{
							reverseLinkedEdge = linkList->at(i).at(j).getDestinationEdge();
							revEdges++;
						}
					}
				}
				if(fwdEdges == 1 && forwardLinkedEdge->hignCoverageAndMatepairFlag == false) // If only one link to the forward edge.
				{
					forwardLinkedEdge->hignCoverageAndMatepairFlag = true;			// Mark it for lower bound 1
					forwardLinkedEdge->getReverseEdge()->hignCoverageAndMatepairFlag = true;	// Mark the reverse edge for lower bound 1.
					cout << "Marking Edge Forward: (" <<forwardLinkedEdge->getSourceRead()->getReadNumber() << "," << forwardLinkedEdge->getDestinationRead()->getReadNumber() << ")" << endl;

				}
				if(revEdges == 1 && reverseLinkedEdge->hignCoverageAndMatepairFlag == false)	// if only one link to the reverse edge.
				{
					reverseLinkedEdge->hignCoverageAndMatepairFlag = true;						// Mark it for lower bound 1
					reverseLinkedEdge->getReverseEdge()->hignCoverageAndMatepairFlag = true; 	// Mark the reverse edge for lower bound 1.
					cout << "Marking Edge Reverse: (" <<reverseLinkedEdge->getSourceRead()->getReadNumber() << "," << reverseLinkedEdge->getDestinationRead()->getReadNumber() << ")" << endl;
				}
			}
		}
	}
	CLOCKSTOP;
}

void MatePairGraph::printMatePairLinkageGraph()
{
	CLOCKSTART;
	for(UINT64 i = 1; i < linkList->size(); i++)
	{
		cout<< "EDGE: " << i << endl;
		cout << "=======================================" << endl;
		for(UINT64 j = 0 ; j < linkList->at(i).size(); j++)
		{
			MatePairLinks currentLink = linkList->at(i).at(j);
			cout << "Edges1: (" <<currentLink.getSourceEdge()->getSourceRead()->getReadNumber() << "," << currentLink.getSourceEdge()->getDestinationRead()->getReadNumber() << ")" << endl;
			cout << "Edge1 ID: " << currentLink.getSourceEdge()->getEdgeID() << endl;
			cout << "Edge1 OverlapOffset: " << currentLink.getSourceEdge()->getOverlapOffset() << endl;
			cout << "Reads in Edge1: " << currentLink.getSourceEdge()->getListOfReads()->size() << endl;
			cout << "Edges2: (" <<currentLink.getDestinationEdge()->getSourceRead()->getReadNumber() << "," << currentLink.getDestinationEdge()->getDestinationRead()->getReadNumber() << ")" << endl;
			cout << "Edge2 ID: " << currentLink.getDestinationEdge()->getEdgeID() << endl;
			cout << "Edge2 OverlapOffset: " << currentLink.getDestinationEdge()->getOverlapOffset() << endl;
			cout << "Reads in Edge2: " << currentLink.getDestinationEdge()->getListOfReads()->size() << endl;
			cout << "Support: " << currentLink.getSupport() << endl;
			cout << "isTransitie: " << currentLink.isTransitive << endl;
			cout << "Average gap distance: " << currentLink.getAverageGapDistance() << endl;
			cout << "Type:" << currentLink.getOrient() << endl << endl;
		}
	}
	CLOCKSTOP;
}




