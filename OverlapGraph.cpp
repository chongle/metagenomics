/*
 * OverlapGraph.cpp
 *
 * Created on: April 22, 2013
 * Author: Md. Bahlul Haider
 */

#include "Common.h"
#include "OverlapGraph.h"


/**********************************************************************************************************************
	Default Constructor
**********************************************************************************************************************/
OverlapGraph::OverlapGraph(void)
{
	// Initialize the variables.
	estimatedGenomeSize = 0;
	numberOfNodes = 0;
	numberOfEdges = 0;
	flowComputed = false;

}



/**********************************************************************************************************************
	Another Constructor. Build the overlap grpah using the hash table.
**********************************************************************************************************************/
OverlapGraph::OverlapGraph(HashTable *ht)
{
	// Initialize the variables.
	estimatedGenomeSize = 0;
	numberOfNodes = 0;
	numberOfEdges = 0;
	flowComputed = false;
//	removeContainedReadsFromHashTable(ht);
	buildOverlapGraphFromHashTable(ht);

}



/**********************************************************************************************************************
	Default destructor.
**********************************************************************************************************************/
OverlapGraph::~OverlapGraph()
{
	// Free the memory used by the overlap graph.
	for(UINT64 i = 0; i < graph->size(); i++)
	{
		for(UINT64 j = 0; j< graph->at(i)->size(); j++)
		{
			delete graph->at(i)->at(j);
		}
		delete graph->at(i);
	}
	delete graph;
}
/**********************************************************************************************************************
	Contained reads removal from hash table
**********************************************************************************************************************/
bool OverlapGraph::removeContainedReadsFromHashTable(HashTable *ht)
{
	this->totaledgenumber = 0;
	filePointer.open("cleantest.fasta");
	CLOCKSTART;
	estimatedGenomeSize = 0;
	numberOfNodes = 0;
	numberOfEdges = 0;
	flowComputed = false;
	hashTable = ht;
	dataSet = ht->getDataset();
	UINT64 counter = 0;


	graph = new vector< vector<Edge *> * >;
	graph->reserve(dataSet->getNumberOfUniqueReads()+1);



	for(UINT64 i = 0; i <= dataSet->getNumberOfUniqueReads(); i++) // Initialization
	{
		vector<Edge *> *newList = new vector<Edge *>;
		graph->push_back(newList);

	}

	markContainedReads();


	filePointer.close();
//YAO	delete hashTable;	// Do not need the hash table any more.

	CLOCKSTOP;
	return true;
}


/**********************************************************************************************************************
	Build the overlap graph from hash table
**********************************************************************************************************************/
bool OverlapGraph::buildOverlapGraphFromHashTable(HashTable *ht)
{
	this->totaledgenumber = 0;
	filePointer.open("omegatemp.txt");

	estimatedGenomeSize = 0;
	numberOfNodes = 0;
	numberOfEdges = 0;
	flowComputed = false;
	hashTable = ht;
	dataSet = ht->getDataset();
	UINT64 counter = 0;
	vector<nodeType> *exploredReads = new vector<nodeType>;
	exploredReads->reserve(dataSet->getNumberOfUniqueReads()+1);

	vector<UINT64> * queue = new vector<UINT64>;
	queue->reserve(dataSet->getNumberOfUniqueReads()+1);

	vector<markType> *markedNodes = new vector<markType>;
	markedNodes->reserve(dataSet->getNumberOfUniqueReads()+1);

	graph = new vector< vector<Edge *> * >;
	graph->reserve(dataSet->getNumberOfUniqueReads()+1);



	for(UINT64 i = 0; i <= dataSet->getNumberOfUniqueReads(); i++) // Initialization
	{
		vector<Edge *> *newList = new vector<Edge *>;
		graph->push_back(newList);
		exploredReads->push_back(UNEXPLORED);
		queue->push_back(0);
		markedNodes->push_back(VACANT);
	}

//YAO	markContainedReads();

//YAO	this->dataSet->readMatePairsFromFile();

	for(UINT64 i = 1; i <= dataSet->getNumberOfUniqueReads(); i++)
	{
		if(exploredReads->at(i) == UNEXPLORED)
		{
			UINT64 start = 0, end = 0; 											// Initialize queue start and end.
			queue->at(end++) = i;
			while(start < end) 													// This loop will explore all connected component starting from read i.
			{
				counter++;
				UINT64 read1 = queue->at(start++);
				if(exploredReads->at(read1) == UNEXPLORED)
				{
					insertAllEdgesOfRead(read1, exploredReads);					// Explore current node.
					exploredReads->at(read1) = EXPLORED;
				}
				if(graph->at(read1)->size() != 0) 								// Read has some edges (required only for the first read when a new queue starts.
				{
					if(exploredReads->at(read1) == EXPLORED) 					// Explore unexplored neighbors first.
					{
						for(UINT64 index1 = 0; index1 < graph->at(read1)->size(); index1++ )
						{
							UINT64 read2 = graph->at(read1)->at(index1)->getDestinationRead()->getReadNumber();
							if(exploredReads->at(read2) == UNEXPLORED) 			// Not explored.
							{
								queue->at(end++) = read2; 						// Put in the queue.
								insertAllEdgesOfRead(read2, exploredReads);
								exploredReads->at(read2) = EXPLORED;
							}
						}
//YAO						markTransitiveEdges(read1, markedNodes); // Mark transitive edges
						exploredReads->at(read1) = EXPLORED_AND_TRANSITIVE_EDGES_MARKED;
					}
/*					if(exploredReads->at(read1) == EXPLORED_AND_TRANSITIVE_EDGES_MARKED)
					{
						for(UINT64 index1 = 0; index1 < graph->at(read1)->size(); index1++) 				// Then explore all neighbour's neighbors
						{
							UINT64 read2 = graph->at(read1)->at(index1)->getDestinationRead()->getReadNumber();
							if(exploredReads->at(read2) == EXPLORED)
							{
								for(UINT64 index2 = 0; index2 < graph->at(read2)->size(); index2++) 		// Explore all neighbors neighbors
								{
									UINT64 read3 = graph->at(read2)->at(index2)->getDestinationRead()->getReadNumber();
									if(exploredReads->at(read3) == UNEXPLORED) 				// Not explored
									{
										queue->at(end++) = read3; 					// Put in the queue
										insertAllEdgesOfRead(read3, exploredReads);
										exploredReads->at(read3) = EXPLORED;
									}
								}
//YAO								markTransitiveEdges(read2, markedNodes); // Mark transitive edge
								exploredReads->at(read2) = EXPLORED_AND_TRANSITIVE_EDGES_MARKED;
							}
						}
//YAO						removeTransitiveEdges(read1); // Remove the transitive edges
					}*/
				}
				if(counter%100000==0)	// Show the progress.
					cout<<"counter: " << setw(10) << counter << " Nodes: " << setw(10) << numberOfNodes << " Edges: " << setw(10) << numberOfEdges/2 << endl;
			}
		}
	}
	cout<<"total counter: " << setw(10) << counter << " Nodes: " << setw(10) << numberOfNodes << " Edges: " << setw(10) << numberOfEdges/2 << endl;
	cout<<"total edge number: "<<this->totaledgenumber<<endl;
	delete exploredReads;
	delete queue;
	delete markedNodes;
	filePointer.close();
//YAO	delete hashTable;	// Do not need the hash table any more.
	//YAO
	//do
	//{
	//	 counter = contractCompositePaths();
	//	 counter += removeDeadEndNodes();
	//} while (counter > 0);
	//////////////////////////////

	return true;
}


/**********************************************************************************************************************
	This function check if a read contains other small reads. If a read is contained in more than one super read
	then it is assigned to the longest such super read.
**********************************************************************************************************************/
void OverlapGraph::markContainedReads(void)
{

//YAO	if(dataSet->longestReadLength == dataSet->shortestReadLength) // If all reads are of same length, then no need to do look for contained reads.
//YAO	{
//YAO		cout << "All reads are of same length. No contained reads." << endl;
//YAO		CLOCKSTOP;
//YAO		return;
//YAO	}
	UINT64 counter = 0;
	for(UINT64 i = 1; i < graph->size(); i++) // For each read
	{
		Read *read1 = dataSet->getReadFromID(i); // Get the read
		string readString = read1->getStringForward(); // Get the forward of the read
		string subString;
		for(UINT64 j = 1; j < read1->getReadLength() - hashTable->getHashStringLength(); j++) // fGr each substring of read1 of length getHashStringLength
		{
			subString = readString.substr(j,hashTable->getHashStringLength()); // Get the substring from read1
			vector<UINT64> * listOfReads=hashTable->getListOfReads(subString); // Search the substring in the hash table
			if(!listOfReads->empty()) // If other reads contain the substring as prefix or suffix
			{
				for(UINT64 k = 0; k < listOfReads->size(); k++) // For each read in the list.
				{
					UINT64 data = listOfReads->at(k); // We used bit operation in the hash table to store read ID and orientation
					Read *read2 = dataSet->getReadFromID(data & 0X3FFFFFFFFFFFFFFF); 	// Least significant 62 bits store the read number.
																						// Most significant 2 bits store the orientation.
																						// Orientation 0 means prefix of forward of the read
																						// Orientation 1 means suffix of forward of the read
																						// Orientation 2 means prefix of reverse of the read
																						// Orientation 3 means prefix of reverse of the read

					if(readString.length() > read2->getStringForward().length() && checkOverlapForContainedRead(read1,read2,(data >> 62),j)) // read1 need to be longer than read2 in order to contain read2
																																			 // Check if the remaining of the strings also match
					{
						if(read2->superReadID == 0) // This is the rist super read found. we store the ID of the super read.
						{
							read2->superReadID = i;
							counter ++;
						}
						else // Already found some previous super read. Now we have another super read.
						{
							if(readString.length() > dataSet->getReadFromID(read2->superReadID)->getReadLength()) // This super read is longer than the previous super read. Update the super read ID.
								read2->superReadID = i;
						}
					}
				}
			}
		}
		if(i%1000000 == 0)
			cout << setw(10) << counter << " contained reads in " << setw(10) << i << " super reads." << endl;
	}

	// Get some statistics
	UINT64 containedReads = 0, nonContainedReads = 0;
	for(UINT64 i = 1; i < graph->size(); i++)
	{
		Read *rr = dataSet->getReadFromID(i);
		if(rr->superReadID == 0) // Count the number of reads that are not contained by some other reads.
		{
			nonContainedReads++;
			std::stringstream sstm;

			sstm<<">"<<rr->getReadNumber()<<endl;
			sstm<<rr->getStringForward()<<endl;
			filePointer<<sstm.str();
		}
		else					// Count the number of reads that are contained by some other read.
			containedReads++;
	}
	cout<< endl << setw(10) << nonContainedReads << " Non-contained reads. (Keep as is)" << endl;
	cout<< setw(10) << containedReads << " contained reads. (Need to change their mate-pair information)" << endl;


}


/**********************************************************************************************************************
	Hash table search found that a proper substring of read1 is a prefix or suffix of read2 or reverse complement of
	read2 (depents on the orient).
	orient 0 means prefix of forward of the read2
	orient 1 means suffix of forward of the read2
	orient 2 means prefix of reverse of the read2
	orient 3 means prefix of reverse of the read2
	We need to check if the remaining of the stings match to see if read2 is contained in read1.
**********************************************************************************************************************/
bool OverlapGraph::checkOverlapForContainedRead(Read *read1, Read *read2, UINT64 orient, UINT64 start)
{
	string string1=read1->getStringForward(); // Get the forward of read1
	UINT64 hashStringLength = hashTable->getHashStringLength(), lengthRemaining1, lengthRemaining2;
	string string2 = (orient == 0 || orient== 1) ? read2->getStringForward() : read2->getStringReverse(); // Get the string in read2 based on the orientation.
	if(orient == 0 || orient == 2)
									// orient 0
									//   >--------MMMMMMMMMMMMMMM*******------> read1      M means match found by hash table
									//            MMMMMMMMMMMMMMM*******>       read2      * means we need to check these characters for match
									//				OR
									// orient 2
									//	 >---*****MMMMMMMMMMMMMMM*******------> read1
									//		      MMMMMMMMMMMMMMM*******<	    Reverese complement of read2
	{
		lengthRemaining1 = string1.length() - start - hashStringLength; 	// This is the remaining of read1
		lengthRemaining2 = string2.length() - hashStringLength; 	// This is the remaining of read2
		if(lengthRemaining1 >= lengthRemaining2)
		{
			return string1.substr(start + hashStringLength, lengthRemaining2) == string2.substr(hashStringLength, lengthRemaining2); // If the remaining of the string match, then read2 is contained in read1
		}
	}
	else							// orient 1
									//   >---*****MMMMMMMMMMMMMMM-------------> read1      M means match found by hash table
									//      >*****MMMMMMMMMMMMMMM       		read2      * means we need to check these characters for match
									//				OR
									// orient 3
									//	 >---*****MMMMMMMMMMMMMMM-------------> read1
									//		<*****MMMMMMMMMMMMMMM				Reverse Complement of Read2
	{
		lengthRemaining1 = start;
		lengthRemaining2 = string2.length() - hashStringLength;
		if(lengthRemaining1 >= lengthRemaining2)
		{
			return string1.substr(start - lengthRemaining2, lengthRemaining2) == string2.substr(0, lengthRemaining2); // If the remaining of the string match, then read2 is contained in read1
		}
	}
	return false;

}



/**********************************************************************************************************************
	Checks if two read overlaps.
	Hash table search found that a proper substring of read1 is a prefix or suffix of read2 or reverse complement of
	read2 (depents on the orient).
	Orientation 0 means prefix of forward of the read2
	Orientation 1 means suffix of forward of the read2
	Orientation 2 means prefix of reverse of the read2
	Orientation 3 means prefix of reverse of the read2
	We need to check if the remaining of the stings match to see if read1 and read2 overlap.
**********************************************************************************************************************/
bool OverlapGraph::checkOverlap(Read *read1, Read *read2, UINT64 orient, UINT64 start)
{
	string string1=read1->getStringForward(); // Get the forward string of read1
	UINT64 hashStringLength = hashTable->getHashStringLength();
	string string2 = (orient == 0 || orient== 1) ? read2->getStringForward() : read2->getStringReverse(); // Get the string from read2 according to orient.
	if(orient == 0 || orient == 2)		// orient 0
										//   >--------MMMMMMMMMMMMMMM*************> 			read1      M means match found by hash table
										//            MMMMMMMMMMMMMMM*************------->      read2      * means we need to check these characters for match
										//				OR
										// orient 2
										//	 >---*****MMMMMMMMMMMMMMM*************> 			read1
										//		      MMMMMMMMMMMMMMM*************-------<	    Reverese complement of read2
	{
		if(string1.length()- start - hashStringLength >= string2.length() - hashStringLength) // The overlap must continue till the end.
			return false;
		return string1.substr(start + hashStringLength, string1.length()-(start + hashStringLength)) == string2.substr(hashStringLength,  string1.length()-(start + hashStringLength)); // If the remaining strings match.
	}
	else								// orient 1
										//   	>********MMMMMMMMMMMMMMM-------------> 			read1      M means match found by hash table
										//  >----********MMMMMMMMMMMMMMM       		    		read2      * means we need to check these characters for match
										//				OR
										// orient 3
										//	 	>********MMMMMMMMMMMMMMM-------------> 			read1
										//	<----********MMMMMMMMMMMMMMM						Reverse Complement of Read2
	{
		if(string2.length()-hashStringLength < start)
			return false;
		return string1.substr(0, start) == string2.substr(string2.length()-hashStringLength-start, start); // If the remaining strings match.
	}
}



/**********************************************************************************************************************
	Insert all edges of a read in the overlap graph
**********************************************************************************************************************/
bool OverlapGraph::insertAllEdgesOfRead(UINT64 readNumber, vector<nodeType> * exploredReads)
{
	Read *read1 = dataSet->getReadFromID(readNumber); 	// Get the current read read1.
	string readString = read1->getStringForward(); 		// Get the forward string of read1.
	string subString;

	//this window [1, len-key-1] is to remove the redundancy, also make sure the minimumoverlap = key -1 is fulfilled
	for(UINT64 j = 1; j < read1->getReadLength()-hashTable->getHashStringLength(); j++) // For each proper substring of length getHashStringLength of read1
	{
		subString = readString.substr(j,hashTable->getHashStringLength());  // Get the proper substring s of read1.
		vector<UINT64> * listOfReads=hashTable->getListOfReads(subString); // Search the string in the hash table.

		if(!listOfReads->empty()) // If there are some reads that contain s as prefix or suffix of the read or their reverse complement
		{
			for(UINT64 k = 0; k < listOfReads->size(); k++) // For each such reads.
			{
				UINT64 data = listOfReads->at(k);			// We used bit operations in the hash table. Most significant 2 bits store orientation and least significant 62 bits store read ID.
				UINT16 overlapOffset;

				Read *read2 = dataSet->getReadFromID(data & 0X3FFFFFFFFFFFFFFF); 	// Least significant 62 bits store the read number.
				if(exploredReads->at(read2->getReadNumber())!= UNEXPLORED) 			// No need to discover the same edge again. All edges of read2 is already inserted in the graph.
						continue;
				if(read1->superReadID == 0 && read2->superReadID == 0 && checkOverlap(read1,read2,(data >> 62),j)) // Both read need to be non contained.
				{
				std::stringstream sstm;

				int orientation = data>>62;
				cout<<read1->getReadNumber() <<" "<< " >>> " << read2->getReadNumber() <<" orientation: "<<orientation<<" position: "<<j<<endl;

				if(read1->getStringForward()<read2->getStringForward())
				sstm<<" orientation: "<<orientation<<" position: "<<j<<"||||"<<read1->getStringForward()<<"||||"<<read2->getStringForward()<<endl;
				else
					sstm<<" orientation: "<<orientation<<" position: "<<j<<"||||"<<read2->getStringForward()<<"||||"<<read1->getStringForward()<<endl;
				filePointer<<sstm.str();
					this->totaledgenumber++;
					switch (data >> 62) // Most significant 2 bit represents  00 - prefix forward, 01 - suffix forward, 10 -  prefix reverse, 11 -  suffix reverse.
					{
						case 0: orientation = 3; overlapOffset = read1->getReadLength() - j; break; 				// 3 = r1>------->r2
						case 1: orientation = 0; overlapOffset = hashTable->getHashStringLength() + j; break; 		// 0 = r1<-------<r2
						case 2: orientation = 2; overlapOffset = read1->getReadLength() - j; break; 				// 2 = r1>-------<r2
						case 3: orientation = 1; overlapOffset = hashTable->getHashStringLength() + j; break; 		// 1 = r2<------->r2
					}
					insertEdge(read1,read2,orientation,read1->getStringForward().length()-overlapOffset); 			// Insert the edge in the graph.
				}
			}
		}
	}
//YAO	if(graph->at(readNumber)->size() != 0)
//YAO		sort(graph->at(readNumber)->begin(),graph->at(readNumber)->end(), compareEdges); // Sort the list of edges of the current node according to the overlap offset (ascending).
	return true;
}



/**********************************************************************************************************************
	Insert an edge in the overlap graph.
**********************************************************************************************************************/
bool OverlapGraph::insertEdge(Edge * edge)
{

	UINT64 ID = edge->getSourceRead()->getReadNumber(); // This is the source read.
	if(graph->at(ID)->empty()) 							// If there is no edge incident to the node
		numberOfNodes++;								// Then a new node is inserted in the graph. Number of nodes increased.
	graph->at(ID)->push_back(edge);						// Insert the edge in the list of edges of ID
		numberOfEdges++;								// Increase the number of edges.
	updateReadLocations(edge);							// If the current edge contains some reads, then we need to update their location information.
	return true;
}



/**********************************************************************************************************************
	Insert an edge in the graph.
**********************************************************************************************************************/
bool OverlapGraph::insertEdge(Read *read1, Read *read2, UINT8 orient, UINT16 overlapOffset)
{
	Edge * edge1 = new Edge(read1,read2,orient,overlapOffset);								// Create a new edge in the graph to insert.
	UINT16 overlapOffsetReverse = read2->getReadLength() + overlapOffset - read1->getReadLength();	// Set the overlap offset accordingly for the reverse edge. Note that read lengths are different.
																						// If read lengths are the same. Then the reverse edge has the same overlap offset.
	Edge * edge2 = new Edge(read2,read1,twinEdgeOrientation(orient),overlapOffsetReverse);		// Create a new edge for the reverses string.

	edge1->setReverseEdge(edge2);		// Set the reverse edge pointer.
	edge2->setReverseEdge(edge1);		// Set the reverse edge pinter.
	insertEdge(edge1);					// Insert the edge in the overlap graph.
	insertEdge(edge2);					// Insert the edge in the overlap graph.
	return true;
}


/**********************************************************************************************************************
	Update location of reads for the new edge.
	The new edge many contain some reads. We need to update the location information of all such read.
**********************************************************************************************************************/
bool OverlapGraph::updateReadLocations(Edge *edge)
{
	UINT64 distance = 0;
	for(UINT64 i = 0; i < edge->getListOfReads()->size(); i++) 				// For each read in the edge
	{
		distance += edge->getListOfOverlapOffsets()->at(i);					// Distance of the read in the edge.
		Read *read = dataSet->getReadFromID(edge->getListOfReads()->at(i));	// Get the read object.
		if(edge->getListOfOrientations()->at(i) == 1)						// Orientation 1 means that the forward string of the read is contained in this edge.
		{
			read->getListOfEdgesForward()->push_back(edge);					// Insert the current edge in the list of forward edges.
			read->getLocationOnEdgeForward()->push_back(distance);			// Also insert the distance within the edge.
			read->getListOfEdgesForward()->resize(read->getListOfEdgesForward()->size());
			read->getLocationOnEdgeForward()->resize(read->getLocationOnEdgeForward()->size());
		}
		else																// Orientation 0 means that the reverser string of the read is contained in this edge.
		{
			read->getListOfEdgesReverse()->push_back(edge);					// Insert the current edge in the list of reverser edges.
			read->getLocationOnEdgeReverse()->push_back(distance);			// Also insert the distance withing the edge.
			read->getListOfEdgesReverse()->resize(read->getListOfEdgesReverse()->size());
			read->getLocationOnEdgeReverse()->resize(read->getLocationOnEdgeReverse()->size());
		}
	}
	return true;
}



/**********************************************************************************************************************
	Orientation of a reverse edge;
	Twin edge of Orientation 0 = <-------< is Orientation 3 = >------->
	Twin edge of Orientation 1 = <-------> is Orientation 1 = <------->
	Twin edge of Orientation 2 = >-------< is Orientation 2 = >-------<
	Twin edge of Orientation 3 = >-------> is Orientation 0 = <-------<
**********************************************************************************************************************/
UINT8 OverlapGraph::twinEdgeOrientation(UINT8 orientation)
{
	UINT8 returnValue;
	if(orientation == 0)
		returnValue = 3;
	else if(orientation == 1)
		returnValue = 1;
	else if(orientation == 2)
		returnValue = 2;
	else if(orientation == 3)
		returnValue = 0;
	else
		MYEXIT("Unsupported edge orientation.")
	return returnValue;
}

