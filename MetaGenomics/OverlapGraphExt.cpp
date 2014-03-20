
#include "Common.h"
#include "OverlapGraph.h"

//#include "CS2/cs2.h"

void saveReadsToFastaFile(ofstream& filePointer, Dataset * dataSet, bool flowComputed, const vector<INT64>& meanOfInsertSizes,
			  const vector<INT64>& sdOfInsertSizes, const vector< vector<Edge *> * > *graph)
{
    UINT64 i, j;//, start;
    string sSeqRead;
    filePointer<<">-2\t"; // "-2" indicates general information
    filePointer<<dataSet->getNumberOfUniqueReads()<<","<<(int)(flowComputed)<<",";
    filePointer<<dataSet->numberOfPairedDatasets<<","<<meanOfInsertSizes.size()<<",(";
    for (i=0; i< meanOfInsertSizes.size();i++)
    {
	if (i>0)
	    filePointer<<",";
	filePointer <<"(" <<meanOfInsertSizes.at(i) <<","<< sdOfInsertSizes.at(i) <<")";
    }
    filePointer<<")"<<endl;
    filePointer<<"NoSeq"<<endl;
    for(i = 1; i < graph->size(); i++)
    {
	filePointer<<">-1\t"; // -1 indicates reads
	Read *r = dataSet->getReadFromID(i);
	filePointer << r->getReadNumber() <<"," << r->superReadID<<","<<r->getFrequency();
	filePointer << "," << r->getMatePairList()->size() << ",(";
	for (j=0; j< r->getMatePairList()->size(); j++ )
	{
	    if (j>0)
		filePointer << ",";
	    filePointer <<"(" <<r->getMatePairList()->at(j).matePairID<<",";
	    filePointer <<(int)(r->getMatePairList()->at(j).matePairOrientation) << ",";
	    filePointer <<(int)(r->getMatePairList()->at(j).datasetNumber) << ")";
	}
	filePointer<<")"<<endl;
	sSeqRead = r->getStringForward();
	filePointer<<sSeqRead<<endl;
//	start = 0;
//	do
//	{
		// print 100 bases in a line
//	    filePointer<<sSeqRead.substr(start, 100)<<endl;
//	    start += 100;
//	} while (start < sSeqRead.length());
    }
    filePointer<<">-3\t"<<endl; // -3 end
    filePointer<<"ENDREAD"<<endl;
}

bool OverlapGraph::saveGraphToFastaFile(string fileName)
{
    UINT64 i, j, k, iSeqId=0, source, destination, start;
    ofstream filePointer;
    string sSeq;
    Edge * e;
    filePointer.open(fileName.c_str());
    if(filePointer == NULL)
	MYEXIT("Unable to open file: "+fileName);

    saveReadsToFastaFile(filePointer, dataSet, flowComputed, meanOfInsertSizes, sdOfInsertSizes, graph);
    for (i=0; i< graph->size(); i++)
    {
	if(!graph->at(i)->empty())
	{
	    for(j = 0; j < graph->at(i)->size(); j++)
	    {
	    	// for each edge
		e = graph->at(i)->at(j);
		source = e->getSourceRead()->getReadNumber();
		destination = e->getDestinationRead()->getReadNumber();
		if(source < destination || (source == destination && e < e->getReverseEdge()))
		{
			// sequence ID start from 0 and increment for each edge
		    filePointer<<">"<<iSeqId;
		    // print the read ID of source read and destination read
		    filePointer<<"\t"<<source<<","<<destination;
		    // print the orientation and offset and flow and coverage of the edge
		    filePointer<<","<<(int)(e->getOrientation())<<","<<e->getOverlapOffset();
		    filePointer<<","<<e->flow<<","<<e->coverageDepth;
		    filePointer<<",(";
		    for(k = 0; k < e->getListOfReads()->size(); k++)
		    {
			// skip the first read, because you don't need to print comma
		    if (k > 0)
			    filePointer<<",";
	    	// go each contained read in this composite edge
		    // print read number, overlap offset and orientation
			filePointer<<"("<<e->getListOfReads()->at(k);
			filePointer<<","<<e->getListOfOverlapOffsets()->at(k);
			filePointer<<","<<(int)(e->getListOfOrientations()->at(k))<<")";
		    }
		    filePointer<<")"<<endl;	
		    // print out the sequence of the edge, but this is not used by the read function
		    sSeq = getStringInEdge(e);
		    start = 0;
		    do
		    {
		    	// print 100 bases in a line
			filePointer<<sSeq.substr(start, 100)<<endl;
			start += 100;
		    } while (start < sSeq.length());
		    iSeqId++;
		}
	    }
	}
    }
    filePointer.close();
    return true;
}

template <typename T> T transferStr(string & inputStr)
{
    T target;
    stringstream stream;
    stream<< inputStr;
    stream>> target;
    return target;
}

// sMajorEdgeInfo = "113277,179601,3,446606,0,0"
// source = 113277, destination = 179601, orientation =3, overlapOffset = 446606, flow =0, coverageDepth =0
void parseMajorEdgeInfo(const string & sMajorEdgeInfo, UINT64 & source, UINT64 & destination,
			UINT64 & orientation, UINT64 & overlapOffset, UINT64 & flow,
			UINT64 & coverageDepth)
{
    size_t iPosBeg=0, iPosSep;
    iPosSep = sMajorEdgeInfo.find(",", iPosBeg);
    string sStrValue;
    sStrValue = sMajorEdgeInfo.substr(iPosBeg, iPosSep-iPosBeg);
    source  = transferStr<UINT64>(sStrValue);
    iPosBeg = iPosSep+1;
    iPosSep = sMajorEdgeInfo.find(",", iPosBeg);
    sStrValue = sMajorEdgeInfo.substr(iPosBeg, iPosSep-iPosBeg);
    destination = transferStr<UINT64>(sStrValue);
    iPosBeg = iPosSep+1;
    iPosSep = sMajorEdgeInfo.find(",", iPosBeg);
    sStrValue = sMajorEdgeInfo.substr(iPosBeg, iPosSep-iPosBeg);
    orientation = transferStr<UINT64>(sStrValue);
    iPosBeg = iPosSep+1;
    iPosSep = sMajorEdgeInfo.find(",", iPosBeg);
    sStrValue = sMajorEdgeInfo.substr(iPosBeg, iPosSep-iPosBeg);
    overlapOffset = transferStr<UINT64>(sStrValue);  
    iPosBeg = iPosSep+1;
    iPosSep = sMajorEdgeInfo.find(",", iPosBeg);
    sStrValue = sMajorEdgeInfo.substr(iPosBeg, iPosSep-iPosBeg);
    flow = transferStr<UINT64>(sStrValue); 
    iPosBeg = iPosSep+1;
    iPosSep = sMajorEdgeInfo.find(",", iPosBeg);
    sStrValue = sMajorEdgeInfo.substr(iPosBeg, iPosSep-iPosBeg);
    coverageDepth = transferStr<UINT64>(sStrValue); 
}

void parseAnElement(const string & sInternalElement, vector<UINT64> *listReads,
		    vector<UINT16> *listOverlapOffsets, vector<UINT8> *listOrientations,
		    UINT64 & length)
{
    string sReadId, sOverlapOffset, sOrientation;
    size_t iPos1, iPos2;
    UINT64 currentReadId;
    UINT16 currentOverlapOffset;
    UINT8  currentOrientation;
    iPos1 = sInternalElement.find(",");
    iPos2 = sInternalElement.find(",", iPos1+1);
    sReadId = sInternalElement.substr(0, iPos1);
    currentReadId = transferStr<UINT64>(sReadId);
    listReads->push_back(currentReadId);
    sOverlapOffset = sInternalElement.substr(iPos1+1, iPos2-iPos1-1);
    currentOverlapOffset = transferStr<UINT16>(sOverlapOffset);
    listOverlapOffsets->push_back(currentOverlapOffset);
    length += currentOverlapOffset;
    sOrientation = sInternalElement.substr(iPos2+1);
    currentOrientation = (UINT8) transferStr<int>(sOrientation);
    listOrientations->push_back(currentOrientation);
}

// length is the length of the edge
void parseCompositionEdgeInfo(const string & sCompositionEdge, vector<UINT64> *listReads, 
			      vector<UINT16> *listOverlapOffsets, vector<UINT8> *listOrientations,
			      UINT64 & length)
{
    string sRestStr, sInternalElement;
    size_t iPos1, iPos2;// iPos1 for "(", iPos2 for ")"
    if (sCompositionEdge != "")
    {
	sRestStr = sCompositionEdge.substr(0);
	iPos1 = sRestStr.find("(");
	while (iPos1 != string::npos)
	{
	    iPos2 = sRestStr.find(")", iPos1+1);
	    sInternalElement = sRestStr.substr(iPos1+1, iPos2-iPos1-1);
	    // populate the vectors in this function
	    parseAnElement(sInternalElement, listReads, listOverlapOffsets, listOrientations, length);
	    sRestStr = sRestStr.substr(iPos2+1);
	    iPos1 = sRestStr.find("(");
	}
    }
}

UINT8 twinEdgeOrientation(UINT8 orientation) // copy from OverLapgraph.cpp
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


void parseEdgeInfo(const string & sEdgeInfo, Dataset * dataSet, Edge *edgeForward, Edge *edgeReverse)
{
	// variables to be populated from sEdgeInfo
    UINT64 source, destination, orientation, overlapOffset, flow, coverageDepth;
    string sMajorEdgeInfo, sCompositionEdge;
    size_t iPos;
    Read *read1, *read2;
    vector<UINT64> *listReads = new vector<UINT64>;
    vector<UINT16> *listOverlapOffsets = new vector<UINT16>;
    vector<UINT8> *listOrientations = new vector<UINT8>;
    vector<UINT64> *listReadsReverse = new vector<UINT64>;
    vector<UINT16> *listOverlapOffsetsReverse = new vector<UINT16>;
    vector<UINT8> *listOrientationsReverse = new vector<UINT8>;
    UINT64 length = 0, size, j, length1, length2, overlapOffsetForward, revereseOverlap;
    UINT64 revereseOverlapOffset;

    // sEdgeInfo = "113277,179601,3,446606,0,0,((194758,3,0),(143828,1,0), ...."
    iPos = sEdgeInfo.find("(");

    // sMajorEdgeInfo = "113277,179601,3,446606,0,0"
    sMajorEdgeInfo = sEdgeInfo.substr(0, iPos-1);
    // sCompositionEdge = "((194758,3,0),(143828,1,0), ...."
    sCompositionEdge = sEdgeInfo.substr(iPos+1, sEdgeInfo.length()-iPos-2);

    // example return:
    // source = 113277, destination = 179601, orientation =3, overlapOffset = 446606, flow =0, coverageDepth =0
    parseMajorEdgeInfo(sMajorEdgeInfo, source, destination, orientation, overlapOffset,
		       flow, coverageDepth);

    // parse the reads and get overlap, orientation and length of the whole edge
    parseCompositionEdgeInfo(sCompositionEdge, listReads, listOverlapOffsets, listOrientations, length);

    // get the pointers for source and destination
    read1 = dataSet->getReadFromID(source);
    read2 = dataSet->getReadFromID(destination);
    size = listReads->size();
    for(j = 0; j < size; j++)
    {
    	// populate reverse edge that goes from destination to source
	listReadsReverse->push_back(listReads->at(size-j-1));
	if(j == 0) // last/first read
	{
	    length1 = read2->getReadLength();
	    // calculate overlap forward on the original edge
	    overlapOffsetForward = overlapOffset - length;
	}
	else // any read within the edge
	{
	    length1 = dataSet->getReadFromID(listReads->at(size-j))->getReadLength();
	    overlapOffsetForward = listOverlapOffsets->at(size-j);
	}
	length2 = dataSet->getReadFromID(listReads->at(size-j-1))->getReadLength();
	revereseOverlap = length1 + overlapOffsetForward - length2;
	listOverlapOffsetsReverse->push_back(revereseOverlap);
	listOrientationsReverse->push_back(!(listOrientations->at(size-j-1)));
    }
    revereseOverlapOffset =  overlapOffset + read2->getReadLength() - read1->getReadLength();
		// make the forward edge.
    edgeForward->makeEdge(read1, read2, orientation, overlapOffset, listReads, listOverlapOffsets,
			  listOrientations);
		// make the reverse edge.
    edgeReverse->makeEdge(read2, read1, twinEdgeOrientation(orientation), revereseOverlapOffset,
			  listReadsReverse, listOverlapOffsetsReverse, listOrientationsReverse);

    edgeForward->setReverseEdge(edgeReverse); // set the reverse edge pointer
    edgeForward->flow = flow;				// Add the flow in the forward edge.
    edgeForward->coverageDepth = coverageDepth; // by yingfeng
    edgeReverse->setReverseEdge(edgeForward); // set the reverse edge pinter.
    edgeReverse->flow = flow;
    edgeReverse->coverageDepth = coverageDepth; //by yingfeng
}

void parseGeneralInfo(const string & sGenernal, Dataset * dataSet, bool & flowComputed, 
		      vector<INT64> & meanOfInsertSizes, vector<INT64> & sdOfInsertSizes)
{
    size_t iPosBeg, iPosEnd, iFirstSep;
    UINT64 numOfUniqueRead,numOfMatePairDataset,tempValue;
    //INT64 tempInt64Value;
    string sRegular, sComposition, sValueStr;
    
    iPosBeg = sGenernal.find("\t", 2)+1;
    iPosEnd = sGenernal.find("(", iPosBeg);
    sRegular= sGenernal.substr(iPosBeg, iPosEnd-iPosBeg-1);
    sComposition = sGenernal.substr(iPosEnd+1);
    iPosEnd = sComposition.rfind(")");
    sComposition = sComposition.substr(0, iPosEnd);
    iPosEnd = sRegular.find(",");
    sValueStr = sRegular.substr(0, iPosEnd);
    numOfUniqueRead = transferStr<UINT64>(sValueStr);
    dataSet->numberOfUniqueReads = numOfUniqueRead;
    iPosBeg = iPosEnd+1;
    iPosEnd = sRegular.find(",", iPosBeg) ;
    sValueStr = sGenernal.substr(iPosBeg, iPosEnd - iPosBeg);
    tempValue = transferStr<UINT64>(sValueStr);
    flowComputed = tempValue;
    iPosBeg = iPosEnd + 1;
    iPosEnd = sGenernal.find(",", iPosBeg);
    sValueStr = sGenernal.substr(iPosBeg, iPosEnd - iPosBeg);
    numOfMatePairDataset = transferStr<UINT64>(sValueStr);
    dataSet->numberOfPairedDatasets = numOfMatePairDataset;
    //sValueStr = sGenernal.substr(iPosEnd+1);
    //tempValue = transferStr<UINT64>(sValueStr);
    iPosBeg = sComposition.find("(");
    while (iPosBeg != string::npos)
    {
	iPosEnd   = sComposition.find(")", iPosBeg + 1);
	iFirstSep = sComposition.find(",", iPosBeg + 1);
	sValueStr = sComposition.substr(iPosBeg+1, iFirstSep - iPosBeg - 1);
	tempValue = transferStr<UINT64>(sValueStr);
	meanOfInsertSizes.push_back(tempValue);
	sValueStr = sComposition.substr(iFirstSep+1, iPosEnd - iFirstSep - 1);
	tempValue = transferStr<UINT64>(sValueStr);
	sdOfInsertSizes.push_back(tempValue);
	iPosBeg = sComposition.find("(", iPosEnd);
    }
}

void parseReadInfo(const string & sReadInfo, const string & sSeqForward, Dataset * dataSet) 
{
    string sReadGeneral, sReadCompoisition, sValueStr;
    size_t iPosBeg, iPosEnd, iFirstSep, iSecondSep;
    UINT64 tempValue, matePairID, orient, datasetNumber;
    Read *r = new Read;
    r->setRead(sSeqForward);
    iPosBeg = sReadInfo.find("\t", 2)+1;
    iPosEnd = sReadInfo.find("(", iPosBeg);
    sReadGeneral = sReadInfo.substr(iPosBeg, iPosEnd-iPosBeg-1);
    iPosBeg = iPosEnd;
    iPosEnd = sReadInfo.rfind(")");
    sReadCompoisition = sReadInfo.substr(iPosBeg+1, iPosEnd - iPosBeg - 1);
    iPosEnd = sReadGeneral.find(",");
    sValueStr = sReadGeneral.substr(0, iPosEnd);
    tempValue = transferStr<UINT64>(sValueStr);
    r->setReadNumber(tempValue);
    iPosBeg = iPosEnd + 1;
    iPosEnd = sReadGeneral.find(",", iPosBeg);
    sValueStr = sReadGeneral.substr(iPosBeg, iPosEnd - iPosBeg);
    tempValue = transferStr<UINT64>(sValueStr);
    r->superReadID = tempValue;
    iPosBeg = iPosEnd + 1;
    iPosEnd = sReadGeneral.find(",", iPosBeg);
    sValueStr = sReadGeneral.substr(iPosBeg, iPosEnd - iPosBeg);
    tempValue = transferStr<UINT64>(sValueStr);    
    r->setFrequency(tempValue);
    iPosBeg = sReadCompoisition.find("(");
    while (iPosBeg != string::npos)
    {
	iPosEnd = sReadCompoisition.find(")", iPosBeg + 1);
	iFirstSep = sReadCompoisition.find(",", iPosBeg + 1);
	iSecondSep= sReadCompoisition.find(",", iFirstSep+1);
	sValueStr = sReadCompoisition.substr(iPosBeg+1, iFirstSep-iPosBeg-1);
	matePairID = transferStr<UINT64>(sValueStr);
	sValueStr = sReadCompoisition.substr(iFirstSep+1, iSecondSep-iFirstSep-1);
	orient = transferStr<UINT64>(sValueStr);
	sValueStr = sReadCompoisition.substr(iSecondSep+1, iPosEnd-iSecondSep-1);
	datasetNumber = transferStr<UINT64>(sValueStr);
	r->addMatePair(matePairID,orient,datasetNumber);
	iPosBeg = sReadCompoisition.find("(", iPosEnd);
    }
    dataSet->addRead(r);
}

bool OverlapGraph::readGraphFromFastaFile(string fileName)
{
    string sCurrentLine, sEdgeInfo, sCurrentReadForward;
    ifstream filePointer;
    filePointer.open(fileName.c_str());
    size_t iPos;
    Edge *edgeForward, *edgeReverse;
    UINT64 i;
    
    if (filePointer == NULL)
	MYEXIT("Unable to open file: "+fileName);
    // create an empty graph
    graph = new vector< vector<Edge *> * >;
    
    while(!filePointer.eof())
    {
	getline(filePointer, sCurrentLine);
	if (sCurrentLine == "")
	    continue;
	if (sCurrentLine.substr(0,4) == ">-3\t")
	{
	    getline(filePointer, sCurrentLine);
	    break;
	}
	if (sCurrentLine.substr(0,4) == ">-2\t")
	{
	    parseGeneralInfo(sCurrentLine, dataSet, flowComputed, this->meanOfInsertSizes, this->sdOfInsertSizes);
	    continue;
	}
	if (sCurrentLine.substr(0,4) == ">-1\t")
	{
	    getline(filePointer, sCurrentReadForward);
	    parseReadInfo(sCurrentLine, sCurrentReadForward, dataSet);
	    continue;
	}
    }
    
    graph->reserve(dataSet->getNumberOfUniqueReads()+1);

    // initiaze the edges in the graph
    for(i = 0; i <= dataSet->getNumberOfUniqueReads(); i++) // initialize the list
    {
	vector<Edge *> *newList = new vector<Edge *>;
	graph->push_back(newList);
    }
    while (! filePointer.eof() )
    {
    	// get a line at a time.
	getline(filePointer, sCurrentLine);
	if ((sCurrentLine != "") && (sCurrentLine.at(0)=='>'))
	{
		// get the header line startin with >
		// e.g. the header line is  ">0  113277,179601,3,446606,0,0,((194758,3,0),(143828,1,0), ...."
	    iPos = sCurrentLine.find("\t", 2);
	    // sEdgeInfo = "113277,179601,3,446606,0,0,((194758,3,0),(143828,1,0), ...."
	    sEdgeInfo = sCurrentLine.substr(iPos+1);
	    edgeForward = new Edge();
	    edgeReverse = new Edge();

	    // dataSet is populated already
	    // this function populate edgeForward, edgeReverse
	    parseEdgeInfo(sEdgeInfo, dataSet, edgeForward, edgeReverse);
	    insertEdge(edgeForward);
	    insertEdge(edgeReverse);
	}
    }
    
    filePointer.close();
    
    return true;
}
