#include "Config.h"
#include "OverlapGraphConstructor.h"
#include "PairedReadMerger.h"
#include "ErrorCorrector.h"


int main(int argc, char **argv)
{

	if(!Config::setConfig(argc, argv)){
		cout << "Error: wrong configurations" << endl;
		return false;
	}

	if (Config::getOperation() == "ConstructGraph"){
		// start graph construction
		OverlapGraphConstructor* graphconstructor = new OverlapGraphConstructor();
		if(!graphconstructor->start())
		{
			cout<< "Failed to run OverlapGraphConstructor start()";
		}
		delete graphconstructor;
	}
	else if (Config::getOperation() == "MergePairedEndReads"){
		// start merge reads
		PairedReadMerger* merger = new PairedReadMerger();
		if(!merger->start())
		{
			cout<< "Failed to run PairedReadMerger start()";
		}
		delete merger;

	}
	else if (Config::getOperation() == "CorrectErrors"){
		// start correct errors
		ErrorCorrector* errorCorrector = new ErrorCorrector();
		if(!errorCorrector->start())
		{
			cout<< "Failed to run ErrorCorrector start()";
		}
		delete errorCorrector;
	}
	else if (Config::getOperation() == "RemoveContainedReads"){
		// start remove container reads
	}
	else{
		cout << "Error: invalid operation type."<< endl;
		return false;
	}
	return true;
}



