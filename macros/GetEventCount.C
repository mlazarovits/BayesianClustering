#include <iostream>
#include "TFile.h"
#include "TTree.h"

void GetEventCount(string flist){
	if(gSystem->AccessPathName(flist.c_str())){ 
		cout << "Error: file " << flist << " doesn't exist." << endl; 
		return; 
	}	
	std::ifstream infile(flist);
	TChain* ch = new TChain("tree/configtree");
	ch->SetTitle(flist.c_str());
	string file;
	cout << "TChaining files in " << flist << "..." << endl;
	while(std::getline(infile,file)){
		if( file[0] == '#' ) continue;
		if(gSystem->AccessPathName(file.c_str())){
			cout << "Skipping file " << file << " - not found." << endl;
			continue;
		}
		//std::cout << "--  adding file: " << file << std::endl;
		ch->Add(file.c_str()); //skip non-recoverable files
	}
	TBranch* b_ntot;
	int ntot;

	TBranch* b_nfiltered;
	int nfiltered;

	TBranch* b_sumEvtWt;
	float sumEvtWt;
	
	ch->SetBranchAddress("nTotEvts",&ntot,&b_ntot);
	ch->SetBranchAddress("nFltrdEvts",&nfiltered,&b_nfiltered);
	ch->SetBranchAddress("sumEvtWgt",&sumEvtWt,&b_sumEvtWt);

	int nentries = ch->GetEntries();
	int totEvts = 0;
	float totEvtWt = 0;
	for(int i = 0; i < nentries; i++){
		ch->GetEntry(i);
		totEvts += ntot;
		totEvtWt += sumEvtWt;
	}
	cout << "Total events ran over " << totEvts << " with total event weight " << totEvtWt << endl;

}
