#include "JetProducer.hh"
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>

int main(int argc, char* argv[]){
	string file = "rootfiles/JetHT_R17_v18_JetHT_AOD_Run2017F_17Nov2017_output140.root";

	TFile* f = TFile::Open(file.c_str());
	
	string oname = "spikeCheck.root";
	TFile* ofile = TFile::Open(oname.c_str(),"RECREATE");

	JetProducer jp_woSpikes(f);
	jp_woSpikes.RejectSpikes(true);
	

	TH1D* nRhs_wSpikes = new TH1D("nRhs_wspikes","nRhs_wspikes",100,0,1000);
	TH1D* nRhs_woSpikes = new TH1D("nRhs_wospikes","nRhs_wospikes",100,0,1000);
	TH1D* nRhs_woOvwSpikesRatio = new TH1D("nRhs_woOvwSpikesRatio","nRhs_woOvwSpikesRatio",50,0.,1.1); 	

	TH2D* rhE_e4e1 = new TH2D("rhE_e4e1","rhE_e4e1;rhE;e4e1",50,0,1e3,50,0,0.12);
	TH2D* swCross_Ecut = new TH2D("swCross_Ecut","swCross_Ecut;swCross;Ecu",50,0,0.12,50,0,0.1);

	TH2D* swCross_rhTime = new TH2D("swCross_rhTime","swCross_rhTime",100,0.2,1.2,120,-60,60);
	TH2D* swCross_rhTime_reject = new TH2D("swCross_rhTime_reject","swCross_rhTime_reject",100,0.2,1.2,120,-60,60);
	TH2D* swCross_rhTime_accept = new TH2D("swCross_rhTime_accept","swCross_rhTime_accept",100,0.2,1.2,120,-60,60);

	TH2D* swCross_rhE_reject = new TH2D("swCross_rhE_reject","swCross_rhE_reject",50,0.2,1.2,120,0,500);
	TH2D* swCross_rhE_accept = new TH2D("swCross_rhE_accept","swCross_rhE_accept",50,0.2,1.2,120,0,500);
	
	TH2D* rhTime_rhE = new TH2D("rhTime_rhE","rhTime_rhE",120,-60,60,120,0,500);
	TH2D* rhTime_rhE_reject = new TH2D("rhTime_rhE_reject","rhTime_rhE_reject",120,-60,60,120,0,500);
	TH2D* rhTime_rhE_accept = new TH2D("rhTime_rhE_accept","rhTime_rhE_accept",120,-60,60,120,0,500);

	TFile* f2 = TFile::Open(file.c_str());
	TTree* tree = (TTree*)f2->Get("tree/llpgtree");
	ReducedBase* base = jp_woSpikes.GetBase();//new ReducedBase(tree);
	int nevts = base->fChain->GetEntries();

	cout << "file " << file << " has " << base->fChain->GetEntries() << " entries" << endl;
	JetProducer jp_wSpikes(f2);
	jp_wSpikes.RejectSpikes(false);

	vector<Jet> rhs_woSpikes, rhs_wSpikes;
	int nrhs;
	int skip = 1;
	for(int i = 0; i < nevts; i+=skip){
		if(i % 100 == 0) cout << "Event #" << i << endl;

		jp_wSpikes.GetRecHits(rhs_wSpikes,i);
		nRhs_wSpikes->Fill((int)rhs_wSpikes.size());
		
		jp_woSpikes.GetRecHits(rhs_woSpikes,i);
		nRhs_woSpikes->Fill((int)rhs_woSpikes.size());

		nRhs_woOvwSpikesRatio->Fill((double)rhs_woSpikes.size() / (double)rhs_wSpikes.size());

		base->GetEntry(i);
		nrhs = base->ECALRecHit_swCross->size();
		for(int r = 0; r < nrhs; r++){
			if(base->ECALRecHit_energy->at(r) < 4) continue;
			rhE_e4e1->Fill(base->ECALRecHit_energy->at(r), 1 - base->ECALRecHit_swCross->at(r));
                        swCross_Ecut->Fill(1 - base->ECALRecHit_swCross->at(r), 0.02*log10(base->ECALRecHit_energy->at(r))+0.02);
			swCross_rhTime->Fill(base->ECALRecHit_swCross->at(r), base->ECALRecHit_time->at(r));
	//cout << "sw " << base->ECALRecHit_swCross->at(r) << " 1 - sw " << 1 - base->ECALRecHit_swCross->at(r) << endl;		
        		rhTime_rhE->Fill(base->ECALRecHit_time->at(r), base->ECALRecHit_energy->at(r));
	                if(1 - base->ECALRecHit_swCross->at(r) < 0.02*log10(base->ECALRecHit_energy->at(r))+0.02){
				swCross_rhTime_reject->Fill(base->ECALRecHit_swCross->at(r), base->ECALRecHit_time->at(r));
				swCross_rhE_reject->Fill(base->ECALRecHit_swCross->at(r), base->ECALRecHit_energy->at(r));
        			rhTime_rhE_reject->Fill(base->ECALRecHit_time->at(r), base->ECALRecHit_energy->at(r));
			}
			else{
				swCross_rhTime_accept->Fill(base->ECALRecHit_swCross->at(r), base->ECALRecHit_time->at(r));
				swCross_rhE_accept->Fill(base->ECALRecHit_swCross->at(r), base->ECALRecHit_energy->at(r));
        			rhTime_rhE_accept->Fill(base->ECALRecHit_time->at(r), base->ECALRecHit_energy->at(r));
			}
		}

		

	}

	ofile->cd();
	nRhs_wSpikes->Write();
	nRhs_woSpikes->Write();
	rhE_e4e1->Write();
	swCross_Ecut->Write();
	nRhs_woOvwSpikesRatio->Write();
	swCross_rhTime->Write();
	swCross_rhTime_reject->Write();
	swCross_rhTime_accept->Write();
	swCross_rhE_reject->Write();
	swCross_rhE_accept->Write();
	rhTime_rhE->Write();
	rhTime_rhE_reject->Write();
	rhTime_rhE_accept->Write();

	cout << "Wrote histograms to " << ofile->GetName() << endl;

	delete ofile;
	return 0;
	
};
