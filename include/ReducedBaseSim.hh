//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu May 22 13:56:29 2025 by ROOT version 6.30/06
// from TTree ReducedBaseSim/ReducedBaseSim
// found on file: simNtuples_testFatJets_ttbar.root
//////////////////////////////////////////////////////////
#ifndef ReducedBaseSim_h
#define ReducedBaseSim_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

using namespace std;

class ReducedBaseSim {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           event;
   vector<double>  *ECALRecHit_energy;
   vector<double>  *ECALRecHit_rhx;
   vector<double>  *ECALRecHit_rhy;
   vector<double>  *ECALRecHit_rhz;
   vector<double>  *ECALRecHit_time;
   vector<double>  *ECALRecHit_eta;
   vector<double>  *ECALRecHit_phi;
   vector<unsigned int> *ECALRecHit_ID;
   Int_t           nRHs;
   Double_t        PV_x;
   Double_t        PV_y;
   Double_t        PV_z;
   vector<double>  *ECALSpike_energy;
   Int_t           nSpikes;
   Int_t           nRecoParticles;
   vector<double>  *Jet_genEta;
   vector<double>  *Jet_genPhi;
   vector<double>  *Jet_genEnergy;
   vector<double>  *Jet_genPt;
   vector<double>  *Jet_genPz;
   vector<double>  *Jet_genMass;
   Int_t           Jet_genNJet;
   vector<vector<int> > *Jet_genConstituentIdxs;
   vector<int>     *Jet_genNConstituents;
   vector<double>  *FatJet_genEta;
   vector<double>  *FatJet_genPhi;
   vector<double>  *FatJet_genEnergy;
   vector<double>  *FatJet_genPt;
   vector<double>  *FatJet_genPz;
   vector<double>  *FatJet_genMass;
   Int_t           FatJet_genNFatJet;
   vector<vector<int> > *FatJet_genConstituentIdxs;
   vector<int>     *FatJet_genNConstituents;
   vector<int>     *Top_decayId;
   vector<double>  *Jet_eta;
   vector<double>  *Jet_phi;
   vector<double>  *Jet_energy;
   vector<double>  *Jet_pt;
   vector<double>  *Jet_mass;
   vector<vector<unsigned int> > *Jet_RhIDs;
   Int_t           Jet_NJet;
   vector<double>  *Track_px;
   vector<double>  *Track_py;
   vector<double>  *Track_pz;
   vector<double>  *Track_eta;
   vector<double>  *Track_phi;
   vector<double>  *genpart_eta;
   vector<double>  *genpart_phi;
   vector<double>  *genpart_energy;
   vector<double>  *genpart_pt;
   vector<double>  *genpart_pz;
   vector<double>  *genpart_mass;
   vector<int>     *genpart_id;
   vector<int>     *genpart_momIdx;
   vector<int>     *genpart_idx;
   Int_t           genpart_ngenpart;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_ECALRecHit_energy;   //!
   TBranch        *b_ECALRecHit_rhx;   //!
   TBranch        *b_ECALRecHit_rhy;   //!
   TBranch        *b_ECALRecHit_rhz;   //!
   TBranch        *b_ECALRecHit_time;   //!
   TBranch        *b_ECALRecHit_eta;   //!
   TBranch        *b_ECALRecHit_phi;   //!
   TBranch        *b_ECALRecHit_ID;   //!
   TBranch        *b_nRHs;   //!
   TBranch        *b_PV_x;   //!
   TBranch        *b_PV_y;   //!
   TBranch        *b_PV_z;   //!
   TBranch        *b_ECALSpike_energy;   //!
   TBranch        *b_nSpikes;   //!
   TBranch        *b_nRecoParticles;   //!
   TBranch        *b_Jet_genEta;   //!
   TBranch        *b_Jet_genPhi;   //!
   TBranch        *b_Jet_genEnergy;   //!
   TBranch        *b_Jet_genPt;   //!
   TBranch        *b_Jet_genPz;   //!
   TBranch        *b_Jet_genMass;   //!
   TBranch        *b_Jet_genNJet;   //!
   TBranch        *b_Jet_genConstituentIdxs;   //!
   TBranch        *b_Jet_genNConstituents;   //!
   TBranch        *b_FatJet_genEta;   //!
   TBranch        *b_FatJet_genPhi;   //!
   TBranch        *b_FatJet_genEnergy;   //!
   TBranch        *b_FatJet_genPt;   //!
   TBranch        *b_FatJet_genPz;   //!
   TBranch        *b_FatJet_genMass;   //!
   TBranch        *b_FatJet_genNFatJet;   //!
   TBranch        *b_FatJet_genConstituentIdxs;   //!
   TBranch        *b_FatJet_genNConstituents;   //!
   TBranch        *b_Top_decayId;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_energy;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_RhIDs;   //!
   TBranch        *b_Jet_NJet;   //!
   TBranch        *b_Track_px;   //!
   TBranch        *b_Track_py;   //!
   TBranch        *b_Track_pz;   //!
   TBranch        *b_Track_eta;   //!
   TBranch        *b_Track_phi;   //!
   TBranch        *b_genpart_eta;   //!
   TBranch        *b_genpart_phi;   //!
   TBranch        *b_genpart_energy;   //!
   TBranch        *b_genpart_pt;   //!
   TBranch        *b_genpart_pz;   //!
   TBranch        *b_genpart_mass;   //!
   TBranch        *b_genpart_id;   //!
   TBranch        *b_genpart_momIdx;   //!
   TBranch        *b_genpart_idx;   //!
   TBranch        *b_genpart_ngenpart;   //!

   ReducedBaseSim(TTree *tree=0);
   virtual ~ReducedBaseSim();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   //virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

//#ifdef ReducedBaseSim_cxx
inline ReducedBaseSim::ReducedBaseSim(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("simNtuples_testFatJets_ttbar.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("simNtuples_testFatJets_ttbar.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("simNtuples_testFatJets_ttbar.root:/tree");
      dir->GetObject("ReducedBaseSim",tree);

   }
   Init(tree);
}

inline ReducedBaseSim::~ReducedBaseSim()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

inline Int_t ReducedBaseSim::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

inline Long64_t ReducedBaseSim::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

inline void ReducedBaseSim::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   ECALRecHit_energy = 0;
   ECALRecHit_rhx = 0;
   ECALRecHit_rhy = 0;
   ECALRecHit_rhz = 0;
   ECALRecHit_time = 0;
   ECALRecHit_eta = 0;
   ECALRecHit_phi = 0;
   ECALRecHit_ID = 0;
   ECALSpike_energy = 0;
   Jet_genEta = 0;
   Jet_genPhi = 0;
   Jet_genEnergy = 0;
   Jet_genPt = 0;
   Jet_genPz = 0;
   Jet_genMass = 0;
   Jet_genConstituentIdxs = 0;
   Jet_genNConstituents = 0;
   FatJet_genEta = 0;
   FatJet_genPhi = 0;
   FatJet_genEnergy = 0;
   FatJet_genPt = 0;
   FatJet_genPz = 0;
   FatJet_genMass = 0;
   FatJet_genConstituentIdxs = 0;
   FatJet_genNConstituents = 0;
   Top_decayId = 0;
   Jet_eta = 0;
   Jet_phi = 0;
   Jet_energy = 0;
   Jet_pt = 0;
   Jet_mass = 0;
   Jet_RhIDs = 0;
   Track_px = 0;
   Track_py = 0;
   Track_pz = 0;
   Track_eta = 0;
   Track_phi = 0;
   genpart_eta = 0;
   genpart_phi = 0;
   genpart_energy = 0;
   genpart_pt = 0;
   genpart_pz = 0;
   genpart_mass = 0;
   genpart_id = 0;
   genpart_momIdx = 0;
   genpart_idx = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("ECALRecHit_energy", &ECALRecHit_energy, &b_ECALRecHit_energy);
   fChain->SetBranchAddress("ECALRecHit_rhx", &ECALRecHit_rhx, &b_ECALRecHit_rhx);
   fChain->SetBranchAddress("ECALRecHit_rhy", &ECALRecHit_rhy, &b_ECALRecHit_rhy);
   fChain->SetBranchAddress("ECALRecHit_rhz", &ECALRecHit_rhz, &b_ECALRecHit_rhz);
   fChain->SetBranchAddress("ECALRecHit_time", &ECALRecHit_time, &b_ECALRecHit_time);
   fChain->SetBranchAddress("ECALRecHit_eta", &ECALRecHit_eta, &b_ECALRecHit_eta);
   fChain->SetBranchAddress("ECALRecHit_phi", &ECALRecHit_phi, &b_ECALRecHit_phi);
   fChain->SetBranchAddress("ECALRecHit_ID", &ECALRecHit_ID, &b_ECALRecHit_ID);
   fChain->SetBranchAddress("nRHs", &nRHs, &b_nRHs);
   fChain->SetBranchAddress("PV_x", &PV_x, &b_PV_x);
   fChain->SetBranchAddress("PV_y", &PV_y, &b_PV_y);
   fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
   fChain->SetBranchAddress("ECALSpike_energy", &ECALSpike_energy, &b_ECALSpike_energy);
   fChain->SetBranchAddress("nSpikes", &nSpikes, &b_nSpikes);
   fChain->SetBranchAddress("nRecoParticles", &nRecoParticles, &b_nRecoParticles);
   fChain->SetBranchAddress("Jet_genEta", &Jet_genEta, &b_Jet_genEta);
   fChain->SetBranchAddress("Jet_genPhi", &Jet_genPhi, &b_Jet_genPhi);
   fChain->SetBranchAddress("Jet_genEnergy", &Jet_genEnergy, &b_Jet_genEnergy);
   fChain->SetBranchAddress("Jet_genPt", &Jet_genPt, &b_Jet_genPt);
   fChain->SetBranchAddress("Jet_genPz", &Jet_genPz, &b_Jet_genPz);
   fChain->SetBranchAddress("Jet_genMass", &Jet_genMass, &b_Jet_genMass);
   fChain->SetBranchAddress("Jet_genNJet", &Jet_genNJet, &b_Jet_genNJet);
   fChain->SetBranchAddress("Jet_genConstituentIdxs", &Jet_genConstituentIdxs, &b_Jet_genConstituentIdxs);
   fChain->SetBranchAddress("Jet_genNConstituents", &Jet_genNConstituents, &b_Jet_genNConstituents);
   fChain->SetBranchAddress("FatJet_genEta", &FatJet_genEta, &b_FatJet_genEta);
   fChain->SetBranchAddress("FatJet_genPhi", &FatJet_genPhi, &b_FatJet_genPhi);
   fChain->SetBranchAddress("FatJet_genEnergy", &FatJet_genEnergy, &b_FatJet_genEnergy);
   fChain->SetBranchAddress("FatJet_genPt", &FatJet_genPt, &b_FatJet_genPt);
   fChain->SetBranchAddress("FatJet_genPz", &FatJet_genPz, &b_FatJet_genPz);
   fChain->SetBranchAddress("FatJet_genMass", &FatJet_genMass, &b_FatJet_genMass);
   fChain->SetBranchAddress("FatJet_genNFatJet", &FatJet_genNFatJet, &b_FatJet_genNFatJet);
   fChain->SetBranchAddress("FatJet_genConstituentIdxs", &FatJet_genConstituentIdxs, &b_FatJet_genConstituentIdxs);
   fChain->SetBranchAddress("FatJet_genNConstituents", &FatJet_genNConstituents, &b_FatJet_genNConstituents);
   fChain->SetBranchAddress("Top_decayId", &Top_decayId, &b_Top_decayId);
   fChain->SetBranchAddress("Jet_eta", &Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_phi", &Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_energy", &Jet_energy, &b_Jet_energy);
   fChain->SetBranchAddress("Jet_pt", &Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_mass", &Jet_mass, &b_Jet_mass);
   fChain->SetBranchAddress("Jet_RhIDs", &Jet_RhIDs, &b_Jet_RhIDs);
   fChain->SetBranchAddress("Jet_NJet", &Jet_NJet, &b_Jet_NJet);
   fChain->SetBranchAddress("Track_px", &Track_px, &b_Track_px);
   fChain->SetBranchAddress("Track_py", &Track_py, &b_Track_py);
   fChain->SetBranchAddress("Track_pz", &Track_pz, &b_Track_pz);
   fChain->SetBranchAddress("Track_eta", &Track_eta, &b_Track_eta);
   fChain->SetBranchAddress("Track_phi", &Track_phi, &b_Track_phi);
   fChain->SetBranchAddress("genpart_eta", &genpart_eta, &b_genpart_eta);
   fChain->SetBranchAddress("genpart_phi", &genpart_phi, &b_genpart_phi);
   fChain->SetBranchAddress("genpart_energy", &genpart_energy, &b_genpart_energy);
   fChain->SetBranchAddress("genpart_pt", &genpart_pt, &b_genpart_pt);
   fChain->SetBranchAddress("genpart_pz", &genpart_pz, &b_genpart_pz);
   fChain->SetBranchAddress("genpart_mass", &genpart_mass, &b_genpart_mass);
   fChain->SetBranchAddress("genpart_id", &genpart_id, &b_genpart_id);
   fChain->SetBranchAddress("genpart_momIdx", &genpart_momIdx, &b_genpart_momIdx);
   fChain->SetBranchAddress("genpart_idx", &genpart_idx, &b_genpart_idx);
   fChain->SetBranchAddress("genpart_ngenpart", &genpart_ngenpart, &b_genpart_ngenpart);
   Notify();
}

inline Bool_t ReducedBaseSim::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

inline void ReducedBaseSim::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

inline Int_t ReducedBaseSim::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
//#endif // #ifdef ReducedBaseSim_cxx
