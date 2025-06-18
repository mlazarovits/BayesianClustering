//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jun 18 12:44:49 2025 by ROOT version 6.30/06
// from TTree ReducedBaseSim/ReducedBaseSim
// found on file: root://cmseos.fnal.gov//store/user/malazaro/SimNtuples/condorSimNtuples_ttbar_defaultv9p9.root
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
   Double_t        PV_t;
   Double_t        ootPV_x;
   Double_t        ootPV_y;
   Double_t        ootPV_z;
   Double_t        ootPV_t;
   vector<double>  *ECALSpike_energy;
   Int_t           nSpikes;
   Int_t           nRecoParticles;
   vector<double>  *AK4Jet_genEta;
   vector<double>  *AK4Jet_genPhi;
   vector<double>  *AK4Jet_genEnergy;
   vector<double>  *AK4Jet_genPt;
   vector<double>  *AK4Jet_genPz;
   vector<double>  *AK4Jet_genMass;
   Int_t           AK4Jet_genNJet;
   vector<vector<int> > *AK4Jet_genConstituentIdxs;
   vector<int>     *AK4Jet_genNConstituents;
   vector<double>  *AK8Jet_genEta;
   vector<double>  *AK8Jet_genPhi;
   vector<double>  *AK8Jet_genEnergy;
   vector<double>  *AK8Jet_genPt;
   vector<double>  *AK8Jet_genPz;
   vector<double>  *AK8Jet_genMass;
   Int_t           AK8Jet_genNJet;
   vector<vector<int> > *AK8Jet_genConstituentIdxs;
   vector<int>     *AK8Jet_genNConstituents;
   vector<double>  *AK15Jet_genEta;
   vector<double>  *AK15Jet_genPhi;
   vector<double>  *AK15Jet_genEnergy;
   vector<double>  *AK15Jet_genPt;
   vector<double>  *AK15Jet_genPz;
   vector<double>  *AK15Jet_genMass;
   Int_t           AK15Jet_genNJet;
   vector<vector<int> > *AK15Jet_genConstituentIdxs;
   vector<int>     *AK15Jet_genNConstituents;
   vector<int>     *Top_decayId;
   vector<double>  *AK4Jet_eta;
   vector<double>  *AK4Jet_phi;
   vector<double>  *AK4Jet_energy;
   vector<double>  *AK4Jet_pt;
   vector<double>  *AK4Jet_mass;
   vector<vector<unsigned int> > *AK4Jet_RhIDs;
   Int_t           AK4Jet_NJet;
   vector<double>  *AK8Jet_eta;
   vector<double>  *AK8Jet_phi;
   vector<double>  *AK8Jet_energy;
   vector<double>  *AK8Jet_pt;
   vector<double>  *AK8Jet_mass;
   vector<vector<unsigned int> > *AK8Jet_RhIDs;
   Int_t           AK8Jet_NJet;
   vector<double>  *AK15Jet_eta;
   vector<double>  *AK15Jet_phi;
   vector<double>  *AK15Jet_energy;
   vector<double>  *AK15Jet_pt;
   vector<double>  *AK15Jet_mass;
   vector<vector<unsigned int> > *AK15Jet_RhIDs;
   Int_t           AK15Jet_NJet;
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
   TBranch        *b_PV_t;   //!
   TBranch        *b_ootPV_x;   //!
   TBranch        *b_ootPV_y;   //!
   TBranch        *b_ootPV_z;   //!
   TBranch        *b_ootPV_t;   //!
   TBranch        *b_ECALSpike_energy;   //!
   TBranch        *b_nSpikes;   //!
   TBranch        *b_nRecoParticles;   //!
   TBranch        *b_AK4Jet_genEta;   //!
   TBranch        *b_AK4Jet_genPhi;   //!
   TBranch        *b_AK4Jet_genEnergy;   //!
   TBranch        *b_AK4Jet_genPt;   //!
   TBranch        *b_AK4Jet_genPz;   //!
   TBranch        *b_AK4Jet_genMass;   //!
   TBranch        *b_AK4Jet_genNJet;   //!
   TBranch        *b_AK4Jet_genConstituentIdxs;   //!
   TBranch        *b_AK4Jet_genNConstituents;   //!
   TBranch        *b_AK8Jet_genEta;   //!
   TBranch        *b_AK8Jet_genPhi;   //!
   TBranch        *b_AK8Jet_genEnergy;   //!
   TBranch        *b_AK8Jet_genPt;   //!
   TBranch        *b_AK8Jet_genPz;   //!
   TBranch        *b_AK8Jet_genMass;   //!
   TBranch        *b_AK8Jet_genNJet;   //!
   TBranch        *b_AK8Jet_genConstituentIdxs;   //!
   TBranch        *b_AK8Jet_genNConstituents;   //!
   TBranch        *b_AK15Jet_genEta;   //!
   TBranch        *b_AK15Jet_genPhi;   //!
   TBranch        *b_AK15Jet_genEnergy;   //!
   TBranch        *b_AK15Jet_genPt;   //!
   TBranch        *b_AK15Jet_genPz;   //!
   TBranch        *b_AK15Jet_genMass;   //!
   TBranch        *b_AK15Jet_genNJet;   //!
   TBranch        *b_AK15Jet_genConstituentIdxs;   //!
   TBranch        *b_AK15Jet_genNConstituents;   //!
   TBranch        *b_Top_decayId;   //!
   TBranch        *b_AK4Jet_eta;   //!
   TBranch        *b_AK4Jet_phi;   //!
   TBranch        *b_AK4Jet_energy;   //!
   TBranch        *b_AK4Jet_pt;   //!
   TBranch        *b_AK4Jet_mass;   //!
   TBranch        *b_AK4Jet_RhIDs;   //!
   TBranch        *b_AK4Jet_NJet;   //!
   TBranch        *b_AK8Jet_eta;   //!
   TBranch        *b_AK8Jet_phi;   //!
   TBranch        *b_AK8Jet_energy;   //!
   TBranch        *b_AK8Jet_pt;   //!
   TBranch        *b_AK8Jet_mass;   //!
   TBranch        *b_AK8Jet_RhIDs;   //!
   TBranch        *b_AK8Jet_NJet;   //!
   TBranch        *b_AK15Jet_eta;   //!
   TBranch        *b_AK15Jet_phi;   //!
   TBranch        *b_AK15Jet_energy;   //!
   TBranch        *b_AK15Jet_pt;   //!
   TBranch        *b_AK15Jet_mass;   //!
   TBranch        *b_AK15Jet_RhIDs;   //!
   TBranch        *b_AK15Jet_NJet;   //!
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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://cmseos.fnal.gov//store/user/malazaro/SimNtuples/condorSimNtuples_ttbar_defaultv9p9.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root://cmseos.fnal.gov//store/user/malazaro/SimNtuples/condorSimNtuples_ttbar_defaultv9p9.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("root://cmseos.fnal.gov//store/user/malazaro/SimNtuples/condorSimNtuples_ttbar_defaultv9p9.root:/tree");
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
   AK4Jet_genEta = 0;
   AK4Jet_genPhi = 0;
   AK4Jet_genEnergy = 0;
   AK4Jet_genPt = 0;
   AK4Jet_genPz = 0;
   AK4Jet_genMass = 0;
   AK4Jet_genConstituentIdxs = 0;
   AK4Jet_genNConstituents = 0;
   AK8Jet_genEta = 0;
   AK8Jet_genPhi = 0;
   AK8Jet_genEnergy = 0;
   AK8Jet_genPt = 0;
   AK8Jet_genPz = 0;
   AK8Jet_genMass = 0;
   AK8Jet_genConstituentIdxs = 0;
   AK8Jet_genNConstituents = 0;
   AK15Jet_genEta = 0;
   AK15Jet_genPhi = 0;
   AK15Jet_genEnergy = 0;
   AK15Jet_genPt = 0;
   AK15Jet_genPz = 0;
   AK15Jet_genMass = 0;
   AK15Jet_genConstituentIdxs = 0;
   AK15Jet_genNConstituents = 0;
   Top_decayId = 0;
   AK4Jet_eta = 0;
   AK4Jet_phi = 0;
   AK4Jet_energy = 0;
   AK4Jet_pt = 0;
   AK4Jet_mass = 0;
   AK4Jet_RhIDs = 0;
   AK8Jet_eta = 0;
   AK8Jet_phi = 0;
   AK8Jet_energy = 0;
   AK8Jet_pt = 0;
   AK8Jet_mass = 0;
   AK8Jet_RhIDs = 0;
   AK15Jet_eta = 0;
   AK15Jet_phi = 0;
   AK15Jet_energy = 0;
   AK15Jet_pt = 0;
   AK15Jet_mass = 0;
   AK15Jet_RhIDs = 0;
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
   fChain->SetBranchAddress("PV_t", &PV_t, &b_PV_t);
   fChain->SetBranchAddress("ootPV_x", &ootPV_x, &b_ootPV_x);
   fChain->SetBranchAddress("ootPV_y", &ootPV_y, &b_ootPV_y);
   fChain->SetBranchAddress("ootPV_z", &ootPV_z, &b_ootPV_z);
   fChain->SetBranchAddress("ootPV_t", &ootPV_t, &b_ootPV_t);
   fChain->SetBranchAddress("ECALSpike_energy", &ECALSpike_energy, &b_ECALSpike_energy);
   fChain->SetBranchAddress("nSpikes", &nSpikes, &b_nSpikes);
   fChain->SetBranchAddress("nRecoParticles", &nRecoParticles, &b_nRecoParticles);
   fChain->SetBranchAddress("AK4Jet_genEta", &AK4Jet_genEta, &b_AK4Jet_genEta);
   fChain->SetBranchAddress("AK4Jet_genPhi", &AK4Jet_genPhi, &b_AK4Jet_genPhi);
   fChain->SetBranchAddress("AK4Jet_genEnergy", &AK4Jet_genEnergy, &b_AK4Jet_genEnergy);
   fChain->SetBranchAddress("AK4Jet_genPt", &AK4Jet_genPt, &b_AK4Jet_genPt);
   fChain->SetBranchAddress("AK4Jet_genPz", &AK4Jet_genPz, &b_AK4Jet_genPz);
   fChain->SetBranchAddress("AK4Jet_genMass", &AK4Jet_genMass, &b_AK4Jet_genMass);
   fChain->SetBranchAddress("AK4Jet_genNJet", &AK4Jet_genNJet, &b_AK4Jet_genNJet);
   fChain->SetBranchAddress("AK4Jet_genConstituentIdxs", &AK4Jet_genConstituentIdxs, &b_AK4Jet_genConstituentIdxs);
   fChain->SetBranchAddress("AK4Jet_genNConstituents", &AK4Jet_genNConstituents, &b_AK4Jet_genNConstituents);
   fChain->SetBranchAddress("AK8Jet_genEta", &AK8Jet_genEta, &b_AK8Jet_genEta);
   fChain->SetBranchAddress("AK8Jet_genPhi", &AK8Jet_genPhi, &b_AK8Jet_genPhi);
   fChain->SetBranchAddress("AK8Jet_genEnergy", &AK8Jet_genEnergy, &b_AK8Jet_genEnergy);
   fChain->SetBranchAddress("AK8Jet_genPt", &AK8Jet_genPt, &b_AK8Jet_genPt);
   fChain->SetBranchAddress("AK8Jet_genPz", &AK8Jet_genPz, &b_AK8Jet_genPz);
   fChain->SetBranchAddress("AK8Jet_genMass", &AK8Jet_genMass, &b_AK8Jet_genMass);
   fChain->SetBranchAddress("AK8Jet_genNJet", &AK8Jet_genNJet, &b_AK8Jet_genNJet);
   fChain->SetBranchAddress("AK8Jet_genConstituentIdxs", &AK8Jet_genConstituentIdxs, &b_AK8Jet_genConstituentIdxs);
   fChain->SetBranchAddress("AK8Jet_genNConstituents", &AK8Jet_genNConstituents, &b_AK8Jet_genNConstituents);
   fChain->SetBranchAddress("AK15Jet_genEta", &AK15Jet_genEta, &b_AK15Jet_genEta);
   fChain->SetBranchAddress("AK15Jet_genPhi", &AK15Jet_genPhi, &b_AK15Jet_genPhi);
   fChain->SetBranchAddress("AK15Jet_genEnergy", &AK15Jet_genEnergy, &b_AK15Jet_genEnergy);
   fChain->SetBranchAddress("AK15Jet_genPt", &AK15Jet_genPt, &b_AK15Jet_genPt);
   fChain->SetBranchAddress("AK15Jet_genPz", &AK15Jet_genPz, &b_AK15Jet_genPz);
   fChain->SetBranchAddress("AK15Jet_genMass", &AK15Jet_genMass, &b_AK15Jet_genMass);
   fChain->SetBranchAddress("AK15Jet_genNJet", &AK15Jet_genNJet, &b_AK15Jet_genNJet);
   fChain->SetBranchAddress("AK15Jet_genConstituentIdxs", &AK15Jet_genConstituentIdxs, &b_AK15Jet_genConstituentIdxs);
   fChain->SetBranchAddress("AK15Jet_genNConstituents", &AK15Jet_genNConstituents, &b_AK15Jet_genNConstituents);
   fChain->SetBranchAddress("Top_decayId", &Top_decayId, &b_Top_decayId);
   fChain->SetBranchAddress("AK4Jet_eta", &AK4Jet_eta, &b_AK4Jet_eta);
   fChain->SetBranchAddress("AK4Jet_phi", &AK4Jet_phi, &b_AK4Jet_phi);
   fChain->SetBranchAddress("AK4Jet_energy", &AK4Jet_energy, &b_AK4Jet_energy);
   fChain->SetBranchAddress("AK4Jet_pt", &AK4Jet_pt, &b_AK4Jet_pt);
   fChain->SetBranchAddress("AK4Jet_mass", &AK4Jet_mass, &b_AK4Jet_mass);
   fChain->SetBranchAddress("AK4Jet_RhIDs", &AK4Jet_RhIDs, &b_AK4Jet_RhIDs);
   fChain->SetBranchAddress("AK4Jet_NJet", &AK4Jet_NJet, &b_AK4Jet_NJet);
   fChain->SetBranchAddress("AK8Jet_eta", &AK8Jet_eta, &b_AK8Jet_eta);
   fChain->SetBranchAddress("AK8Jet_phi", &AK8Jet_phi, &b_AK8Jet_phi);
   fChain->SetBranchAddress("AK8Jet_energy", &AK8Jet_energy, &b_AK8Jet_energy);
   fChain->SetBranchAddress("AK8Jet_pt", &AK8Jet_pt, &b_AK8Jet_pt);
   fChain->SetBranchAddress("AK8Jet_mass", &AK8Jet_mass, &b_AK8Jet_mass);
   fChain->SetBranchAddress("AK8Jet_RhIDs", &AK8Jet_RhIDs, &b_AK8Jet_RhIDs);
   fChain->SetBranchAddress("AK8Jet_NJet", &AK8Jet_NJet, &b_AK8Jet_NJet);
   fChain->SetBranchAddress("AK15Jet_eta", &AK15Jet_eta, &b_AK15Jet_eta);
   fChain->SetBranchAddress("AK15Jet_phi", &AK15Jet_phi, &b_AK15Jet_phi);
   fChain->SetBranchAddress("AK15Jet_energy", &AK15Jet_energy, &b_AK15Jet_energy);
   fChain->SetBranchAddress("AK15Jet_pt", &AK15Jet_pt, &b_AK15Jet_pt);
   fChain->SetBranchAddress("AK15Jet_mass", &AK15Jet_mass, &b_AK15Jet_mass);
   fChain->SetBranchAddress("AK15Jet_RhIDs", &AK15Jet_RhIDs, &b_AK15Jet_RhIDs);
   fChain->SetBranchAddress("AK15Jet_NJet", &AK15Jet_NJet, &b_AK15Jet_NJet);
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
