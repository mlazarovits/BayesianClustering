//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Mar 12 12:43:37 2024 by ROOT version 6.28/06
// from TTree ReducedBaseSim/ReducedBaseSim
// found on file: rootfiles/simNtuples_ttbar.root
//////////////////////////////////////////////////////////
#ifndef ReducedBaseSim_hh
#define ReducedBaseSim_hh

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
   vector<double>  *Jet_genMass;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_ECALRecHit_energy;   //!
   TBranch        *b_ECALRecHit_rhx;   //!
   TBranch        *b_ECALRecHit_rhy;   //!
   TBranch        *b_ECALRecHit_rhz;   //!
   TBranch        *b_ECALRecHit_time;   //!
   TBranch        *b_ECALRecHit_eta;   //!
   TBranch        *b_ECALRecHit_phi;   //!
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
   TBranch        *b_Jet_genMass;   //!

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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("rootfiles/simNtuples_ttbar.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("rootfiles/simNtuples_ttbar.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("rootfiles/simNtuples_ttbar.root:/tree");
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
   ECALSpike_energy = 0;
   Jet_genEta = 0;
   Jet_genPhi = 0;
   Jet_genEnergy = 0;
   Jet_genPt = 0;
   Jet_genMass = 0;
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
   fChain->SetBranchAddress("Jet_genMass", &Jet_genMass, &b_Jet_genMass);
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
