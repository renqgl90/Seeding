//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jun 26 04:22:09 2015 by ROOT version 5.34/30
// from TTree ClusterMCHitAndTrackStudy/Events
// found on file: SciFi-Tuple-Debug_1000.root
//////////////////////////////////////////////////////////

#ifndef PerfectTracking_ClusterStudy_h
#define PerfectTracking_ClusterStudy_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class PerfectTracking_ClusterStudy : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Int_t           numberMCHitToCluster;
   UInt_t          ClusterCharge;
   UInt_t          ClusterSize;
   Double_t        ClusterFraction;
   UInt_t          ClusterChannelID;
   UInt_t          ClusterChannelIDSipmCell;
   UInt_t          ClusterChannelIDSipmID;
   Bool_t          ClusterChannelIDMat;
   UInt_t          ClusterChannelIDModule;
   UInt_t          ClusterChannelLayer;
   UInt_t          ClusterChannelQuarter;
   Bool_t          isX;
   Bool_t          isU;
   Bool_t          isV;
   Bool_t          isT1;
   Bool_t          isT2;
   Bool_t          isT3;
   Int_t           layer;
   Int_t           MCParticlePID;
   Bool_t          MCParticleIsLong;
   Bool_t          MCParticleIsSeed;
   Double_t        MCParticleP;
   Double_t        MCParticlePt;
   Double_t        MCParticleGamma;
   Double_t        MCParticleBeta;
   Double_t        MCParticleVirtualMass;
   Double_t        MCParticleCharge;
   Double_t        MCParticlePseudoRapidity;
   Bool_t          MCParticleAccT;
   Bool_t          MCParticleAccTT;
   UInt_t          zone;
   Double_t        MCHit_X;
   Double_t        MCHit_Y;
   Double_t        MCHit_Z;
   Float_t         PrHit_XatYEq0;
   Float_t         PrHit_ZatYEq0;
   Float_t         PrHit_XatYMCHit;
   Float_t         PrHit_ZatYMCHit;
   Double_t        XResidual;
   Double_t        ZResidual;
   Float_t         Hit_dxDy;
   Float_t         Hit_werr;
   Float_t         Hit_w;
   Float_t         Hit_coord;
   Bool_t          Hit_isUsed;
   Float_t         Hit_yMin;
   Float_t         Hit_yMax;
   Int_t           Hit_Zone;
   Float_t         Hit_dzDy_manually;
   Double_t        TrackChi2NDOFSeed;
   Int_t           TrackNDoFSeed;
   Double_t        TrackChi2NDoFFwd;
   Int_t           TrackNDoFFwd;
   Bool_t          assocSeed;
   Bool_t          isUsedBySeed;
   Bool_t          assocFwd;
   Bool_t          isUsedByFwd;

   // List of branches
   TBranch        *b_numberMCHitToCluster;   //!
   TBranch        *b_ClusterCharge;   //!
   TBranch        *b_ClusterSize;   //!
   TBranch        *b_ClusterFraction;   //!
   TBranch        *b_ClusterChannelID;   //!
   TBranch        *b_ClusterChannelIDSipmCell;   //!
   TBranch        *b_ClusterChannelIDSipmID;   //!
   TBranch        *b_ClusterChannelIDMat;   //!
   TBranch        *b_ClusterChannelIDModule;   //!
   TBranch        *b_ClusterChannelLayer;   //!
   TBranch        *b_ClusterChannelQuarter;   //!
   TBranch        *b_isX;   //!
   TBranch        *b_isU;   //!
   TBranch        *b_isV;   //!
   TBranch        *b_isT1;   //!
   TBranch        *b_isT2;   //!
   TBranch        *b_isT3;   //!
   TBranch        *b_layer;   //!
   TBranch        *b_MCParticlePID;   //!
   TBranch        *b_MCParticleIsLong;   //!
   TBranch        *b_MCParticleIsSeed;   //!
   TBranch        *b_MCParticleP;   //!
   TBranch        *b_MCParticlePt;   //!
   TBranch        *b_MCParticleGamma;   //!
   TBranch        *b_MCParticleBeta;   //!
   TBranch        *b_MCParticleVirtualMass;   //!
   TBranch        *b_MCParticleCharge;   //!
   TBranch        *b_MCParticlePseudoRapidity;   //!
   TBranch        *b_MCParticleAccT;   //!
   TBranch        *b_MCParticleAccTT;   //!
   TBranch        *b_zone;   //!
   TBranch        *b_MCHit_X;   //!
   TBranch        *b_MCHit_Y;   //!
   TBranch        *b_MCHit_Z;   //!
   TBranch        *b_PrHit_XatYEq0;   //!
   TBranch        *b_PrHit_ZatYEq0;   //!
   TBranch        *b_PrHit_XatYMCHit;   //!
   TBranch        *b_PrHit_ZatYMCHit;   //!
   TBranch        *b_XResidual;   //!
   TBranch        *b_ZResidual;   //!
   TBranch        *b_Hit_dxDy;   //!
   TBranch        *b_Hit_werr;   //!
   TBranch        *b_Hit_w;   //!
   TBranch        *b_Hit_coord;   //!
   TBranch        *b_Hit_isUsed;   //!
   TBranch        *b_Hit_yMin;   //!
   TBranch        *b_Hit_yMax;   //!
   TBranch        *b_Hit_Zone;   //!
   TBranch        *b_Hit_dzDy_manually;   //!
   TBranch        *b_TrackChi2NDOFSeed;   //!
   TBranch        *b_TrackNDoFSeed;   //!
   TBranch        *b_TrackChi2NDoFFwd;   //!
   TBranch        *b_TrackNDoFFwd;   //!
   TBranch        *b_assocSeed;   //!
   TBranch        *b_isUsedBySeed;   //!
   TBranch        *b_assocFwd;   //!
   TBranch        *b_isUsedByFwd;   //!

   PerfectTracking_ClusterStudy(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~PerfectTracking_ClusterStudy() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(PerfectTracking_ClusterStudy,0);
};

#endif

#ifdef PerfectTracking_ClusterStudy_cxx
void PerfectTracking_ClusterStudy::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("numberMCHitToCluster", &numberMCHitToCluster, &b_numberMCHitToCluster);
   fChain->SetBranchAddress("ClusterCharge", &ClusterCharge, &b_ClusterCharge);
   fChain->SetBranchAddress("ClusterSize", &ClusterSize, &b_ClusterSize);
   fChain->SetBranchAddress("ClusterFraction", &ClusterFraction, &b_ClusterFraction);
   fChain->SetBranchAddress("ClusterChannelID", &ClusterChannelID, &b_ClusterChannelID);
   fChain->SetBranchAddress("ClusterChannelIDSipmCell", &ClusterChannelIDSipmCell, &b_ClusterChannelIDSipmCell);
   fChain->SetBranchAddress("ClusterChannelIDSipmID", &ClusterChannelIDSipmID, &b_ClusterChannelIDSipmID);
   fChain->SetBranchAddress("ClusterChannelIDMat", &ClusterChannelIDMat, &b_ClusterChannelIDMat);
   fChain->SetBranchAddress("ClusterChannelIDModule", &ClusterChannelIDModule, &b_ClusterChannelIDModule);
   fChain->SetBranchAddress("ClusterChannelLayer", &ClusterChannelLayer, &b_ClusterChannelLayer);
   fChain->SetBranchAddress("ClusterChannelQuarter", &ClusterChannelQuarter, &b_ClusterChannelQuarter);
   fChain->SetBranchAddress("isX", &isX, &b_isX);
   fChain->SetBranchAddress("isU", &isU, &b_isU);
   fChain->SetBranchAddress("isV", &isV, &b_isV);
   fChain->SetBranchAddress("isT1", &isT1, &b_isT1);
   fChain->SetBranchAddress("isT2", &isT2, &b_isT2);
   fChain->SetBranchAddress("isT3", &isT3, &b_isT3);
   fChain->SetBranchAddress("layer", &layer, &b_layer);
   fChain->SetBranchAddress("MCParticlePID", &MCParticlePID, &b_MCParticlePID);
   fChain->SetBranchAddress("MCParticleIsLong", &MCParticleIsLong, &b_MCParticleIsLong);
   fChain->SetBranchAddress("MCParticleIsSeed", &MCParticleIsSeed, &b_MCParticleIsSeed);
   fChain->SetBranchAddress("MCParticleP", &MCParticleP, &b_MCParticleP);
   fChain->SetBranchAddress("MCParticlePt", &MCParticlePt, &b_MCParticlePt);
   fChain->SetBranchAddress("MCParticleGamma", &MCParticleGamma, &b_MCParticleGamma);
   fChain->SetBranchAddress("MCParticleBeta", &MCParticleBeta, &b_MCParticleBeta);
   fChain->SetBranchAddress("MCParticleVirtualMass", &MCParticleVirtualMass, &b_MCParticleVirtualMass);
   fChain->SetBranchAddress("MCParticleCharge", &MCParticleCharge, &b_MCParticleCharge);
   fChain->SetBranchAddress("MCParticlePseudoRapidity", &MCParticlePseudoRapidity, &b_MCParticlePseudoRapidity);
   fChain->SetBranchAddress("MCParticleAccT", &MCParticleAccT, &b_MCParticleAccT);
   fChain->SetBranchAddress("MCParticleAccTT", &MCParticleAccTT, &b_MCParticleAccTT);
   fChain->SetBranchAddress("zone", &zone, &b_zone);
   fChain->SetBranchAddress("MCHit_X", &MCHit_X, &b_MCHit_X);
   fChain->SetBranchAddress("MCHit_Y", &MCHit_Y, &b_MCHit_Y);
   fChain->SetBranchAddress("MCHit_Z", &MCHit_Z, &b_MCHit_Z);
   fChain->SetBranchAddress("PrHit_XatYEq0", &PrHit_XatYEq0, &b_PrHit_XatYEq0);
   fChain->SetBranchAddress("PrHit_ZatYEq0", &PrHit_ZatYEq0, &b_PrHit_ZatYEq0);
   fChain->SetBranchAddress("PrHit_XatYMCHit", &PrHit_XatYMCHit, &b_PrHit_XatYMCHit);
   fChain->SetBranchAddress("PrHit_ZatYMCHit", &PrHit_ZatYMCHit, &b_PrHit_ZatYMCHit);
   fChain->SetBranchAddress("XResidual", &XResidual, &b_XResidual);
   fChain->SetBranchAddress("ZResidual", &ZResidual, &b_ZResidual);
   fChain->SetBranchAddress("Hit_dxDy", &Hit_dxDy, &b_Hit_dxDy);
   fChain->SetBranchAddress("Hit_werr", &Hit_werr, &b_Hit_werr);
   fChain->SetBranchAddress("Hit_w", &Hit_w, &b_Hit_w);
   fChain->SetBranchAddress("Hit_coord", &Hit_coord, &b_Hit_coord);
   fChain->SetBranchAddress("Hit_isUsed", &Hit_isUsed, &b_Hit_isUsed);
   fChain->SetBranchAddress("Hit_yMin", &Hit_yMin, &b_Hit_yMin);
   fChain->SetBranchAddress("Hit_yMax", &Hit_yMax, &b_Hit_yMax);
   fChain->SetBranchAddress("Hit_Zone", &Hit_Zone, &b_Hit_Zone);
   fChain->SetBranchAddress("Hit_dzDy_manually", &Hit_dzDy_manually, &b_Hit_dzDy_manually);
   fChain->SetBranchAddress("TrackChi2NDOFSeed", &TrackChi2NDOFSeed, &b_TrackChi2NDOFSeed);
   fChain->SetBranchAddress("TrackNDoFSeed", &TrackNDoFSeed, &b_TrackNDoFSeed);
   fChain->SetBranchAddress("TrackChi2NDoFFwd", &TrackChi2NDoFFwd, &b_TrackChi2NDoFFwd);
   fChain->SetBranchAddress("TrackNDoFFwd", &TrackNDoFFwd, &b_TrackNDoFFwd);
   fChain->SetBranchAddress("assocSeed", &assocSeed, &b_assocSeed);
   fChain->SetBranchAddress("isUsedBySeed", &isUsedBySeed, &b_isUsedBySeed);
   fChain->SetBranchAddress("assocFwd", &assocFwd, &b_assocFwd);
   fChain->SetBranchAddress("isUsedByFwd", &isUsedByFwd, &b_isUsedByFwd);
}

Bool_t PerfectTracking_ClusterStudy::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef PerfectTracking_ClusterStudy_cxx
