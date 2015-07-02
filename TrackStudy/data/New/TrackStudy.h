//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun 30 00:16:17 2015 by ROOT version 6.02/08
// from TTree TrackStudy/Events
// found on file: SciFi-Tuple-Debug_1000.root
//////////////////////////////////////////////////////////

#ifndef TrackStudy_h
#define TrackStudy_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.

class TrackStudy : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        event;
   Double_t        run;
   Double_t        nPV;
   Int_t           FiredLayers;
   Float_t         FiredLayers_Counter[24];   //[FiredLayers]
   Int_t           CheatedSeeding_NHits;
   Int_t           MC_ass;
   Float_t         MCHit_ty[100];   //[MC_ass]
   Float_t         MCHit_tx[100];   //[MC_ass]
   Float_t         MCHit_p[100];   //[MC_ass]
   Float_t         MCHit_pathlength[100];   //[MC_ass]
   Float_t         MCHit_Assoc_X[100];   //[MC_ass]
   Float_t         MCHit_Assoc_Y[100];   //[MC_ass]
   Float_t         MCHit_Assoc_Z[100];   //[MC_ass]
   ULong64_t       N_MCHit_Assoc;
   Int_t           PrHit;
   Float_t         PrHit_LHCbID[100];   //[PrHit]
   Float_t         PrHit_Xat0[100];   //[PrHit]
   Float_t         PrHit_Zat0[100];   //[PrHit]
   Float_t         PrHit_dxDy[100];   //[PrHit]
   Float_t         PrHit_planeCode[100];   //[PrHit]
   Float_t         PrHit_zone[100];   //[PrHit]
   Float_t         PrHit_isX[100];   //[PrHit]
   Float_t         PrHit_yMin[100];   //[PrHit]
   Float_t         PrHit_yMax[100];   //[PrHit]
   Float_t         PrHit_w2[100];   //[PrHit]
   Float_t         PrHit_dzDy[100];   //[PrHit]
   Int_t           ChID;
   Float_t         ChID_Fraction[100];   //[ChID]
   Int_t           SipmCell;
   Float_t         ChID_SipmCell[100];   //[SipmCell]
   Int_t           Sip;
   Float_t         ChID_ID[100];   //[ChID]
   ULong64_t       N_PrHit_Assoc;
   Int_t           nStereo_Upper;
   Int_t           nStereo_Lower;
   Int_t           nTotal_noMultiple;
   Bool_t          Reconstructible_PrHit;
   Bool_t          Reconstructed_10Hits_6and4;
   Bool_t          Reconstructed_9Hits_5and4;
   Bool_t          Reconstructed_9Hits_3and6;
   Int_t           LowerX_noMultiple;
   Int_t           UpperX_noMultiple;
   Int_t           Stereo_noMultiple;
   Int_t           nT1_X_noMultiple;
   Int_t           nT1_UV_noMultiple;
   Int_t           nT2_X_noMultiple;
   Int_t           nT2_UV_noMultiple;
   Int_t           nT3_X_noMultiple;
   Int_t           nT3_UV_noMultiple;
   Int_t           nU;
   Int_t           nV;
   Int_t           nT1;
   Int_t           nT2;
   Int_t           nT3;
   Int_t           nUV;
   Int_t           nXUp;
   Int_t           nxDown;
   Int_t           MC;
   Float_t         MC_Hit_X[100];   //[MC]
   Float_t         MC_Hit_Y[100];   //[MC]
   Float_t         MC_Hit_Z[100];   //[MC]
   Float_t         MC_Hit_P[100];   //[MC]
   Float_t         MC_Hit_PathLenght[100];   //[MC]
   Float_t         MC_Hit_Energy[100];   //[MC]
   Float_t         MC_Hit_Particle_P[100];   //[MC]
   Float_t         MC_Hit_dxdz[100];   //[MC]
   Float_t         MC_Hit_dydz[100];   //[MC]
   Float_t         MC_time[100];   //[MC]
   ULong64_t       Number_MCHit_size;
   Bool_t          MC_Hit_hasDuplicate;
   Bool_t          fullInfo;
   Bool_t          isSeed;
   Bool_t          isLong;
   Bool_t          isDown;
   Bool_t          over5;
   Bool_t          trigger;
   Bool_t          isInVelo;
   Bool_t          isInUT;
   Bool_t          isElectron;
   Bool_t          accT;
   Bool_t          accTT;
   Bool_t          OVTXunder100mm;
   Double_t        pseudoRapidity;
   Bool_t          Eta_in25;
   Double_t        eta;
   Double_t        P;
   Double_t        Pt;
   Double_t        MC_px;
   Double_t        MC_py;
   Double_t        MC_pz;
   Double_t        MC_Ovtx_x;
   Double_t        MC_Ovtx_y;
   Double_t        MC_Ovtx_z;
   Double_t        MC_Charge;
   Bool_t          strange_Long;
   Bool_t          strange_Long_more5;
   Bool_t          hasT;
   Bool_t          isLong_more5;
   Bool_t          fromB;
   Bool_t          isLong_fromB;
   Bool_t          isLong_fromB_more5;
   Bool_t          strange_UT_T;
   Bool_t          strange_UT_T_more5;
   Bool_t          strange_UT_T_noVelo;
   Bool_t          strange_UT_T_noVelo_more5;
   Bool_t          strange_fromDB_UT_T;
   Bool_t          strange_fromDB_UT_T_noVelo;
   Bool_t          strange_fromDB_UT_T_noVelo_more5;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_run;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_FiredLayers;   //!
   TBranch        *b_FiredLayers_Counter;   //!
   TBranch        *b_CheatedSeeding_NHits;   //!
   TBranch        *b_MC_ass;   //!
   TBranch        *b_MCHit_ty;   //!
   TBranch        *b_MCHit_tx;   //!
   TBranch        *b_MCHit_p;   //!
   TBranch        *b_MCHit_pathlength;   //!
   TBranch        *b_MCHit_Assoc_X;   //!
   TBranch        *b_MCHit_Assoc_Y;   //!
   TBranch        *b_MCHit_Assoc_Z;   //!
   TBranch        *b_N_MCHit_Assoc;   //!
   TBranch        *b_PrHit;   //!
   TBranch        *b_PrHit_LHCbID;   //!
   TBranch        *b_PrHit_Xat0;   //!
   TBranch        *b_PrHit_Zat0;   //!
   TBranch        *b_PrHit_dxDy;   //!
   TBranch        *b_PrHit_planeCode;   //!
   TBranch        *b_PrHit_zone;   //!
   TBranch        *b_PrHit_isX;   //!
   TBranch        *b_PrHit_yMin;   //!
   TBranch        *b_PrHit_yMax;   //!
   TBranch        *b_PrHit_w2;   //!
   TBranch        *b_PrHit_dzDy;   //!
   TBranch        *b_ChID;   //!
   TBranch        *b_ChID_Fraction;   //!
   TBranch        *b_SipmCell;   //!
   TBranch        *b_ChID_SipmCell;   //!
   TBranch        *b_Sip;   //!
   TBranch        *b_ChID_ID;   //!
   TBranch        *b_N_PrHit_Assoc;   //!
   TBranch        *b_nStereo_Upper;   //!
   TBranch        *b_nStereo_Lower;   //!
   TBranch        *b_nTotal_noMultiple;   //!
   TBranch        *b_Reconstructible_PrHit;   //!
   TBranch        *b_Reconstructed_10Hits_6and4;   //!
   TBranch        *b_Reconstructed_9Hits_5and4;   //!
   TBranch        *b_Reconstructed_9Hits_3and6;   //!
   TBranch        *b_LowerX_noMultiple;   //!
   TBranch        *b_UpperX_noMultiple;   //!
   TBranch        *b_Stereo_noMultiple;   //!
   TBranch        *b_nT1_X_noMultiple;   //!
   TBranch        *b_nT1_UV_noMultiple;   //!
   TBranch        *b_nT2_X_noMultiple;   //!
   TBranch        *b_nT2_UV_noMultiple;   //!
   TBranch        *b_nT3_X_noMultiple;   //!
   TBranch        *b_nT3_UV_noMultiple;   //!
   TBranch        *b_nU;   //!
   TBranch        *b_nV;   //!
   TBranch        *b_nT1;   //!
   TBranch        *b_nT2;   //!
   TBranch        *b_nT3;   //!
   TBranch        *b_nUV;   //!
   TBranch        *b_nXUp;   //!
   TBranch        *b_nxDown;   //!
   TBranch        *b_MC;   //!
   TBranch        *b_MC_Hit_X;   //!
   TBranch        *b_MC_Hit_Y;   //!
   TBranch        *b_MC_Hit_Z;   //!
   TBranch        *b_MC_Hit_P;   //!
   TBranch        *b_MC_Hit_PathLenght;   //!
   TBranch        *b_MC_Hit_Energy;   //!
   TBranch        *b_MC_Hit_Particle_P;   //!
   TBranch        *b_MC_Hit_dxdz;   //!
   TBranch        *b_MC_Hit_dydz;   //!
   TBranch        *b_MC_time;   //!
   TBranch        *b_Number_MCHit_size;   //!
   TBranch        *b_MC_Hit_hasDuplicate;   //!
   TBranch        *b_fullInfo;   //!
   TBranch        *b_isSeed;   //!
   TBranch        *b_isLong;   //!
   TBranch        *b_isDown;   //!
   TBranch        *b_over5;   //!
   TBranch        *b_trigger;   //!
   TBranch        *b_isInVelo;   //!
   TBranch        *b_isInUT;   //!
   TBranch        *b_isElectron;   //!
   TBranch        *b_accT;   //!
   TBranch        *b_accTT;   //!
   TBranch        *b_OVTXunder100mm;   //!
   TBranch        *b_pseudoRapidity;   //!
   TBranch        *b_Eta_in25;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_P;   //!
   TBranch        *b_Pt;   //!
   TBranch        *b_MC_px;   //!
   TBranch        *b_MC_py;   //!
   TBranch        *b_MC_pz;   //!
   TBranch        *b_MC_Ovtx_x;   //!
   TBranch        *b_MC_Ovtx_y;   //!
   TBranch        *b_MC_Ovtx_z;   //!
   TBranch        *b_MC_Charge;   //!
   TBranch        *b_strange_Long;   //!
   TBranch        *b_strange_Long_more5;   //!
   TBranch        *b_hasT;   //!
   TBranch        *b_isLong_more5;   //!
   TBranch        *b_fromB;   //!
   TBranch        *b_isLong_fromB;   //!
   TBranch        *b_isLong_fromB_more5;   //!
   TBranch        *b_strange_UT_T;   //!
   TBranch        *b_strange_UT_T_more5;   //!
   TBranch        *b_strange_UT_T_noVelo;   //!
   TBranch        *b_strange_UT_T_noVelo_more5;   //!
   TBranch        *b_strange_fromDB_UT_T;   //!
   TBranch        *b_strange_fromDB_UT_T_noVelo;   //!
   TBranch        *b_strange_fromDB_UT_T_noVelo_more5;   //!

   TrackStudy(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~TrackStudy() { }
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

   ClassDef(TrackStudy,0);
};

#endif

#ifdef TrackStudy_cxx
void TrackStudy::Init(TTree *tree)
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

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("FiredLayers", &FiredLayers, &b_FiredLayers);
   fChain->SetBranchAddress("FiredLayers_Counter", FiredLayers_Counter, &b_FiredLayers_Counter);
   fChain->SetBranchAddress("CheatedSeeding_NHits", &CheatedSeeding_NHits, &b_CheatedSeeding_NHits);
   fChain->SetBranchAddress("MC_ass", &MC_ass, &b_MC_ass);
   fChain->SetBranchAddress("MCHit_ty", MCHit_ty, &b_MCHit_ty);
   fChain->SetBranchAddress("MCHit_tx", MCHit_tx, &b_MCHit_tx);
   fChain->SetBranchAddress("MCHit_p", MCHit_p, &b_MCHit_p);
   fChain->SetBranchAddress("MCHit_pathlength", MCHit_pathlength, &b_MCHit_pathlength);
   fChain->SetBranchAddress("MCHit_Assoc_X", MCHit_Assoc_X, &b_MCHit_Assoc_X);
   fChain->SetBranchAddress("MCHit_Assoc_Y", MCHit_Assoc_Y, &b_MCHit_Assoc_Y);
   fChain->SetBranchAddress("MCHit_Assoc_Z", MCHit_Assoc_Z, &b_MCHit_Assoc_Z);
   fChain->SetBranchAddress("N_MCHit_Assoc", &N_MCHit_Assoc, &b_N_MCHit_Assoc);
   fChain->SetBranchAddress("PrHit", &PrHit, &b_PrHit);
   fChain->SetBranchAddress("PrHit_LHCbID", PrHit_LHCbID, &b_PrHit_LHCbID);
   fChain->SetBranchAddress("PrHit_Xat0", PrHit_Xat0, &b_PrHit_Xat0);
   fChain->SetBranchAddress("PrHit_Zat0", PrHit_Zat0, &b_PrHit_Zat0);
   fChain->SetBranchAddress("PrHit_dxDy", PrHit_dxDy, &b_PrHit_dxDy);
   fChain->SetBranchAddress("PrHit_planeCode", PrHit_planeCode, &b_PrHit_planeCode);
   fChain->SetBranchAddress("PrHit_zone", PrHit_zone, &b_PrHit_zone);
   fChain->SetBranchAddress("PrHit_isX", PrHit_isX, &b_PrHit_isX);
   fChain->SetBranchAddress("PrHit_yMin", PrHit_yMin, &b_PrHit_yMin);
   fChain->SetBranchAddress("PrHit_yMax", PrHit_yMax, &b_PrHit_yMax);
   fChain->SetBranchAddress("PrHit_w2", PrHit_w2, &b_PrHit_w2);
   fChain->SetBranchAddress("PrHit_dzDy", PrHit_dzDy, &b_PrHit_dzDy);
   fChain->SetBranchAddress("ChID", &ChID, &b_ChID);
   fChain->SetBranchAddress("ChID_Fraction", ChID_Fraction, &b_ChID_Fraction);
   fChain->SetBranchAddress("SipmCell", &SipmCell, &b_SipmCell);
   fChain->SetBranchAddress("ChID_SipmCell", ChID_SipmCell, &b_ChID_SipmCell);
   fChain->SetBranchAddress("Sip", &Sip, &b_Sip);
   fChain->SetBranchAddress("ChID_ID", ChID_ID, &b_ChID_ID);
   fChain->SetBranchAddress("N_PrHit_Assoc", &N_PrHit_Assoc, &b_N_PrHit_Assoc);
   fChain->SetBranchAddress("nStereo_Upper", &nStereo_Upper, &b_nStereo_Upper);
   fChain->SetBranchAddress("nStereo_Lower", &nStereo_Lower, &b_nStereo_Lower);
   fChain->SetBranchAddress("nTotal_noMultiple", &nTotal_noMultiple, &b_nTotal_noMultiple);
   fChain->SetBranchAddress("Reconstructible_PrHit", &Reconstructible_PrHit, &b_Reconstructible_PrHit);
   fChain->SetBranchAddress("Reconstructed_10Hits_6and4", &Reconstructed_10Hits_6and4, &b_Reconstructed_10Hits_6and4);
   fChain->SetBranchAddress("Reconstructed_9Hits_5and4", &Reconstructed_9Hits_5and4, &b_Reconstructed_9Hits_5and4);
   fChain->SetBranchAddress("Reconstructed_9Hits_3and6", &Reconstructed_9Hits_3and6, &b_Reconstructed_9Hits_3and6);
   fChain->SetBranchAddress("LowerX_noMultiple", &LowerX_noMultiple, &b_LowerX_noMultiple);
   fChain->SetBranchAddress("UpperX_noMultiple", &UpperX_noMultiple, &b_UpperX_noMultiple);
   fChain->SetBranchAddress("Stereo_noMultiple", &Stereo_noMultiple, &b_Stereo_noMultiple);
   fChain->SetBranchAddress("nT1_X_noMultiple", &nT1_X_noMultiple, &b_nT1_X_noMultiple);
   fChain->SetBranchAddress("nT1_UV_noMultiple", &nT1_UV_noMultiple, &b_nT1_UV_noMultiple);
   fChain->SetBranchAddress("nT2_X_noMultiple", &nT2_X_noMultiple, &b_nT2_X_noMultiple);
   fChain->SetBranchAddress("nT2_UV_noMultiple", &nT2_UV_noMultiple, &b_nT2_UV_noMultiple);
   fChain->SetBranchAddress("nT3_X_noMultiple", &nT3_X_noMultiple, &b_nT3_X_noMultiple);
   fChain->SetBranchAddress("nT3_UV_noMultiple", &nT3_UV_noMultiple, &b_nT3_UV_noMultiple);
   fChain->SetBranchAddress("nU", &nU, &b_nU);
   fChain->SetBranchAddress("nV", &nV, &b_nV);
   fChain->SetBranchAddress("nT1", &nT1, &b_nT1);
   fChain->SetBranchAddress("nT2", &nT2, &b_nT2);
   fChain->SetBranchAddress("nT3", &nT3, &b_nT3);
   fChain->SetBranchAddress("nUV", &nUV, &b_nUV);
   fChain->SetBranchAddress("nXUp", &nXUp, &b_nXUp);
   fChain->SetBranchAddress("nxDown", &nxDown, &b_nxDown);
   fChain->SetBranchAddress("MC", &MC, &b_MC);
   fChain->SetBranchAddress("MC_Hit_X", MC_Hit_X, &b_MC_Hit_X);
   fChain->SetBranchAddress("MC_Hit_Y", MC_Hit_Y, &b_MC_Hit_Y);
   fChain->SetBranchAddress("MC_Hit_Z", MC_Hit_Z, &b_MC_Hit_Z);
   fChain->SetBranchAddress("MC_Hit_P", MC_Hit_P, &b_MC_Hit_P);
   fChain->SetBranchAddress("MC_Hit_PathLenght", MC_Hit_PathLenght, &b_MC_Hit_PathLenght);
   fChain->SetBranchAddress("MC_Hit_Energy", MC_Hit_Energy, &b_MC_Hit_Energy);
   fChain->SetBranchAddress("MC_Hit_Particle_P", MC_Hit_Particle_P, &b_MC_Hit_Particle_P);
   fChain->SetBranchAddress("MC_Hit_dxdz", MC_Hit_dxdz, &b_MC_Hit_dxdz);
   fChain->SetBranchAddress("MC_Hit_dydz", MC_Hit_dydz, &b_MC_Hit_dydz);
   fChain->SetBranchAddress("MC_time", MC_time, &b_MC_time);
   fChain->SetBranchAddress("Number_MCHit_size", &Number_MCHit_size, &b_Number_MCHit_size);
   fChain->SetBranchAddress("MC_Hit_hasDuplicate", &MC_Hit_hasDuplicate, &b_MC_Hit_hasDuplicate);
   fChain->SetBranchAddress("fullInfo", &fullInfo, &b_fullInfo);
   fChain->SetBranchAddress("isSeed", &isSeed, &b_isSeed);
   fChain->SetBranchAddress("isLong", &isLong, &b_isLong);
   fChain->SetBranchAddress("isDown", &isDown, &b_isDown);
   fChain->SetBranchAddress("over5", &over5, &b_over5);
   fChain->SetBranchAddress("trigger", &trigger, &b_trigger);
   fChain->SetBranchAddress("isInVelo", &isInVelo, &b_isInVelo);
   fChain->SetBranchAddress("isInUT", &isInUT, &b_isInUT);
   fChain->SetBranchAddress("isElectron", &isElectron, &b_isElectron);
   fChain->SetBranchAddress("accT", &accT, &b_accT);
   fChain->SetBranchAddress("accTT", &accTT, &b_accTT);
   fChain->SetBranchAddress("OVTXunder100mm", &OVTXunder100mm, &b_OVTXunder100mm);
   fChain->SetBranchAddress("pseudoRapidity", &pseudoRapidity, &b_pseudoRapidity);
   fChain->SetBranchAddress("Eta_in25", &Eta_in25, &b_Eta_in25);
   fChain->SetBranchAddress("eta", &eta, &b_eta);
   fChain->SetBranchAddress("P", &P, &b_P);
   fChain->SetBranchAddress("Pt", &Pt, &b_Pt);
   fChain->SetBranchAddress("MC_px", &MC_px, &b_MC_px);
   fChain->SetBranchAddress("MC_py", &MC_py, &b_MC_py);
   fChain->SetBranchAddress("MC_pz", &MC_pz, &b_MC_pz);
   fChain->SetBranchAddress("MC_Ovtx_x", &MC_Ovtx_x, &b_MC_Ovtx_x);
   fChain->SetBranchAddress("MC_Ovtx_y", &MC_Ovtx_y, &b_MC_Ovtx_y);
   fChain->SetBranchAddress("MC_Ovtx_z", &MC_Ovtx_z, &b_MC_Ovtx_z);
   fChain->SetBranchAddress("MC_Charge", &MC_Charge, &b_MC_Charge);
   fChain->SetBranchAddress("strange_Long", &strange_Long, &b_strange_Long);
   fChain->SetBranchAddress("strange_Long_more5", &strange_Long_more5, &b_strange_Long_more5);
   fChain->SetBranchAddress("hasT", &hasT, &b_hasT);
   fChain->SetBranchAddress("isLong_more5", &isLong_more5, &b_isLong_more5);
   fChain->SetBranchAddress("fromB", &fromB, &b_fromB);
   fChain->SetBranchAddress("isLong_fromB", &isLong_fromB, &b_isLong_fromB);
   fChain->SetBranchAddress("isLong_fromB_more5", &isLong_fromB_more5, &b_isLong_fromB_more5);
   fChain->SetBranchAddress("strange_UT_T", &strange_UT_T, &b_strange_UT_T);
   fChain->SetBranchAddress("strange_UT_T_more5", &strange_UT_T_more5, &b_strange_UT_T_more5);
   fChain->SetBranchAddress("strange_UT_T_noVelo", &strange_UT_T_noVelo, &b_strange_UT_T_noVelo);
   fChain->SetBranchAddress("strange_UT_T_noVelo_more5", &strange_UT_T_noVelo_more5, &b_strange_UT_T_noVelo_more5);
   fChain->SetBranchAddress("strange_fromDB_UT_T", &strange_fromDB_UT_T, &b_strange_fromDB_UT_T);
   fChain->SetBranchAddress("strange_fromDB_UT_T_noVelo", &strange_fromDB_UT_T_noVelo, &b_strange_fromDB_UT_T_noVelo);
   fChain->SetBranchAddress("strange_fromDB_UT_T_noVelo_more5", &strange_fromDB_UT_T_noVelo_more5, &b_strange_fromDB_UT_T_noVelo_more5);
}
void TrackStudy::fitXProjection(PrSeedTrack *track);
Bool_t TrackStudy::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef TrackStudy_cxx
