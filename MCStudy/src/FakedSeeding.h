//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Mar  5 18:45:43 2015 by ROOT version 5.34/23
// from TTree XLayers/Branch are filled with X Layers Info for tracks with 12 Hits
// found on file: SelectedSearchWindows.root
//////////////////////////////////////////////////////////

#ifndef FakedSeeding_h
#define FakedSeeding_h
#include <TROOT.h>
#include "Track.h"
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
//#include "Track.h"
//#include "LinParFit.h"

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class FakedSeeding : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
  
   Double_t        X1st;
   Double_t        X2nd;
   Double_t        X3rd;
   Double_t        X4th;
   Double_t        X5th;
   Double_t        X6th;
   Double_t        Z1st;
   Double_t        Z2nd;
   Double_t        Z3rd;
   Double_t        Z4th;
   Double_t        Z5th;
   Double_t        Z6th;
   Double_t        Y1st;
   Double_t        Y2nd;
   Double_t        Y3rd;
   Double_t        Y4th;
   Double_t        Y5th;
   Double_t        Y6th;
   Double_t        Ux1st;
   Double_t        Ux2nd;
   Double_t        Ux3rd;
   Double_t        Uz1st;
   Double_t        Uz2nd;
   Double_t        Uz3rd;
   Double_t        Uy1st;
   Double_t        Uy2nd;
   Double_t        Uy3rd;
   Double_t        Vx1st;
   Double_t        Vx2nd;
   Double_t        Vx3rd;
   Double_t        Vz1st;
   Double_t        Vz2nd;
   Double_t        Vz3rd;
   Double_t        Vy1st;
   Double_t        Vy2nd;
   Double_t        Vy3rd;
   Double_t        Pz;
   Double_t        Px;
   Double_t        P;
   Double_t        Q;
   Double_t        DeltaXLastFirst;
   Double_t        Pol3_a;
   Double_t        Pol3_b;
   Double_t        Pol3_c;
   Double_t        Pol3_d;
   Double_t        Pol3_e;
   Double_t        OVTX_Z;
   Int_t           PID;

   // List of branches
   TBranch        *b_X1st;   //!
   TBranch        *b_X2nd;   //!
   TBranch        *b_X3rd;   //!
   TBranch        *b_X4th;   //!
   TBranch        *b_X5th;   //!
   TBranch        *b_X6th;   //!
   TBranch        *b_Z1st;   //!
   TBranch        *b_Z2nd;   //!
   TBranch        *b_Z3rd;   //!
   TBranch        *b_Z4th;   //!
   TBranch        *b_Z5th;   //!
   TBranch        *b_Z6th;   //!
   TBranch        *b_Y1st;   //!
   TBranch        *b_Y2nd;   //!
   TBranch        *b_Y3rd;   //!
   TBranch        *b_Y4th;   //!
   TBranch        *b_Y5th;   //!
   TBranch        *b_Y6th;   //!
   TBranch        *b_Ux1st;   //!
   TBranch        *b_Ux2nd;   //!
   TBranch        *b_Ux3rd;   //!
   TBranch        *b_Uz1st;   //!
   TBranch        *b_Uz2nd;   //!
   TBranch        *b_Uz3rd;   //!
   TBranch        *b_Uy1st;   //!
   TBranch        *b_Uy2nd;   //!
   TBranch        *b_Uy3rd;   //!
   TBranch        *b_Vx1st;   //!
   TBranch        *b_Vx2nd;   //!
   TBranch        *b_Vx3rd;   //!
   TBranch        *b_Vz1st;   //!
   TBranch        *b_Vz2nd;   //!
   TBranch        *b_Vz3rd;   //!
   TBranch        *b_Vy1st;   //!
   TBranch        *b_Vy2nd;   //!
   TBranch        *b_Vy3rd;   //!
   TBranch        *b_Pz;   //!
   TBranch        *b_Px;   //!
   TBranch        *b_P;   //!
   TBranch        *b_Q;   //!
   TBranch        *b_DeltaXLastFirst;   //!
   TBranch        *b_Pol3_a;   //!
   TBranch        *b_Pol3_b;   //!
   TBranch        *b_Pol3_c;   //!
   TBranch        *b_Pol3_d;   //!
   TBranch        *b_Pol3_e;   //!
   TBranch        *b_OVTX_Z;   //!
   TBranch        *b_PID;   //!


   FakedSeeding(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~FakedSeeding() { }
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
   virtual void    SolveParabola(Hit * hit1,Hit * hit2,Hit * hit3,double &a, double &b,double &c);
   virtual void    SolveParabola2(Hit * hit1,Hit * hit2,Hit * hit3,double &a, double &b,double &c);
   virtual double Extrapol(double Z,double a,double b,double c);//adz2
   virtual double Extrapol2(double Z1,double a1,double b1,double c1);
   bool fitXProjection(PrSeedTrack * track);
   //virtual void   FITTA(PrSeedTrack * track);
   
 private:
   int m_Total;
   int m_L0Selected;
   int m_L1Selected;
   int m_L2_Par1_Selected;
   int m_L2_Par2_Selected;
   int m_L2_Selection_ParCubic;
   float m_zReference;
   float m_dRatio;
   //double m_sigmaSmearing;
   bool doFit;
   //ClassDef(FakedSeeding,0);
};

#endif

#ifdef FakedSeeding_cxx
void FakedSeeding::Init(TTree *tree)
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
   fChain->SetBranchAddress("X1st", &X1st, &b_X1st);
   fChain->SetBranchAddress("X2nd", &X2nd, &b_X2nd);
   fChain->SetBranchAddress("X3rd", &X3rd, &b_X3rd);
   fChain->SetBranchAddress("X4th", &X4th, &b_X4th);
   fChain->SetBranchAddress("X5th", &X5th, &b_X5th);
   fChain->SetBranchAddress("X6th", &X6th, &b_X6th);
   fChain->SetBranchAddress("Z1st", &Z1st, &b_Z1st);
   fChain->SetBranchAddress("Z2nd", &Z2nd, &b_Z2nd);
   fChain->SetBranchAddress("Z3rd", &Z3rd, &b_Z3rd);
   fChain->SetBranchAddress("Z4th", &Z4th, &b_Z4th);
   fChain->SetBranchAddress("Z5th", &Z5th, &b_Z5th);
   fChain->SetBranchAddress("Z6th", &Z6th, &b_Z6th);
   fChain->SetBranchAddress("Y1st", &Y1st, &b_Y1st);
   fChain->SetBranchAddress("Y2nd", &Y2nd, &b_Y2nd);
   fChain->SetBranchAddress("Y3rd", &Y3rd, &b_Y3rd);
   fChain->SetBranchAddress("Y4th", &Y4th, &b_Y4th);
   fChain->SetBranchAddress("Y5th", &Y5th, &b_Y5th);
   fChain->SetBranchAddress("Y6th", &Y6th, &b_Y6th);
   fChain->SetBranchAddress("Ux1st", &Ux1st, &b_Ux1st);
   fChain->SetBranchAddress("Ux2nd", &Ux2nd, &b_Ux2nd);
   fChain->SetBranchAddress("Ux3rd", &Ux3rd, &b_Ux3rd);
   fChain->SetBranchAddress("Uz1st", &Uz1st, &b_Uz1st);
   fChain->SetBranchAddress("Uz2nd", &Uz2nd, &b_Uz2nd);
   fChain->SetBranchAddress("Uz3rd", &Uz3rd, &b_Uz3rd);
   fChain->SetBranchAddress("Uy1st", &Uy1st, &b_Uy1st);
   fChain->SetBranchAddress("Uy2nd", &Uy2nd, &b_Uy2nd);
   fChain->SetBranchAddress("Uy3rd", &Uy3rd, &b_Uy3rd);
   fChain->SetBranchAddress("Vx1st", &Vx1st, &b_Vx1st);
   fChain->SetBranchAddress("Vx2nd", &Vx2nd, &b_Vx2nd);
   fChain->SetBranchAddress("Vx3rd", &Vx3rd, &b_Vx3rd);
   fChain->SetBranchAddress("Vz1st", &Vz1st, &b_Vz1st);
   fChain->SetBranchAddress("Vz2nd", &Vz2nd, &b_Vz2nd);
   fChain->SetBranchAddress("Vz3rd", &Vz3rd, &b_Vz3rd);
   fChain->SetBranchAddress("Vy1st", &Vy1st, &b_Vy1st);
   fChain->SetBranchAddress("Vy2nd", &Vy2nd, &b_Vy2nd);
   fChain->SetBranchAddress("Vy3rd", &Vy3rd, &b_Vy3rd);
   fChain->SetBranchAddress("Pz", &Pz, &b_Pz);
   fChain->SetBranchAddress("Px", &Px, &b_Px);
   fChain->SetBranchAddress("P", &P, &b_P);
   fChain->SetBranchAddress("Q", &Q, &b_Q);
   fChain->SetBranchAddress("DeltaXLastFirst", &DeltaXLastFirst, &b_DeltaXLastFirst);
   fChain->SetBranchAddress("Pol3_a", &Pol3_a, &b_Pol3_a);
   fChain->SetBranchAddress("Pol3_b", &Pol3_b, &b_Pol3_b);
   fChain->SetBranchAddress("Pol3_c", &Pol3_c, &b_Pol3_c);
   fChain->SetBranchAddress("Pol3_d", &Pol3_d, &b_Pol3_d);
   fChain->SetBranchAddress("Pol3_e", &Pol3_e, &b_Pol3_e);
   fChain->SetBranchAddress("OVTX_Z", &OVTX_Z, &b_OVTX_Z);
   fChain->SetBranchAddress("PID", &PID, &b_PID);
}

Bool_t FakedSeeding::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef FakedSeeding_cxx
