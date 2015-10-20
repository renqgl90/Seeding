#define TrackStudy_cxx
// The class definition in TrackStudy.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("TrackStudy.C")
// Root > T->Process("TrackStudy.C","some options")
// Root > T->Process("TrackStudy.C+")
//
#include <functional>
#include "TrackStudy.h"
#include <TH2.h>
//#include <TStyle.h>
#include "PatHit.h"
#include "Track.h"
#include "FTCluster.h"
#include "MCHit.h"
#include <iostream>
#include <cmath>
#include <cerrno>
#include <cfenv>
#include <cstring>
#include <vector>
#include <string>
#include <iostream>
#include <TROOT.h>
#include <TMath.h>

//-------------------Declaration of the Ntuple and TFile to be filled with the Chi2 Fits on MCHits
char NameFile[100]="Analysed_Track.root";
TFile*file = new TFile(NameFile,"RECREATE");
TTree * t1 = new TTree("treee","Track_Fit");
//Parabolic XZ
//Track n Hits

Int_t Track_Zone;

//------------------------------------------Track Fit Business with MCHits
//------------ Parabolic XZ
Double_t Chi2_ParabolaXZ;
Double_t ax_par;
Double_t bx_par;
Double_t cx_par;
//------------ Cubic XZ
Double_t Chi2_CubicXZ;
Double_t ax_cub;
Double_t bx_cub;
Double_t cx_cub;
Double_t dx_cub;

Int_t worstPlane;

Double_t X_T3;
Double_t X_T1;
Double_t Coord_T3;
Double_t Coord_T1;
Double_t Z_T3;
Double_t Z_T1;
Double_t Y_T3 ;
Double_t Y_T1 ;
//----------- Line Y
Double_t Chi2_LineY;
Double_t ay_line;
Double_t by_line;

//----------- Parabolic Y
Double_t Chi2_ParabolaY;
Double_t ay_par;
Double_t by_par;
Double_t cy_par;

//-----------Track Flags

Double_t X0Back;

Double_t Momentum;
Double_t Px;
Double_t Py;
Double_t Track_eta;
Double_t Pt;
Double_t Track_Phi;
Double_t Track_Pt;
Int_t MCParticleiD;
Bool_t IsElectron;
Bool_t ISLONG;

//-------------- Fit XZ with correction of dRatio

//Fit with parametrised dRatio
Double_t Chi2_XZDRATIO;
Double_t Chi2_XZDRATIOFixed;


//------------- Fit XZ with a parabola but with real PrHit : no clone removal !
Double_t ax_XStepHitFit;
Double_t bx_XStepHitFit;
Double_t cx_XStepHitFit;
Double_t Chi2_XStepFit;
Int_t NXHits;
Double_t Chi2_perDoF_XStepFit;
Int_t NClustersX;
Int_t NSplittedX;
Int_t NSplittedUV;
Int_t NPv;
Int_t NUV;


//Coord hough transformation
Double_t minY;
Double_t maxY;
Double_t minCoord;
Double_t maxCoord;
Double_t avgCoord;

TH1D * Chi2CloneX = new TH1D("Chi2ClonesX","Chi2DoFClonesX;Counts;Chi2",100,0,300);
TH1D * Chi2NoCloneX = new TH1D("Chi2Clone","Chi2",100,0,300);

//NHIts study
//Glancing of Z position
TH1D *zMCHitClones = new TH1D("zMCHitClones","zMCHitClones;z[mm];Counts",1000,7600,9500);
TH2D *zMCHitClonevsP = new TH2D("zMCHitClones_Vs_P","zMCHit Clones vs P;z[mm];P[MeV]",200,7600,9500,200,1000,12000);
TH2D *xvsyClones = new TH2D("xVsy_Clones","x vs y Clones;x[mm];y[mm]",300,-3000,3000,300,-2500,2500);
//TH3D *MCHitCloneXYZ = new TH3D("X Y Z MCHitClone","")
//Pair UV Study
TH1D * DeltaYT1 = new TH1D("DeltaYT1","DeltaYT1",400,-100,100);
TH2D * DeltaYT1vsyT1 = new TH2D("DeltaYT1vsyT1","DeltaYT1vsyT1",400,-3000,3000,100,-200,200);
TH1D * DeltaXCombinedT1 = new TH1D("DeltaXCombinedT1","DeltaXCombinedT1",100,-20,20);


TH1D * DeltaRemaining_Case0 = new TH1D("DeltaRemaing_Case0","DeltaRemaining_Case0",200,-5,5);



Double_t DeltaXBestCombT1;
Double_t DeltaXBestCombCorrSlopeY;
Double_t TyPairT1;
Double_t YPairT1;
Double_t XPairT1;
Double_t DeltaUT2Projected;
Double_t DeltaVT2Projected;
Double_t YUT2Projected;
Double_t YVT2Projected;
Double_t OVTX_Z;
Bool_t ISDOWN;
Bool_t PhysicsInterest;
Bool_t Case2Accepted;

Double_t Chi2X_Case2;

Int_t nRemSeed1_Case0;
Int_t nRemSeed2_Case0;
Int_t nRemSeed1_Case1;
Int_t nRemSeed2_Case1;
Int_t nRemSeed1_Case2;
Int_t nRemSeed2_Case2;
//XATZ STUDY

Double_t DeltaSeed1_Case0[3];
Double_t DeltaSeed2_Case0[3];
Double_t DeltaSeed1_Case1[3];
Double_t DeltaSeed2_Case1[3];
Double_t DeltaSeed1_Case2[3];
Double_t DeltaSeed2_Case2[3];

                       


//Case 0
Bool_t Case0Accepted;
Double_t Chi2X_Case0;
Double_t Chi2X_Case0Cub;

int NX_Case0;
int NUV_Case0;
Double_t ax_FitHitCase0_Cub;
Double_t bx_FitHitCase0_Cub;
Double_t cx_FitHitCase0_Cub;
Double_t dx_FitHitCase0_Cub;

Double_t ax_FitHitCase0;
Double_t bx_FitHitCase0;
Double_t cx_FitHitCase0;
Double_t ax_fullFitCase0;
Double_t bx_fullFitCase0;
Double_t cx_fullFitCase0;
Double_t ay_fullFitCase0;
Double_t by_fullFitCase0;
Double_t Chi2_fullFitCase0;
Double_t MaxChi2_fullFitCase0;
int NHits_fullFitCase0;
Double_t y0_fullFitCase0;
Double_t tx_inf_Case0;
Double_t x0_Case0;
Double_t x_inf_Case0;
Double_t Del_x_inf_Case0;
Double_t tx_pickedcombination_Case0;
Double_t DelxProjectedSeed1_Case0;
Double_t DelxProjectedSeed2_Case0;
Double_t MaxChi2Hit_Case0;
Double_t ax_Case0_Seed1;
Double_t bx_Case0_Seed1;
Double_t cx_Case0_Seed1;
Double_t ax_Case0_Seed2;
Double_t bx_Case0_Seed2;
Double_t cx_Case0_Seed2;


Double_t FitLineY_ay_case0;
Double_t FitLineY_by_case0;
Double_t FitLineY_Chi2_case0;
Double_t FitLineY_Chi2DoF_case0;
Double_t FitParY_Chi2_case0;
Double_t FitParY_Chi2DoF_case0;


Double_t FitLineY_ay_case1;
Double_t FitLineY_by_case1;
Double_t FitLineY_Chi2_case1;
Double_t FitLineY_Chi2DoF_case1;
Double_t FitParY_Chi2_case1;
Double_t FitParY_Chi2DoF_case1;
Double_t FitLineY_ay_case2;
Double_t FitLineY_by_case2;
Double_t FitLineY_Chi2_case2;
Double_t FitLineY_Chi2DoF_case2;
Double_t FitParY_Chi2_case2;
Double_t FitParY_Chi2DoF_case2;
Int_t WorstPlane;
Int_t nbHitsInWorst;
//Case 1
Bool_t Case1Accepted;
Double_t Chi2X_Case1;
int NX_Case1;
int NUV_Case1;
Double_t ax_FitHitCase1;
Double_t bx_FitHitCase1;
Double_t cx_FitHitCase1;
Double_t ax_fullFitCase1;
Double_t bx_fullFitCase1;
Double_t cx_fullFitCase1;
Double_t ay_fullFitCase1;
Double_t by_fullFitCase1;
Double_t Chi2_fullFitCase1;
Double_t MaxChi2_fullFitCase1;
int NHits_fullFitCase1;
Double_t y0_fullFitCase1;
Double_t tx_inf_Case1;
Double_t x0_Case1;
Double_t x_inf_Case1;
Double_t Del_x_inf_Case1;
Double_t tx_pickedcombination_Case1;
Double_t DelxProjectedSeed1_Case1;
Double_t DelxProjectedSeed2_Case1;
Double_t MaxChi2Hit_Case1;
Double_t ax_Case1_Seed1;
Double_t bx_Case1_Seed1;
Double_t cx_Case1_Seed1;
Double_t ax_Case1_Seed2;
Double_t bx_Case1_Seed2;
Double_t cx_Case1_Seed2;

//Case 2
int NX_Case2;
int NUV_Case2;
Double_t ax_FitHitCase2 ;
Double_t bx_FitHitCase2;
Double_t cx_FitHitCase2;
Double_t ax_fullFitCase2;
Double_t bx_fullFitCase2;
Double_t cx_fullFitCase2;
Double_t ay_fullFitCase2;
Double_t by_fullFitCase2;
Double_t Chi2_fullFitCase2;
Double_t MaxChi2_fullFitCase2;
int NHits_fullFitCase2;
Double_t y0_fullFitCase2;
Double_t tx_inf_Case2;
Double_t x0_Case2;
Double_t x_inf_Case2;
Double_t Del_x_inf_Case2;
Double_t tx_pickedcombination_Case2;
Double_t DelxProjectedSeed1_Case2;
Double_t DelxProjectedSeed2_Case2;
Double_t MaxChi2Hit_Case2;
Double_t ax_Case2_Seed1;
Double_t bx_Case2_Seed1;
Double_t cx_Case2_Seed1;
Double_t ax_Case2_Seed2;
Double_t bx_Case2_Seed2;
Double_t cx_Case2_Seed2;



void TrackStudy::Begin(TTree * /*tree*/)
{
   //Y segment fit with XZ fit
  //case 0
  t1->Branch("Track_Zone",&Track_Zone,"Track_Zone/I");
  t1->Branch("FitLineY_ay_case0",& FitLineY_ay_case0,"FitLineY_ay_case0/D");
   t1->Branch("FitLineY_by_case0",& FitLineY_by_case0,"FitLineY_by_case0/D");
   t1->Branch("FitLineY_Chi2_case0",& FitLineY_Chi2_case0,"FitLineY_Chi2_case0/D");
   t1->Branch("FitLineY_Chi2DoF_case0",& FitLineY_Chi2DoF_case0,"FitLineY_Chi2DoF_case0/D");
   t1->Branch("BestUV_NHits_case0",& NUV_Case0, "BestUV_NHits_case0/I");

   t1->Branch("FitLineY_ay_case1",& FitLineY_ay_case1,"FitLineY_ay_case1/D");
   t1->Branch("FitLineY_by_case1",& FitLineY_by_case1,"FitLineY_by_case1/D");
   t1->Branch("FitLineY_Chi2_case1",& FitLineY_Chi2_case1,"FitLineY_Chi2_case1/D");
   t1->Branch("FitLineY_Chi2DoF_case1",& FitLineY_Chi2DoF_case1,"FitLineY_Chi2DoF_case1/D");
   t1->Branch("BestUV_NHits_case1",& NUV_Case1, "BestUV_NHits_case1/I");

   t1->Branch("FitLineY_ay_case2",& FitLineY_ay_case2,"FitLineY_ay_case2/D");
   t1->Branch("FitLineY_by_case2",& FitLineY_by_case2,"FitLineY_by_case2/D");
   t1->Branch("FitLineY_Chi2_case2",& FitLineY_Chi2_case2,"FitLineY_Chi2_case2/D");
   t1->Branch("FitLineY_Chi2DoF_case2",& FitLineY_Chi2DoF_case2,"FitLineY_Chi2DoF_case2/D");
   t1->Branch("BestUV_NHits_case2",& NUV_Case2, "BestUV_NHits_case2/I");

   //Full fit case 0
   t1->Branch("ax_fullFitCase0",&ax_fullFitCase0,"ax_fullFitCase0/D");
   t1->Branch("bx_fullFitCase0",&bx_fullFitCase0,"bx_fullFitCase0/D");
   t1->Branch("cx_fullFitCase0",&cx_fullFitCase0,"cx_fullFitCase0/D");
   t1->Branch("ay_fullFitCase0",&ay_fullFitCase0,"ay_fullFitCase0/D");
   t1->Branch("by_fullFitCase0",&by_fullFitCase0,"by_fullFitCase0/D");
   t1->Branch("Chi2_fullFitCase0",&Chi2_fullFitCase0,"Chi2_fullFitCase0/D");
   t1->Branch("MaxChi2_fullFitCase0",&MaxChi2_fullFitCase0,"MaxChi2_fullFitCase0/D");
   t1->Branch("NHits_fullFitCase0",&NHits_fullFitCase0,"NHits_fullFitCase0/I");
   t1->Branch("y0_fullFitCase0",&y0_fullFitCase0,"y0_fullFitCase0/D");


   //Full fit case 1

   t1->Branch("ax_fullFitCase1",&ax_fullFitCase1,"ax_fullFitCase1/D");
   t1->Branch("bx_fullFitCase1",&bx_fullFitCase1,"bx_fullFitCase1/D");
   t1->Branch("cx_fullFitCase1",&cx_fullFitCase1,"cx_fullFitCase1/D");
   t1->Branch("ay_fullFitCase1",&ay_fullFitCase1,"ay_fullFitCase1/D");
   t1->Branch("by_fullFitCase1",&by_fullFitCase1,"by_fullFitCase1/D");
   t1->Branch("Chi2_fullFitCase1",&Chi2_fullFitCase1,"Chi2_fullFitCase1/D");
   t1->Branch("MaxChi2_fullFitCase1",&MaxChi2_fullFitCase1,"MaxChi2_fullFitCase1/D");
   t1->Branch("NHits_fullFitCase1",&NHits_fullFitCase1,"NHits_fullFitCase1/I");
   t1->Branch("y0_fullFitCase1",&y0_fullFitCase1,"y0_fullFitCase1/D");

   //Full fit case 2
   t1->Branch("ax_fullFitCase2",&ax_fullFitCase2,"ax_fullFitCase2/D");
   t1->Branch("bx_fullFitCase2",&bx_fullFitCase2,"bx_fullFitCase2/D");
   t1->Branch("cx_fullFitCase2",&cx_fullFitCase2,"cx_fullFitCase2/D");
   t1->Branch("ay_fullFitCase2",&ay_fullFitCase2,"ay_fullFitCase2/D");
   t1->Branch("by_fullFitCase2",&by_fullFitCase2,"by_fullFitCase2/D");
   t1->Branch("Chi2_fullFitCase2",&Chi2_fullFitCase2,"Chi2_fullFitCase2/D");
   t1->Branch("MaxChi2_fullFitCase2",&MaxChi2_fullFitCase2,"MaxChi2_fullFitCase2/D");
   t1->Branch("NHits_fullFitCase2",&NHits_fullFitCase2,"NHits_fullFitCase2/I");
   t1->Branch("y0_fullFitCase2",&y0_fullFitCase2,"y0_fullFitCase2/D");

   // Fit XZ Case 0 ,1 ,2
   t1->Branch("ax_FitHitCase0_Cub",&ax_FitHitCase0_Cub,"ax_FitHitCase0_Cub/D");
   t1->Branch("bx_FitHitCase0_Cub",&bx_FitHitCase0_Cub,"bx_FitHitCase0_Cub/D");
   t1->Branch("cx_FitHitCase0_Cub",&cx_FitHitCase0_Cub,"cx_FitHitCase0_Cub/D");
   t1->Branch("dx_FitHitCase0_Cub",&dx_FitHitCase0_Cub,"dx_FitHitCase0_Cub/D");

   t1->Branch("ax_FitHitCase0",&ax_FitHitCase0,"ax_FitHitCase0/D");
   t1->Branch("bx_FitHitCase0",&bx_FitHitCase0,"bx_FitHitCase0/D");
   t1->Branch("cx_FitHitCase0",&cx_FitHitCase0,"cx_FitHitCase0/D");
   t1->Branch("ax_FitHitCase1",&ax_FitHitCase1,"ax_FitHitCase1/D");
   t1->Branch("bx_FitHitCase1",&bx_FitHitCase1,"bx_FitHitCase1/D");
   t1->Branch("cx_FitHitCase1",&cx_FitHitCase1,"cx_FitHitCase1/D");
   t1->Branch("ax_FitHitCase2",&ax_FitHitCase2,"ax_FitHitCase2/D");
   t1->Branch("bx_FitHitCase2",&bx_FitHitCase2,"bx_FitHitCase2/D");
   t1->Branch("cx_FitHitCase2",&cx_FitHitCase2,"cx_FitHitCase2/D");

//Solve Parabola Seed 1 /2 Case 0,1,2
   t1->Branch("MaxChi2Hit_Case0",&MaxChi2Hit_Case0,"MaxChi2Hit_Case0/D");
   t1->Branch("MaxChi2Hit_Case1",&MaxChi2Hit_Case1,"MaxChi2Hit_Case1/D");
   t1->Branch("MaxChi2Hit_Case2",&MaxChi2Hit_Case2,"MaxChi2Hit_Case2/D");
   t1->Branch("ax_Case0_Seed1",&ax_Case0_Seed1,"ax_Case0_Seed1/D");
   t1->Branch("bx_Case0_Seed1",&bx_Case0_Seed1,"bx_Case0_Seed1/D");
   t1->Branch("cx_Case0_Seed1",&cx_Case0_Seed1,"cx_Case_Seed1/D");
   t1->Branch("ax_Case0_Seed2",&ax_Case0_Seed2,"ax_Case0_Seed2/D");
   t1->Branch("bx_Case0_Seed2",&bx_Case0_Seed2,"bx_Case0_Seed2/D");
   t1->Branch("cx_Case0_Seed2",&cx_Case0_Seed2,"cx_Case0_Seed2/D");

   t1->Branch("WorstPlane",&WorstPlane,"WorstPlane/I");
   t1->Branch("NbHitsInWorst",&nbHitsInWorst,"NbHitsInWorst/I");
   t1->Branch("ax_Case1_Seed1",&ax_Case1_Seed1,"ax_Case1_Seed1/D");
   t1->Branch("bx_Case1_Seed1",&bx_Case1_Seed1,"bx_Case1_Seed1/D");
   t1->Branch("cx_Case1_Seed1",&cx_Case1_Seed1,"cx_Case1_Seed1/D");
   t1->Branch("ax_Case1_Seed2",&ax_Case1_Seed2,"ax_Case1_Seed2/D");
   t1->Branch("bx_Case1_Seed2",&bx_Case1_Seed2,"bx_Case1_Seed2/D");
   t1->Branch("cx_Case1_Seed2",&cx_Case1_Seed2,"cx_Case1_Seed2/D");
   t1->Branch("ax_Case2_Seed1",&ax_Case2_Seed1,"ax_Case2_Seed1/D");
   t1->Branch("bx_Case2_Seed1",&bx_Case2_Seed1,"bx_Case2_Seed1/D");
   t1->Branch("cx_Case2_Seed1",&cx_Case2_Seed1,"cx_Case2_Seed1/D");

   t1->Branch("ax_Case2_Seed2",&ax_Case2_Seed2,"ax_Case2_Seed2/D");
   t1->Branch("bx_Case2_Seed2",&bx_Case2_Seed2,"bx_Case2_Seed2/D");
   t1->Branch("cx_Case2_Seed2",&cx_Case2_Seed2,"cx_Case2_Seed2/D");

   t1->Branch("Chi2X_Case0",&Chi2X_Case0,"Chi2X_Case0/D");
   //t1->Branch("Chi2X_Case0Cub",&Chi2X_Case0Cub,"Chi2X_Case0Cub/D");

   
   t1->Branch("Chi2X_Case1",&Chi2X_Case1,"Chi2X_Case1/D");
   t1->Branch("Chi2X_Case2",&Chi2X_Case2,"Chi2X_Case2/D");
   t1->Branch("tx_inf_Case0",&tx_inf_Case0,"tx_inf_Case0Case0/D");
   t1->Branch("tx_inf_Case1",&tx_inf_Case1,"tx_inf_Case0Case1/D");
   t1->Branch("tx_inf_Case2",&tx_inf_Case2,"tx_inf_Case0Case2/D");
   t1->Branch("x0_Case0",&x0_Case0,"x0_Case0/D");
   t1->Branch("x0_Case1",&x0_Case1,"x0_Case1/D");
   t1->Branch("x0_Case2",&x0_Case2,"x0_Case2/D");
   t1->Branch("x_inf_Case0",&x_inf_Case0,"x_inf_Case0/D");
   t1->Branch("x_inf_Case1",&x_inf_Case1,"x_inf_Case1/D");
   t1->Branch("x_inf_Case2",&x_inf_Case2,"x_inf_Case2/D");
   t1->Branch("Del_x_inf_Case0",&Del_x_inf_Case0,"Del_x_inf_Case0/D");
   t1->Branch("Del_x_inf_Case1",&Del_x_inf_Case1,"Del_x_inf_Case1/D");
   t1->Branch("Del_x_inf_Case2",&Del_x_inf_Case2,"Del_x_inf_Case2/D");
   t1->Branch("tx_pickedcombination_Case0",&tx_pickedcombination_Case0,"tx_pickedcombination_Case0/D");
   t1->Branch("tx_pickedcombination_Case1",&tx_pickedcombination_Case1,"tx_pickedcombination_Case1/D");
   t1->Branch("tx_pickedcombination_Case2",&tx_pickedcombination_Case2,"tx_pickedcombination_Case2/D");

   t1->Branch("DelxProjectedSeed1_Case0",&DelxProjectedSeed1_Case0,"DelxProjectedSeed1_Case0/D");
   t1->Branch("DelxProjectedSeed1_Case1",&DelxProjectedSeed1_Case1,"DelxProjectedSeed1_Case1/D");
   t1->Branch("DelxProjectedSeed1_Case2",&DelxProjectedSeed1_Case2,"DelxProjectedSeed1_Case2/D");
   t1->Branch("DelxProjectedSeed2_Case0",&DelxProjectedSeed2_Case0,"DelxProjectedSeed2_Case0/D");
   t1->Branch("DelxProjectedSeed2_Case1",&DelxProjectedSeed2_Case1,"DelxProjectedSeed2_Case1/D");
   t1->Branch("DelxProjectedSeed2_Case2",&DelxProjectedSeed2_Case2,"DelxProjectedSeed2_Case2/D");

   t1->Branch("NX_Case0",&NX_Case0,"NX_Case0/I");
   t1->Branch("NX_Case1",&NX_Case1,"NX_Case1/I");
   t1->Branch("NX_Case2",&NX_Case2,"NX_Case2/I");


   t1->Branch("Case0Accepted",&Case0Accepted,"Case0Accepted/O");
   t1->Branch("Case1Accepted",&Case1Accepted,"Case1Accepted/O");
   t1->Branch("Case2Accepted",&Case2Accepted,"Case2Accepted/O");
   t1->Branch("PhysicsInterest",&PhysicsInterest,"PhysicsInterest/O");
   t1->Branch("isDown",&ISDOWN,"isDown/O");
   t1->Branch("X0Back",&X0Back,"X0Back/D");
   t1->Branch("ISLONG",&ISLONG,"ISLONG/O");
   t1->Branch("OVTX_Z",&OVTX_Z,"OVTX_Z/D");
   t1->Branch("minCoord",&minCoord,"minCoord/D");
   t1->Branch("maxCoord",&maxCoord,"maxCoord/D");
   t1->Branch("avgCoord",&avgCoord,"avgCoord/D");
   t1->Branch("minY",&minY,"minY/D");
   t1->Branch("maxY",&maxY,"maxY/D");
   t1->Branch("NUV",&NUV,"NUV/I");
   //t1->Branch("nCombT1",nCombT1,"nCombT1/D");
   //t1->Branch("DeltaYT1",)
   t1->Branch("DeltaXT1BestCorrSlopeTy",&DeltaXBestCombCorrSlopeY,"DeltaXT1BestCorrSlopeTy/D");

   t1->Branch("YUT2Projected",&YUT2Projected,"YUT2Projected/D");
   t1->Branch("YVT2Projected",&YVT2Projected,"YVT2Projected/D");
   t1->Branch("DeltaUT2Projected",&DeltaUT2Projected,"DeltaUT2Projected/D");
   t1->Branch("DeltaVT2Projected",&DeltaVT2Projected,"DeltaVT2Projected/D");
   t1->Branch("XPairT1",&XPairT1,"XPairT1/D");
   t1->Branch("YPairT1",&YPairT1,"YPairT1/D");
   t1->Branch("TYPairT1",&TyPairT1,"YPairT1/D");


   t1->Branch("DeltaXBestCombinationT1",&DeltaXBestCombT1,"DeltaXBestCombT1/D");
   t1->Branch("Chi2_XZDRATIOFixed",&Chi2_XZDRATIOFixed,"Chi2_XZDRATIOFixed/D");
   t1->Branch("Chi2_XZDRATIO",&Chi2_XZDRATIO,"Chi2_XZDRATIO/D");
   t1->Branch("Chi2_XStepFit",&Chi2_XStepFit,"Chi2_XStepFit/D");
   t1->Branch("Chi2_perDoF_XStepFit",&Chi2_perDoF_XStepFit,"Chi2_perDoF_XStepFit/d");
   t1->Branch("ax_XStepHitFit",&ax_XStepHitFit,"ax_XStepHitFit/D");
   t1->Branch("bx_XStepHitFit",&bx_XStepHitFit,"bx_XStepHitFit/D");
   t1->Branch("cx_XStepHitFit",&cx_XStepHitFit,"cx_XStepHitFit/D");
   t1->Branch("NHits_XStepHitFit",&NClustersX,"NHits_XStepHitFit/I");
   t1->Branch("NSplittedX",&NSplittedX,"NSplittedX/I");
   t1->Branch("NSplittedUV",&NSplittedUV,"NSplittedUV/I");

   //TH1D * Hist_MCHitContribDuplicate = new TH1D("
   t1->Branch("Chi2_ParabolaXZ",&Chi2_ParabolaXZ,"Chi2_ParabolaXZ/D");
   t1->Branch("ax_parXZ",&ax_par,"ax_parXZ/D");
   t1->Branch("bx_parXZ",&bx_par,"bx_parXZ/D");
   t1->Branch("cx_parXZ",&cx_par,"cx_parXZ/D");
   t1->Branch("Chi2_CubicXZ",&Chi2_CubicXZ,"Chi2_CubicXZ/D");
   t1->Branch("ax_cubXZ",&ax_cub,"ax_cubXZ/D");
   t1->Branch("bx_cubXZ",&bx_cub,"bx_cubXZ/D");
   t1->Branch("cx_cubXZ",&cx_cub,"cx_cubXZ/D");
   t1->Branch("dx_cubXZ",&dx_cub,"dx_cubXZ/D");
   t1->Branch("Z_T3",&Z_T3,"Z_T3/D");
   t1->Branch("Z_T1",&Z_T1,"Z_T1/D");    
   t1->Branch("Y_T3",&Y_T3,"Y_T3/D");
   t1->Branch("Y_T1",&Y_T1,"Y_T1/D");
   t1->Branch("X_T3",&X_T3,"X_T3/D");
   t1->Branch("X_T1",&X_T1,"X_T1/D");
   t1->Branch("Coord_T3",&Coord_T3,"Coord_T3/D");
   t1->Branch("Coord_T1",&Coord_T1,"Coord_T1/D");
   
   t1->Branch("Chi2_LineY",&Chi2_LineY,"Chi2_LineY/D");
   t1->Branch("ay_line",&ay_line,"ay_line/D");
   t1->Branch("by_line",&by_line,"by_line/D");
   t1->Branch("Chi2_ParabolaY",&Chi2_ParabolaY,"Chi2_ParabolaY/D");
   t1->Branch("ay_par",&ay_par,"ay_par/D");
   t1->Branch("by_par",&by_par,"by_par/D");
   t1->Branch("cy_par",&cy_par,"cy_par/D");
   t1->Branch("Track_P",&Momentum,"Track_P/D");
   t1->Branch("Track_Px",&Px,"Track_Px/D");
   t1->Branch("Track_Py",&Py,"Track_Py/D");
   t1->Branch("Track_eta",&Track_eta,"Track_Px/D");
   t1->Branch("Track_Pt",&Track_Pt,"Track_Pt/D");
   t1->Branch("Track_Phi",&Track_Phi,"Track_Phi/D");
   t1->Branch("nPV",&nPV,"nPV/I");
   t1->Branch("isElectron",&IsElectron,"isElectron/O");
   t1->Branch("MCParticleID",&MCParticleiD,"MCParticleID/I");


   //Remaining XATZ
   t1->Branch("nSeed1Case0Rem",&nRemSeed1_Case0,"nSeed1Case0Rem/I");
   t1->Branch("nSeed2Case0Rem",&nRemSeed2_Case0,"nSeed2Case0Rem/I");
   t1->Branch("DelSeed1_Case0Rem", &DeltaSeed1_Case0,"DelSeed1_Case0Rem[3]/D");
   t1->Branch("DelSeed2_Case0Rem", &DeltaSeed2_Case0,"DelSeed2_Case0Rem[3]/D");
   
   t1->Branch("nSeed1Case1Rem",&nRemSeed1_Case1,"nSeed1Case1Rem/I");
   t1->Branch("nSeed2Case1Rem",&nRemSeed2_Case1,"nSeed2Case1Rem/I");
   t1->Branch("DelSeed1_Case1Rem", &DeltaSeed1_Case1,"DelSeed1_Case1Rem[3]/D");
   t1->Branch("DelSeed2_Case1Rem", &DeltaSeed2_Case1,"DelSeed2_Case1Rem[3]/D");
   
   t1->Branch("nSeed1Case2Rem",&nRemSeed1_Case2,"nSeed1Case2Rem/I");
   t1->Branch("nSeed2Case2Rem",&nRemSeed2_Case2,"nSeed2Case2Rem/I");
   t1->Branch("DelSeed1_Case2Rem", &DeltaSeed1_Case2,"DelSeed1_Case2Rem[3]/D");
   t1->Branch("DelSeed2_Case2Rem", &DeltaSeed2_Case2,"DelSeed2_Case2Rem[3]/D");
   TString option = GetOption();
   m_dRatio0 = -0.000262;
   //Counters For Efficiencies on Long tracks
   m_long_NGeantHit=0;
   m_long_NPrHit=0;
   m_long_NMCHitIntoCluster=0;
   m_long_NLayer_Geant=0;
   m_long_NLayer_PrHit=0;
   m_long_NLayer_MCHitInCluster=0;



   //XZ Study
   m_physical=0;
   m_physicalCase0Intrinsic=0;
   m_physicalCase0IntrinsicMore4=0;
   m_physicalCase0IntrinsicMore5=0;

   m_physicalCase0IntrinsicWithduplicates=0;
   m_physicalCase0IntrinsicWithduplicatesMore4=0;
   m_physicalCase0IntrinsicWithduplicatesMore5=0;

   m_physicalCase1Intrinsic=0;
   m_physicalCase1IntrinsicMore4=0;
   m_physicalCase1IntrinsicMore5=0;

   m_physicalCase1IntrinsicWithduplicates=0;
   m_physicalCase1IntrinsicWithduplicatesMore4=0;
   m_physicalCase1IntrinsicWithduplicatesMore5=0;

   m_physicalCase2Intrinsic=0;
   m_physicalCase2IntrinsicMore4=0;
   m_physicalCase2IntrinsicMore5=0;

   m_physicalCase2IntrinsicWithduplicates=0;
   m_physicalCase2IntrinsicWithduplicatesMore4=0;
   m_physicalCase2IntrinsicWithduplicatesMore5=0;
   m_nbTrackCloneX=0;
   m_nbTrackCloneUV=0;
   m_nbTrackCloneUVandX=0;
   //zReference Position to be used in the Fits
   m_zReference = 8520.;
}

void TrackStudy::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t TrackStudy::Process(Long64_t entry, Bool_t debug)
{

   // std::cout<"Processing Entry = "<<entry<<std::endl;
   fChain->GetTree()->GetEntry(entry);
   // std::cout<"Processing Entry = "<<entry<<std::endl;
   //Bool_t debug = true;
   if(!(PrHit>1)) return kTRUE;
   if(debug) std::cout<<"Ciao"<<std::endl;
   // HERE DEBUG AND PRINT AT VIDEO OF HIT CONTENT
   Double_t PCut =0;
   ISLONG= isLong;
   ISDOWN = isDown;
   IsElectron = isElectron;
   MCParticleiD = MCParticleID;
   OVTX_Z = MC_Ovtx_z;

   PhysicsInterest= (fromDecay || fromPrimaryVertex) && (isLong || isDown ) && P>PCut;
   if(PhysicsInterest) m_physical++;
   if(!PhysicsInterest || std::fabs(MCParticleID==11)) return kTRUE;
   //Count number of Clusters / MCHits from Geant // MCHits ending into Cluster
   CountHits(debug);

   //Study Fit Chi2 values etc on the MCHits
   FitHits(debug);

   Int_t nX=0;
   Int_t zone = 0;
   //NClonesX=0;
   if(nXUp>=4)  zone = 1;
   if(nxDown>=4) zone = 0;
   //Generate the tracks : Only X Projection + UV Segment + Full Track
   PrSeedTrack xProj = PrSeedTrack(zone, m_zReference);
   PrSeedTrack UVSegment = PrSeedTrack(zone, m_zReference);
   PrSeedTrack FullTrack = PrSeedTrack(zone, m_zReference);
   for(int i =0;i<PrHit;i++){
      PatHit hit = PatHit();
      hit.setHit(PrHit_Xat0[i],PrHit_Zat0[i],PrHit_dxDy[i],PrHit_dzDy[i],std::sqrt(PrHit_w2[i]),PrHit_yMin[i],PrHit_yMax[i],PrHit_zone[i],PrHit_planeCode[i],PrHit_isX[i],PrHit_LHCbID[i]);
      FTCluster cluster;
      cluster.setCluster(ChID_Fraction[i],PrHit_Size[i],ChID_Charge[i],ChID_SipmID[i],ChID_SipmCell[i],ChID_Module[i],ChID_Layer[i],ChID_Mat[i],ChID_Quarter[i]);
      hit.SetCluster(cluster);
      FullTrack.addHit(hit);
      if(hit.isX()) xProj.addHit(hit);
      if(!hit.isX()) UVSegment.addHit(hit);
   }
   fitXProjection(xProj,debug);
   ///xProj.PrintHits();
   //xProj.PrintTrack();
   X0Back =(Double_t) (xProj.ax() - m_zReference*xProj.bx()+ 2.458e8*xProj.cx());
   //std::cout<<"X0Back = " <<X0Back<<std::endl;

   xProj.setChi2X(xProj.Chi2());
   ax_XStepHitFit=xProj.ax();
   bx_XStepHitFit=xProj.bx();
   cx_XStepHitFit=xProj.cx();
   Chi2_XStepFit=xProj.Chi2();
   NClustersX = xProj.hits().size();
   Chi2_perDoF_XStepFit = xProj.Chi2()/((double)xProj.hits().size()-3.);
   CountClones(FullTrack,debug);

   StereoSearch(xProj,UVSegment,debug);
   //Case 0 first - 0 last =11;

   Case0Accepted =  XZStudyCase0(xProj,UVSegment,debug);
   // Case 1 first -2 last =11;
   Case1Accepted = XZStudyCase1(xProj,UVSegment,debug);
   //Case 2 first -1 last =8;
   Case2Accepted = XZStudyCase2(xProj,UVSegment,debug);

   t1->Fill();
   if(PrHit<9) return kTRUE;
   return kTRUE;
}
void TrackStudy::solveParabola2(const PatHit hit1,const PatHit hit2,const PatHit hit3,double& a1, double& b1,double& c1){
   const double z1 = hit1.z() - m_zReference;
   const double z2 = hit2.z() - m_zReference;
   const double z3 = hit3.z() - m_zReference;
   const double x1 = hit1.x();
   const double x2 = hit2.x();
   const double x3 = hit3.x();
   //const double e = m_dRatio;
   //double m_dRatio0 = -0.000246;
   const double corrZ1 = 1+m_dRatio0*z1;
   const double corrZ2 = 1+m_dRatio0*z2;
   const double corrZ3 = 1+m_dRatio0*z3;
   const double det = (z1*z1)*corrZ1*z2 + z1*(z3*z3)*corrZ3 + (z2*z2)*corrZ2*z3 - z2*(z3*z3)*corrZ3 - z1*(z2*z2)*corrZ2 - z3*(z1*z1)*corrZ1;
   if( std::fabs(det) < 1e-8 )
   {
      a1 = 0.0;
      b1 = 0.0;
      c1 = 0.0;
      return;
   }
   const double det1 = (x1)*z2 + z1*(x3) + (x2)*z3 - z2*(x3) - z1*(x2) - z3*(x1);
   const double det2 = (z1*z1)*corrZ1*x2 + x1*(z3*z3)*corrZ3 + (z2*z2)*corrZ2*x3 - x2*(z3*z3)*corrZ3 - x1*(z2*z2)*corrZ2 - x3*(z1*z1)*corrZ1;
   const double det3 = (z1*z1)*corrZ1*z2*x3 + z1*(z3*z3)*corrZ3*x2 + (z2*z2)*corrZ2*z3*x1 - z2*(z3*z3)*corrZ3*x1 - z1*(z2*z2)*corrZ2*x3 - z3*(z1*z1)*corrZ1*x2;
   a1 = det1/det;
   b1 = det2/det;
   c1 = det3/det;
}


Bool_t TrackStudy::fitXProjectionCubic(PrSeedTrack & track, Bool_t debug)
{
float mat[10];
float rhs[4];
std::vector<PatHit> Hits = track.hits();
track.setParameters2(0,0,0,0,0,0);
for(int loop =0; 3>loop; ++loop){
      std::fill(mat,mat+10,0.);
      std::fill(rhs,rhs+4,0.);
      for( int i=0;i<Hits.size();i++){
         const float w = Hits[i].w2();
         if(debug) std::cout<<"Hit w2"<<Hits[i].w2()<<std::endl;
         const float dz = Hits[i].z(0.)-m_zReference;
         //float deta = dz*dz*(1.+m_dRatio0*dz);
         //if(m_usedRatio){
         //deta = dz*dz*(1-dRatio*dz);
         //}
         float dist = track.distance2(Hits[i]);
         if(debug) std::cout<<"i"<<setw(20)<<i<<setw(20)<<"Hit X"<<setw(20)<<Hits[i].x(0)<<setw(20)<<"Track x"<<setw(20)<<track.x(Hits[i].z(0.))<<setw(20)<<"Distance"<< track.distance2(Hits[i])<<std::endl;
         mat[0]+= w     ;
         mat[1]+= w * dz;
         mat[2]+= w * dz * dz;
         mat[3]+= w * dz * dz;
         mat[4]+= w * dz * dz * dz;
         mat[5]+= w * dz * dz * dz * dz ;
         mat[6]+= w * dz * dz * dz;
         mat[7]+= w * dz * dz * dz * dz ;
         mat[8]+= w * dz * dz * dz * dz * dz;
         mat[9]+= w * dz * dz * dz * dz * dz * dz;
         //right hand side
         rhs[0]+= w * dist;
         rhs[1]+= w * dist * dz;
         rhs[2]+= w * dist * dz*dz;
         rhs[3]+= w * dist * dz*dz*dz;
      }

      ROOT::Math::CholeskyDecomp<float,4> decomp(mat);
      if(!decomp){
      std::cout<<"Failed to decompose matrix"<<std::endl;
      return false;
      }
      decomp.Solve(rhs);
      track.updateParameters2(rhs[0],rhs[1],rhs[2],0.,0.,rhs[3]);
   }
   Float_t chi2 = track.Chi2();
   //std::cout<<"Chi2 Cubic = "<<chi2<<std::endl;
   //std::cout<<"ax"<<setw(20)<<"bx"<<setw(20)<<"cx"<<setw(20)<<"dx"<<setw(20)<<std::endl;
   //std::cout<<track.ax_cub()<<setw(20)<<track.bx_cub()<<setw(20)<<track.cx_cub()<<setw(12)<<track.dx_cub()<<std::endl;
   //std::cout<<"DRATIO =  "<<track.dx_cub()/track.cx_cub();
   return true;

}




Bool_t TrackStudy::fitXProjection(PrSeedTrack & track, Bool_t debug)
{
   float mat[6];
   float rhs[3];
   std::vector<PatHit> Hits = track.hits();

   if(debug) std::cout<<"Hits Size"<<Hits.size()<<std::endl;
   if(debug) std::cout<<"Hits Size from track"<<track.hits().size();
   float dRatio=m_dRatio0;
   track.setdRatio(dRatio);
   for(int loop = 0;3>loop;++loop){
      std::fill(mat,mat+6,0.);
      std::fill(rhs,rhs+3,0.);
      if(debug) std::cout<<"******************* Fit Loop "<<loop<<"*********************"<<std::endl;
      for( int i=0;i<Hits.size();i++){
         const float w = Hits[i].w2();
         if(debug) std::cout<<"Hit w2"<<Hits[i].w2()<<std::endl;
         const float dz = Hits[i].z(0.)-m_zReference;
         float deta = dz*dz*(1.+m_dRatio0*dz);
         //if(m_usedRatio){
         //deta = dz*dz*(1-dRatio*dz);
         //}
         float dist = track.distance(Hits[i]);
         if(debug) std::cout<<"i"<<setw(20)<<i<<setw(20)<<"Hit X"<<setw(20)<<Hits[i].x(0)<<setw(20)<<"Track x"<<setw(20)<<track.x(Hits[i].z(0.))<<setw(20)<<"Distance"<< track.distance(Hits[i])<<std::endl;
         mat[0]+= w     ;
         mat[1]+= w * dz;
         mat[2]+= w * dz * dz;
         mat[3]+= w * deta;
         mat[4]+= w * dz * deta;
         mat[5]+= w * deta * deta;
         //right hand side
         rhs[0]+= w * dist;
         rhs[1]+= w * dist * dz;
         rhs[2]+= w * dist * deta;
      }

      ROOT::Math::CholeskyDecomp<float,3> decomp(mat);
      if(!decomp){
         std::cout<<"Failed to decompose matrix"<<std::endl;
         return false;
      }
      decomp.Solve(rhs);
      track.updateParameters(rhs[0],rhs[1],rhs[2],0.,0.,0.);
   }
   Float_t chi2 = track.Chi2();
   track.setChi2X(chi2);
   track.setNDOFX(track.hits().size());
   return true;
}


void TrackStudy::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void TrackStudy::Terminate()
{
   file->Write();
   std::cout<<"Long Tracks"<<std::endl;
   std::cout<<"Number of Geant MCHit                             "<<setw(20)<<m_long_NGeantHit<<std::endl;
   std::cout<<"Number of different Layers with Geant MCHit       "<<setw(20)<<m_long_NLayer_Geant<<std::endl;

   std::cout<<"Number of MCHit into Cluster                      "<<setw(20)<<m_long_NMCHitIntoCluster<<std::endl;
   std::cout<<"Number of different Layers with MCHit into Cluster"<<setw(20)<<m_long_NLayer_MCHitInCluster<<std::endl;

   std::cout<<"Number of Clusters \t                             " <<m_long_NPrHit<<std::endl;
   std::cout<<"Number of diffrent Layers Clusters \t             " <<m_long_NLayer_PrHit<<std::endl;
   zMCHitClones->Draw();
   zMCHitClonevsP->Draw("colz");


   std::cout<<"Physical Tracks =                                 "<<setw(20)<<m_physical<<std::endl;
   std::cout<<"**********************"<<setw(20)<<"Case0"<<setw(20)<<"**********************"<<std::endl;
   std::cout<<"Physical Tracks (Case0 intrinsic)=                "<<setw(20)<<m_physicalCase0Intrinsic<<std::endl;
   std::cout<<"Physical Tracks (Case0 intrinsic) duplicates      "<<setw(20)<<m_physicalCase0IntrinsicWithduplicates<<std::endl;
   std::cout<<"Physical Tracks (Case0 intrinsic) >=4             "<<setw(20)<<m_physicalCase0IntrinsicMore4<<std::endl;
   std::cout<<"Physical Tracks (Case0 intrinsic) >=4 duplicates  "<<setw(20)<<m_physicalCase0IntrinsicWithduplicatesMore4<<std::endl;
   std::cout<<"Physical Tracks (Case0 intrinsic) >=5             "<<setw(20)<<m_physicalCase0IntrinsicMore5<<std::endl;
   std::cout<<"Physical Tracks (Case0 intrinsic) duplicates      "<<setw(20)<<m_physicalCase0IntrinsicWithduplicatesMore5<<std::endl;
   std::cout<<"**********************"<<setw(20)<<"Case1"<<setw(20)<<"**********************"<<std::endl;

   std::cout<<"Physical Tracks (Case1 intrinsic)=                "<<setw(20)<<m_physicalCase1Intrinsic<<std::endl;
   std::cout<<"Physical Tracks (Case1 intrinsic) duplicates      "<<setw(20)<<m_physicalCase1IntrinsicWithduplicates<<std::endl;
   std::cout<<"Physical Tracks (Case1 intrinsic) >=4             "<<setw(20)<<m_physicalCase1IntrinsicMore4<<std::endl;
   std::cout<<"Physical Tracks (Case1 intrinsic) >=4 duplicates  "<<setw(20)<<m_physicalCase1IntrinsicWithduplicatesMore4<<std::endl;
   std::cout<<"Physical Tracks (Case1 intrinsic) >=5             "<<setw(20)<<m_physicalCase1IntrinsicMore5<<std::endl;
   std::cout<<"Physical Tracks (Case1 intrinsic) duplicates      "<<setw(20)<<m_physicalCase1IntrinsicWithduplicatesMore5<<std::endl;
   std::cout<<"**********************"<<setw(20)<<"Case2"<<setw(20)<<"**********************"<<std::endl;

   std::cout<<"Physical Tracks (Case2 intrinsic)=                "<<setw(20)<<m_physicalCase2Intrinsic<<std::endl;
   std::cout<<"Physical Tracks (Case2 intrinsic) duplicates      "<<setw(20)<<m_physicalCase2IntrinsicWithduplicates<<std::endl;
   std::cout<<"Physical Tracks (Case2 intrinsic) >=4             "<<setw(20)<<m_physicalCase2IntrinsicMore4<<std::endl;
   std::cout<<"Physical Tracks (Case2 intrinsic) >=4 duplicates  "<<setw(20)<<m_physicalCase2IntrinsicWithduplicatesMore4<<std::endl;
   std::cout<<"Physical Tracks (Case2 intrinsic) >=5=            "<<setw(20)<<m_physicalCase2IntrinsicMore5<<std::endl;
   std::cout<<"Physical Tracks (Case2 intrinsic) duplicates      "<<setw(20)<<m_physicalCase2IntrinsicWithduplicatesMore5<<std::endl;
   std::cout<<"NbTrack with Splitted Cluster on X layer                            "<<setw(20)<<m_nbTrackCloneX<<std::endl;
   std::cout<<"NbTrack with Splitted Cluster on UV layer                          "<<setw(20)<<m_nbTrackCloneUV<<std::endl;
   std::cout<<"NbTrack with Splitted Cluster X and UV                            "<<setw(20)<<m_nbTrackCloneUVandX<<std::endl;
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}
void TrackStudy::CountHits(Bool_t debug){

   m_long_NGeantHit+=MC;
   m_long_NPrHit+=PrHit;
   m_long_NMCHitIntoCluster+=MC_ass;
   //Count different Hit Per Layer;
   //Count different MCHit Geant PerLayer;
   for(Int_t i =0;i<MC; i++){
      if(std::fabs(MC_Hit_Z[i]-MC_Hit_Z[i-1]) <40 && i>0) continue;
      m_long_NLayer_Geant+=1;
   }
   for(Int_t i =0;i<PrHit; i++){
      if(std::fabs(PrHit_Zat0[i]-PrHit_Zat0[i-1]) <40 && i>0) continue;
      m_long_NLayer_PrHit+=1;
   }
   for(Int_t i=0;i<MC_ass; i++){
      if(std::fabs(MCHit_Assoc_Z[i]-MCHit_Assoc_Z[i-1]) <40 && i>0) continue;
      m_long_NLayer_MCHitInCluster+=1;
   }
}

void TrackStudy::FitHits(Bool_t debug){

   std::vector<MCHit> track_MC(100); //Vector of MCHit going into cluster
   std::vector<MCHit> Particle_GeantHit(100);//Vector of MCHit even not going into clusters
   for (Int_t i = 0 ;  i< CheatedSeeding_NHits; i++){
      PatHit hitdeb = PatHit();
      hitdeb.setHit(PrHit_Xat0[i],PrHit_Zat0[i],PrHit_dxDy[i],PrHit_dzDy[i],std::sqrt(PrHit_w2[i]),PrHit_yMin[i],PrHit_yMax[i],PrHit_zone[i],PrHit_planeCode[i],PrHit_isX[i],PrHit_LHCbID[i]);
      FTCluster cluster;
      cluster.setCluster(ChID_Fraction[i],PrHit_Size[i],ChID_Charge[i],ChID_SipmID[i],ChID_SipmCell[i],ChID_Module[i],ChID_Layer[i],ChID_Mat[i],ChID_Quarter[i]);
      hitdeb.SetCluster(cluster);
   }
   for (Int_t i = 0;  i< MC_ass; i++) {
      MCHit mcHit =  MCHit();
      if(i>0 && std::fabs(MCHit_Assoc_Z[i]-MCHit_Assoc_Z[i-1])<40) {
         zMCHitClones->Fill(MCHit_Assoc_Z[i]);
         zMCHitClonevsP->Fill(MCHit_Assoc_Z[i],P);
         xvsyClones->Fill(MCHit_Assoc_X[i],MCHit_Assoc_Y[i]);
         if(debug) std::cout<<"MCHit i = "<<i<<"\t Z "<<setw(20)<< MCHit_Assoc_Z[i]<<"\t X"<<setw(20)<<MCHit_Assoc_X[i]<<"\t Y "<<setw(20)<<MCHit_Assoc_Y[i]<<"\tpathlenght "<<setw(20)<<MCHit_pathlength[i]<<""<<setw(20)<<std::endl;
         if(debug) std::cout<<"MCHit i -1 = "<<i<<"\t Z "<<setw(20)<< MCHit_Assoc_Z[i-1]<<"\t X"<<setw(20)<<MCHit_Assoc_X[i]<<"\t Y "<<setw(20)<<MCHit_Assoc_Y[i-1]<<"\t pathlenght   "<<setw(20)<<MCHit_pathlength[i-1]<<""<<setw(20)<<std::endl;
         continue; //Merge them or keep the best one?
      }
      mcHit.setMCHit(MCHit_Assoc_X[i], MCHit_Assoc_Y[i],MCHit_Assoc_Z[i],MCHit_tx[i],MCHit_ty[i],MCHit_p[i],MCHit_pathlength[i],P,MC_px,MC_py,MC_pz);
      track_MC.push_back(mcHit);
   }
   Track_eta =eta;
   Track_Phi = MC_Phi;
   Track_Pt = MC_Pt;
   NPv = nPV;
   Float_t ErrorX = 0.100; //Fixed error for the MCHits
   //How good is our track model for the full Fit, i.e., y(z) = a + b *z; x(z)
   //Fit a Parabola on X and straight line on y on top of the MCHits : should be preselect removing the glancing?
   double solution_xz_par[3];
   std::fill(solution_xz_par,solution_xz_par+3,0.);
   std::vector<double> dsolution_xz_par;
   //Fit a Cubic on X
   double solution_xz_cubic[4];
   std::fill(solution_xz_cubic,solution_xz_cubic+4,0.);
   std::vector<double> dsolution_xz_cubic;
   //Fit a Line on Y
   double solution_yz_line[2];
   std::fill(solution_yz_line,solution_yz_line+2,0.);
   std::vector<double> dsolution_yz_line;
   //Fit a Parabola in Y
   double solution_yz_par[3];
   std::fill(solution_yz_par,solution_yz_par+3,0.);
   std::vector<double> dsolution_yz_par;
   double dz=0.;
   double resXZ_par=0.;
   double resYZ_lin=0.;
   double resXZ_cub=0.;
   double resYZ_par=0.;
   for(int j= 0; j<9;j++)
   {
      LinParFit<double> fit_CubicXZ(4);
      LinParFit<double> fit_parabolaXZ(3);
      LinParFit<double> fit_parabolaY(3);
      LinParFit<double> fit_LineY(2);
      //double ErrorX = gRandom->Gaus(0,100);
      //double ErrorX =
      for (int i = 0; i<MC_ass; i++){
         dz=(double)((double)MCHit_Assoc_Z[i]-m_zReference);
         resXZ_par=(double)((double)MCHit_Assoc_X[i]- (solution_xz_par[0]+ solution_xz_par[1]*dz+ solution_xz_par[2]*dz*dz ));
         resXZ_cub=(double)((double)MCHit_Assoc_X[i]- (solution_xz_cubic[0]+ solution_xz_cubic[1]*dz+ solution_xz_cubic[2]*dz*dz + solution_xz_cubic[3]*dz*dz*dz));
         resYZ_lin=(double)((double)MCHit_Assoc_Y[i]- (solution_yz_line[0]+solution_yz_line[1]*dz));
         resYZ_par=(double)((double)MCHit_Assoc_Y[i]- (solution_yz_par[0]+solution_yz_par[1]*dz+solution_yz_par[2]*dz*dz));
         fit_parabolaXZ.accumulate(resXZ_par,ErrorX,dz);
         fit_LineY.accumulate(resYZ_lin,ErrorX/(TMath::Sin(5./TMath::Pi()*180)),dz);
         fit_parabolaY.accumulate(resYZ_par,1.14,dz);
         fit_CubicXZ.accumulate(resXZ_cub,0.100,dz);
      }
      bool ok = true;
      if(!fit_parabolaXZ.solve()){
         if(debug) std::cout<<"Not Able to Fit XZ with a parabola for a MCHit list of size "<<MC_ass<<"\t and Momentum "<<P<<std::endl;
         ok=false;
      }
      if(!fit_LineY.solve())
      {
         if(debug) std::cout<<"Not Able to Fit Y with a Line for a MCHit list of size"<<MC_ass<<"\t and Momentum "<<P<<std::endl;
         ok=false;
      }
      if(!fit_CubicXZ.solve()){
         if(debug) std::cout<<"Not Able to Fit Cubic XZ for a MCHit list of size"<<MC_ass<<"\t and Momentum "<<P<<std::endl;
         ok=false;
      }
      if( !fit_parabolaY.solve() ){
         ok=false;
      }
      //if(!ok) return ;
      dsolution_yz_line = fit_LineY.solution();
      dsolution_xz_par = fit_parabolaXZ.solution();
      dsolution_xz_cubic=fit_CubicXZ.solution();
      dsolution_yz_par=fit_parabolaY.solution();
      //Parabola Fit XZ
      solution_xz_par[0]+=dsolution_xz_par[0];
      solution_xz_par[1]+=dsolution_xz_par[1];
      solution_xz_par[2]+=dsolution_xz_par[2];

      solution_yz_line[0]+=dsolution_yz_line[0];
      solution_yz_line[1]+=dsolution_yz_line[1];

      solution_xz_cubic[0]+=dsolution_xz_cubic[0];
      solution_xz_cubic[1]+=dsolution_xz_cubic[1];
      solution_xz_cubic[2]+=dsolution_xz_cubic[2];
      solution_xz_cubic[3]+=dsolution_xz_cubic[3];

      solution_yz_par[0]+=dsolution_yz_par[0];
      solution_yz_par[1]+=dsolution_yz_par[1];
      solution_yz_par[2]+=dsolution_yz_par[2];
   }

   ax_cub=solution_xz_cubic[0];
   bx_cub=solution_xz_cubic[1];
   cx_cub=solution_xz_cubic[2];
   dx_cub=solution_xz_cubic[3];

   ax_par=solution_xz_par[0];
   bx_par=solution_xz_par[1];
   cx_par=solution_xz_par[2];

   ay_line=solution_yz_line[0];
   by_line=solution_yz_line[1];

   ay_par=solution_yz_par[0];
   by_par=solution_yz_par[1];
   cy_par=solution_yz_par[2];
   Chi2_ParabolaXZ=0.;
   Chi2_CubicXZ=0.;
   Chi2_LineY=0.;
   Chi2_ParabolaY=0.;
   for(int i= 0; i<MC_ass;i++){
      dz=(double)((double)MCHit_Assoc_Z[i]-m_zReference);
      Chi2_ParabolaXZ+= std::pow( ((double)MCHit_Assoc_X[i]- (solution_xz_par[0]+ solution_xz_par[1]*dz+ solution_xz_par[2]*dz*dz))/0.100,2)/(MC_ass-2);
      Chi2_CubicXZ+= std::pow( ((double)MCHit_Assoc_X[i]- (solution_xz_cubic[0]+ solution_xz_cubic[1]*dz+ solution_xz_cubic[2]*dz*dz +solution_xz_cubic[3]*dz*dz*dz))/0.100,2)/(MC_ass-4);
      Chi2_LineY+=std::pow( ((double)MCHit_Assoc_Y[i]-(solution_yz_line[0]+solution_yz_line[1]*dz))/1.14,2)/(MC_ass-2);
      Chi2_ParabolaY+=std::pow( ((double)MCHit_Assoc_Y[i]-(solution_yz_par[0]+solution_yz_par[1]*dz+solution_yz_par[2]*dz*dz))/1.14,2)/(MC_ass-3);
   }
   Double_t Radius = std::sqrt(ay_line*ay_line*(std::abs(ay_line)/1000.) + ax_par*ax_par*(std::abs(ax_par)/2000.));
   Double_t dRatio = 0.000262 + 1.66e-08*Radius+1.03e-11*Radius*Radius;

   //Will Fit the X projection , so only using the X Hits
   double mat[6];
   double rhs[3];
   double solution[3];
   solution[0]=0.;
   solution[1]=0.;
   solution[2]=0.;
   Bool_t Fit = false;
   for(int loop = 0;3>loop;++loop){
      std::fill(mat,mat+6,0.);
      std::fill(rhs,rhs+3,0.);
      for( int i = 0;i<MC_ass;i++){
         const double err = std::pow(1./0.100 ,2);
         const double dz = MCHit_Assoc_Z[i]-m_zReference;
         double deta  = dz*dz*(1- dRatio*dz);
         double dist = MCHit_Assoc_X[i]-(solution[0]+solution[1]*dz+solution[2]*dz*dz*(1-dRatio*dz));
         mat[0]+= err     ;
         mat[1]+= err * dz;
         mat[2]+= err * dz * dz;
         mat[3]+= err * deta;
         mat[4]+= err * dz * deta;
         mat[5]+= err * deta * deta;
         //right hand side
         rhs[0]+= err * dist;
         rhs[1]+= err * dist * dz;
         rhs[2]+= err * dist * deta;
      }

      ROOT::Math::CholeskyDecomp<double,3> decomp(mat);
      if(!decomp){
         //  if(debug) std::cout<<"Failed to decompose matrix"<<std::endl;
         Fit = false;
      }
      if(decomp) Fit = true;
      decomp.Solve(rhs);
      solution[0]+=rhs[0];
      solution[1]+=rhs[1];
      solution[2]+=rhs[2];
   }

   Chi2_XZDRATIO = 0;
   for( int i = 0;i<MC_ass;i++){
      float dz = MCHit_Assoc_Z[i]-m_zReference;
      float dist = MCHit_Assoc_X[i]-(solution[0]+solution[1]*dz+solution[2]*dz*dz*(1-dRatio*dz));
      Chi2_XZDRATIO+= (dist*dist)/(std::pow(0.100,2));
   }
   Chi2_XZDRATIO = Chi2_XZDRATIO/(MC_ass-4);
   Momentum  = P;
   double mat1[6];
   double rhs1[3];
   double solution1[3];
   solution1[0]=0.;
   solution1[1]=0.;
   solution1[2]=0.;
   Bool_t Fit1 = false;
   dRatio = 0.000262;
   for(int loop = 0;3>loop;++loop){
      std::fill(mat1,mat1+6,0.);
      std::fill(rhs1,rhs1+3,0.);
      for( int i = 0;i<MC_ass;i++){
         const double err1 = std::pow(1./0.100 ,2);
         const double dz1 = MCHit_Assoc_Z[i]-m_zReference;
         double deta1  = dz1*dz1*(1- dRatio*dz1);
         double dist1 = MCHit_Assoc_X[i]-(solution1[0]+solution1[1]*dz1+solution1[2]*dz1*dz1*(1-dRatio*dz1));
         mat1[0]+= err1     ;
         mat1[1]+= err1 * dz1;
         mat1[2]+= err1 * dz1 * dz1;
         mat1[3]+= err1 * deta1;
         mat1[4]+= err1 * dz1 * deta1;
         mat1[5]+= err1 * deta1 * deta1;
         //right hand side
         rhs1[0]+= err1 * dist1;
         rhs1[1]+= err1 * dist1 * dz1;
         rhs1[2]+= err1 * dist1 * deta1;
      }

      ROOT::Math::CholeskyDecomp<double,3> decomp(mat1);
      if(!decomp){
         //  if(debug) std::cout<<"Failed to decompose matrix"<<std::endl;
         Fit1 = false;
      }
      if(decomp) Fit1 = true;
      decomp.Solve(rhs1);
      solution1[0]+=rhs1[0];
      solution1[1]+=rhs1[1];
      solution1[2]+=rhs1[2];
   }
   if(Fit){
      if(debug) std::cout<<"ax"<<setw(20)<<"bx"<<setw(20)<<"cx"<<setw(20)<<std::endl;
      if(debug) std::cout<<solution[0]<<setw(20)<<solution[1]<<setw(20)<<solution[2]<<std::endl;
   }
   Chi2_XZDRATIOFixed = 0;
   for( int i = 0;i<MC_ass;i++){
      float dz1 = MCHit_Assoc_Z[i]-m_zReference;
      float dist1 = MCHit_Assoc_X[i]-(solution1[0]+solution1[1]*dz1+solution1[2]*dz1*dz1*(1-dRatio*dz1));
      Chi2_XZDRATIOFixed+= (dist1*dist1)/(std::pow(0.100,2));
   }
   Chi2_XZDRATIOFixed = Chi2_XZDRATIOFixed/(MC_ass-4);
   if(debug) std::cout<<"Chi2 XZDRATIOFIXED "<<setw(20)<<Chi2_XZDRATIOFixed<<std::endl;
   if(debug) std::cout<<"Chi2 = "<<Chi2_XZDRATIO;
}


void TrackStudy::CountClones(PrSeedTrack fullTrack,Bool_t debug){
   fullTrack.sortbyZ();
   int clX=0;
   int clUV=0;
   for(int i = 0; i<fullTrack.hits().size();i++){
      PatHit hit = fullTrack.hits()[i];
      if(hit.isX() && i>0 && hit.planeCode()==fullTrack.hits()[i-1].planeCode() &&
      hit.id()!=fullTrack.hits()[i-1].id()) {
         clX++;
      }
      if(!hit.isX() && i>0 && hit.planeCode()==fullTrack.hits()[i-1].planeCode() &&
      hit.id()!=fullTrack.hits()[i-1].id()){
         clUV++;
      }
   }
   NSplittedX = clX;
   NSplittedUV = clUV;
   if(clX>0){
      m_nbTrackCloneX++;
   }
   if(clUV>0){
      m_nbTrackCloneUV++;
   }
   if(clX>0 && clUV>0){
      m_nbTrackCloneUVandX++;
   }
}

bool TrackStudy::XZStudyCase0(PrSeedTrack& xProj, PrSeedTrack& UVSegment,Bool_t debug){

   NUV_Case0 = -10;
   ax_Case0_Seed1 = -99.;
   bx_Case0_Seed1 = -99.;
   cx_Case0_Seed1 = -99.;
   ax_Case0_Seed2 = -9999.;
   bx_Case0_Seed2 = -99.;
   cx_Case0_Seed2 = -99.;
   Chi2X_Case0 = -99.;
   Chi2X_Case0Cub = -99.;
   NX_Case0 = -1;
   tx_inf_Case0 = -99.;
   x0_Case0 = -99.;
   x_inf_Case0 = -99.;
   Del_x_inf_Case0 = -99.;
   tx_pickedcombination_Case0 = -99.;
   DelxProjectedSeed2_Case0 = -99.;
   DelxProjectedSeed1_Case0 = -99.;
   MaxChi2Hit_Case0 = -99.;
   //Full Fit
   ax_FitHitCase0 = -99.;
   bx_FitHitCase0 = -99.;
   cx_FitHitCase0 = -99.;

   ax_FitHitCase0_Cub = -99.;
   bx_FitHitCase0_Cub = -99.;
   cx_FitHitCase0_Cub = -99.;
   dx_FitHitCase0_Cub = -99;

   WorstPlane = -99;
   nbHitsInWorst = -99;
   ax_fullFitCase0 = -9999.;
   bx_fullFitCase0= -99.;
   cx_fullFitCase0 = -99.;
   ay_fullFitCase0 = -9999.;
   by_fullFitCase0 = -99.;
   Chi2_fullFitCase0 =  -99.;
   MaxChi2_fullFitCase0 = -99.;
   NHits_fullFitCase0= -1;
   y0_fullFitCase0 = -9999.;
   Y_T1 = -999999.;
   Y_T3 = +999999.;
   X_T1 = -999999.;
   X_T3 = +999999.;
   Z_T1 = -999999.;
   Z_T3 = -999999.;
   Coord_T1 = -999999.;
   Coord_T3 = +999999.;
   FitLineY_ay_case0 = -999;
   FitLineY_by_case0 = -999;
   FitLineY_Chi2_case0= -999;
   FitLineY_Chi2DoF_case0= -999;
   FitParY_Chi2_case0= -999;
   FitParY_Chi2DoF_case0= -999;


   // if(debug) std::cout<<"Start XZ study"<<std::endl;
   //extract the T1X hit
   //std::vector<PatHit>  XZHit;
   std::vector<PatHit> Case0First ;
   Case0First.reserve(10);


   std::vector<PatHit> Case0Last ;
   Case0Last.reserve(10);




   std::vector<PatHit> ParabolaSeedHits1;
   ParabolaSeedHits1.reserve(10);
   std::vector<PatHit> ParabolaSeedHits2;
   ParabolaSeedHits2.reserve(10);
   std::vector<PatHit> RemainingCase0Seed1;
   RemainingCase0Seed1.reserve(10);
   std::vector<PatHit> RemainingCase0Seed2;
   RemainingCase0Seed2.reserve(10);
   //std::vector<PatHit> RemainingCase1Seed1(20);
   //std::vector<PatHit> RemainingCase1Seed2(20);
   //std::vector<PatHit> RemainingCase2Seed1(20);
   //std::vector<PatHit> RemainingCase2Seed2(20);

   std::vector<int> remainingCase0Seed1N(3,0);
   std::vector<int> remainingCase0Seed2N(3,0);
   //std::fill(remainingCase0Seed1N,)
   //std::vector<PatHit> test;
   FTCluster cluster;
   Int_t zone = 0;
   //NClonesX=0;
   if(nXUp>=4)  zone = 1;
   if(nxDown>=4) zone = 0;
   PrSeedTrack xHits = PrSeedTrack(zone, m_zReference);
   PrSeedTrack UV = PrSeedTrack(zone,m_zReference);
   //if(PrHit<=8) return kTRUE;
   //xHits.sortbyZ();
   //bool cloneX =false;
   //bool cloneUV = false;

   for( int i = 0; i< PrHit;i++)
   {
      PatHit hit = PatHit();
      hit.setHit(PrHit_Xat0[i],PrHit_Zat0[i],PrHit_dxDy[i],PrHit_dzDy[i],std::sqrt(PrHit_w2[i]),PrHit_yMin[i],PrHit_yMax[i],PrHit_zone[i],PrHit_planeCode[i],PrHit_isX[i],PrHit_LHCbID[i]);
      cluster.setCluster(ChID_Fraction[i],PrHit_Size[i],ChID_Charge[i],ChID_SipmID[i],ChID_SipmCell[i],ChID_Module[i],ChID_Layer[i],ChID_Mat[i],ChID_Quarter[i]);
      hit.SetCluster(cluster);
      if(!PrHit_isX[i])
      {
         UV.addHit(hit);
      }
      if(PrHit_isX[i])
      {
         xHits.addHit(hit);

         if(hit.planeCode() ==4)
         {
            if(debug) std::cout<<"ParabolaSeed Hits 4 loaded"<<std::endl;
            ParabolaSeedHits1.push_back(hit);
         }
         if(hit.planeCode() ==7)
         {
            if(debug) std::cout<<"ParabolaSeed Hits 7 loaded"<<std::endl;
            ParabolaSeedHits2.push_back(hit);
         }
         if(hit.planeCode() == 0){ //T1-1X first in case 0 and case 2
            Case0First.push_back(hit);
            // Case2First.push_back(hit);
            if(debug) std::cout<<"ParabolaSeed Hits 0 loaded"<<std::endl;
         }
         if(hit.planeCode() == 3){ //T1-2X first in case 1
            // Case1First.push_back(hit);
            if(debug) std::cout<<"ParabolaSeed Hits 3 loaded"<<std::endl;
         }
         if(hit.planeCode() == 11){ //T3-2X last in case 0 and case 1
            Case0Last.push_back(hit);
            // Case1Last.push_back(hit);
            if(debug) std::cout<<"ParabolaSeed Hits 11 loaded"<<std::endl;
         }
         //Case 0 Seed 1
         if( hit.planeCode() == 3 || hit.planeCode() ==7 || hit.planeCode() == 8) {
            if(debug) std::cout<<"Plane 3,4,8 Hits loaded"<<std::endl;

            if(hit.planeCode() == 3) remainingCase0Seed1N[0]++;
            if(hit.planeCode() == 7) remainingCase0Seed1N[1]++;
            if(hit.planeCode() == 8) remainingCase0Seed1N[2]++;
            if(debug) std::cout<<"Planes"<<setw(20)<<hit.planeCode()<<"\t Hits done"<<std::endl;
            RemainingCase0Seed1.push_back(hit);
         }

         //Case 0 Seed 2
         if(  hit.planeCode() == 3 || hit.planeCode() ==4 || hit.planeCode() == 8){
            if(debug) std::cout<<"Planes"<<setw(20)<<hit.planeCode()<<"\t Hits done"<<std::endl;
            if(hit.planeCode() ==3) remainingCase0Seed2N[0]++;
            if(hit.planeCode() ==4) remainingCase0Seed2N[1]++;
            if(hit.planeCode() ==8) remainingCase0Seed2N[2]++;
            RemainingCase0Seed2.push_back(hit);
         }
      }
      if(debug) std::cout<<"Done"<<std::endl;
   }
   bool PrelimSele_Case0 = false;
   if(Case0First.size()!=0 && Case0Last.size()!=0 && (ParabolaSeedHits1.size()!=0 || ParabolaSeedHits2.size()!=0)) PrelimSele_Case0 = true;


   
   //std::vector<double> DeltaZ;
   if(PrelimSele_Case0)
   {
      m_physicalCase0Intrinsic++;
      PrSeedTracks xProjections;
      for(int i = 0; i<Case0First.size();i++)
      {
         if(debug) std::cout<<"First Hit "<<std::endl;
         //Case0First[i].PrintHit();
         if(debug) std::cout<<"Loop first Layer"<<std::endl;
         for(int j = 0;j<Case0Last.size();j++)
         {
           for(int k=0; k< ParabolaSeedHits1.size();k ++)
           {
             if(debug) std::cout<<"Parabola Seed Hit 1"<<std::endl;
             //ParabolaSeedHits1[k].PrintHit();
             if(debug) std::cout<<"Loop ParabolaSeed 1"<<std::endl;
             //PrSeedTrack1 xProjection();
             PrSeedTrack xProjection = PrSeedTrack(0,m_zReference);
             xProjection.addHit(ParabolaSeedHits1[k]);
             xProjection.addHit(Case0First[i]);
             xProjection.addHit(Case0Last[j]);
             //xProjection.PrintTrack();
             double a=0.;
             double b=0.;
             double c=0.;
             solveParabola2(Case0First[i],ParabolaSeedHits1[k],Case0Last[j],a,b,c );
             ax_Case0_Seed1 = a;
             bx_Case0_Seed1 = b;
             cx_Case0_Seed1 = c;
             if(debug) std::cout<<"a"<<setw(20)<<a<<std::endl;
             if(debug) std::cout<<"b"<<setw(20)<<b<<std::endl;
             if(debug) std::cout<<"c"<<setw(20)<<c<<std::endl;
               //std::vector<PatHit> bestRemaining(3);
             int best1_kk;
             int best2_kk;
             int best3_kk;
             // PatHit hit1_1;
             // PatHit hit2_1;
             // PatHit hit3_1;
             std::vector<double> minVal(3);
             minVal[0]=1000; minVal[1]=1000; minVal[2]=1000;
             std::vector<double> xAtz;
             
             for(int kk =0;kk<remainingCase0Seed1N[0] ; kk++){
               if(debug) std::cout<<"Remaining Seed 1 -1 "<<std::endl;
               //RemainingCase0Seed1[kk].PrintHit();
               double dz = RemainingCase0Seed1[kk].z() -m_zReference;
               double xAtZ = a*dz*dz*(1.+m_dRatio0*dz) + b*dz +c;
               double delta = std::fabs(RemainingCase0Seed1[kk].x(0)-xAtZ);
               if(debug) std::cout<<"delta Remaining 1 -1 "<<setw(20)<<delta<<std::endl;
               xAtz.push_back(xAtZ);
               if(std::fabs(RemainingCase0Seed1[kk].x(0)-xAtZ)<minVal[0]){
                 minVal[0]=std::fabs(delta);
                    best1_kk=kk;
                    //hit1_1 = RemainingCase0Seed1[kk];
               }
             }
             if(minVal[0]<999){
               xProjection.addHit(RemainingCase0Seed1[best1_kk]);
               DeltaRemaining_Case0->Fill(RemainingCase0Seed1[best1_kk].x(0.)-xAtz[best1_kk]);
               // DeltaRestantiCase0_Seed1.push_back( RemainingCase0Seed1[best1_kk].x(0.)-xAtz[best1_kk]); 
             }
             xAtz.clear();
             for(int kk =remainingCase0Seed1N[0];kk< ( remainingCase0Seed1N[0]+remainingCase0Seed1N[1]) ; kk++){
               double dz = RemainingCase0Seed1[kk].z() - m_zReference;
               double xAtZ = a*dz*dz*(1.+m_dRatio0*dz) + b*dz +c;
               xAtz.push_back(xAtZ);
               double delta = std::fabs(RemainingCase0Seed1[kk].x(0)-xAtZ);
               if(debug) std::cout<<"delta Remaining 1 -2 "<<setw(20)<<delta<<std::endl;
                  DeltaRemaining_Case0->Fill(RemainingCase0Seed1[kk].x(0.)-xAtZ);
                  
                  if(std::fabs(delta)<minVal[1]){
                    minVal[1]=std::fabs(delta);
                    best2_kk=kk;
                    //hit2_1 = RemainingCase0Seed1[kk];
                  }
             }
             if(minVal[1]<999){
               xProjection.addHit(RemainingCase0Seed1[best2_kk]);
               DeltaRemaining_Case0->Fill(RemainingCase0Seed1[best2_kk].x(0.)-xAtz[best2_kk]);
               // DeltaRestantiCase0_Seed1.push_back( RemainingCase0Seed1[best2_kk].x(0.)-xAtz[best2_kk]);
               }
             xAtz.clear();
             for( int kk= (remainingCase0Seed1N[0]+remainingCase0Seed1N[1]);kk<( remainingCase0Seed1N[0]+remainingCase0Seed1N[1]+remainingCase0Seed1N[2] ); kk++){
               double dz = RemainingCase0Seed1[kk].z() -m_zReference;
               double xAtZ = a*dz*dz *(1.+m_dRatio0*dz)+ b*dz +c;
               xAtz.push_back(xAtZ);
               double delta = std::fabs(RemainingCase0Seed1[kk].x(0)-xAtZ);
               if(debug) std::cout<<"delta Remaining 1 -3"<<setw(20)<<delta<<std::endl;
               
               if(std::fabs(RemainingCase0Seed1[kk].x(0)-xAtZ)<minVal[2]){
                 minVal[2]=std::fabs(delta);
                 best3_kk=kk;
                 DeltaRemaining_Case0->Fill(RemainingCase0Seed1[kk].x(0.)-xAtZ);
                 //hit3_1 = RemainingCase0Seed1[kk];
               }
             }
             if(minVal[2]<999){
               xProjection.addHit(RemainingCase0Seed1[best3_kk]);
               DeltaRemaining_Case0->Fill(RemainingCase0Seed1[best3_kk].x(0.)-xAtz[best3_kk]);
               // DeltaRestantiCase0_Seed1.push_back( RemainingCase0Seed1[best3_kk].x(0.)-xAtz[best1_kk]);
             }
             if(debug) std::cout<<"Pushing Back track Case Seed 1"<<std::endl;
             //xProjection.PrintHits();
             if(debug) std::cout<<"Case Seed 1"<<std::endl;
             
             xProjections.push_back(xProjection);
             //T2
           }
           for(int k=0; k< ParabolaSeedHits2.size();k ++)
           {
             if(debug) std::cout<<"Parabola Seed Hit 1"<<std::endl;
             //ParabolaSeedHits2[k].PrintHit();
             //PrSeedTrack1 xProjection();
             PrSeedTrack xProjection = PrSeedTrack(0,m_zReference);
             xProjection.addHit(Case0First[i]);
             xProjection.addHit(ParabolaSeedHits2[k]);
             xProjection.addHit(Case0Last[j]);
             
             double a=0;
             double b=0;
             double c=0;
             solveParabola2(Case0First[i],ParabolaSeedHits2[k],Case0Last[j],a,b,c );
             ax_Case0_Seed2 = a;
             bx_Case0_Seed2 = b;
             cx_Case0_Seed2 = c;
             
             if(debug) std::cout<<"a "<<setw(20)<<a<<std::endl;
             if(debug) std::cout<<"b "<<setw(20)<<b<<std::endl;
             if(debug) std::cout<<"c "<<setw(20)<<c<<std::endl;
             int best1_kk;
             int best2_kk;
             int best3_kk;
             // PatHit hit1_1;
             // PatHit hit2_1;
               // PatHit hit3_1;
             std::vector<double> minVal(3);
             // std::vector<double> DeltaRestantiCase0_Seed2;
             // DeltaRestantiCase0_Seed2.reserve(3);
             minVal[0]=1000; minVal[1]=1000; minVal[2]=1000;
             for(int kk =0;kk<remainingCase0Seed2N[0] ; kk++){
               double dz = RemainingCase0Seed2[kk].z() -m_zReference;
               //double m_dRatio0 = -0.000246;
               double xAtZ = a*dz*dz*(1.+m_dRatio0*dz) + b*dz +c;
               double delta = std::fabs(RemainingCase0Seed2[kk].x(0)-xAtZ);
               if(debug) std::cout<<"delta 2-1 "<<delta<<std::endl;
               if(std::fabs(RemainingCase0Seed2[kk].x(0)-xAtZ)<minVal[0]){
                 minVal[0]=std::fabs(delta);
                 best1_kk = kk;
                 //hit1_1 = RemainingCase0Seed2[kk];
               }
             }
             if(minVal[0]<999){
               xProjection.addHit(RemainingCase0Seed2[best1_kk]);
               // DeltaRestantiCase0_Seed2.push_back( RemainingCase0Seed2[best1_kk].x(0.)-xAtZ[best1_kk]);
             }
             for(int kk =remainingCase0Seed2N[0];kk<(remainingCase0Seed2N[0]+remainingCase0Seed2N[1]) ; kk++){
               double dz = RemainingCase0Seed1[kk].z() -m_zReference;
               double xAtZ = a*dz*dz*(1.+m_dRatio0*dz) + b*dz +c;
               double delta = std::fabs(RemainingCase0Seed2[kk].x(0)-xAtZ);
               if(debug) std::cout<<"Delta 2-2 "<<setw(20)<<delta<<std::endl;
               if(std::fabs(delta)<minVal[1]){
                 minVal[1]=std::fabs(delta);
                 best2_kk=kk;
                 // hit2_1 = RemainingCase0Seed2[kk];
               }
             }
             if(minVal[1]<999){
               xProjection.addHit(RemainingCase0Seed2[best2_kk]);
               // DeltaRestantiCase0_Seed2.push_back( RemainingCase0Seed2[best2_kk].x(0.)-xAtZ[best2_kk]);
             }
             for(int  kk=(remainingCase0Seed2N[0]+remainingCase0Seed2N[1]) ;kk<(remainingCase0Seed2N[0]+remainingCase0Seed2N[1]+remainingCase0Seed2N[2]); kk++){
               double dz = RemainingCase0Seed2[kk].z() -m_zReference;
               // double m_dRatio0 = -0.000246;
               double xAtZ = a*dz*dz *(1.+m_dRatio0*dz) + b*dz +c;

               double delta = std::fabs(RemainingCase0Seed2[kk].x(0)-xAtZ);
               if(debug) std::cout<<"Delta 2-3"<<setw(20)<<delta<<std::endl;
               
               if(std::fabs(RemainingCase0Seed2[kk].x(0)-xAtZ)<minVal[2]){
                 minVal[2]=std::fabs(delta);
                 best3_kk=kk;
                 //hit3_1 = RemainingCase0Seed2[kk];
               }
             }
             if(minVal[2]<999){
               xProjection.addHit(RemainingCase0Seed2[best3_kk]);
               // DeltaRestantiCase0_Seed2.push_back( RemainingCase0Seed2[best3_kk].x(0.)-xAtZ[best3_kk]);
               //                  DeltaRemaining_Case0->Fill(RemainingCase0Seed1[kk].x(0.)-xAtZ);
               
             }
             if(debug) std::cout<<"Pushing Back track Case Seed 1"<<std::endl;
             //xProjection.PrintHits();
             if(debug) std::cout<<"Case Seed 1"<<std::endl;
             
             xProjections.push_back(xProjection);
           }
           
         }
      }
      if(debug) std::cout<<"I found xCandidates"<<xProjections.size()<<std::endl;
      for(int i =0 ; i< xProjections.size() ; i++){
         PrSeedTrack track= xProjections[i];
         track.sortbyZ();
         // PrSeedTrack track = xProjections[0];
         if(debug) std::cout<<"Number of XZ Projections before erase"<<xProjections.size()<<std::endl;
         if(xProjections.size()>1)
         {
           xProjections.erase(std::unique(xProjections.begin(),xProjections.end(),[](PrSeedTrack a, PrSeedTrack b)->bool{
                 if(a.hits().size() != b.hits().size()) return false;
                 a.sortbyZ();
                 b.sortbyZ();
                 bool flag = true;
                 for(int i =0;i<a.hits().size();i++){
                   if( (a.hits()[i].x() != b.hits()[i].x()) && a.hits()[i].z()!=a.hits()[i].z()) {  flag = false; break;}
                 }
                 return flag;
               }),xProjections.end());
         }
      }
      if(xProjections.size()==0){std::cout<<"I didn't find any xProjection, error case 0"<<std::endl; xProj.PrintHits();}
      if(xProjections.size()>1) m_physicalCase0IntrinsicWithduplicates++;
      if(xProjections.size()>1 && xProjections[0].hits().size()>=4) m_physicalCase0IntrinsicWithduplicatesMore4++;
      if(xProjections.size()>1 && xProjections[0].hits().size()>=5) m_physicalCase0IntrinsicWithduplicatesMore5++;
      //if(xProjections.size()>0){
      if(xProjections[0].hits().size()>=4) m_physicalCase0IntrinsicMore4++;
      if(xProjections[0].hits().size()>=5) m_physicalCase0IntrinsicMore5++;

      //New

      if(xProjections.size()>1)
      {
         std::vector<std::pair<PrSeedTrack,double>> trackChi2;
         for(int i = 0;i<xProjections.size();i++){
            fitXProjection(xProjections[i],debug);
         }
         std::sort(xProjections.begin(),xProjections.end(),[](PrSeedTrack val1, PrSeedTrack val2)->bool{return val1.chi2XnDoF()<val2.chi2XnDoF();});
      }
      PrSeedTrack bestXZ = xProjections[0];
      bool OK = false;
      bool OK2 = false;
      OK = fitXProjection(bestXZ,debug);
      OK2 = fitXProjectionCubic(bestXZ,debug);
      if( OK != OK2){
        std::cout<<"OK Parabola = "<<OK<<"\n OK2 Cubic = "<<OK<<std::endl;
      }
      if(OK){
        ax_FitHitCase0 = bestXZ.ax();
        bx_FitHitCase0 = bestXZ.bx();
        cx_FitHitCase0 = bestXZ.cx();
      }
      if(OK2){
        ax_FitHitCase0_Cub = bestXZ.ax_cub();
        bx_FitHitCase0_Cub = bestXZ.bx_cub();
        cx_FitHitCase0_Cub = bestXZ.cx_cub();
        dx_FitHitCase0_Cub = bestXZ.dx_cub();
      }
      if(!OK){
         bestXZ.PrintHits();
         std::cout<<"Case 0:  Chi2 ="<<setw(20)<<bestXZ.chi2X()<<std::endl;
         std::cout<<"         NDOF ="<<setw(20)<<bestXZ.ndofX()<<std::endl;
         std::cout<<"     Chi2NdoF ="<<setw(20)<<bestXZ.chi2XnDoF()<<std::endl;
      }

      //Case 0 first Plane = 0; Last = 11;
      PatHit first = bestXZ.hitInPlane(0);
      // first.PrintHit();
      PatHit last = bestXZ.hitInPlane(11);
      //last.PrintHit();
      Track_Zone = first.zone();
      if(first.zone() != last.zone()) {
         std::cout<<"Killed Track because Up->Down or vice versa"<<std::endl;
         Track_Zone = -99;
         PrelimSele_Case0 = false;
      }
      tx_inf_Case0 = first.x(0.)/first.z(0.);
      //if(NX_Case0 <4){ bestXZ.PrintHits(); }
      NX_Case0= bestXZ.hits().size();
      if(NX_Case0 <4){ Case0Accepted = false; }
      Chi2X_Case0 =bestXZ.chi2X();
      MaxChi2Hit_Case0 = bestXZ.maxChi2Hit();
      tx_inf_Case0 =first.x(0.)/first.z(0.);
      tx_pickedcombination_Case0 = (last.x(0.) - first.x(0.))/(last.z(0.)-first.z(0.));
      x_inf_Case0 =first.x(0.) + tx_inf_Case0 * ( last.z(0.) - first.z(0.));
      Del_x_inf_Case0 = last.x(0.)-x_inf_Case0;
      tx_pickedcombination_Case0 = (last.x(0.) - first.x(0.))/(last.z(0.)-first.z(0.));
      x0_Case0 = first.x(0.) - tx_pickedcombination_Case0 *(first.z(0.) );
      if(bestXZ.hasHitInPlane(4)){
         PatHit parSeed1 = bestXZ.hitInPlane(4);
         Double_t xProectedSeed1 = first.x(0.)+tx_pickedcombination_Case0*(parSeed1.z()-first.z(0.));
         DelxProjectedSeed1_Case0 = parSeed1.x(0.)-xProectedSeed1;
      }
      if(bestXZ.hasHitInPlane(7)){
        PatHit parSeed2 =  bestXZ.hitInPlane(7);
        Double_t xProectedSeed2 = first.x(0.)+tx_pickedcombination_Case0*(parSeed2.z()-first.z(0.));
        DelxProjectedSeed2_Case0 = parSeed2.x(0.)-xProectedSeed2;
      }
      
      int Case = 0;
      Remaining(bestXZ,debug, Case);
      // PatHit parSeed1 = bestXZ.hitInPlane(4);
      // PatHit parSeed2 = bestXZ.hitInPlane(7);
      // std::vector<PatHit> RemainingSeed1;
      // std::vector<PatHit> RemainingSeed2;
      //PatHit parSeed;
      
      for(int i = 0; i<PrHit; i++){
         PatHit hit = PatHit();
         hit.setHit(PrHit_Xat0[i],PrHit_Zat0[i],PrHit_dxDy[i],PrHit_dzDy[i],std::sqrt(PrHit_w2[i]),PrHit_yMin[i],PrHit_yMax[i],PrHit_zone[i],PrHit_planeCode[i],PrHit_isX[i],PrHit_LHCbID[i]);
         FTCluster cluster;
         cluster.setCluster(ChID_Fraction[i],PrHit_Size[i],ChID_Charge[i],ChID_SipmID[i],ChID_SipmCell[i],ChID_Module[i],ChID_Layer[i],ChID_Mat[i],ChID_Quarter[i]);
         hit.SetCluster(cluster);
         if(!hit.isX()) hit.setCoord( -(hit.x(0.)-xProj.x(hit.z(0.)))/(hit.dxDy() * hit.z(0)));
         if(!hit.isX())bestXZ.addHit(hit);
      }
      bool ok =  FitSimultaneouslyXY(bestXZ, UVSegment,debug);

      int nT1U = bestXZ.nHitsInPlane( 1 );
      int nT1V = bestXZ.nHitsInPlane( 2 );
      int nT2U = bestXZ.nHitsInPlane( 5 );
      int nT2V = bestXZ.nHitsInPlane( 6 );
      int nT3U = bestXZ.nHitsInPlane( 9 );
      int nT3V = bestXZ.nHitsInPlane( 10 );

      std::vector<PatHit> hitT1U = bestXZ.hitsinPlane(1);
      std::vector<PatHit> hitT1V = bestXZ.hitsinPlane(2);
      std::vector<PatHit> hitT2U = bestXZ.hitsinPlane(5);
      std::vector<PatHit> hitT2V = bestXZ.hitsinPlane(6);
      std::vector<PatHit> hitT3U = bestXZ.hitsinPlane(9);
      std::vector<PatHit> hitT3V = bestXZ.hitsinPlane(10);

      std::vector<PatHit> HitsUVT1 = hitT1U;
      HitsUVT1.insert(HitsUVT1.end(),hitT1V.begin(),hitT1V.end());

      std::vector<PatHit> HitsUVT2 = hitT2U;
      HitsUVT2.insert(HitsUVT2.end(),hitT2V.begin(),hitT2V.end());

      std::vector<PatHit> HitsUVT3 = hitT3U;
      HitsUVT3.insert(HitsUVT3.end(),hitT3V.begin(),hitT3V.end());

      //if(nT1U>1){std::cout<<"Has duplicate in T1U"<<std::endl; bestXZ.PrintHits(); bestXZ.removeAllHitInPlane(1); std::cout<<"removed duplicates in plane 1"<<std::endl; bestXZ.PrintHits();}
      PrSeedTracks bestXZplusUV;
      if(nT1U+nT1V ==0) {PrelimSele_Case0 = false; std::cout<<"NO Hits in U + V stereo T1!"<<std::endl;}
      if(nT2U+nT2V ==0) {PrelimSele_Case0 = false; std::cout<<"NO Hits in U + V stereo T2!"<<std::endl;}
      if(nT3U+nT3V ==0) {PrelimSele_Case0 = false; std::cout<<"NO Hits in U + V stereo T3!"<<std::endl;}
      bool FindBest = false;
      PrSeedTrack beforeFindBest = bestXZ;

      if( (nT1U+nT1V)*(nT2U+nT2V)*(nT3U+nT3V) != 0){ //reconstructibility
         if(nT1U>1){bestXZ.removeAllHitInPlane(1); FindBest= true;}
         if(nT1V>1){bestXZ.removeAllHitInPlane(2); FindBest= true;}
         if(nT2U>1){bestXZ.removeAllHitInPlane(5); FindBest = true;}
         if(nT2V>1){bestXZ.removeAllHitInPlane(6); FindBest = true;}
         if(nT3U>1){bestXZ.removeAllHitInPlane(9); FindBest = true;}
         if(nT3V>1){bestXZ.removeAllHitInPlane(10); FindBest = true;}
         PrSeedTracks xzplussingleUVcombinations;
         bestXZ.sortbyZ();
      }
      double minChi2 = bestXZ.Chi2()+1e3;
      PrSeedTrack BestCombination = PrSeedTrack(0, m_zReference);
      if(FindBest){ // if you enter here you have removed in one of the layer all the hits
         for(int i = 0 ; i< HitsUVT1.size();i++){
            if(bestXZ.nHitsInPlane(HitsUVT1[i].planeCode()) == 0){
               bestXZ.addHit(HitsUVT1[i]);
            }
            if(bestXZ.nHitsInPlane(HitsUVT1[i].planeCode()) >0 ){
               bestXZ.substituteHitinSameLayer( HitsUVT1[i]); //will be dummy if already present 1 single hit, if not it really makes the susbtitution
            }
            for(int j = 0 ; j<HitsUVT2.size();j++){
               if(bestXZ.nHitsInPlane(HitsUVT2[j].planeCode())==0){
                  bestXZ.addHit(HitsUVT2[j]);
               }
               if(bestXZ.nHitsInPlane(HitsUVT2[j].planeCode())>0){
                  bestXZ.substituteHitinSameLayer( HitsUVT2[j]);
               }
               for( int k = 0; k< HitsUVT3.size(); k++){
                  if(bestXZ.nHitsInPlane(HitsUVT3[k].planeCode())==0){
                     bestXZ.addHit(HitsUVT3[k]);
                  }
                  if(bestXZ.nHitsInPlane(HitsUVT3[k].planeCode())>0){
                     bestXZ.substituteHitinSameLayer(HitsUVT3[k]);

                     bestXZ.sortbyZ();
                     //bestXZ.PrintTrack();
                     bool Ok =  FitSimultaneouslyXY(bestXZ, UVSegment, false);
                     if(bestXZ.Chi2()< minChi2){
                        minChi2 = bestXZ.Chi2();
                        BestCombination = bestXZ;
                     }
                  }
               }
            }
         }
         if(debug){
         //need to find best combination of the UV layers hit
            if(nT1U>1)
            {  std::cout<<"Had duplicate in T1U"<<std::endl;}
            if(nT1V>1)
            {  std::cout<<"Had duplicate in T1V"<<std::endl;}
            if(nT1V>1)
            {  std::cout<<"Had duplicate in T1V"<<std::endl;}
            if(nT2U>1)
            {  std::cout<<"Had duplicate in T2U"<<std::endl;}
            if(nT2V>1)
            {  std::cout<<"Had duplicate in T2V"<<std::endl;}
            if(nT3U>1)
            {  std::cout<<"Had duplicate in T3U"<<std::endl;}
            if(nT3V>1)
            {  std::cout<<"Had duplicate in T3V"<<std::endl;}
         }
         // std::cout<<"\n \n ===================Before Finding the best track==================="<<std::endl;
         // beforeFindBest.PrintHits();
         // beforeFindBest.PrintTrack();
         // std::cout<<"\n \n ===========================Best Combination=====================\n\n"<<std::endl;
         // BestCombination.PrintHits();
         BestCombination.PrintTrack();
         // std::cout<<"\n \n"<<std::endl;
         // for(layer1){
         bestXZ=BestCombination;
      }
      PrSeedTrack BestUVSegment = PrSeedTrack(bestXZ.zone(), m_zReference);
      for(int i =0; i<bestXZ.hits().size(); i++){
         if( !bestXZ.hits()[i].isX()){
            BestUVSegment.addHit(bestXZ.hits()[i]);
         }
      }
      if(debug) {
         std::cout<<"Best UV Segment Case0"<<std::endl;
         BestUVSegment.PrintHits();
      }

      //Fit a single Line for It.
      BestUVSegment.setParameters( ax_FitHitCase0, bx_FitHitCase0, cx_FitHitCase0 ,0,0);
      FitLine(BestUVSegment,debug, 0);
      NUV_Case0 = (int)BestUVSegment.hits().size();
      
      //std::cout<<"ax = "<<bestXZ.ax();
      ax_fullFitCase0 = (Double_t )bestXZ.ax();
      bx_fullFitCase0= (Double_t )bestXZ.bx();
      cx_fullFitCase0 =  (Double_t )bestXZ.cx();
      ay_fullFitCase0 =  (Double_t )bestXZ.ay();
      by_fullFitCase0 =  (Double_t )bestXZ.by();
      Chi2_fullFitCase0 =  (Double_t )bestXZ.Chi2();
      MaxChi2_fullFitCase0 =  (Double_t )bestXZ.maxChi2Hit();
      NHits_fullFitCase0= (int)bestXZ.hits().size();
      y0_fullFitCase0 =  (Double_t )bestXZ.y(0.);
      if(debug) std::cout<<"Number of XZ Projections after erase"<<xProjections.size()<<std::endl;
      bestXZ.sortbyZ();
      
      WorstPlane = bestXZ.worstPlane();
      nbHitsInWorst = 0;
      for(int i = 0;i<PrHit;i++){
        if(PrHit_planeCode[i]==WorstPlane) nbHitsInWorst++;
      }
      if(nbHitsInWorst>1){
        bestXZ.PrintHits();
        std::cout<<"Nb Hits in worst Plane is"<<setw(20)<<nbHitsInWorst<<setw(20)<<" and worst plane = "<<WorstPlane<<std::endl;
      }
      
      if(bestXZ.hasHitInPlane(1)){
        PatHit FirstY = bestXZ.hitInPlane(1);
        X_T1 = FirstY.x(0.);
        Y_T1 = FirstY.coord()*FirstY.z(0.);
        Coord_T1 = FirstY.coord();
        Z_T1 = FirstY.z();
        //std::cout<<FirstY.coord()<<std::endl;
      }
      if(!bestXZ.hasHitInPlane(1)&& bestXZ.hasHitInPlane(2)){
        PatHit FirstY = bestXZ.hitInPlane(2);
        X_T1 = FirstY.x(0.);
        Y_T1 = FirstY.coord()*FirstY.z(0.);
        Coord_T1 = FirstY.coord();
        Z_T1 = FirstY.z();
      }
      if(bestXZ.hasHitInPlane(10)){
        PatHit LastY = bestXZ.hitInPlane(10);
        X_T3 = LastY.x(0.);
        Y_T3 =LastY.coord()*LastY.z(0.);
        Coord_T3 = LastY.coord();
        Z_T3 = LastY.z(0.);
      }
      if(!bestXZ.hasHitInPlane(10) && bestXZ.hasHitInPlane(9)){
        PatHit LastY = bestXZ.hitInPlane(9);
        Y_T3 = LastY.coord()*LastY.z(0.);
        Coord_T3 = LastY.coord();
        Z_T3 = LastY.z(0.);
        X_T3 = LastY.x(0.);
      }
      //BestUVSegment.sortbyZ();
      //BestUVSegment.PrintHits();
   }
   //End Case 0 Prelim Sele
   return PrelimSele_Case0;
}

void TrackStudy::Remaining(PrSeedTrack & bestXZ, Bool_t debug, int Case){
  int firstplane = -99;
  int lastplane = -99;
  int seed_1 = 4;
  int seed_2 = 7;
  if(Case ==0){
    firstplane = 0;
    lastplane = 11;
  }
  if(Case==1){
    firstplane = 3;
    lastplane = 11;
  }
  if(Case==2){
    firstplane = 0;
    lastplane = 8;
  }
  PatHit first = bestXZ.hitInPlane(firstplane);
  if( first.planeCode() == -1){
    std::cout<<" No Hit in First Plane "<<std::endl;
  }
  PatHit last = bestXZ.hitInPlane(lastplane);
  if( last.planeCode() == -1){
    std::cout<<"No Hit in Last Plane"<<std::endl;
  }
  PatHit seed1 = bestXZ.hitInPlane(seed_1);
  PatHit seed2 = bestXZ.hitInPlane(seed_2);
  std::vector<PatHit> remainingSeed1;
  std::vector<PatHit> remainingSeed2;
  bool hasSeed1 = true;
  bool hasSeed2 = true;
  if( seed1.planeCode() == -1){
    hasSeed1 = false;
    //std::cout<<"No Hit in Plane 4"<<std::endl;
  }
  if( seed2.planeCode() == -1){
    hasSeed2 = false;
  }
  for( int i = 0; i<bestXZ.hits().size(); i++){
    PatHit hit = bestXZ.hits()[i];
    if( hit.planeCode() != firstplane &&
        hit.planeCode() != lastplane &&
        hit.planeCode() != seed_1){
      if(hasSeed1){ // You have an hit in plane seed 1
        remainingSeed1.push_back(hit);
      }
    }
    if( hit.planeCode() != firstplane &&
        hit.planeCode() != lastplane &&
        hit.planeCode() != seed_2){
      if(hasSeed2){
        remainingSeed2.push_back(hit);
      }
    }
  }
  double a_1 = 0; double a_2 = 0;
  double b_1 = 0; double b_2 = 0;
  double c_1 = 0; double c_2 = 0;
  std::vector<double> DeltaSeed1;
  std::vector<double> DeltaSeed2;
  if(Case ==0){
    nRemSeed1_Case0 = 0;
    DeltaSeed1_Case0[0] =-99.;
    DeltaSeed1_Case0[1] =-99.;
    DeltaSeed1_Case0[2] =-99.;
    DeltaSeed2_Case0[0] =-99.;
    DeltaSeed2_Case0[1] =-99.;
    DeltaSeed2_Case0[2] =-99.;
  }
  if(Case ==1){
    nRemSeed1_Case1= 0;
    DeltaSeed1_Case1[0] =-99.;
    DeltaSeed1_Case1[1] =-99.;
    DeltaSeed1_Case1[2] =-99.;
    DeltaSeed2_Case1[0] =-99.;
    DeltaSeed2_Case1[1] =-99.;
    DeltaSeed2_Case1[2] =-99.;
  }
  if(Case ==2){
    nRemSeed1_Case2 = 0;
    DeltaSeed1_Case2[0] = -99.;
    DeltaSeed1_Case2[1] = -99.;
    DeltaSeed1_Case2[2] = -99.;
    DeltaSeed2_Case2[0] =-99.;
    DeltaSeed2_Case2[1] =-99.;
    DeltaSeed2_Case2[2] =-99.;
  }
  
  if(hasSeed1){
    if(Case ==0) nRemSeed1_Case0 = (Int_t)remainingSeed1.size();
    if(Case ==1) nRemSeed1_Case1 = (Int_t)remainingSeed1.size();
    if(Case ==2) nRemSeed1_Case2 = (Int_t)remainingSeed1.size();

    solveParabola2( first, seed1 , last, a_1 , b_1, c_1);
    for( int i =0 ; i< remainingSeed1.size();i++){
      //remainingSeed1[i].PrintHit();
      double dz = remainingSeed1[i].z(0.)-m_zReference;
      double xatz = a_1*dz*dz* (1.+m_dRatio0*dz) + b_1*dz + c_1;
      double delta = remainingSeed1[i].x(0.)-xatz;
      //DeltaSeed1.push_back(delta);
      if(Case ==0) DeltaSeed1_Case0[i] = delta;
      if(Case ==1) DeltaSeed1_Case1[i] = delta;
      if(Case ==2) DeltaSeed1_Case2[i] = delta;
    }
  }
  if(Case==0) nRemSeed2_Case0 = 0;
  if(hasSeed2){
    if(Case ==0) nRemSeed2_Case0 = (Int_t)remainingSeed2.size();
    if(Case ==1) nRemSeed2_Case1 = (Int_t)remainingSeed2.size();
    if(Case ==2) nRemSeed2_Case2 = (Int_t)remainingSeed2.size();

    solveParabola2( first,seed2, last, a_2,b_2,c_2);
    for( int i =0; i<remainingSeed2.size(); i++){
      double dz = remainingSeed2[i].z(0.)-m_zReference;
      double xatz = a_2*dz*dz*(1.+m_dRatio0*dz) + b_2*dz + c_2;
      double delta = remainingSeed2[i].x(0.)-xatz;
      //DeltaSeed2.push_back(delta);
      if(Case==0) DeltaSeed2_Case0[i] = delta;
      if(Case==1) DeltaSeed2_Case1[i] = delta;
      if(Case==2) DeltaSeed2_Case2[i] = delta;
    }
  }
}
bool TrackStudy::XZStudyCase1( PrSeedTrack& xProj, PrSeedTrack& UVSegment , Bool_t debug){
  NUV_Case1 = -10;
  ax_Case1_Seed1 = -99.;
  bx_Case1_Seed1 = -99.;
  cx_Case1_Seed1 = -99.;
  ax_Case1_Seed1 = -99.;
  bx_Case1_Seed1 = -99.;
  cx_Case1_Seed1 = -99.;

   Chi2X_Case1 = -99.;
   NX_Case1 = -1;
   tx_inf_Case1 = -99.;
   x0_Case1 = -99.;
   x_inf_Case1 = -99.;
   Del_x_inf_Case1 = -99.;
   tx_pickedcombination_Case1 = -99.;
   DelxProjectedSeed2_Case1 = -99.;
   DelxProjectedSeed1_Case1 = -99.;
   MaxChi2Hit_Case1= -99;

   ax_FitHitCase1 = -99.;
   bx_FitHitCase1 = -99.;
   cx_FitHitCase1 = -99.;

   ax_fullFitCase1 = -9999.;
   bx_fullFitCase1= -99.;
   cx_fullFitCase1 = -99.;
   ay_fullFitCase1 = -9999.;
   by_fullFitCase1 = -99.;
   Chi2_fullFitCase1 =  -99.;
   MaxChi2_fullFitCase1 = -99.;
   NHits_fullFitCase1 = -1;
   y0_fullFitCase1 = -9999.;



   std::vector<PatHit> Case1First ;
   Case1First.reserve(10);
   std::vector<PatHit> Case1Last  ;
   Case1Last.reserve(10);
   std::vector<PatHit> ParabolaSeedHits1;
   ParabolaSeedHits1.reserve(10);
   std::vector<PatHit> ParabolaSeedHits2;
   ParabolaSeedHits2.reserve(10);
   std::vector<PatHit> RemainingCase1Seed1;
   RemainingCase1Seed1.reserve(10);
   std::vector<PatHit> RemainingCase1Seed2;
   RemainingCase1Seed2.reserve(10);
   std::vector<int> remainingCase1Seed1N(3,0);
   std::vector<int> remainingCase1Seed2N(3,0);




   for(int i =0; i<xProj.hits().size();i++){
      PatHit hit = xProj.hits()[i];
      if(hit.planeCode()==3){
         Case1First.push_back(hit);
      }
      if(hit.planeCode()==11){
         Case1Last.push_back(hit);
      }
      if(hit.planeCode()==4){
         ParabolaSeedHits1.push_back(hit);
      }
      if(hit.planeCode()==7){
         ParabolaSeedHits2.push_back(hit);
      }
      if( hit.planeCode() == 0 || hit.planeCode() ==7 || hit.planeCode() == 8) {
         if(debug) std::cout<<"Plane 3,4,8 Hits loaded"<<std::endl;

         if(hit.planeCode() == 0) remainingCase1Seed1N[0]++;
         if(hit.planeCode() == 7) remainingCase1Seed1N[1]++;
         if(hit.planeCode() == 8) remainingCase1Seed1N[2]++;
         RemainingCase1Seed1.push_back(hit);
      }
      if(  hit.planeCode() == 0 || hit.planeCode() ==4 || hit.planeCode() == 8){
         if(debug) std::cout<<"Planes"<<setw(20)<<hit.planeCode()<<"\t Hits done"<<std::endl;
         if(hit.planeCode() ==0) remainingCase1Seed2N[0]++;
         if(hit.planeCode() ==4) remainingCase1Seed2N[1]++;
         if(hit.planeCode() ==8) remainingCase1Seed2N[2]++;
         RemainingCase1Seed2.push_back(hit);
      }
   }
   bool PrelimSele_Case1 = false;
   if(Case1First.size()!=0 && Case1Last.size()!=0 && (ParabolaSeedHits1.size()!=0 || ParabolaSeedHits2.size()!=0)){
      PrelimSele_Case1 = true;
   }

   if(PrelimSele_Case1){
      m_physicalCase1Intrinsic++;
      PrSeedTracks xProjections;
      for(int i = 0; i<Case1First.size();i++)
      {
         if(debug) std::cout<<"First Hit "<<std::endl;
         //Case1First[i].PrintHit();
         if(debug) std::cout<<"Loop first Layer"<<std::endl;
         for(int j = 0;j<Case1Last.size();j++)
         {
            if(debug) std::cout<<"Second Hit"<<std::endl;
            if(debug) std::cout<<"Loop Last Layer"<<std::endl;

            if(debug) std::cout<<"Loop Last "<<std::endl;
            if(debug) std::cout<<"parabolaSeedHits1 sizee   "<<ParabolaSeedHits1.size()<<std::endl;
            //FInd Given the parabolaSeedHit1 the best combination
            for(int k=0; k< ParabolaSeedHits1.size();k ++)
            {
               if(debug) std::cout<<"Parabola Seed Hit 1"<<std::endl;
               //ParabolaSeedHits1[k].PrintHit();
               if(debug) std::cout<<"Loop ParabolaSeed 1"<<std::endl;
               //PrSeedTrack1 xProjection();
               PrSeedTrack xProjection = PrSeedTrack(0,m_zReference);
               xProjection.addHit(ParabolaSeedHits1[k]);
               xProjection.addHit(Case1First[i]);
               xProjection.addHit(Case1Last[j]);
               //xProjection.PrintTrack();
               double a=0.;
               double b=0.;
               double c=0.;
               solveParabola2(Case1First[i],ParabolaSeedHits1[k],Case1Last[j],a,b,c );
               ax_Case1_Seed1 = a;
               bx_Case1_Seed1 = b;
               cx_Case1_Seed1 = c;
               if(debug) std::cout<<"a"<<setw(20)<<a<<std::endl;
               if(debug) std::cout<<"b"<<setw(20)<<b<<std::endl;
               if(debug) std::cout<<"c"<<setw(20)<<c<<std::endl;

               //std::vector<PatHit> bestRemaining(3);
               int best1_kk;
               int best2_kk;
               int best3_kk;
               // PatHit hit1_1;
               // PatHit hit2_1;
               // PatHit hit3_1;
               std::vector<double> minVal(3);
               // std::vector<double> DeltaRestantiCase1_Seed1;
               // std::vector<double> DeltaRestantiCase1_Seed2;
               // DeltaRestantiCase1_Seed1.reserve(3);
               // DeltaRestantiCase1_Seed2.reserve(3);
               minVal[0]=1000; minVal[1]=1000; minVal[2]=1000;
               for(int kk =0;kk<remainingCase1Seed1N[0] ; kk++){
                  if(debug) std::cout<<"Remaining Seed 1 -1 "<<std::endl;
                  //RemainingCase1Seed1[kk].PrintHit();
                  double dz = RemainingCase1Seed1[kk].z() -m_zReference;
                  double xAtZ = a*dz*dz*(1.+m_dRatio0*dz) + b*dz +c;
                  double delta = std::fabs(RemainingCase1Seed1[kk].x(0)-xAtZ);
                  if(debug) std::cout<<"delta Remaining 1 -1 "<<setw(20)<<delta<<std::endl;
                  if(std::fabs(RemainingCase1Seed1[kk].x(0)-xAtZ)<minVal[0]){
                     minVal[0]=std::fabs(delta);
                     best1_kk=kk;
                     //hit1_1 = RemainingCase1Seed1[kk];
                  }
               }
               if(minVal[0]<999){
                  xProjection.addHit(RemainingCase1Seed1[best1_kk]);
                  // DeltaRestantiCase1_Seed1.push_back( RemainingCase1Seed1[best1_kk].x(0.)-xAtz[best1_kk]);
               }
               for(int kk =remainingCase1Seed1N[0];kk< ( remainingCase1Seed1N[0]+remainingCase1Seed1N[1]) ; kk++){
                  double dz = RemainingCase1Seed1[kk].z() - m_zReference;
                  double xAtZ = a*dz*dz*(1.+m_dRatio0*dz) + b*dz +c;
                  double delta = std::fabs(RemainingCase1Seed1[kk].x(0)-xAtZ);
                  if(debug) std::cout<<"delta Remaining 1 -2 "<<setw(20)<<delta<<std::endl;

                  if(std::fabs(delta)<minVal[1]){
                     minVal[1]=std::fabs(delta);
                     best2_kk=kk;
                     // DeltaRestantiCase1_Seed1.push_back( RemainingCase1Seed1[best2_kk].x(0.)-xAtz[best2_kk]);
                     //hit2_1 = RemainingCase1Seed1[kk];
                  }
               }
               if(minVal[1]<999){
                  xProjection.addHit(RemainingCase1Seed1[best2_kk]);
               }
               for( int kk= (remainingCase1Seed1N[0]+remainingCase1Seed1N[1]);kk<( remainingCase1Seed1N[0]+remainingCase1Seed1N[1]+remainingCase1Seed1N[2] ); kk++){
                  double dz = RemainingCase1Seed1[kk].z() -m_zReference;
                  double xAtZ = a*dz*dz *(1.+m_dRatio0*dz)+ b*dz +c;
                  double delta = std::fabs(RemainingCase1Seed1[kk].x(0)-xAtZ);
                  if(debug) std::cout<<"delta Remaining 1 -3"<<setw(20)<<delta<<std::endl;

                  if(std::fabs(RemainingCase1Seed1[kk].x(0)-xAtZ)<minVal[2]){
                     minVal[2]=std::fabs(delta);
                     best3_kk=kk;
                     //hit3_1 = RemainingCase1Seed1[kk];
                  }
               }
               if(minVal[2]<999){
                  xProjection.addHit(RemainingCase1Seed1[best3_kk]);
                  // DeltaRestantiCase1_Seed1.push_back( RemainingCase1Seed1[best3_kk].x(0.)-xAtz[best3_kk]); 
               }
               if(debug) std::cout<<"Pushing Back track Case Seed 1"<<std::endl;
               //xProjection.PrintHits();
               if(debug) std::cout<<"Case Seed 1"<<std::endl;

               xProjections.push_back(xProjection);
               //T2
            }
            for(int k=0; k< ParabolaSeedHits2.size();k ++)
            {
               if(debug) std::cout<<"Parabola Seed Hit 1"<<std::endl;
               //ParabolaSeedHits2[k].PrintHit();
               //PrSeedTrack1 xProjection();
               PrSeedTrack xProjection = PrSeedTrack(0,m_zReference);
               xProjection.addHit(Case1First[i]);
               xProjection.addHit(ParabolaSeedHits2[k]);
               xProjection.addHit(Case1Last[j]);

               double a=0;
               double b=0;
               double c=0;
               solveParabola2(Case1First[i],ParabolaSeedHits2[k],Case1Last[j],a,b,c );
               ax_Case1_Seed2 = a;
               bx_Case1_Seed2 = b;
               cx_Case1_Seed2 = c;
               if(debug) std::cout<<"a "<<setw(20)<<a<<std::endl;
               if(debug) std::cout<<"b "<<setw(20)<<b<<std::endl;
               if(debug) std::cout<<"c "<<setw(20)<<c<<std::endl;
               int best1_kk;
               int best2_kk;
               int best3_kk;
               // PatHit hit1_1;
               // PatHit hit2_1;
               // PatHit hit3_1;
               std::vector<double> minVal(3);
               minVal[0]=1000; minVal[1]=1000; minVal[2]=1000;
               for(int kk =0;kk<remainingCase1Seed2N[0] ; kk++){
                 double dz = RemainingCase1Seed2[kk].z() -m_zReference;
                 //double m_dRatio0 = -0.000246;
                 double xAtZ = a*dz*dz*(1.+m_dRatio0*dz) + b*dz +c;
                 double delta = std::fabs(RemainingCase1Seed2[kk].x(0)-xAtZ);
                 if(debug) std::cout<<"delta 2-1 "<<delta<<std::endl;
                 if(std::fabs(RemainingCase1Seed2[kk].x(0)-xAtZ)<minVal[0]){
                   minVal[0]=std::fabs(delta);
                   best1_kk = kk;
                   //hit1_1 = RemainingCase1Seed2[kk];
                 }
               }
               if(minVal[0]<999){
                  xProjection.addHit(RemainingCase1Seed2[best1_kk]);
               }
               for(int kk =remainingCase1Seed2N[0];kk<(remainingCase1Seed2N[0]+remainingCase1Seed2N[1]) ; kk++){
                  double dz = RemainingCase1Seed1[kk].z() -m_zReference;
                  double xAtZ = a*dz*dz*(1.+m_dRatio0*dz) + b*dz +c;
                  double delta = std::fabs(RemainingCase1Seed2[kk].x(0)-xAtZ);
                  if(debug) std::cout<<"Delta 2-2 "<<setw(20)<<delta<<std::endl;
                  if(std::fabs(delta)<minVal[1]){
                     minVal[1]=std::fabs(delta);
                     best2_kk=kk;
                     // hit2_1 = RemainingCase1Seed2[kk];
                  }
               }
               if(minVal[1]<999){
                  xProjection.addHit(RemainingCase1Seed2[best2_kk]);
               }
               for(int  kk=(remainingCase1Seed2N[0]+remainingCase1Seed2N[1]) ;kk<(remainingCase1Seed2N[0]+remainingCase1Seed2N[1]+remainingCase1Seed2N[2]); kk++){
                  double dz = RemainingCase1Seed2[kk].z() -m_zReference;
                  // double m_dRatio0 = -0.000246;
                  double xAtZ = a*dz*dz *(1.+m_dRatio0*dz) + b*dz +c;

                  double delta = std::fabs(RemainingCase1Seed2[kk].x(0)-xAtZ);
                  if(debug) std::cout<<"Delta 2-3"<<setw(20)<<delta<<std::endl;

                  if(std::fabs(RemainingCase1Seed2[kk].x(0)-xAtZ)<minVal[2]){
                     minVal[2]=std::fabs(delta);
                     best3_kk=kk;
                     //hit3_1 = RemainingCase1Seed2[kk];
                  }
               }
               if(minVal[2]<999){
                  xProjection.addHit(RemainingCase1Seed2[best3_kk]);
                  // DeltaRemaining_Case0
               }
               if(debug) std::cout<<"Pushing Back track Case Seed 1"<<std::endl;
               //xProjection.PrintHits();
               if(debug) std::cout<<"Case Seed 1"<<std::endl;

               xProjections.push_back(xProjection);

            }

         }
      }
      for(int i =0 ; i< xProjections.size() ; i++){
         PrSeedTrack track= xProjections[i];
         track.sortbyZ();
         // PrSeedTrack track = xProjections[0];
         if(debug) std::cout<<"Number of XZ Projections before erase"<<xProjections.size()<<std::endl;
         if(xProjections.size()>1)
         {
           xProjections.erase(std::unique(xProjections.begin(),xProjections.end(),[](PrSeedTrack a, PrSeedTrack b)->bool{
                 if(a.hits().size() != b.hits().size()) return false;
                 a.sortbyZ();
                 b.sortbyZ();
                 bool flag = true;
                 for(int i =0;i<a.hits().size();i++){
                   if( (a.hits()[i].x() != b.hits()[i].x()) && a.hits()[i].z()!=a.hits()[i].z()) {  flag = false; break;}
                 }
                 return flag;
               }),xProjections.end());
         }
      }
      if(xProjections.size()==0){std::cout<<"I didn't find any xProjection, error case 1"<<std::endl; xProj.PrintHits();}
      if(xProjections.size()>1) m_physicalCase1IntrinsicWithduplicates++;
      if(xProjections.size()>1 && xProjections[0].hits().size()>=4) m_physicalCase1IntrinsicWithduplicatesMore4++;
      if(xProjections.size()>1 && xProjections[0].hits().size()>=5) m_physicalCase1IntrinsicWithduplicatesMore5++;
      //if(xProjections.size()>0){
      if(xProjections[0].hits().size()>=4) m_physicalCase1IntrinsicMore4++;
      if(xProjections[0].hits().size()>=5) m_physicalCase1IntrinsicMore5++;
      //}

      if(xProjections.size()>1)
      {
         std::vector<std::pair<PrSeedTrack,double>> trackChi2;
         for(int i = 0;i<xProjections.size();i++){
            fitXProjection(xProjections[i],debug);
         }
         std::sort(xProjections.begin(),xProjections.end(),[](PrSeedTrack val1, PrSeedTrack val2)->bool{return val1.chi2XnDoF()<val2.chi2XnDoF();});
      }
      PrSeedTrack bestXZ = xProjections[0];
      bool OK = false;
      bestXZ.sortbyZ();
      OK = fitXProjection(bestXZ,debug);
      // if(bestXZ.chi2X()/bestXZ.ndofX()>4. && Momentum>5000){
      //    std::cout<<"NbSplitted"<<setw(15)<<NSplittedX<<setw(20)<<"P"<<setw(15)<<Momentum<<std::endl;
      //    bestXZ.PrintHits();
      // }
      ax_FitHitCase1 = bestXZ.ax();
      bx_FitHitCase1 = bestXZ.bx();
      cx_FitHitCase1 = bestXZ.cx();

      if(!OK){
         bestXZ.PrintHits();
         std::cout<<"Case 0:  Chi2 ="<<setw(20)<<bestXZ.chi2X()<<std::endl;
         std::cout<<"         NDOF ="<<setw(20)<<bestXZ.ndofX()<<std::endl;
         std::cout<<"     Chi2NdoF ="<<setw(20)<<bestXZ.chi2XnDoF()<<std::endl;
      }


      //Case 0 first Plane = 0; Last = 11;
      PatHit first = bestXZ.hitInPlane(3);
      // first.PrintHit();
      PatHit last = bestXZ.hitInPlane(11);
      //last.PrintHit();
      Track_Zone = first.zone();
      if(first.zone() != last.zone()) {
         std::cout<<"Killed Track because Up->Down or vice versa"<<std::endl;
         PrelimSele_Case1 = false;
         Track_Zone =-99;
      }
      tx_inf_Case1 = first.x(0.)/first.z(0.);
      NX_Case1= bestXZ.hits().size();
      if(NX_Case1 <4){Case1Accepted = false;}
      Chi2X_Case1 =bestXZ.chi2X();
      MaxChi2Hit_Case1 = bestXZ.maxChi2Hit();
      tx_inf_Case1 =first.x(0.)/first.z(0.);
      // last_case1 = last.x(0.);
      // first_case1 = first.x(0.);

      tx_pickedcombination_Case1 = (last.x(0.) - first.x(0.))/(last.z(0.)-first.z(0.));
      x_inf_Case1 = first.x(0.)+ tx_inf_Case1 * ( last.z(0.) - first.z(0.));
      Del_x_inf_Case1 = last.x(0.)-x_inf_Case1;
      tx_pickedcombination_Case1 = (last.x(0.) - first.x(0.))/(last.z(0.)-first.z(0.));
      x0_Case1 = first.x(0.) - tx_pickedcombination_Case1 *(first.z(0.) );
      if(bestXZ.hasHitInPlane(4)){
         PatHit parSeed1 = bestXZ.hitInPlane(4);
         Double_t xProectedSeed1 = first.x(0.)+tx_pickedcombination_Case1*( parSeed1.z(0.)-first.z(0.) );
         DelxProjectedSeed1_Case1 = parSeed1.x(0.)-xProectedSeed1;
      }
      if(bestXZ.hasHitInPlane(7)){
         PatHit parSeed2 =  bestXZ.hitInPlane(7);
         Double_t xProectedSeed2 = first.x(0.)+tx_pickedcombination_Case1*(parSeed2.z()-first.z(0.));
         DelxProjectedSeed2_Case1 = parSeed2.x(0.)-xProectedSeed2;
      }


      int Case = 1;
      Remaining(bestXZ,debug,Case);
      
      for(int i = 0; i<PrHit; i++){
         PatHit hit = PatHit();
         hit.setHit(PrHit_Xat0[i],PrHit_Zat0[i],PrHit_dxDy[i],PrHit_dzDy[i],std::sqrt(PrHit_w2[i]),PrHit_yMin[i],PrHit_yMax[i],PrHit_zone[i],PrHit_planeCode[i],PrHit_isX[i],PrHit_LHCbID[i]);
         FTCluster cluster;
         cluster.setCluster(ChID_Fraction[i],PrHit_Size[i],ChID_Charge[i],ChID_SipmID[i],ChID_SipmCell[i],ChID_Module[i],ChID_Layer[i],ChID_Mat[i],ChID_Quarter[i]);
         hit.SetCluster(cluster);
         if(!hit.isX()) hit.setCoord( -(hit.x(0.)-xProj.x(hit.z(0.)))/(hit.dxDy() * hit.z(0)));
         if(!hit.isX())bestXZ.addHit(hit);
      }
      bool ok =  FitSimultaneouslyXY(bestXZ, UVSegment,debug);

      int nT1U = bestXZ.nHitsInPlane( 1 );
      int nT1V = bestXZ.nHitsInPlane( 2 );
      int nT2U = bestXZ.nHitsInPlane( 5 );
      int nT2V = bestXZ.nHitsInPlane( 6 );
      int nT3U = bestXZ.nHitsInPlane( 9 );
      int nT3V = bestXZ.nHitsInPlane( 10 );

      std::vector<PatHit> hitT1U = bestXZ.hitsinPlane(1);
      std::vector<PatHit> hitT1V = bestXZ.hitsinPlane(2);
      std::vector<PatHit> hitT2U = bestXZ.hitsinPlane(5);
      std::vector<PatHit> hitT2V = bestXZ.hitsinPlane(6);
      std::vector<PatHit> hitT3U = bestXZ.hitsinPlane(9);
      std::vector<PatHit> hitT3V = bestXZ.hitsinPlane(10);

      std::vector<PatHit> HitsUVT1 = hitT1U;
      HitsUVT1.insert(HitsUVT1.end(),hitT1V.begin(),hitT1V.end());

      std::vector<PatHit> HitsUVT2 = hitT2U;
      HitsUVT2.insert(HitsUVT2.end(),hitT2V.begin(),hitT2V.end());

      std::vector<PatHit> HitsUVT3 = hitT3U;
      HitsUVT3.insert(HitsUVT3.end(),hitT3V.begin(),hitT3V.end());

      //if(nT1U>1){std::cout<<"Has duplicate in T1U"<<std::endl; bestXZ.PrintHits(); bestXZ.removeAllHitInPlane(1); std::cout<<"removed duplicates in plane 1"<<std::endl; bestXZ.PrintHits();}
      PrSeedTracks bestXZplusUV;
      if(nT1U+nT1V ==0) {PrelimSele_Case1 = false; std::cout<<"NO Hits in U + V stereo T1!"<<std::endl;}
      if(nT2U+nT2V ==0) {PrelimSele_Case1 = false; std::cout<<"NO Hits in U + V stereo T2!"<<std::endl;}
      if(nT3U+nT3V ==0) {PrelimSele_Case1 = false; std::cout<<"NO Hits in U + V stereo T3!"<<std::endl;}
      bool FindBest = false;
      PrSeedTrack beforeFindBest = bestXZ;

      if( (nT1U+nT1V)*(nT2U+nT2V)*(nT3U+nT3V) != 0){ //reconstructibility
         if(nT1U>1){bestXZ.removeAllHitInPlane(1); FindBest= true;}
         if(nT1V>1){bestXZ.removeAllHitInPlane(2); FindBest= true;}
         if(nT2U>1){bestXZ.removeAllHitInPlane(5); FindBest = true;}
         if(nT2V>1){bestXZ.removeAllHitInPlane(6); FindBest = true;}
         if(nT3U>1){bestXZ.removeAllHitInPlane(9); FindBest = true;}
         if(nT3V>1){bestXZ.removeAllHitInPlane(10); FindBest = true;}
         PrSeedTracks xzplussingleUVcombinations;
         bestXZ.sortbyZ();
      }
      double minChi2 = bestXZ.Chi2()+1e3;
      PrSeedTrack BestCombination = PrSeedTrack(0, m_zReference);
      if(FindBest){ // if you enter here you have removed in one of the layer all the hits
         for(int i = 0 ; i< HitsUVT1.size();i++){
            if(bestXZ.nHitsInPlane(HitsUVT1[i].planeCode()) == 0){
               bestXZ.addHit(HitsUVT1[i]);
            }
            if(bestXZ.nHitsInPlane(HitsUVT1[i].planeCode()) >0 ){
               bestXZ.substituteHitinSameLayer( HitsUVT1[i]); //will be dummy if already present 1 single hit, if not it really makes the susbtitution
            }
            for(int j = 0 ; j<HitsUVT2.size();j++){
               if(bestXZ.nHitsInPlane(HitsUVT2[j].planeCode())==0){
                  bestXZ.addHit(HitsUVT2[j]);
               }
               if(bestXZ.nHitsInPlane(HitsUVT2[j].planeCode())>0){
                  bestXZ.substituteHitinSameLayer( HitsUVT2[j]);
               }
               for( int k = 0; k< HitsUVT3.size(); k++){
                  if(bestXZ.nHitsInPlane(HitsUVT3[k].planeCode())==0){
                     bestXZ.addHit(HitsUVT3[k]);
                  }
                  if(bestXZ.nHitsInPlane(HitsUVT3[k].planeCode())>0){
                     bestXZ.substituteHitinSameLayer(HitsUVT3[k]);

                     bestXZ.sortbyZ();
                     //bestXZ.PrintTrack();
                     bool Ok =  FitSimultaneouslyXY(bestXZ, UVSegment, false);
                     if(bestXZ.Chi2()< minChi2){
                        minChi2 = bestXZ.Chi2();
                        BestCombination = bestXZ;
                     }
                  }
               }
            }
         }
         if(debug){
         //need to find best combination of the UV layers hit
            if(nT1U>1)
            {  std::cout<<"Had duplicate in T1U"<<std::endl;}
            if(nT1V>1)
            {  std::cout<<"Had duplicate in T1V"<<std::endl;}
            if(nT1V>1)
            {  std::cout<<"Had duplicate in T1V"<<std::endl;}
            if(nT2U>1)
            {  std::cout<<"Had duplicate in T2U"<<std::endl;}
            if(nT2V>1)
            {  std::cout<<"Had duplicate in T2V"<<std::endl;}
            if(nT3U>1)
            {  std::cout<<"Had duplicate in T3U"<<std::endl;}
            if(nT3V>1)
            {  std::cout<<"Had duplicate in T3V"<<std::endl;}
         }
         // std::cout<<"\n \n ===================Before Finding the best track==================="<<std::endl;
         // beforeFindBest.PrintHits();
         // beforeFindBest.PrintTrack();
         // std::cout<<"\n \n ===========================Best Combination=====================\n\n"<<std::endl;
         // BestCombination.PrintHits();
         // BestCombination.PrintTrack();
         // std::cout<<"\n \n"<<std::endl;
         // for(layer1){
         bestXZ=BestCombination;
      }
      PrSeedTrack BestUVSegment = PrSeedTrack(bestXZ.zone(), m_zReference);
      for(int i =0; i<bestXZ.hits().size(); i++){
        //add only UVHits to bestXZ
         if( !bestXZ.hits()[i].isX()){
           BestUVSegment.addHit(bestXZ.hits()[i]);
         }
      }
      if(debug){
         std::cout<<"Best UV Segment Case 1"<<std::endl;
         BestUVSegment.PrintHits();
      }
      BestUVSegment.setParameters(ax_FitHitCase1,bx_FitHitCase1,cx_FitHitCase1,0,0);
      FitLine(BestUVSegment,debug,1);
      NUV_Case1 = (int)BestUVSegment.hits().size();
      //std::cout<<"ax = "<<bestXZ.ax();
      ax_fullFitCase1 = (Double_t )bestXZ.ax();
      bx_fullFitCase1= (Double_t )bestXZ.bx();
      cx_fullFitCase1 =  (Double_t )bestXZ.cx();
      ay_fullFitCase1 =  (Double_t )bestXZ.ay();
      by_fullFitCase1 =  (Double_t )bestXZ.by();
      Chi2_fullFitCase1 =  (Double_t )bestXZ.Chi2();
      MaxChi2_fullFitCase1 =  (Double_t )bestXZ.maxChi2Hit();
      NHits_fullFitCase1= (int)bestXZ.hits().size();
      y0_fullFitCase1 =  (Double_t )bestXZ.y(0.);
      if(debug) std::cout<<"Number of XZ Projections after erase"<<xProjections.size()<<std::endl;



   }
   return PrelimSele_Case1;
}

bool TrackStudy::XZStudyCase2(PrSeedTrack& xProj, PrSeedTrack& UVSegment ,Bool_t debug){
  NUV_Case2 = -10;
  ax_Case2_Seed1 = -99.;
  bx_Case2_Seed1 = -99.;
  cx_Case2_Seed1 = -99.;
  ax_Case2_Seed2 = -99.;
  bx_Case2_Seed2 = -99.;
  cx_Case2_Seed2 = -99.;
  Chi2X_Case2 = -99.;
  NX_Case2 = -1;
  tx_inf_Case2 = -99.;
  x0_Case2 = -99.;
  x_inf_Case2 = -99.;
  Del_x_inf_Case2 = -99.;
  tx_pickedcombination_Case2 = -99.;
  DelxProjectedSeed2_Case2 = -99.;
  DelxProjectedSeed1_Case2 = -99.;
  MaxChi2Hit_Case2 = -99;


  ax_FitHitCase2 = -99.;
  bx_FitHitCase2 = -99.;
  cx_FitHitCase2 = -99.;

  ax_fullFitCase2 = -9999.;
  bx_fullFitCase2= -99.;
  cx_fullFitCase2 = -99.;
  ay_fullFitCase2 = -9999.;
  by_fullFitCase2 = -99.;
  Chi2_fullFitCase2 =  -99.;
  MaxChi2_fullFitCase2 = -99.;
  NHits_fullFitCase2 = -1;
  y0_fullFitCase2 = -9999.;



  std::vector<PatHit> Case2First ;
   Case2First.reserve(10);
   std::vector<PatHit> Case2Last  ;
   Case2Last.reserve(10);
   std::vector<PatHit> ParabolaSeedHits1;
   ParabolaSeedHits1.reserve(10);
   std::vector<PatHit> ParabolaSeedHits2;
   ParabolaSeedHits2.reserve(10);
   std::vector<PatHit> RemainingCase2Seed1;
   RemainingCase2Seed1.reserve(10);
   std::vector<PatHit> RemainingCase2Seed2;
   RemainingCase2Seed2.reserve(10);
   std::vector<int> remainingCase2Seed1N(3,0);
   std::vector<int> remainingCase2Seed2N(3,0);

   for(int i =0; i<xProj.hits().size();i++){
      PatHit hit = xProj.hits()[i];
      if(hit.planeCode()==0){
         Case2First.push_back(hit);
      }
      if(hit.planeCode()==8){
         Case2Last.push_back(hit);
      }
      if(hit.planeCode()==4){
         ParabolaSeedHits1.push_back(hit);
      }
      if(hit.planeCode()==7){
         ParabolaSeedHits2.push_back(hit);
      }
      if( hit.planeCode() == 3 || hit.planeCode() ==7 || hit.planeCode() == 11) {
         if(debug) std::cout<<"Plane 3,4,8 Hits loaded"<<std::endl;

         if(hit.planeCode() == 3) remainingCase2Seed1N[0]++;
         if(hit.planeCode() == 7) remainingCase2Seed1N[1]++;
         if(hit.planeCode() == 11) remainingCase2Seed1N[2]++;
         RemainingCase2Seed1.push_back(hit);
      }
      if(  hit.planeCode() == 3 || hit.planeCode() ==4 || hit.planeCode() == 11){
         if(debug) std::cout<<"Planes"<<setw(20)<<hit.planeCode()<<"\t Hits done"<<std::endl;
         if(hit.planeCode() ==3) remainingCase2Seed2N[0]++;
         if(hit.planeCode() ==4) remainingCase2Seed2N[1]++;
         if(hit.planeCode() ==11) remainingCase2Seed2N[2]++;
         RemainingCase2Seed2.push_back(hit);
      }
   }
   bool PrelimSele_Case2 = false;
   if(Case2First.size()!=0 && Case2Last.size()!=0 && (ParabolaSeedHits1.size()!=0 || ParabolaSeedHits2.size()!=0)){
      PrelimSele_Case2 = true;
   }


   if(PrelimSele_Case2){
      m_physicalCase2Intrinsic++;
      PrSeedTracks xProjections;
      for(int i = 0; i<Case2First.size();i++)
      {
         if(debug) std::cout<<"First Hit "<<std::endl;
         //Case2First[i].PrintHit();
         if(debug) std::cout<<"Loop first Layer"<<std::endl;
         for(int j = 0;j<Case2Last.size();j++)
         {
            if(debug) std::cout<<"Second Hit"<<std::endl;
            //Case2Last[j].PrintHit();
            if(debug) std::cout<<"Loop Last Layer"<<std::endl;
            // double deltaz = Case2Last[j].z()-Case2First[i].z();
            // double tx_pickedcombination = ( Case2Last[j].x(0.)-Case2First[i].x(0.))  /( Case2Last[j].z()-Case2First[i].z());
            // double xFirstProjected = Case2First[i].x(0.) + Case2First[i].x(0.)/Case2First[i].z(0.)*( deltaz);
            // double txinf = Case2First[i].x(0.)/Case2First[i].z(0.);
            //PrSeedTrack xProjection1(0, m_zReference);
            //PatHit fHit = Case2First[i];

            //PatHit lHit = Case2Last[j];
            if(debug) std::cout<<"Loop Last "<<std::endl;
            if(debug) std::cout<<"parabolaSeedHits1 sizee   "<<ParabolaSeedHits1.size()<<std::endl;
            //FInd Given the parabolaSeedHit1 the best combination
            for(int k=0; k< ParabolaSeedHits1.size();k ++)
            {
               if(debug) std::cout<<"Parabola Seed Hit 1"<<std::endl;
               //ParabolaSeedHits1[k].PrintHit();
               if(debug) std::cout<<"Loop ParabolaSeed 1"<<std::endl;
               //PrSeedTrack1 xProjection();
               PrSeedTrack xProjection = PrSeedTrack(0,m_zReference);
               xProjection.addHit(ParabolaSeedHits1[k]);
               xProjection.addHit(Case2First[i]);
               xProjection.addHit(Case2Last[j]);
               //xProjection.PrintTrack();
               double a=0.;
               double b=0.;
               double c=0.;
               solveParabola2(Case2First[i],ParabolaSeedHits1[k],Case2Last[j],a,b,c );
               ax_Case2_Seed1 = a;
               bx_Case2_Seed1 = b;
               cx_Case2_Seed1 = c;
               if(debug) std::cout<<"a"<<setw(20)<<a<<std::endl;
               if(debug) std::cout<<"b"<<setw(20)<<b<<std::endl;
               if(debug) std::cout<<"c"<<setw(20)<<c<<std::endl;

               //std::vector<PatHit> bestRemaining(3);
               int best1_kk;
               int best2_kk;
               int best3_kk;

               std::vector<double> minVal(3);
               minVal[0]=1000; minVal[1]=1000; minVal[2]=1000;
               for(int kk =0;kk<remainingCase2Seed1N[0] ; kk++){
                  if(debug) std::cout<<"Remaining Seed 1 -1 "<<std::endl;
                  //RemainingCase2Seed1[kk].PrintHit();
                  double dz = RemainingCase2Seed1[kk].z() -m_zReference;
                  double xAtZ = a*dz*dz*(1.+m_dRatio0*dz) + b*dz +c;
                  double delta = std::fabs(RemainingCase2Seed1[kk].x(0)-xAtZ);
                  if(debug) std::cout<<"delta Remaining 1 -1 "<<setw(20)<<delta<<std::endl;
                  if(std::fabs(RemainingCase2Seed1[kk].x(0)-xAtZ)<minVal[0]){
                     minVal[0]=std::fabs(delta);
                     best1_kk=kk;
                     //hit1_1 = RemainingCase2Seed1[kk];
                  }
               }
               if(minVal[0]<999){
                  xProjection.addHit(RemainingCase2Seed1[best1_kk]);
               }
               for(int kk =remainingCase2Seed1N[0];kk< ( remainingCase2Seed1N[0]+remainingCase2Seed1N[1]) ; kk++){
                  double dz = RemainingCase2Seed1[kk].z() - m_zReference;
                  double xAtZ = a*dz*dz*(1.+m_dRatio0*dz) + b*dz +c;
                  double delta = std::fabs(RemainingCase2Seed1[kk].x(0)-xAtZ);
                  if(debug) std::cout<<"delta Remaining 1 -2 "<<setw(20)<<delta<<std::endl;

                  if(std::fabs(delta)<minVal[1]){
                     minVal[1]=std::fabs(delta);
                     best2_kk=kk;
                     //hit2_1 = RemainingCase2Seed1[kk];
                  }
               }
               if(minVal[1]<999){
                  xProjection.addHit(RemainingCase2Seed1[best2_kk]);
               }
               for( int kk= (remainingCase2Seed1N[0]+remainingCase2Seed1N[1]);kk<( remainingCase2Seed1N[0]+remainingCase2Seed1N[1]+remainingCase2Seed1N[2] ); kk++){
                  double dz = RemainingCase2Seed1[kk].z() -m_zReference;
                  double xAtZ = a*dz*dz *(1.+m_dRatio0*dz)+ b*dz +c;
                  double delta = std::fabs(RemainingCase2Seed1[kk].x(0)-xAtZ);
                  if(debug) std::cout<<"delta Remaining 1 -3"<<setw(20)<<delta<<std::endl;

                  if(std::fabs(RemainingCase2Seed1[kk].x(0)-xAtZ)<minVal[2]){
                     minVal[2]=std::fabs(delta);
                     best3_kk=kk;
                     //hit3_1 = RemainingCase2Seed1[kk];
                  }
               }
               if(minVal[2]<999){
                  xProjection.addHit(RemainingCase2Seed1[best3_kk]);
               }
               if(debug) std::cout<<"Pushing Back track Case Seed 1"<<std::endl;
               //xProjection.PrintHits();
               if(debug) std::cout<<"Case Seed 1"<<std::endl;

               xProjections.push_back(xProjection);
               //T2
            }
            for(int k=0; k< ParabolaSeedHits2.size();k ++)
            {
               if(debug) std::cout<<"Parabola Seed Hit 1"<<std::endl;
               //ParabolaSeedHits2[k].PrintHit();
               //PrSeedTrack1 xProjection();
               PrSeedTrack xProjection = PrSeedTrack(0,m_zReference);
               xProjection.addHit(Case2First[i]);
               xProjection.addHit(ParabolaSeedHits2[k]);
               xProjection.addHit(Case2Last[j]);

               double a=0;
               double b=0;
               double c=0;
               solveParabola2(Case2First[i],ParabolaSeedHits2[k],Case2Last[j],a,b,c );
               ax_Case2_Seed2 = a;
               bx_Case2_Seed2 = b;
               cx_Case2_Seed2 = c;
               if(debug) std::cout<<"a "<<setw(20)<<a<<std::endl;
               if(debug) std::cout<<"b "<<setw(20)<<b<<std::endl;
               if(debug) std::cout<<"c "<<setw(20)<<c<<std::endl;
               int best1_kk;
               int best2_kk;
               int best3_kk;
               // PatHit hit1_1;
               // PatHit hit2_1;
               // PatHit hit3_1;
               std::vector<double> minVal(3);
               minVal[0]=1000; minVal[1]=1000; minVal[2]=1000;
               for(int kk =0;kk<remainingCase2Seed2N[0] ; kk++){
                  double dz = RemainingCase2Seed2[kk].z() -m_zReference;
                  //double m_dRatio0 = -0.000246;
                  double xAtZ = a*dz*dz*(1.+m_dRatio0*dz) + b*dz +c;
                  double delta = std::fabs(RemainingCase2Seed2[kk].x(0)-xAtZ);
                  if(debug) std::cout<<"delta 2-1 "<<delta<<std::endl;
                  if(std::fabs(RemainingCase2Seed2[kk].x(0)-xAtZ)<minVal[0]){
                     minVal[0]=std::fabs(delta);
                     best1_kk = kk;
                     //hit1_1 = RemainingCase2Seed2[kk];
                  }
               }
               if(minVal[0]<999){
                  xProjection.addHit(RemainingCase2Seed2[best1_kk]);
               }
               for(int kk =remainingCase2Seed2N[0];kk<(remainingCase2Seed2N[0]+remainingCase2Seed2N[1]) ; kk++){
                  double dz = RemainingCase2Seed1[kk].z() -m_zReference;
                  double xAtZ = a*dz*dz*(1.+m_dRatio0*dz) + b*dz +c;
                  double delta = std::fabs(RemainingCase2Seed2[kk].x(0)-xAtZ);
                  if(debug) std::cout<<"Delta 2-2 "<<setw(20)<<delta<<std::endl;
                  if(std::fabs(delta)<minVal[1]){
                     minVal[1]=std::fabs(delta);
                     best2_kk=kk;
                     // hit2_1 = RemainingCase2Seed2[kk];
                  }
               }
               if(minVal[1]<999){
                  xProjection.addHit(RemainingCase2Seed2[best2_kk]);
               }
               for(int  kk=(remainingCase2Seed2N[0]+remainingCase2Seed2N[1]) ;kk<(remainingCase2Seed2N[0]+remainingCase2Seed2N[1]+remainingCase2Seed2N[2]); kk++){
                  double dz = RemainingCase2Seed2[kk].z() -m_zReference;
                  // double m_dRatio0 = -0.000246;
                  double xAtZ = a*dz*dz *(1.+m_dRatio0*dz) + b*dz +c;

                  double delta = std::fabs(RemainingCase2Seed2[kk].x(0)-xAtZ);
                  if(debug) std::cout<<"Delta 2-3"<<setw(20)<<delta<<std::endl;

                  if(std::fabs(RemainingCase2Seed2[kk].x(0)-xAtZ)<minVal[2]){
                     minVal[2]=std::fabs(delta);
                     best3_kk=kk;
                     //hit3_1 = RemainingCase2Seed2[kk];
                  }
               }
               if(minVal[2]<999){
                  xProjection.addHit(RemainingCase2Seed2[best3_kk]);
               }
               if(debug) std::cout<<"Pushing Back track Case Seed 1"<<std::endl;
               //xProjection.PrintHits();
               if(debug) std::cout<<"Case Seed 1"<<std::endl;

               xProjections.push_back(xProjection);

            }

         }
      }
      for(int i =0 ; i< xProjections.size() ; i++){
         PrSeedTrack track= xProjections[i];
         track.sortbyZ();
         // PrSeedTrack track = xProjections[0];
         if(debug) std::cout<<"Number of XZ Projections before erase"<<xProjections.size()<<std::endl;
         if(xProjections.size()>1)
         {
            xProjections.erase(std::unique(xProjections.begin(),xProjections.end(),[](PrSeedTrack a, PrSeedTrack b)->bool{
               if(a.hits().size() != b.hits().size()) return false;
               a.sortbyZ();
               b.sortbyZ();
               bool flag = true;
               for(int i =0;i<a.hits().size();i++){
                  if( (a.hits()[i].x() != b.hits()[i].x()) && a.hits()[i].z()!=a.hits()[i].z()) {  flag = false; break;}
               }
               return flag;
            }),xProjections.end());
         }
      }
      if(xProjections.size()==0){std::cout<<"I didn't find any xProjection, error case 1"<<std::endl; xProj.PrintHits();}
      if(xProjections.size()>1) m_physicalCase2IntrinsicWithduplicates++;
      if(xProjections.size()>1 && xProjections[0].hits().size()>=4) m_physicalCase2IntrinsicWithduplicatesMore4++;
      if(xProjections.size()>1 && xProjections[0].hits().size()>=5) m_physicalCase2IntrinsicWithduplicatesMore5++;
      //if(xProjections.size()>0){
      if(xProjections[0].hits().size()>=4) m_physicalCase2IntrinsicMore4++;
      if(xProjections[0].hits().size()>=5) m_physicalCase2IntrinsicMore5++;
      //}


      if(xProjections.size()>1)
      {
         std::vector<std::pair<PrSeedTrack,double>> trackChi2;
         for(int i = 0;i<xProjections.size();i++){
            fitXProjection(xProjections[i],debug);
         }
         std::sort(xProjections.begin(),xProjections.end(),[](PrSeedTrack val1, PrSeedTrack val2)->bool{return val1.chi2XnDoF()<val2.chi2XnDoF();});
      }
      PrSeedTrack bestXZ = xProjections[0];
      bool OK = false;
      OK = fitXProjection(bestXZ,debug);
      if(!OK){
         bestXZ.PrintHits();
         std::cout<<"Case 0:  Chi2 ="<<setw(20)<<bestXZ.chi2X()<<std::endl;
         std::cout<<"         NDOF ="<<setw(20)<<bestXZ.ndofX()<<std::endl;
         std::cout<<"     Chi2NdoF ="<<setw(20)<<bestXZ.chi2XnDoF()<<std::endl;
      }

      ax_FitHitCase2 = bestXZ.ax();
      bx_FitHitCase2 = bestXZ.bx();
      cx_FitHitCase2 = bestXZ.cx();

      //Case 0 first Plane = 0; Last = 11;
      PatHit first = bestXZ.hitInPlane(0);
      // first.PrintHit();
      PatHit last = bestXZ.hitInPlane(8);
      //last.PrintHit();
      Track_Zone =first.zone();
      if(first.zone() != last.zone()) {
        std::cout<<"Killed Track because Up->Down or vice versa"<<std::endl;
        Track_Zone = -99;
        PrelimSele_Case2 = false;
      }
      tx_inf_Case2 = first.x(0.)/first.z(0.);
      NX_Case2= bestXZ.hits().size();
      if(NX_Case2 <4 ){Case2Accepted = false;}
      Chi2X_Case2 =bestXZ.chi2X();
      MaxChi2Hit_Case2 = bestXZ.maxChi2Hit();
      tx_inf_Case2 =first.x(0.)/first.z(0.);
      tx_pickedcombination_Case2 = (last.x(0.) - first.x(0.))/(last.z(0.)-first.z(0.));
      x_inf_Case2 = tx_inf_Case2 * ( last.z(0.) - first.z(0.)) +first.x(0.);
      Del_x_inf_Case2 =  last.x(0.)-x_inf_Case2;
      tx_pickedcombination_Case2 = (last.x(0.) - first.x(0.))/(last.z(0.)-first.z(0.));
      x0_Case2 = first.x(0) - tx_pickedcombination_Case2 *(first.z(0.) );
      if(bestXZ.hasHitInPlane(4)){
         PatHit parSeed1 = bestXZ.hitInPlane(4);
         Double_t xProectedSeed1 = first.x(0.)+tx_pickedcombination_Case2*(parSeed1.z()-first.z(0.));
         DelxProjectedSeed1_Case2 = parSeed1.x(0.)-xProectedSeed1;
      }
      if(bestXZ.hasHitInPlane(7)){
         PatHit parSeed2 =  bestXZ.hitInPlane(7);
         Double_t xProectedSeed2 = first.x(0.)+tx_pickedcombination_Case2*(parSeed2.z()-first.z(0.));
         DelxProjectedSeed2_Case2 = parSeed2.x(0.)-xProectedSeed2;
      }
      
      int Case = 2;
      Remaining(bestXZ,debug,Case);
      
      for(int i = 0; i<PrHit; i++){
         PatHit hit = PatHit();
         hit.setHit(PrHit_Xat0[i],PrHit_Zat0[i],PrHit_dxDy[i],PrHit_dzDy[i],std::sqrt(PrHit_w2[i]),PrHit_yMin[i],PrHit_yMax[i],PrHit_zone[i],PrHit_planeCode[i],PrHit_isX[i],PrHit_LHCbID[i]);
         FTCluster cluster;
         cluster.setCluster(ChID_Fraction[i],PrHit_Size[i],ChID_Charge[i],ChID_SipmID[i],ChID_SipmCell[i],ChID_Module[i],ChID_Layer[i],ChID_Mat[i],ChID_Quarter[i]);
         hit.SetCluster(cluster);
         if(!hit.isX()) hit.setCoord(- (hit.x(0.)-xProj.x(hit.z(0.)))/(hit.dxDy() * hit.z(0)));
         if(!hit.isX())bestXZ.addHit(hit);

      }
      bool ok =  FitSimultaneouslyXY(bestXZ, UVSegment,debug);

      int nT1U = bestXZ.nHitsInPlane( 1 );
      int nT1V = bestXZ.nHitsInPlane( 2 );
      int nT2U = bestXZ.nHitsInPlane( 5 );
      int nT2V = bestXZ.nHitsInPlane( 6 );
      int nT3U = bestXZ.nHitsInPlane( 9 );
      int nT3V = bestXZ.nHitsInPlane( 10 );

      std::vector<PatHit> hitT1U = bestXZ.hitsinPlane(1);
      std::vector<PatHit> hitT1V = bestXZ.hitsinPlane(2);
      std::vector<PatHit> hitT2U = bestXZ.hitsinPlane(5);
      std::vector<PatHit> hitT2V = bestXZ.hitsinPlane(6);
      std::vector<PatHit> hitT3U = bestXZ.hitsinPlane(9);
      std::vector<PatHit> hitT3V = bestXZ.hitsinPlane(10);

      std::vector<PatHit> HitsUVT1 = hitT1U;
      HitsUVT1.insert(HitsUVT1.end(),hitT1V.begin(),hitT1V.end());

      std::vector<PatHit> HitsUVT2 = hitT2U;
      HitsUVT2.insert(HitsUVT2.end(),hitT2V.begin(),hitT2V.end());

      std::vector<PatHit> HitsUVT3 = hitT3U;
      HitsUVT3.insert(HitsUVT3.end(),hitT3V.begin(),hitT3V.end());

      //if(nT1U>1){std::cout<<"Has duplicate in T1U"<<std::endl; bestXZ.PrintHits(); bestXZ.removeAllHitInPlane(1); std::cout<<"removed duplicates in plane 1"<<std::endl; bestXZ.PrintHits();}
      PrSeedTracks bestXZplusUV;
      if(nT1U+nT1V ==0) {PrelimSele_Case2 = false; std::cout<<"NO Hits in U + V stereo T1!"<<std::endl;}
      if(nT2U+nT2V ==0) {PrelimSele_Case2 = false; std::cout<<"NO Hits in U + V stereo T2!"<<std::endl;}
      if(nT3U+nT3V ==0) {PrelimSele_Case2 = false; std::cout<<"NO Hits in U + V stereo T3!"<<std::endl;}
      bool FindBest = false;
      PrSeedTrack beforeFindBest = bestXZ;

      if( (nT1U+nT1V)*(nT2U+nT2V)*(nT3U+nT3V) != 0){ //reconstructibility
         if(nT1U>1){bestXZ.removeAllHitInPlane(1); FindBest= true;}
         if(nT1V>1){bestXZ.removeAllHitInPlane(2); FindBest= true;}
         if(nT2U>1){bestXZ.removeAllHitInPlane(5); FindBest = true;}
         if(nT2V>1){bestXZ.removeAllHitInPlane(6); FindBest = true;}
         if(nT3U>1){bestXZ.removeAllHitInPlane(9); FindBest = true;}
         if(nT3V>1){bestXZ.removeAllHitInPlane(10); FindBest = true;}
         PrSeedTracks xzplussingleUVcombinations;
         bestXZ.sortbyZ();
      }
      double minChi2 = bestXZ.Chi2()+1e3;
      PrSeedTrack BestCombination = PrSeedTrack(0, m_zReference);
      if(FindBest){ // if you enter here you have removed in one of the layer all the hits
         for(int i = 0 ; i< HitsUVT1.size();i++){
            if(bestXZ.nHitsInPlane(HitsUVT1[i].planeCode()) == 0){
               bestXZ.addHit(HitsUVT1[i]);
            }
            if(bestXZ.nHitsInPlane(HitsUVT1[i].planeCode()) >0 ){
               bestXZ.substituteHitinSameLayer( HitsUVT1[i]); //will be dummy if already present 1 single hit, if not it really makes the susbtitution
            }
            for(int j = 0 ; j<HitsUVT2.size();j++){
               if(bestXZ.nHitsInPlane(HitsUVT2[j].planeCode())==0){
                  bestXZ.addHit(HitsUVT2[j]);
               }
               if(bestXZ.nHitsInPlane(HitsUVT2[j].planeCode())>0){
                  bestXZ.substituteHitinSameLayer( HitsUVT2[j]);
               }
               for( int k = 0; k< HitsUVT3.size(); k++){
                  if(bestXZ.nHitsInPlane(HitsUVT3[k].planeCode())==0){
                     bestXZ.addHit(HitsUVT3[k]);
                  }
                  if(bestXZ.nHitsInPlane(HitsUVT3[k].planeCode())>0){
                     bestXZ.substituteHitinSameLayer(HitsUVT3[k]);

                     bestXZ.sortbyZ();
                     //bestXZ.PrintTrack();
                     bool Ok =  FitSimultaneouslyXY(bestXZ, UVSegment, false);
                     if(bestXZ.Chi2()< minChi2){
                        minChi2 = bestXZ.Chi2();
                        BestCombination = bestXZ;
                        //BestCombination
                     }
                  }
               }
            }
         }
         //need to find best combination of the UV layers hit
         if(debug){
            if(nT1U>1)
            {  std::cout<<"Had duplicate in T1U"<<std::endl;}
            if(nT1V>1)
            {  std::cout<<"Had duplicate in T1V"<<std::endl;}
            if(nT1V>1)
            {  std::cout<<"Had duplicate in T1V"<<std::endl;}
            if(nT2U>1)
            {  std::cout<<"Had duplicate in T2U"<<std::endl;}
            if(nT2V>1)
            {  std::cout<<"Had duplicate in T2V"<<std::endl;}
            if(nT3U>1)
            {  std::cout<<"Had duplicate in T3U"<<std::endl;}
            if(nT3V>1)
            {  std::cout<<"Had duplicate in T3V"<<std::endl;}
         }
         // std::cout<<"\n \n ===================Before Finding the best track==================="<<std::endl;
         // beforeFindBest.PrintHits();
         // beforeFindBest.PrintTrack();
         // std::cout<<"\n \n ===========================Best Combination=====================\n\n"<<std::endl;
         // BestCombination.PrintHits();
         // BestCombination.PrintTrack();
         // std::cout<<"\n \n"<<std::endl;
         // for(layer1){
         bestXZ=BestCombination;
      }
      PrSeedTrack BestUVSegment = PrSeedTrack(bestXZ.zone(), m_zReference);
      for(int i =0; i<bestXZ.hits().size(); i++){
         if( !bestXZ.hits()[i].isX()){
            BestUVSegment.addHit(bestXZ.hits()[i]);
         }
      }
      if(debug){
         std::cout<<"Best UV Segment Case 2"<<std::endl;
         BestUVSegment.PrintHits();
      }
      BestUVSegment.setParameters( ax_FitHitCase2, bx_FitHitCase2, cx_FitHitCase2,0,0);
      FitLine(BestUVSegment,debug,2);
      NUV_Case2 = (int) BestUVSegment.hits().size();
      //std::cout<<"ax = "<<bestXZ.ax();
      ax_fullFitCase2 = (Double_t )bestXZ.ax();
      bx_fullFitCase2= (Double_t )bestXZ.bx();
      cx_fullFitCase2 =  (Double_t )bestXZ.cx();
      ay_fullFitCase2 =  (Double_t )bestXZ.ay();
      by_fullFitCase2 =  (Double_t )bestXZ.by();
      Chi2_fullFitCase2 =  (Double_t )bestXZ.Chi2();
      MaxChi2_fullFitCase2 =  (Double_t )bestXZ.maxChi2Hit();
      NHits_fullFitCase2= (int)bestXZ.hits().size();
      y0_fullFitCase2 =  (Double_t )bestXZ.y(0.);
      if(debug) std::cout<<"Number of XZ Projections after erase"<<xProjections.size()<<std::endl;


   }

   return PrelimSele_Case2;

}

void TrackStudy::StereoSearch(PrSeedTrack xProj,PrSeedTrack UVSegment,Bool_t debug){
   fitXProjection(xProj,debug);
   std::vector<PatHit> HitsUVT1;
   std::vector<PatHit> HitsUVT2;
   std::vector<PatHit> HitsUVT3;
   std::vector<PatHit> AllUV;
   for(int i = 0 ; i< UVSegment.hits().size();i++) {
      PatHit hit = UVSegment.hits()[i];
      hit.setCoord(- (hit.x(0.)-xProj.x(hit.z(0.)))/(hit.dxDy() * hit.z(0)));
      AllUV.push_back(hit);
      if(!hit.isX()){
         if(hit.planeCode() <4){
            HitsUVT1.push_back(hit);
         }
         if(hit.planeCode() > 3 && hit.planeCode() < 8){
            HitsUVT2.push_back(hit);
         }
         if(hit.planeCode() >8){
            HitsUVT3.push_back(hit);
         }
      }
   }

   std::vector<Double_t> DeltaYCombT1;
   std::vector<Double_t> DeltaXCombT1;
   int nT1U = std::count_if( HitsUVT1.begin(), HitsUVT1.end(), [](PatHit hit){return hit.planeCode()==1 ;});
   int nT1V = std::count_if( HitsUVT1.begin(), HitsUVT1.end(), [](PatHit hit){return hit.planeCode()==2 ;});
   int nT2U = std::count_if( HitsUVT2.begin(), HitsUVT2.end(), [](PatHit hit){return hit.planeCode()==5 ;});
   int nT2V = std::count_if( HitsUVT2.begin(), HitsUVT2.end(), [](PatHit hit){return hit.planeCode()==6 ;});
   int nT3U = std::count_if( HitsUVT3.begin(), HitsUVT3.end(), [](PatHit hit){return hit.planeCode()==9 ;});
   int nT3V = std::count_if( HitsUVT3.begin(), HitsUVT3.end(), [](PatHit hit){return hit.planeCode()==10 ;});
   std::pair<PatHit ,PatHit> bestHitPair;
   if(nT1U > 0 && nT1V >0)
   {
      Double_t BestDeltaCom = 100.;
      if(debug) std::cout<<"In T1 i found possible pairs from a track which has part =     "<<xProj.zone()<<std::endl;
      for(Int_t i = 0; i< nT1U;i++ )
      {
         Double_t yT1U = (HitsUVT1[i].x(0.) - xProj.x(HitsUVT1[i].z(0.)))/HitsUVT1[i].dxDy();
         Double_t yT1VProj = yT1U + yT1U/( HitsUVT1[i].z(0.) - HitsUVT1[nT1U].z(0.));
         for(Int_t j = nT1U; j< nT1U+nT1V; j++)
         {
            //nCombT1++;
            Double_t yT1V = (HitsUVT1[j].x(0.)-xProj.x(HitsUVT1[j].z(0.)))/HitsUVT1[j].dxDy();
            Double_t XCombined = (HitsUVT1[i].x(0.)+HitsUVT1[j].x(0.))/2.;
            Double_t XPredicted = xProj.x( (HitsUVT1[i].z(0.)+ HitsUVT1[j].z(0.) )/2.);

            DeltaYT1->Fill(yT1V-yT1VProj);
            DeltaYT1vsyT1->Fill(yT1U,yT1V-yT1VProj);
            if(debug) std::cout<<"U Hit:   PlaneCode"<<setw(20)<<HitsUVT1[i].planeCode()<<setw(20)<<"dxDy"<<HitsUVT1[i].dxDy()<<setw(20)<<"X at 0"<<HitsUVT1[i].x(0.)<<setw(20)<<"Z at 0"<<setw(20)<<HitsUVT1[i].z(0.)<<setw(20)<<"Yexp"<<setw(20)<<yT1U<<std::endl;
            if(debug) std::cout<<"V Hit:   PlaneCode"<<setw(20)<<HitsUVT1[j].planeCode()<<setw(20)<<"dxDy"<<HitsUVT1[j].dxDy()<<setw(20)<<"X at 0"<<HitsUVT1[j].x(0.)<<setw(20)<<"Z at 0"<<setw(20)<<HitsUVT1[j].z(0.)<<setw(20)<<"Yexp"<<setw(20)<<yT1V<<std::endl;

            DeltaXCombinedT1->Fill(XCombined-XPredicted);
            DeltaYCombT1.push_back(yT1V-yT1VProj);
            DeltaXCombT1.push_back(XCombined-XPredicted);
            if(std::abs(XCombined-XPredicted)<BestDeltaCom)
            {
               BestDeltaCom = std::abs(XCombined-XPredicted);
               bestHitPair = std::make_pair(HitsUVT1[i], HitsUVT1[j]);
            }
            //if(std::abs(XCombined-XPredicted + 2*))
         }
         //if(debug) std::cout<<"X Combination"<<setw(20)<<XCombined<<std::endl;
         //if(debug) std::cout<<"X Predicted  "<<setw(20)<<xProj.x( (HitsUVT1[i].z(0.)+ HitsUVT1[j].z(0.) )/2.)<<std::endl;
      }
   }
   DeltaXBestCombCorrSlopeY =-999;
   DeltaXBestCombT1 = -999;
   TyPairT1 = -999.;
   YPairT1 = -999.;
   XPairT1 = -999.;
   YPairT1 = -999.;
   YUT2Projected = -999.;
   YVT2Projected = -999.;
   DeltaUT2Projected = -999.;
   DeltaVT2Projected = -999.;
   if(debug) std::cout<<"\n \n \n"<<setw(20)<<"New Track"<<std::endl;
   //xProj.PrintHits();
   if(nT1U>0 && nT1V>0){
      if(debug) std::cout<<"Sorted Combinations"<<std::endl;
      std::sort(DeltaXCombT1.begin(),DeltaXCombT1.end(),[](Double_t Value1, Double_t Value2)->bool{return std::abs(Value1)<std::abs(Value2);});
      for(Int_t i = 0 ; i<DeltaXCombT1.size();i++)
      {
         if(DeltaXCombT1.size()>1)
         {
            if(debug) std::cout<<"i"<<setw(20)<<"Delta X Best Combination"<<setw(20)<<DeltaXCombT1[i]<<std::endl;
         }
      }
      DeltaXBestCombT1 = DeltaXCombT1[0];
      if(debug) std::cout<<"DeltaBestCombT1   "<<setw(20)<<DeltaXBestCombT1<<std::endl;
      if(debug) std::cout<<"Compute Delta Best"<<setw(20)<<(bestHitPair.first.x(0.)+bestHitPair.second.x(0.))/2 - xProj.x( (bestHitPair.first.z(0.) + bestHitPair.second.z(0.)) /2.)<<std::endl;
      PatHit bestUT1 = bestHitPair.first;
      PatHit bestVT1 = bestHitPair.second;
      Double_t zMid = (bestUT1.z() + bestVT1.z())/2. ;
      Double_t xMidFromHit = (bestUT1.x() + bestVT1.x()) /2. ;//+ xProj.xSlope(bestUT1.z())*(zMid - bestUT1.z());

      DeltaXBestCombCorrSlopeY = xMidFromHit - xProj.x(zMid) + (zMid-bestUT1.z())*bestUT1.dxDy()*TyPairT1;
      Double_t yT1U = (bestUT1.x(0.)-xProj.x(bestUT1.z(0.)))/bestUT1.dxDy();
      Double_t yT1V = (bestVT1.x(0.)-xProj.x(bestVT1.z(0.)))/bestVT1.dxDy();

      Double_t yAvg = (yT1U+yT1V)/2.;
      Double_t TY = (yT1V-yT1U)/(bestVT1.z(0.) - bestUT1.z(0.));
      TyPairT1 = TY;
      DeltaXBestCombCorrSlopeY = xMidFromHit - xProj.x(zMid) + (zMid-bestUT1.z())*bestUT1.dxDy()*TyPairT1;
      Double_t xMidFromHitCorr = (bestUT1.x() + bestVT1.x())/2. + TY/2.*bestUT1.dxDy()*(zMid-bestUT1.z());
      XPairT1 = (bestUT1.x()+bestVT1.x())/2.;
      YPairT1 = yAvg;
      Double_t zAvg = (bestUT1.z(0.)+bestVT1.z(0.))/2.;
      DeltaUT2Projected = -999.;
      if(nT2U>0)
      {
         Double_t bestDeltaUT2Projected = 200.;
         for(Int_t i =0 ; i<nT2U;i++){
            Double_t yUT2Exp = yAvg + TY*(HitsUVT2[i].z(0.) - zAvg);
            Double_t uExp = yUT2Exp*HitsUVT2[i].dxDy() + xProj.x(HitsUVT2[i].z(0.));
            // if(debug) std::cout<<"T2  -- U -- i ="<<setw(20)<<i<<setw(20)<<"Delta U T2 expected best combination in T1"<<setw(20)<<uExp - HitsUVT2[i].x(0.)<<std::endl;
            Double_t DeltaUT2 = uExp - xProj.x(HitsUVT2[i].z(0.));
            if(DeltaUT2<bestDeltaUT2Projected){
               DeltaUT2Projected = DeltaUT2;
            }//DeltaYT1vsyT1->Fill(yT1U,yT1V-yT1VProj);
         }
      }
      if(nT2V>0)
      {
         for(Int_t i = nT2U ; i<nT2U+nT2V;i++)
         {
            Double_t yVT2Exp = yAvg + TY*(HitsUVT2[i].z(0.) - zAvg);
            Double_t vExp = yVT2Exp*HitsUVT2[i].dxDy() + xProj.x(HitsUVT2[i].z(0.));
            //DeltaYT1vsyT1->Fill(yT1U,yT1V-yT1VProj);
            // if(debug) std::cout<<"T2  -- V -- i"<<setw(20)<<i<<"Delta V T2 expected best combination in T1"<<setw(20)<<vExp - HitsUVT2[i].x(0.)<<std::endl;
            Double_t bestDeltaVTProjected =vExp - HitsUVT2[i].x(0.) ;
         }
      }
      Double_t DeltaXT2PairsCombined =-999.;
      if(nT2V>0 && nT2U>0){
         if(debug) std::cout<<"I found possible pairs in UV T2"<<std::endl;
         Double_t BestDeltaCom = 100.;
         Int_t pairN = 0;
         for(Int_t i = 0; i< nT2U;i++ )
         {
            for(Int_t j=nT2U; j< nT2U+ nT2V;j++)
            {
               pairN++;
               Double_t xExpT2Middle = (HitsUVT2[i].x(0.) + HitsUVT2[j].x(0.))/2. ;
               Double_t DeltaT2UVX = xExpT2Middle- xProj.x( (HitsUVT2[i].z(0.)+HitsUVT2[j].z(0.))/2. );
               // if(debug) std::cout <<"Pair"<<setw(20)<<pairN<<setw(20)<<"Delta In T2 after Projection" <<setw(20)<<DeltaT2UVX<<std::endl;
               DeltaXT2PairsCombined = DeltaT2UVX;
            }
         }
      }
   }

   //processXZ(xProj);

   //Study the UV adding of hits with the typical hough transformation
   //std::vector<PatHit> AllUVHit = AllUV.hits();
   minCoord =AllUV[0].coord();
   maxCoord =AllUV[0].coord();
   avgCoord =0;
   minY = 100000.;
   maxY = -100000.;
   NUV=0;
   for(int i = 0 ; i< AllUV.size();i++){
     double y = AllUV[i].coord()*AllUV[i].z();
     if(y<minY){
       minY=y;}
     if(y>maxY){
       maxY=y;}

      if(AllUV[i].coord() < minCoord) minCoord = AllUV[i].coord();
      if(AllUV[i].coord() > maxCoord) maxCoord = AllUV[i].coord();
      avgCoord+=AllUV[i].coord()/AllUV.size();
      NUV++;
   }
   //fitXProjection(xProj,debug);
   //X0Back =(Double_t) (xProj.ax() - m_zReference*xProj.bx()+ 2.458e8*xProj.cx());
   //std::cout<<"X0Back = "<<X0Back<<std::endl;
}

bool TrackStudy::FitLine( PrSeedTrack &UVSegment, Bool_t debug , Int_t Case ){

   LinParFit<double> fit_LineY(2);
   double solution_yz_line[2];
   double solution_yz_par[3];
   std::fill(solution_yz_par, solution_yz_par+3,0.);
   std::fill(solution_yz_line, solution_yz_line+2,0.);
   std::vector<double> dsolution_yz_par;
   std::vector<double> dsolution_yz_line;
   if(debug) UVSegment.PrintTrack();
   double dz = 0.;
   double resYZ_lin = 0.;
   double resYZ_par = 0.;
   for (int j = 0; j < 9;j ++) {
      //dz = 0.;
      //resYZ_lin = 0.;
      LinParFit<double> fit_LineY(2);
      LinParFit<double> fit_ParY(3);
      for(int i = 0 ; i<UVSegment.hits().size(); i++){
         PatHit hit = UVSegment.hits()[i];
         double dz= (double) (UVSegment.hits()[i].z() - m_zReference);
         double y = (double) (hit.x(0)-UVSegment.x(hit.z()))/hit.dxDy();
         resYZ_lin = y - (solution_yz_line[0] + dz* solution_yz_line[1]) ;
         resYZ_par = y - (solution_yz_par[0] + dz* solution_yz_par[1]+ dz*dz*solution_yz_par[2]);
         double erry = std::fabs(1./(hit.w()*hit.dxDy()));
         fit_LineY.accumulate(resYZ_lin, erry , dz );
         fit_ParY.accumulate(resYZ_par,erry,dz);
      }
      if(!fit_LineY.solve()){
         return false;
         if(debug) std::cout<<"Not able to fit the Line in Y"<<std::endl;
      }
      dsolution_yz_par = fit_ParY.solution();
      dsolution_yz_line = fit_LineY.solution();
      solution_yz_line[0]+=dsolution_yz_line[0];
      solution_yz_line[1]+=dsolution_yz_line[1];

      solution_yz_par[0]+=dsolution_yz_par[0];
      solution_yz_par[1]+=dsolution_yz_par[1];
      solution_yz_par[2]+=dsolution_yz_par[2];
   }
   if(debug) std::cout<<"solution ay" <<setw(15)<<solution_yz_line[0]<<std::endl;
   if(debug) std::cout<<"solution by" <<setw(15)<<solution_yz_line[1]<<std::endl;
   double Chi2_ParY = 0.;
   double Chi2_LineBestUVSegm = 0.;
   for(int i =0; i<UVSegment.hits().size(); i++){
      PatHit hit = UVSegment.hits()[i];
      if(debug) std::cout<<"Err X "<<1./hit.w()<<std::endl;
      double y = ( hit.x() - UVSegment.x(hit.z()))/ hit.dxDy();
      if(debug) std::cout<<"y = "<<setw(15)<<y<<std::endl;
      double err = std::abs( 1./hit.w() * 1./hit.dxDy());
      if(debug) std::cout<<"y pred = "<<setw(15)<<solution_yz_line[0]+solution_yz_line[1]*(hit.z()-m_zReference)<<std::endl;
      if(debug) std::cout<<"Err y ="<<setw(15)<<err<<std::endl;
      if(debug) std::cout<<"Contribution = "<< y - (solution_yz_line[0] + solution_yz_line[1]*(hit.z()-m_zReference))<<std::endl;
      Chi2_LineBestUVSegm+= std::pow(y - (solution_yz_line[0] + solution_yz_line[1]*(hit.z()-m_zReference)) ,2 )/ std::pow(err,2);
      Chi2_ParY+= std::pow( y - solution_yz_par[0]- (hit.z()-m_zReference)* solution_yz_par[1]- std::pow(hit.z()-m_zReference,2)* solution_yz_par[2],2)/std::pow(err,2);
   }
   if(debug) std::cout<<"Chi2 Line in Y "<<setw(15)<<Chi2_LineBestUVSegm<<std::endl;
   if(debug) std::cout<<"Chi2 Par in Y"<<setw(15) <<Chi2_ParY<<std::endl;
   if(debug) std::cout<<"Chi2PerDof Line in Y"<<setw(15)<<Chi2_LineBestUVSegm/(UVSegment.hits().size()-2.)<<std::endl;
   if(debug) std::cout<<"Chi2PerDof Par in Y"<<setw(15)<<Chi2_ParY/(UVSegment.hits().size()-3.)<<std::endl;


   if(Case ==0 ){
     FitLineY_ay_case0 = solution_yz_line[0]-m_zReference*solution_yz_line[1];
     FitLineY_by_case0 = solution_yz_line[1];
     FitLineY_Chi2_case0 = Chi2_LineBestUVSegm;
     FitLineY_Chi2DoF_case0 = Chi2_LineBestUVSegm/(UVSegment.hits().size()- 2. );
     FitParY_Chi2_case0 = Chi2_ParY;
     FitParY_Chi2DoF_case0 = Chi2_ParY/(UVSegment.hits().size()-3. );

   }
   if(Case ==1 ){
     FitLineY_ay_case1 = solution_yz_line[0]-m_zReference*solution_yz_line[1];
     FitLineY_by_case1 = solution_yz_line[1];
     FitLineY_Chi2_case1 = Chi2_LineBestUVSegm;
     FitLineY_Chi2DoF_case1 = Chi2_LineBestUVSegm/(UVSegment.hits().size()- 2. );
     FitParY_Chi2_case1 = Chi2_ParY;
     FitParY_Chi2DoF_case1 = Chi2_ParY/(UVSegment.hits().size()-3. );

   }
   if(Case ==2 ){
     FitLineY_ay_case2 = solution_yz_line[0]-m_zReference*solution_yz_line[1];
     FitLineY_by_case2 = solution_yz_line[1];
     FitLineY_Chi2_case2 = Chi2_LineBestUVSegm;
     FitLineY_Chi2DoF_case2 = Chi2_LineBestUVSegm/(UVSegment.hits().size()- 2. );
     FitParY_Chi2_case2 = Chi2_ParY;
     FitParY_Chi2DoF_case2 = Chi2_ParY/(UVSegment.hits().size()-3. );
   }

   return true;
}
bool TrackStudy::FitSimultaneouslyXY(PrSeedTrack& track,PrSeedTrack& UVSegment,Bool_t debug)
{
   // PrSeedTrack track = XSegment;
   // track.PrintTrack();
   // for(int i =0; i< UVSegment.hits().size(); i++){
   //   PatHit hitUV = UVSegment.hits()[i];
   //   track.addHit(hitUV);
   // }
   //track.setParameters(0,0,0,0,0);
   track.sortbyZ();
   //track.PrintTrack();
   float mat[15];
   float rhs[5];
   int nHitsX =0;
   int nHitsStereo=0;
   float m_dRatio = -2.622e-4;
   track.setdRatio(m_dRatio);
   float zRef = m_zReference;
   for(int loop = 0; 10> loop ; ++loop ){
      //std::cout<<"Loop\t"<<loop<<std::endl;
      if(loop ==1 ){
         float RadiusPosition = std::sqrt( (track.ax()*track.ax()*std::fabs(track.ax())/2000.) +
         (track.y(zRef)*track.y(zRef)*std::fabs(track.y(zRef))/1000.));
         float RadiusSlopes = std::sqrt( track.bx()*track.bx()*std::fabs(track.bx())/0.3 +
         track.by()*track.by()*std::fabs(track.by())/0.1   );
         float dRatioPos = - (2.622e-4 +1.943e-8*RadiusPosition + 1.08e-11*RadiusPosition*RadiusPosition);
         float dRatioSlopes = - ( 2.6098e-4+ 6.31e-5*RadiusSlopes  -0.000156778*RadiusSlopes*RadiusSlopes + 0.000134126*RadiusSlopes*RadiusSlopes*RadiusSlopes);
         track.setdRatio(dRatioPos);
      }
      std::fill(mat,mat+15,0.);
      std::fill(rhs,rhs+5,0.);
      float z0 = m_zReference;
      for ( int i =0; i< track.hits().size(); i++ )
      {
         PatHit hit = track.hits()[i];
         if(hit.isX())
         {
           nHitsX++;
         }else{  nHitsStereo++;}
         float yOnTrack = 0.;
         if(loop >0) yOnTrack = hit.yOnTrack(track.y(0.), track.by());
         float w = hit.w2();
         float dxdy = hit.dxDy();
         //      if(debug)std::cout<<" dxDy "<<dxdy<<std::endl;
         float dz =(hit.z(yOnTrack) - m_zReference);
         //  float deta = 0.;
         float dRatio = m_dRatio;
         if(loop >0) dRatio = track.dRatio();
         float deta = dz*dz*(1+m_dRatio*dz);
         float wdz = w * dz;
         float eta = dz * dz * (1. + dz * m_dRatio);
         float weta = w * deta;
         float wdxdy = w * dxdy;
         float wdxdydz = wdxdy * dz;
         float dist = track.distance( hit );
         //Fill Matrix
         mat[0] += w;
         mat[1] += wdz; mat[2] += wdz * dz;
         mat[3] += weta; mat[4] += weta * dz; mat[5] += weta * eta;
         mat[6] -= wdxdy;mat[7] -= wdxdydz;   mat[8] -= wdxdy * deta;  mat[9] += wdxdy * dxdy;
         mat[10] -= wdxdydz; mat[11] -= wdxdydz * dz;  mat[12] -= wdxdydz * deta;  mat[13] += wdxdydz * dxdy; mat[14] += wdxdydz * dz * dxdy;
         // fill right hand side
         rhs[0] += w * dist;
         rhs[1] += wdz * dist;
         rhs[2] += weta * dist;
         rhs[3] -= wdxdy * dist;
         rhs[4] -= wdxdydz * dist;
      }//Loop over Hits to fill the matrix
      // decompose matrix, protect against numerical troubles
      if(nHitsX < 4 || nHitsStereo < 4) return false;
      ROOT::Math::CholeskyDecomp<float, 5> decomp(mat);
      if (!decomp) return false;
      decomp.Solve(rhs);
      //float yAtZRef = rhs[3];
      //rhs[4]*=1.e3;
      //rhs[2]*=1.e6;
      //rhs[1]*=1.e3;
      rhs[3] -= rhs[4] * z0; // ???? this should be here only if the track y part is ay + b_y*(z-zref)
      track.updateParameters(rhs[0],rhs[1],rhs[2],rhs[3],rhs[4],0.);
      // if (std::abs(rhs[0]) > 1e4 || std::abs(rhs[1]) > 5. ||
      //     std::abs(rhs[2]) > 1e-3 || std::abs(rhs[3]) > 1e4 ||
   }
   track.setChi2Full(track.Chi2());
   track.setNDOFFull(track.hits().size()-5);
   //track.PrintTrack();
   //std::cout<<"End fit"<<std::endl;
   return false;
}
