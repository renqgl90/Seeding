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
Int_t NClonesX;
Int_t NPv;
Int_t NUV;



//Coord hough transformation
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
void TrackStudy::Begin(TTree * /*tree*/)
{
   t1->Branch("PhysicsInterest",&PhysicsInterest,"PhysicsInterest/B");
   t1->Branch("isDown",&ISDOWN,"isDown/B");
   t1->Branch("X0Back",&X0Back,"X0Back/D");
   t1->Branch("ISLONG",&ISLONG,"ISLONG/B");
   t1->Branch("OVTX_Z",&OVTX_Z,"OVTX_Z/D");
   t1->Branch("minCoord",&minCoord,"minCoord/D");
   t1->Branch("maxCoord",&maxCoord,"maxCoord/D");
   t1->Branch("avgCoord",&avgCoord,"avgCoord/D");
   t1->Branch("NUV",&NUV,"NUV/D");
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
   t1->Branch("NClonesX",&NClonesX,"NClonesX/I");
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
   t1->Branch("isElectron",&IsElectron,"isElectron/B");
   t1->Branch("MCParticleID",&MCParticleiD,"MCParticleID/I");


   TString option = GetOption();
   m_dRatio0 = -0.000246;
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

Bool_t TrackStudy::Process(Long64_t entry)
{
   fChain->GetTree()->GetEntry(entry);
   Bool_t debug = false;
   if(!(PrHit>1)) return kTRUE;
   ////std::cout<<"Ciao"<<std::endl;
   // HERE DEBUG AND PRINT AT VIDEO OF HIT CONTENT
   Double_t PCut =2000;
   //Counters of Hits
   //Fill for the TTree

   ISLONG= isLong;
   ISDOWN = isDown;
   IsElectron = isElectron;
   MCParticleiD = MCParticleID;
   OVTX_Z = MC_Ovtx_z;
   //Define our True Sample
   PhysicsInterest= (fromDecay || fromPrimaryVertex) && (isLong || isDown ) ;
   if(!PhysicsInterest || std::fabs(MCParticleID==11)) return kTRUE;
   CountHits();
   FitHits();
   Int_t nX=0;
   Int_t zone = 0;
   NClonesX=0;
   if(nXUp>=4)  zone = 1;
   if(nxDown>=4) zone = 0;

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
   CountClones(FullTrack);
   StereoSearch(xProj,UVSegment);
  //  std::cout<<"FullTrack"<<std::endl;
  //  FullTrack.PrintHits();
  //  std::cout<<"OnlyX"<<std::endl;
  //  xProj.PrintHits();
  //  std::cout<<"UVSegment"<<std::endl;
  //  UVSegment.PrintHits();
   std::vector<PatHit> Hits;
   std::vector<PatHit> HitsUVT1;
   std::vector<PatHit> HitsUVT2;
   std::vector<PatHit> HitsUVT3;
   std::vector<PatHit> AllUV;
   //NClonesUV=0;
   for(Int_t i=0;i<PrHit;i++){
      PatHit hit = PatHit();
      hit.setHit(PrHit_Xat0[i],PrHit_Zat0[i],PrHit_dxDy[i],PrHit_dzDy[i],std::sqrt(PrHit_w2[i]),PrHit_yMin[i],PrHit_yMax[i],PrHit_zone[i],PrHit_planeCode[i],PrHit_isX[i],PrHit_LHCbID[i]);
      FTCluster cluster;
      cluster.setCluster(ChID_Fraction[i],PrHit_Size[i],ChID_Charge[i],ChID_SipmID[i],ChID_SipmCell[i],ChID_Module[i],ChID_Layer[i],ChID_Mat[i],ChID_Quarter[i]);
      hit.SetCluster(cluster);
      //Select only Hits with charge ?
      if(PrHit_isX[i])
      {
         nX++;
         if(i>0 && (PrHit_Zat0[i] - PrHit_Zat0[i-1])<10 ){ //here you have a copy
          NClonesX++;
      }

         Hits.push_back(hit);
      }
      if(!PrHit_isX[i])
      {
         if(PrHit_planeCode[i] < 4 )
         {

            HitsUVT1.push_back(hit);
         }
         if(PrHit_planeCode[i]>3 && PrHit_planeCode[i] <8){
            HitsUVT2.push_back(hit);
         }
         if(PrHit_planeCode[i]>8){

            HitsUVT3.push_back(hit);
         }
      }
   }
   for(Int_t i=0;i<HitsUVT1.size();i++){
      ////std::cout<<i<<setw(20)<<sqrt(1./HitsUVT1[i].w2())<<setw(20)<<HitsUVT1[i].x(0.)<<setw(20)<<HitsUVT1[i].z(0.)<<setw(20)<<HitsUVT1[i].planeCode()<<setw(20)<<HitsUVT1[i].zone()<<setw(20)<<HitsUVT1[i].dxDy()<<setw(20)<<HitsUVT1[i].dzDy()<<setw(20)<<HitsUVT1[i].cluster().charge()<<setw(20)<<HitsUVT1[i].cluster().fraction()<<setw(20)<< HitsUVT1[i].cluster().size() <<setw(20)<<HitsUVT1[i].cluster().charge()<<std::endl;
   }
   PrSeedTrack xProje = PrSeedTrack(zone,m_zReference,Hits);

   if(debug) xProje.PrintHits();

   //xProj.PrintHits();
   fitXProjection(xProje);
   X0Back = xProj.ax() - m_zReference*xProj.bx()+ 2.458e8*xProj.cx();
   for(Int_t i=0;i<PrHit;i++){
      if(!PrHit_isX[i])
      {

         PatHit HitUV = PatHit();
         HitUV.setHit(PrHit_Xat0[i],PrHit_Zat0[i],PrHit_dxDy[i],PrHit_dzDy[i],std::sqrt(PrHit_w2[i]),PrHit_yMin[i],PrHit_yMax[i],PrHit_zone[i],PrHit_planeCode[i],PrHit_isX[i],PrHit_LHCbID[i]);
         FTCluster cluster;
         cluster.setCluster(ChID_Fraction[i],PrHit_Size[i],ChID_Charge[i],ChID_SipmID[i],ChID_SipmCell[i],ChID_Module[i],ChID_Layer[i],ChID_Mat[i],ChID_Quarter[i]);
         HitUV.SetCluster(cluster);
         HitUV.setCoord( (PrHit_Xat0[i]-xProj.x(PrHit_Zat0[i]))/(PrHit_dxDy[i]*PrHit_Zat0[i]));
         AllUV.push_back(HitUV);
      }
   }



   if(NClonesX>0 && !isElectron && Eta_in25 && isLong && P>5000)
   {

      //xProj.PrintTrack();
   }
   //cout<<"ax"<<setw(10)<<"bx"<<setw(10)<<"cx \n"<<endl;
   //cout<<xProj.ax()<<setw(10)<<xProj.bx()<<setw(10)<<xProj.cx()<<endl;
   if(debug){
      //xProj.PrintTrack();
   }
   //Here the study for the Fitting on PrHit
   ////std::cout<<"Chi2"<<setw(20)<<"Chi2 per DoF"<<std::endl;
   ////std::cout<<xProj.Chi2()<<setw(20)<<xProj.Chi2()/((double)xProj.hits().size()-3);

   ax_XStepHitFit=xProj.ax();
   bx_XStepHitFit=xProj.bx();
   cx_XStepHitFit=xProj.cx();
   Chi2_XStepFit=xProj.Chi2();
   NClustersX = xProj.hits().size();
   Chi2_perDoF_XStepFit = xProj.Chi2()/((double)xProj.hits().size()-3.);
   XZStudy(xProj);
   //int num_items1 = std::count_if(v.begin(), v.end(), [](int i) {return i % 3 == 0;});
   //nCombT1=0;
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
      ////std::cout<<"In T1 i found possible pairs from a track which has part =     "<<xProj.zone()<<std::endl;
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
            ////std::cout<<"U Hit:   PlaneCode"<<setw(20)<<HitsUVT1[i].planeCode()<<setw(20)<<"dxDy"<<HitsUVT1[i].dxDy()<<setw(20)<<"X at 0"<<HitsUVT1[i].x(0.)<<setw(20)<<"Z at 0"<<setw(20)<<HitsUVT1[i].z(0.)<<setw(20)<<"Yexp"<<setw(20)<<yT1U<<std::endl;
            ////std::cout<<"V Hit:   PlaneCode"<<setw(20)<<HitsUVT1[j].planeCode()<<setw(20)<<"dxDy"<<HitsUVT1[j].dxDy()<<setw(20)<<"X at 0"<<HitsUVT1[j].x(0.)<<setw(20)<<"Z at 0"<<setw(20)<<HitsUVT1[j].z(0.)<<setw(20)<<"Yexp"<<setw(20)<<yT1V<<std::endl;

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
         //////std::cout<<"X Combination"<<setw(20)<<XCombined<<std::endl;
         //////std::cout<<"X Predicted  "<<setw(20)<<xProj.x( (HitsUVT1[i].z(0.)+ HitsUVT1[j].z(0.) )/2.)<<std::endl;
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
   ////std::cout<<"\n \n \n"<<setw(20)<<"New Track"<<std::endl;
   //xProj.PrintHits();
   if(nT1U>0 && nT1V>0)
   {
      ////std::cout<<"Sorted Combinations"<<std::endl;
      std::sort(DeltaXCombT1.begin(),DeltaXCombT1.end(),[](Double_t Value1, Double_t Value2)->bool{return std::abs(Value1)<std::abs(Value2);});
      for(Int_t i = 0 ; i<DeltaXCombT1.size();i++)
      {
         if(DeltaXCombT1.size()>1)
         {
            ////std::cout<<"i"<<setw(20)<<"Delta X Best Combination"<<setw(20)<<DeltaXCombT1[i]<<std::endl;
         }
      }
      DeltaXBestCombT1 = DeltaXCombT1[0];
      ////std::cout<<"DeltaBestCombT1   "<<setw(20)<<DeltaXBestCombT1<<std::endl;
      ////std::cout<<"Compute Delta Best"<<setw(20)<<(bestHitPair.first.x(0.)+bestHitPair.second.x(0.))/2 - xProj.x( (bestHitPair.first.z(0.) + bestHitPair.second.z(0.)) /2.)<<std::endl;
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
            // ////std::cout<<"T2  -- U -- i ="<<setw(20)<<i<<setw(20)<<"Delta U T2 expected best combination in T1"<<setw(20)<<uExp - HitsUVT2[i].x(0.)<<std::endl;
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
            // ////std::cout<<"T2  -- V -- i"<<setw(20)<<i<<"Delta V T2 expected best combination in T1"<<setw(20)<<vExp - HitsUVT2[i].x(0.)<<std::endl;
            Double_t bestDeltaVTProjected =vExp - HitsUVT2[i].x(0.) ;
         }
      }
      Double_t DeltaXT2PairsCombined =-999.;
      if(nT2V>0 && nT2U>0){
         ////std::cout<<"I found possible pairs in UV T2"<<std::endl;
         Double_t BestDeltaCom = 100.;
         Int_t pairN = 0;
         for(Int_t i = 0; i< nT2U;i++ )
         {
            for(Int_t j=nT2U; j< nT2U+ nT2V;j++)
            {
               pairN++;
               Double_t xExpT2Middle = (HitsUVT2[i].x(0.) + HitsUVT2[j].x(0.))/2. ;
               Double_t DeltaT2UVX = xExpT2Middle- xProj.x( (HitsUVT2[i].z(0.)+HitsUVT2[j].z(0.))/2. );
               // ////std::cout <<"Pair"<<setw(20)<<pairN<<setw(20)<<"Delta In T2 after Projection" <<setw(20)<<DeltaT2UVX<<std::endl;
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
   NUV=0;
   for(int i = 0 ; i< AllUV.size();i++){
      if(AllUV[i].coord() < minCoord) minCoord = AllUV[i].coord();
      if(AllUV[i].coord() > maxCoord) maxCoord = AllUV[i].coord();
      avgCoord+=AllUV[i].coord()/AllUV.size();
      NUV++;
   }

   //if()

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

Bool_t TrackStudy::fitXProjection(PrSeedTrack & track)
{
   float mat[6];
   float rhs[3];
   std::vector<PatHit> Hits = track.hits();

   ////std::cout<<"Hits Size"<<Hits.size()<<std::endl;
   ////std::cout<<"Hits Size from track"<<track.hits().size();
   float dRatio=0;
   for(int loop = 0;3>loop;++loop){

      std::fill(mat,mat+6,0.);
      std::fill(rhs,rhs+3,0.);
      ////std::cout<<"******************* Fit Loop "<<loop<<"*********************"<<std::endl;
      for( int i=0;i<Hits.size();i++){
         const float w = Hits[i].w2();
         ////std::cout<<"Hit w2"<<Hits[i].w2()<<std::endl;
         const float dz = Hits[i].z(0.)-m_zReference;
         float deta = dz*dz;
         //if(m_usedRatio){
         //deta = dz*dz*(1-dRatio*dz);
         //}
         float dist = track.distance(Hits[i]);
         ////std::cout<<"i"<<setw(20)<<i<<setw(20)<<"Hit X"<<setw(20)<<Hits[i].x(0)<<setw(20)<<"Track x"<<setw(20)<<track.x(Hits[i].z(0.))<<setw(20)<<"Distance"<< track.distance(Hits[i])<<std::endl;
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
         ////std::cout<<"Failed to decompose matrix"<<std::endl;
         return false;
      }
      decomp.Solve(rhs);
      track.updateParameters(rhs[0],rhs[1],rhs[2],0.,0.);
   }
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
   std::cout<<"Number of Geant MCHit \t"<<m_long_NGeantHit<<std::endl;
   std::cout<<"Number of different Layers with Geant MCHit \t"<<m_long_NLayer_Geant<<std::endl;

   std::cout<<"Number of MCHit into Cluster \t"<<m_long_NMCHitIntoCluster<<std::endl;
   std::cout<<"Number of different Layers with MCHit into Cluster\t"<<m_long_NLayer_MCHitInCluster<<std::endl;

   std::cout<<"Number of Clusters \t" <<m_long_NPrHit<<std::endl;
   std::cout<<"Number of diffrent Layers Clusters \t" <<m_long_NLayer_PrHit<<std::endl;
   zMCHitClones->Draw();
   zMCHitClonevsP->Draw("colz");


   std::cout<<"Physical Tracks = "<<setw(20)<<m_physical<<std::endl;
   std::cout<<"Physical Tracks (Case0 intrinsic)= "<<setw(20)<<m_physicalCase0Intrinsic<<std::endl;
   std::cout<<"Physical Tracks (Case0 intrinsic) duplicates"<<setw(20)<<m_physicalCase0IntrinsicWithduplicates<<std::endl;
   std::cout<<"Physical Tracks (Case0 intrinsic) >=4= "<<setw(20)<<m_physicalCase0IntrinsicMore4<<std::endl;
   std::cout<<"Physical Tracks (Case0 intrinsic) >=4 duplicates"<<setw(20)<<m_physicalCase0IntrinsicWithduplicatesMore4<<std::endl;
   std::cout<<"Physical Tracks (Case0 intrinsic) >=5= "<<setw(20)<<m_physicalCase0IntrinsicMore5<<std::endl;
   std::cout<<"Physical Tracks (Case0 intrinsic) duplicates"<<setw(20)<<m_physicalCase0IntrinsicWithduplicatesMore5<<std::endl;



   std::cout<<"NbTrackClones X \t"<<m_nbTrackCloneX<<std::endl;
   std::cout<<"NbTrackClones UV \t"<<m_nbTrackCloneUV<<std::endl;
   std::cout<<"NbTrackClones X and UV \t"<<m_nbTrackCloneUVandX<<std::endl;
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}







void TrackStudy::CountHits(){

  if( (isLong || isDown) && (std::abs(MCParticleID)!=11) ){
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
     //Count different MCHitIntoCluster Per Layer;
  }


}


void TrackStudy::FitHits(){

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
         ////std::cout<<"MCHit i = "<<i<<"\t Z \t"<< MCHit_Assoc_Z[i]<<"\t X\t"<<MCHit_Assoc_X[i]<<"\t Y \t"<<MCHit_Assoc_Y[i]<<"\tpathlenght \t"<<MCHit_pathlength[i]<<"\t"<<std::endl;
         ////std::cout<<"MCHit i -1 = "<<i<<"\t Z \t"<< MCHit_Assoc_Z[i-1]<<"\t X\t"<<MCHit_Assoc_X[i]<<"\t Y \t"<<MCHit_Assoc_Y[i-1]<<"\t pathlenght   \t"<<MCHit_pathlength[i-1]<<"\t"<<std::endl;
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
         ////std::cout<<"Not Able to Fit XZ with a parabola for a MCHit list of size "<<MC_ass<<"\t and Momentum "<<P<<std::endl;
         ok=false;
      }
      if(!fit_LineY.solve())
      {
         ////std::cout<<"Not Able to Fit Y with a Line for a MCHit list of size"<<MC_ass<<"\t and Momentum "<<P<<std::endl;
         ok=false;
      }
      if(!fit_CubicXZ.solve()){
         ////std::cout<<"Not Able to Fit Cubic XZ for a MCHit list of size"<<MC_ass<<"\t and Momentum "<<P<<std::endl;
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
         //  ////std::cout<<"Failed to decompose matrix"<<std::endl;
         Fit = false;
      }
      if(decomp) Fit = true;
      decomp.Solve(rhs);
      solution[0]+=rhs[0];
      solution[1]+=rhs[1];
      solution[2]+=rhs[2];
   }
   if(Fit){
      ////std::cout<<"ax"<<setw(20)<<"bx"<<setw(20)<<"cx"<<setw(20)<<std::endl;
      ////std::cout<<solution[0]<<setw(20)<<solution[1]<<setw(20)<<solution[2]<<std::endl;
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
         //  ////std::cout<<"Failed to decompose matrix"<<std::endl;
         Fit1 = false;
      }
      if(decomp) Fit1 = true;
      decomp.Solve(rhs1);
      solution1[0]+=rhs1[0];
      solution1[1]+=rhs1[1];
      solution1[2]+=rhs1[2];
   }
   if(Fit){
      ////std::cout<<"ax"<<setw(20)<<"bx"<<setw(20)<<"cx"<<setw(20)<<std::endl;
      ////std::cout<<solution[0]<<setw(20)<<solution[1]<<setw(20)<<solution[2]<<std::endl;
   }
   Chi2_XZDRATIOFixed = 0;
   for( int i = 0;i<MC_ass;i++){
      float dz1 = MCHit_Assoc_Z[i]-m_zReference;
      float dist1 = MCHit_Assoc_X[i]-(solution1[0]+solution1[1]*dz1+solution1[2]*dz1*dz1*(1-dRatio*dz1));
      Chi2_XZDRATIOFixed+= (dist1*dist1)/(std::pow(0.100,2));
   }
   Chi2_XZDRATIOFixed = Chi2_XZDRATIOFixed/(MC_ass-4);
   ////std::cout<<"Chi2 XZDRATIOFIXED \t"<<Chi2_XZDRATIOFixed<<std::endl;
   ////std::cout<<"Chi2 = "<<Chi2_XZDRATIO;
}

void TrackStudy::XZStudy(PrSeedTrack xProj){
   m_physical++;
   // ////std::cout<<"Start XZ study"<<std::endl;
   //extract the T1X hit
   //std::vector<PatHit>  XZHit;
   std::vector<PatHit> Case0First ;
   Case0First.reserve(10);
   std::vector<PatHit> Case1First ;
   Case1First.reserve(10);

   std::vector<PatHit> Case2First ;
   Case2First.reserve(10);

   std::vector<PatHit> Case0Last ;
   Case0Last.reserve(10);

   std::vector<PatHit> Case1Last  ;
   Case1Last.reserve(10);

   std::vector<PatHit> Case2Last  ;
   Case2Last.reserve(10);
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
    NClonesX=0;
    if(nXUp>=4)  zone = 1;
    if(nxDown>=4) zone = 0;
   PrSeedTrack xHits = PrSeedTrack(zone, m_zReference);
   PrSeedTrack UV = PrSeedTrack(zone,m_zReference);
   //if(PrHit<=8) return kTRUE;
   //xHits.sortbyZ();
   //bool cloneX =false;
   //bool cloneUV = false;
   int clUV =0;
   int clX = 0;
   for( int i = 0; i< PrHit;i++)
   {
      PatHit hit = PatHit();
      hit.setHit(PrHit_Xat0[i],PrHit_Zat0[i],PrHit_dxDy[i],PrHit_dzDy[i],std::sqrt(PrHit_w2[i]),PrHit_yMin[i],PrHit_yMax[i],PrHit_zone[i],PrHit_planeCode[i],PrHit_isX[i],PrHit_LHCbID[i]);
      cluster.setCluster(ChID_Fraction[i],PrHit_Size[i],ChID_Charge[i],ChID_SipmID[i],ChID_SipmCell[i],ChID_Module[i],ChID_Layer[i],ChID_Mat[i],ChID_Quarter[i]);
      hit.SetCluster(cluster);
      if(!PrHit_isX[i])
      {
         if( (i>0) &&  (PrHit_planeCode[i]==PrHit_planeCode[i-1]) && (PrHit_LHCbID[i]!=PrHit_LHCbID[i-1]) )
         {
            clUV++;
            ////std::cout<<"hasCloneUV"<<std::endl;
         }
         UV.addHit(hit);
      }
      if(PrHit_isX[i])
      {
         xHits.addHit(hit);
         if( (i>0) &&  (PrHit_planeCode[i]==PrHit_planeCode[i-1] ) && (PrHit_LHCbID[i]!=PrHit_LHCbID[i-1]) ){
            //cloneX=true;
            clX++;
         }
         //////std::cout<<"HAS X LAYER CLONE"<<std::endl;
         //PatHit hit = PatHit();
         //////std::cout<<"Filling Vectors"<<std::endl;
         //////std::cout<<"Get Hit "<<i<<std::endl;
         //////std::cout<<"PlaneCode"<<hit.planeCode()<<std::endl;
         if(hit.planeCode() ==4)
         {
            ////std::cout<<"ParabolaSeed Hits 4 loaded"<<std::endl;
            ParabolaSeedHits1.push_back(hit);
         }
         if(hit.planeCode() ==7)
         {
            ////std::cout<<"ParabolaSeed Hits 7 loaded"<<std::endl;
            ParabolaSeedHits2.push_back(hit);
         }
         if(hit.planeCode() == 0){ //T1-1X first in case 0 and case 2
            Case0First.push_back(hit);
            // Case2First.push_back(hit);
            ////std::cout<<"ParabolaSeed Hits 0 loaded"<<std::endl;
         }
         if(hit.planeCode() == 3){ //T1-2X first in case 1
            // Case1First.push_back(hit);
            ////std::cout<<"ParabolaSeed Hits 3 loaded"<<std::endl;
         }
         if(hit.planeCode() == 11){ //T3-2X last in case 0 and case 1
            Case0Last.push_back(hit);
            // Case1Last.push_back(hit);
            ////std::cout<<"ParabolaSeed Hits 11 loaded"<<std::endl;
         }

         if(hit.planeCode() == 8){
            // Case2Last.push_back(hit);
            ////std::cout<<"Palane 8 Hits loaded"<<std::endl;
         }

         // int planeCode = hit.planeCode();
         // ////std::cout<<"Get PlaneCode"<<std::endl;

         //Case 0 Seed 1
         if( hit.planeCode() == 3 || hit.planeCode() ==7 || hit.planeCode() == 8) {
            ////std::cout<<"Plane 3,4,8 Hits loaded"<<std::endl;

            if(hit.planeCode() == 3) remainingCase0Seed1N[0]++;
            if(hit.planeCode() == 7) remainingCase0Seed1N[1]++;
            if(hit.planeCode() == 8) remainingCase0Seed1N[2]++;
            ////std::cout<<"Planes\t"<<hit.planeCode()<<"\t Hits done"<<std::endl;

            RemainingCase0Seed1.push_back(hit);
         }

         //Case 0 Seed 2
         if(  hit.planeCode() == 3 || hit.planeCode() ==4 || hit.planeCode() == 8){
            ////std::cout<<"Planes\t"<<hit.planeCode()<<"\t Hits done"<<std::endl;
            if(hit.planeCode() ==3) remainingCase0Seed2N[0]++;
            if(hit.planeCode() ==4) remainingCase0Seed2N[1]++;
            if(hit.planeCode() ==8) remainingCase0Seed2N[2]++;
            RemainingCase0Seed2.push_back(hit);
         }
      }
      ////std::cout<<"Done"<<std::endl;
   }
   bool PrelimSele_Case0 = false;
   bool PrelimSele_Case1 = false;
   bool PrelimSele_Case2 = false;
   if(Case0First.size()!=0 && Case0Last.size()!=0 && (ParabolaSeedHits1.size()!=0 || ParabolaSeedHits2.size()!=0)) PrelimSele_Case0 = true;
   if(Case1First.size()!=0 && Case1Last.size()!=0) PrelimSele_Case1 = true;
   if(Case2First.size()!=0 && Case2Last.size()!=0) PrelimSele_Case2 = true;
   if(clUV>0)
   {
      //UV.PrintHits();
      m_nbTrackCloneUV++;
   }
   if(clX>0)
   {
      //xHits.PrintHits();
      m_nbTrackCloneX++;
   }
   if(clX>0 && clUV>0){
      m_nbTrackCloneUVandX++;
   }
   //std::vector<double> DeltaZ;
   if(PrelimSele_Case0)
   {
      m_physicalCase0Intrinsic++;
      PrSeedTracks xProjections;
      for(int i = 0; i<Case0First.size();i++)
      {
         //std::cout<<"First Hit "<<std::endl;
         //Case0First[i].PrintHit();
         ////std::cout<<"Loop first Layer"<<std::endl;
         for(int j = 0;j<Case0Last.size();j++)
         {
            //std::cout<<"Second Hit"<<std::endl;
            //Case0Last[j].PrintHit();
            ////std::cout<<"Loop Last Layer"<<std::endl;
            // double deltaz = Case0Last[j].z()-Case0First[i].z();
            // double tx_pickedcombination = ( Case0Last[j].x(0.)-Case0First[i].x(0.))  /( Case0Last[j].z()-Case0First[i].z());
            // double xFirstProjected = Case0First[i].x(0.) + Case0First[i].x(0.)/Case0First[i].z(0.)*( deltaz);
            // double txinf = Case0First[i].x(0.)/Case0First[i].z(0.);
            //PrSeedTrack xProjection1(0, m_zReference);
            //PatHit fHit = Case0First[i];

            //PatHit lHit = Case0Last[j];
            ////std::cout<<"Loop Last "<<std::endl;
            ////std::cout<<"parabolaSeedHits1 sizee   "<<ParabolaSeedHits1.size()<<std::endl;
            //FInd Given the parabolaSeedHit1 the best combination
            for(int k=0; k< ParabolaSeedHits1.size();k ++)
            {
               //std::cout<<"Parabola Seed Hit 1"<<std::endl;
               //ParabolaSeedHits1[k].PrintHit();
               ////std::cout<<"Loop ParabolaSeed 1"<<std::endl;
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
               ////std::cout<<"a\t"<<a<<std::endl;
               ////std::cout<<"b\t"<<b<<std::endl;
               ////std::cout<<"c\t"<<c<<std::endl;

               //std::vector<PatHit> bestRemaining(3);
               int best1_kk;
               int best2_kk;
               int best3_kk;
               // PatHit hit1_1;
               // PatHit hit2_1;
               // PatHit hit3_1;
               std::vector<double> minVal(3);
               minVal[0]=1000; minVal[1]=1000; minVal[2]=1000;
               for(int kk =0;kk<remainingCase0Seed1N[0] ; kk++){
                  //std::cout<<"Remaining Seed 1 -1 "<<std::endl;
                  //RemainingCase0Seed1[kk].PrintHit();
                  double dz = RemainingCase0Seed1[kk].z() -m_zReference;
                  double xAtZ = a*dz*dz*(1.+m_dRatio0*dz) + b*dz +c;
                  double delta = std::fabs(RemainingCase0Seed1[kk].x(0)-xAtZ);
                  ////std::cout<<"delta Remaining 1 -1 \t"<<delta<<std::endl;
                  if(std::fabs(RemainingCase0Seed1[kk].x(0)-xAtZ)<minVal[0]){
                     minVal[0]=std::fabs(delta);
                     best1_kk=kk;
                     //hit1_1 = RemainingCase0Seed1[kk];
                  }
               }
               if(minVal[0]<999){
                  xProjection.addHit(RemainingCase0Seed1[best1_kk]);
               }
               for(int kk =remainingCase0Seed1N[0];kk< ( remainingCase0Seed1N[0]+remainingCase0Seed1N[1]) ; kk++){
                  double dz = RemainingCase0Seed1[kk].z() - m_zReference;
                  double xAtZ = a*dz*dz*(1.+m_dRatio0*dz) + b*dz +c;
                  double delta = std::fabs(RemainingCase0Seed1[kk].x(0)-xAtZ);
                  ////std::cout<<"delta Remaining 1 -2 \t"<<delta<<std::endl;

                  if(std::fabs(delta)<minVal[1]){
                     minVal[1]=std::fabs(delta);
                     best2_kk=kk;
                     //hit2_1 = RemainingCase0Seed1[kk];
                  }
               }
               if(minVal[1]<999){
                  xProjection.addHit(RemainingCase0Seed1[best2_kk]);
               }
               for( int kk= (remainingCase0Seed1N[0]+remainingCase0Seed1N[1]);kk<( remainingCase0Seed1N[0]+remainingCase0Seed1N[1]+remainingCase0Seed1N[2] ); kk++){
                  double dz = RemainingCase0Seed1[kk].z() -m_zReference;
                  double xAtZ = a*dz*dz *(1.+m_dRatio0*dz)+ b*dz +c;
                  double delta = std::fabs(RemainingCase0Seed1[kk].x(0)-xAtZ);
                  ////std::cout<<"delta Remaining 1 -3\t"<<delta<<std::endl;

                  if(std::fabs(RemainingCase0Seed1[kk].x(0)-xAtZ)<minVal[2]){
                     minVal[2]=std::fabs(delta);
                     best3_kk=kk;
                     //hit3_1 = RemainingCase0Seed1[kk];
                  }
               }
               if(minVal[2]<999){
                  xProjection.addHit(RemainingCase0Seed1[best3_kk]);
               }
               //std::cout<<"Pushing Back track Case Seed 1"<<std::endl;
               //xProjection.PrintHits();
               //std::cout<<"Case Seed 1"<<std::endl;

               xProjections.push_back(xProjection);
               //T2
            }
            for(int k=0; k< ParabolaSeedHits2.size();k ++)
            {
               //std::cout<<"Parabola Seed Hit 1"<<std::endl;
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
               //std::cout<<"a \t"<<a<<std::endl;
               //std::cout<<"b \t"<<b<<std::endl;
               //std::cout<<"c \t"<<c<<std::endl;
               int best1_kk;
               int best2_kk;
               int best3_kk;
               // PatHit hit1_1;
               // PatHit hit2_1;
               // PatHit hit3_1;
               std::vector<double> minVal(3);
               minVal[0]=1000; minVal[1]=1000; minVal[2]=1000;
               for(int kk =0;kk<remainingCase0Seed2N[0] ; kk++){
                  double dz = RemainingCase0Seed2[kk].z() -m_zReference;
                  //double m_dRatio0 = -0.000246;
                  double xAtZ = a*dz*dz*(1.+m_dRatio0*dz) + b*dz +c;
                  double delta = std::fabs(RemainingCase0Seed2[kk].x(0)-xAtZ);
                  //std::cout<<"delta 2-1 "<<delta<<std::endl;
                  if(std::fabs(RemainingCase0Seed2[kk].x(0)-xAtZ)<minVal[0]){
                     minVal[0]=std::fabs(delta);
                     best1_kk = kk;
                     //hit1_1 = RemainingCase0Seed2[kk];
                  }
               }
               if(minVal[0]<999){
                  xProjection.addHit(RemainingCase0Seed2[best1_kk]);
               }
               for(int kk =remainingCase0Seed2N[0];kk<(remainingCase0Seed2N[0]+remainingCase0Seed2N[1]) ; kk++){
                  double dz = RemainingCase0Seed1[kk].z() -m_zReference;
                  double xAtZ = a*dz*dz*(1.+m_dRatio0*dz) + b*dz +c;
                  double delta = std::fabs(RemainingCase0Seed2[kk].x(0)-xAtZ);
                  //std::cout<<"Delta 2-2 \t"<<delta<<std::endl;
                  if(std::fabs(delta)<minVal[1]){
                     minVal[1]=std::fabs(delta);
                     best2_kk=kk;
                     // hit2_1 = RemainingCase0Seed2[kk];
                  }
               }
               if(minVal[1]<999){
                  xProjection.addHit(RemainingCase0Seed2[best2_kk]);
               }
               for(int  kk=(remainingCase0Seed2N[0]+remainingCase0Seed2N[1]) ;kk<(remainingCase0Seed2N[0]+remainingCase0Seed2N[1]+remainingCase0Seed2N[2]); kk++){
                  double dz = RemainingCase0Seed2[kk].z() -m_zReference;
                  // double m_dRatio0 = -0.000246;
                  double xAtZ = a*dz*dz *(1.+m_dRatio0*dz) + b*dz +c;

                  double delta = std::fabs(RemainingCase0Seed2[kk].x(0)-xAtZ);
                  //std::cout<<"Delta 2-3\t"<<delta<<std::endl;

                  if(std::fabs(RemainingCase0Seed2[kk].x(0)-xAtZ)<minVal[2]){
                     minVal[2]=std::fabs(delta);
                     best3_kk=kk;
                     //hit3_1 = RemainingCase0Seed2[kk];
                  }
               }
               if(minVal[2]<999){
                  xProjection.addHit(RemainingCase0Seed2[best3_kk]);
               }
               //std::cout<<"Pushing Back track Case Seed 1"<<std::endl;
               //xProjection.PrintHits();
               //std::cout<<"Case Seed 1"<<std::endl;

               xProjections.push_back(xProjection);
            }

         }
      }
      //std::cout<<"I found xCandidates"<<xProjections.size()<<std::endl;
      for(int i =0 ; i< xProjections.size() ; i++){
         PrSeedTrack track= xProjections[i];
         track.sortbyZ();
         if(clX>0){
            ////std::cout<<"Track had Clones in X"<<std::endl;
            //track.PrintHits();
         }
         // PrSeedTrack track = xProjections[0];
         ////std::cout<<"Number of XZ Projections before erase"<<xProjections.size()<<std::endl;
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
      if(xProjections.size()==0){std::cout<<"I didn't find any xProjection"<<std::endl; xProj.PrintHits();}
      if(xProjections.size()>1) m_physicalCase0IntrinsicWithduplicates++;
      if(xProjections.size()>1 && xProjections[0].hits().size()>=4) m_physicalCase0IntrinsicWithduplicatesMore4++;
      if(xProjections.size()>1 && xProjections[0].hits().size()>=5) m_physicalCase0IntrinsicWithduplicatesMore5++;
      //if(xProjections.size()>0){
      if(xProjections[0].hits().size()>=4) m_physicalCase0IntrinsicMore4++;
      if(xProjections[0].hits().size()>=5) m_physicalCase0IntrinsicMore5++;
      //}


      ////std::cout<<"Number of XZ Projections after erase"<<xProjections.size()<<std::endl;

   }//End Case 0 Prelim Sele
}

// void AddStereoStudy(){
//
//
//
// }
// void AddStereoStudy(PrSeedTrack xProj){
//
//
// }
