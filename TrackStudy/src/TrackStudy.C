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
#include <TStyle.h>
#include "LinParFit.h"
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
Double_t Momentum;
Double_t Px;
Double_t Py;
Double_t Track_eta;
Double_t Pt;
Double_t Track_Phi;
Double_t Track_Pt;

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
  Double_t DeltaXBestCombCorrSlopeT1;
  Double_t TyPairT1;
  Double_t YPairT1;
  Double_t XPairT1;
  Double_t DeltaUT2Projected;
  Double_t DeltaVT2Projected;
  Double_t YUT2Projected;
  Double_t YVT2Projected;
void TrackStudy::Begin(TTree * /*tree*/)
{
  //t1->Branch("nCombT1",nCombT1,"nCombT1/D");
  //t1->Branch("DeltaYT1",)
  t1->Branch("DeltaXT1BestCorrSlopeT1",&DeltaXBestCombCorrSlopeT1,"DeltaXT1BestCorrSlopeT1/D");

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
  TString option = GetOption();

  //Counters For Efficiencies on Long tracks
  m_long_NGeantHit=0;
  m_long_NPrHit=0;
  m_long_NMCHitIntoCluster=0;
  m_long_NLayer_Geant=0;
  m_long_NLayer_PrHit=0;
  m_long_NLayer_MCHitInCluster=0;



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
  if(!(PrHit>0)) return kTRUE;
  //std::cout<<"Ciao"<<std::endl;
  // HERE DEBUG AND PRINT AT VIDEO OF HIT CONTENT
  Double_t PCut =2000;
  //Counters of Hits
  if(isLong && P>PCut && !isElectron){
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
  if(isElectron) return kTRUE; //No Electrons
  //if(!Eta_in25) return kTRUE; //Eta Cut
  //if(P>5000) return kTRUE;
  //if(!isSeed) return kTRUE;
  //if(P<2000) return kTRUE; //Study done for P>2000 !electron & LongTracks
  if(MC_Ovtx_z>4000) return kTRUE; //Study done for Tracks not generated in interactions (naively)
  //std::vector<PatHit> track(100); //Create the Vector of Hit to be analysed in terms of hit content (remove some of them? (like duplicates in layers?))
  std::vector<MCHit> track_MC(100); //Vector of MCHit going into cluster
  std::vector<MCHit> Particle_GeantHit(100);
  //Create a vector of hits to be placed in the PrSeedTrack to be analysed: check if it's the case to remove clones

  if(debug)
  {
    if(CheatedSeeding_NHits>0)
    {
      std::cout<<"i \t Xat0 \t Zat0 \t planeCode \t zone \t dxDy \t dzDy"<<"\t Fraction"<<"\t Size"<<"\t Charge"<<std::endl;
    }
    for (Int_t i = 0 ;  i< CheatedSeeding_NHits; i++){
      PatHit hitdeb = PatHit();
      hitdeb.setHit(PrHit_Xat0[i],PrHit_Zat0[i],PrHit_dxDy[i],PrHit_dzDy[i],std::sqrt(PrHit_w2[i]),PrHit_yMin[i],PrHit_yMax[i],PrHit_zone[i],PrHit_planeCode[i],PrHit_isX[i],PrHit_LHCbID[i]);
      //hit.setCluser(); Need the corresponding branches in the ntuple
      std::cout<<i<<"\t"<<hitdeb.x(0.)<<"\t"<<hitdeb.z(0.)<<"\t"<<hitdeb.planeCode()<<"\t"<<hitdeb.cluster().charge()<<"\t"<<hitdeb.cluster().fraction()<<"\t"<<hitdeb.cluster().size()<<"\t"<<hitdeb.cluster().charge()<<std::endl;
    }
  }

  std::cout.precision(4);
  for (Int_t i = 0;  i< MC_ass; i++) {
    MCHit mcHit =  MCHit();
    if(i>0 && std::fabs(MCHit_Assoc_Z[i]-MCHit_Assoc_Z[i-1])<40) {
      zMCHitClones->Fill(MCHit_Assoc_Z[i]);
      zMCHitClonevsP->Fill(MCHit_Assoc_Z[i],P);
      xvsyClones->Fill(MCHit_Assoc_X[i],MCHit_Assoc_Y[i]);
      //std::cout<<"MCHit i = "<<i<<"\t Z \t"<< MCHit_Assoc_Z[i]<<"\t X\t"<<MCHit_Assoc_X[i]<<"\t Y \t"<<MCHit_Assoc_Y[i]<<"\tpathlenght \t"<<MCHit_pathlength[i]<<"\t"<<std::endl;
      //std::cout<<"MCHit i -1 = "<<i<<"\t Z \t"<< MCHit_Assoc_Z[i-1]<<"\t X\t"<<MCHit_Assoc_X[i]<<"\t Y \t"<<MCHit_Assoc_Y[i-1]<<"\t pathlenght \t"<<MCHit_pathlength[i-1]<<"\t"<<std::endl;
      continue; //Merge them or keep the best one?
    }
    mcHit.setMCHit(MCHit_Assoc_X[i], MCHit_Assoc_Y[i],MCHit_Assoc_Z[i],MCHit_tx[i],MCHit_ty[i],MCHit_p[i],MCHit_pathlength[i],P,MC_px,MC_py,MC_pz);
    track_MC.push_back(mcHit);
  }
  //Track_Phi=Phi
  //MCHit mcHit_Geant = MCHIt();
  //  std::cout<<"MCHit from Geant with size "<<MC<<"\n"
  //<<"i \t \t X \t \t Y \t \t Z \t \t P \t \t PathLenght \t ParticleP \t time  "<<std::endl;}

  //std::vector<MCHit> CloneHits_OnTrack;
  //for(Int_t i=0; i<MC; i++){
  //  if(MC>15){
  //  std::cout<<i<<"\t \t"<<MC_Hit_X[i]<<"\t \t"<<MC_Hit_Y[i]<<"\t \t"<<MC_Hit_Z[i]<<"\t \t"<<MC_Hit_P[i]<<"\t \t"<<MC_Hit_PathLenght[i]<<"\t \t"<<MC_Hit_Particle_P[i]<<"\t \t"<<MC_Hit_Energy[i]<<std::endl;
  //}
  //}

  //CheckMCHITS(MCHits); //should remove duplicates of MCHits in the list?
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
  //if(MC_ass<14){
  bool doFit=false;
  if(doFit){
    for(int j= 0; j<9;j++){
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
      if(!fit_parabolaXZ.solve()){ std::cout<<"Not Able to Fit XZ with a parabola for a MCHit list of size "<<MC_ass<<"\t and Momentum "<<P<<std::endl; ok=false;}
      if(!fit_LineY.solve()){std::cout<<"Not Able to Fit Y with a Line for a MCHit list of size"<<MC_ass<<"\t and Momentum "<<P<<std::endl; ok=false;;}
      if(!fit_CubicXZ.solve()){std::cout<<"Not Able to Fit Cubic XZ for a MCHit list of size"<<MC_ass<<"\t and Momentum "<<P<<std::endl; ok=false;}
      if(!fit_parabolaY.solve()){std::cout<<"Not Able to Fit Parabolic YZ for a MCHIt list of size"<<MC_ass<<"\t and Momentum "<<P<<std::endl; ok=false;}
      if(!ok) return kTRUE;
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
    //  std::cout<<"Failed to decompose matrix"<<std::endl;
      Fit = false;
    }
    if(decomp) Fit = true;
    decomp.Solve(rhs);
    solution[0]+=rhs[0];
    solution[1]+=rhs[1];
    solution[2]+=rhs[2];
  }
  //if(!Fit) std::cout<<"Fit Failed"<<std::endl;
  if(Fit){
  //std::cout<<"ax"<<setw(20)<<"bx"<<setw(20)<<"cx"<<setw(20)<<std::endl;
  //std::cout<<solution[0]<<setw(20)<<solution[1]<<setw(20)<<solution[2]<<std::endl;
  }
  Chi2_XZDRATIO = 0;
  for( int i = 0;i<MC_ass;i++){
    float dz = MCHit_Assoc_Z[i]-m_zReference;
    float dist = MCHit_Assoc_X[i]-(solution[0]+solution[1]*dz+solution[2]*dz*dz*(1-dRatio*dz));
    Chi2_XZDRATIO+= (dist*dist)/(std::pow(0.100,2));
  }
  Chi2_XZDRATIO = Chi2_XZDRATIO/(MC_ass-4);
  //std::cout<<"Chi2 = "<<Chi2_XZDRATIO;
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
    //  std::cout<<"Failed to decompose matrix"<<std::endl;
      Fit1 = false;
    }
    if(decomp) Fit1 = true;
    decomp.Solve(rhs1);
    solution1[0]+=rhs1[0];
    solution1[1]+=rhs1[1];
    solution1[2]+=rhs1[2];
  }
  //if(!Fit) std::cout<<"Fit Failed"<<std::endl;
  if(Fit){
  //std::cout<<"ax"<<setw(20)<<"bx"<<setw(20)<<"cx"<<setw(20)<<std::endl;
  //std::cout<<solution[0]<<setw(20)<<solution[1]<<setw(20)<<solution[2]<<std::endl;
  }
  Chi2_XZDRATIOFixed = 0;
  for( int i = 0;i<MC_ass;i++){
    float dz1 = MCHit_Assoc_Z[i]-m_zReference;
    float dist1 = MCHit_Assoc_X[i]-(solution1[0]+solution1[1]*dz1+solution1[2]*dz1*dz1*(1-dRatio*dz1));
    Chi2_XZDRATIOFixed+= (dist1*dist1)/(std::pow(0.100,2));
  }
  Chi2_XZDRATIOFixed = Chi2_XZDRATIOFixed/(MC_ass-4);
  //std::cout<<"Chi2 XZDRATIOFIXED \t"<<Chi2_XZDRATIOFixed<<std::endl;
  //std::cout<<"Chi2 = "<<Chi2_XZDRATIO;


  //Study the ySlope business
  //For each track let's plot the distribution of the beta
  //I need to fit the X projection Before...load the PrHits X position

  //FitXProjection(track);
  //Direct Track Study for Fit
  Int_t nX=0;
  Int_t zone = 0;
  NClonesX=0;
  if(nXUp>=4)  zone = 1;
  if(nxDown>=4) zone = 0;
  std::vector<PatHit> Hits;
  std::vector<PatHit> HitsUVT1;
  std::vector<PatHit> HitsUVT2;
  std::vector<PatHit> HitsUVT3;
  //NClonesUV=0;
  for(Int_t i=0;i<PrHit;i++){
    if(PrHit_isX[i])
    {
      nX++;
      if(i>0 && (PrHit_Zat0[i] - PrHit_Zat0[i-1])<10 ){ //here you have a copy
        NClonesX++;
      }
      PatHit hit = PatHit();
      hit.setHit(PrHit_Xat0[i],PrHit_Zat0[i],PrHit_dxDy[i],PrHit_dzDy[i],std::sqrt(PrHit_w2[i]),PrHit_yMin[i],PrHit_yMax[i],PrHit_zone[i],PrHit_planeCode[i],PrHit_isX[i],PrHit_LHCbID[i]);
      //Select only Hits with charge ?
      Hits.push_back(hit);
    }
    if(!PrHit_isX[i])
    {
      //std::cout<<"i"<<setw(20)<<"1/w"<<setw(20)<<"Xat0"<<setw(20)<<"Zat0"<<setw(20)<<"planeCode"<<setw(20)<<"zone"<<setw(20)<<"dxDy"<<setw(20)<<"dzDy"<<setw(20)<<"Fraction"<<setw(20)<<"Size"<<setw(20)<<"Charge \n"<<std::scientific;
      if(PrHit_planeCode[i] < 4 )
      {
        PatHit hitT1UV = PatHit();
        hitT1UV.setHit(PrHit_Xat0[i],PrHit_Zat0[i],PrHit_dxDy[i],PrHit_dzDy[i],std::sqrt(PrHit_w2[i]),PrHit_yMin[i],PrHit_yMax[i],PrHit_zone[i],PrHit_planeCode[i],PrHit_isX[i],PrHit_LHCbID[i]);
        HitsUVT1.push_back(hitT1UV);
        //std::cout<<i<<setw(20)<<sqrt(1./HitsUVT1[i].w2())<<setw(20)<<HitsUVT1[i].x(0.)<<setw(20)<<HitsUVT1[i].z(0.)<<setw(20)<<HitsUVT1[i].planeCode()<<setw(20)<<HitsUVT1[i].zone()<<setw(20)<<HitsUVT1[i].dxDy()<<setw(20)<<HitsUVT1[i].dzDy()<<setw(20)<<HitsUVT1[i].cluster().charge()<<setw(20)<<HitsUVT1[i].cluster().fraction()<<setw(20)<< HitsUVT1[i].cluster().size() <<setw(20)<<HitsUVT1[i].cluster().charge()<<std::endl;
      }
      if(PrHit_planeCode[i]>3 && PrHit_planeCode[i] <8){
        PatHit hitT2UV = PatHit();
        hitT2UV.setHit(PrHit_Xat0[i],PrHit_Zat0[i],PrHit_dxDy[i],PrHit_dzDy[i],std::sqrt(PrHit_w2[i]),PrHit_yMin[i],PrHit_yMax[i],PrHit_zone[i],PrHit_planeCode[i],PrHit_isX[i],PrHit_LHCbID[i]);
        HitsUVT2.push_back(hitT2UV);
      }
      if(PrHit_planeCode[i]>8){
        PatHit hitT3UV = PatHit();
        hitT3UV.setHit(PrHit_Xat0[i],PrHit_Zat0[i],PrHit_dxDy[i],PrHit_dzDy[i],std::sqrt(PrHit_w2[i]),PrHit_yMin[i],PrHit_yMax[i],PrHit_zone[i],PrHit_planeCode[i],PrHit_isX[i],PrHit_LHCbID[i]);
        HitsUVT3.push_back(hitT3UV);
      }
    }
  }
  for(Int_t i=0;i<HitsUVT1.size();i++){
  //std::cout<<i<<setw(20)<<sqrt(1./HitsUVT1[i].w2())<<setw(20)<<HitsUVT1[i].x(0.)<<setw(20)<<HitsUVT1[i].z(0.)<<setw(20)<<HitsUVT1[i].planeCode()<<setw(20)<<HitsUVT1[i].zone()<<setw(20)<<HitsUVT1[i].dxDy()<<setw(20)<<HitsUVT1[i].dzDy()<<setw(20)<<HitsUVT1[i].cluster().charge()<<setw(20)<<HitsUVT1[i].cluster().fraction()<<setw(20)<< HitsUVT1[i].cluster().size() <<setw(20)<<HitsUVT1[i].cluster().charge()<<std::endl;
  }
  PrSeedTrack xProj = PrSeedTrack(zone,m_zReference,Hits);

  if(debug) xProj.PrintHits();

  //xProj.PrintHits();
  fitXProjection(xProj);
  if(NClonesX>0 && !isElectron && Eta_in25 && isLong && P>5000){
    //xProj.PrintTrack();
  }
  if(debug){
    //cout<<"ax"<<setw(10)<<"bx"<<setw(10)<<"cx \n"<<endl;
    //cout<<xProj.ax()<<setw(10)<<xProj.bx()<<setw(10)<<xProj.cx()<<endl;
    if(debug){
    xProj.PrintTrack();}
    //Here the study for the Fitting on PrHit
    //std::cout<<"Chi2"<<setw(20)<<"Chi2 per DoF"<<std::endl;
    //std::cout<<xProj.Chi2()<<setw(20)<<xProj.Chi2()/((double)xProj.hits().size()-3);
  }
  ax_XStepHitFit=xProj.ax();
  bx_XStepHitFit=xProj.bx();
  cx_XStepHitFit=xProj.cx();
  Chi2_XStepFit=xProj.Chi2();
  NClustersX = xProj.hits().size();
  Chi2_perDoF_XStepFit = xProj.Chi2()/((double)xProj.hits().size()-3.);

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
    //std::cout<<"In T1 i found possible pairs from a track which has part =     "<<xProj.zone()<<std::endl;
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
        //std::cout<<"U Hit:   PlaneCode"<<setw(20)<<HitsUVT1[i].planeCode()<<setw(20)<<"dxDy"<<HitsUVT1[i].dxDy()<<setw(20)<<"X at 0"<<HitsUVT1[i].x(0.)<<setw(20)<<"Z at 0"<<setw(20)<<HitsUVT1[i].z(0.)<<setw(20)<<"Yexp"<<setw(20)<<yT1U<<std::endl;
        //std::cout<<"V Hit:   PlaneCode"<<setw(20)<<HitsUVT1[j].planeCode()<<setw(20)<<"dxDy"<<HitsUVT1[j].dxDy()<<setw(20)<<"X at 0"<<HitsUVT1[j].x(0.)<<setw(20)<<"Z at 0"<<setw(20)<<HitsUVT1[j].z(0.)<<setw(20)<<"Yexp"<<setw(20)<<yT1V<<std::endl;

        DeltaXCombinedT1->Fill(XCombined-XPredicted);
        DeltaYCombT1.push_back(yT1V-yT1VProj);
        DeltaXCombT1.push_back(XCombined-XPredicted);
        if(std::abs(XCombined-XPredicted)<BestDeltaCom)
        {
        BestDeltaCom = std::abs(XCombined-XPredicted);
        bestHitPair = std::make_pair(HitsUVT1[i], HitsUVT1[j]);
        }
      }
    //std::cout<<"X Combination"<<setw(20)<<XCombined<<std::endl;
    //std::cout<<"X Predicted  "<<setw(20)<<xProj.x( (HitsUVT1[i].z(0.)+ HitsUVT1[j].z(0.) )/2.)<<std::endl;
    }
  }
  DeltaXBestCombCorrSlopeT1 =-999;
  DeltaXBestCombT1 = -999;
  TyPairT1 = -999.;
  YPairT1 = -999.;
  XPairT1 = -999.;
  YPairT1 = -999.;
  YUT2Projected = -999.;
  YVT2Projected = -999.;
  DeltaUT2Projected = -999.;
  DeltaVT2Projected = -999.;
  std::cout<<"\n \n \n"<<setw(20)<<"New Track"<<std::endl;
  if(nT1U>0 && nT1V>0)
  {
    //std::cout<<"Sorted Combinations"<<std::endl;
    std::sort(DeltaXCombT1.begin(),DeltaXCombT1.end(),[](Double_t Value1, Double_t Value2)->bool{return std::abs(Value1)<std::abs(Value2);});
    for(Int_t i = 0 ; i<DeltaXCombT1.size();i++)
    {
      if(DeltaXCombT1.size()>1)
      {
        //std::cout<<"i"<<setw(20)<<"Delta X Best Combination"<<setw(20)<<DeltaXCombT1[i]<<std::endl;
      }
    }
    DeltaXBestCombT1 = DeltaXCombT1[0];
    //std::cout<<"DeltaBestCombT1   "<<setw(20)<<DeltaXBestCombT1<<std::endl;
    //std::cout<<"Compute Delta Best"<<setw(20)<<(bestHitPair.first.x(0.)+bestHitPair.second.x(0.))/2 - xProj.x( (bestHitPair.first.z(0.) + bestHitPair.second.z(0.)) /2.)<<std::endl;
    PatHit bestUT1 = bestHitPair.first;
    PatHit bestVT1 = bestHitPair.second;
    Double_t zMid = (bestUT1.z() + bestVT1.z())/2. ;
    Double_t xMidFromHit = (bestUT1.x() + bestVT1.x()) /2. ;//+ xProj.xSlope(bestUT1.z())*(zMid - bestUT1.z());
    DeltaXBestCombCorrSlopeT1 = xMidFromHit - xProj.x(zMid);
    Double_t yT1U = (bestUT1.x(0.)-xProj.x(bestUT1.z(0.)))/bestUT1.dxDy();
    Double_t yT1V = (bestVT1.x(0.)-xProj.x(bestVT1.z(0.)))/bestVT1.dxDy();
    Double_t yAvg = (yT1U+yT1V)/2.;
    Double_t TY = (yT1V-yT1U)/(bestVT1.z(0.) - bestUT1.z(0.));
    TyPairT1 = TY;
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
        std::cout<<"T2  -- U -- i ="<<setw(20)<<i<<setw(20)<<"Delta U T2 expected best combination in T1"<<setw(20)<<uExp - HitsUVT2[i].x(0.)<<std::endl;
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
        std::cout<<"T2  -- V -- i"<<setw(20)<<i<<"Delta V T2 expected best combination in T1"<<setw(20)<<vExp - HitsUVT2[i].x(0.)<<std::endl;
        Double_t bestDeltaVTProjected =vExp - HitsUVT2[i].x(0.) ;
      }
    }
    DeltaXT2PairsCombined =-999.;
    if(nT2V>0 && nT2U>0){
      std::cout<<"I found possible pairs in UV T2"<<std::endl;
      Double_t BestDeltaCom = 100.;
      Int_t pairN = 0;
      for(Int_t i = 0; i< nT2U;i++ )
      {
        for(Int_t j=nT2U; j< nT2U+ nT2V;j++)
        {
          pairN++;
          Double_t xExpT2Middle = (HitsUVT2[i].x(0.) + HitsUVT2[j].x(0.))/2. ;
          Double_t DeltaT2UVX = xExpT2Middle- xProj.x( (HitsUVT2[i].z(0.)+HitsUVT2[j].z(0.))/2. );
          std::cout <<"Pair"<<setw(20)<<pairN<<setw(20)<<"Delta In T2 after Projection" <<setw(20)<<DeltaT2UVX<<std::endl;
          DeltaXT2PairsCombined = DeltaT2UVX;
        }
      }
    }
  }
  //if()

  t1->Fill();
  return kTRUE;
}


Bool_t TrackStudy::fitXProjection(PrSeedTrack & track)
{
  float mat[6];
  float rhs[3];
  std::vector<PatHit> Hits = track.hits();

  //std::cout<<"Hits Size"<<Hits.size()<<std::endl;
  //std::cout<<"Hits Size from track"<<track.hits().size();
  float dRatio=0;
  for(int loop = 0;3>loop;++loop){

    std::fill(mat,mat+6,0.);
    std::fill(rhs,rhs+3,0.);
    //std::cout<<"******************* Fit Loop "<<loop<<"*********************"<<std::endl;
    for( int i=0;i<Hits.size();i++){
      const float w = Hits[i].w2();
      //std::cout<<"Hit w2"<<Hits[i].w2()<<std::endl;
      const float dz = Hits[i].z(0.)-m_zReference;
      float deta = dz*dz;
      //if(m_usedRatio){
      //deta = dz*dz*(1-dRatio*dz);
      //}
      float dist = track.distance(Hits[i]);
      //std::cout<<"i"<<setw(20)<<i<<setw(20)<<"Hit X"<<setw(20)<<Hits[i].x(0)<<setw(20)<<"Track x"<<setw(20)<<track.x(Hits[i].z(0.))<<setw(20)<<"Distance"<< track.distance(Hits[i])<<std::endl;
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
  //std::cout<<"Long Tracks"<<std::endl;
  //std::cout<<"Number of Geant MCHit \t"<<m_long_NGeantHit<<std::endl;
  //std::cout<<"Number of different Layers with Geant MCHit \t"<<m_long_NLayer_Geant<<std::endl;

  //std::cout<<"Number of MCHit into Cluster \t"<<m_long_NMCHitIntoCluster<<std::endl;
  //std::cout<<"Number of different Layers with MCHit into Cluster\t"<<m_long_NLayer_MCHitInCluster<<std::endl;

  //std::cout<<"Number of Clusters \t" <<m_long_NPrHit<<std::endl;
  //std::cout<<"Number of diffrent Layers Clusters \t" <<m_long_NLayer_PrHit<<std::endl;
  zMCHitClones->Draw();
  zMCHitClonevsP->Draw("colz");
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

}
// void TrackStudy::WhichPlane(MCHit &hit){
// Int_t zone=0;
// if(hit->y()>0) zone = 1;
// if(hit->z() > 7840. && hit->z() < 7870.) hit->setLayer(0,zone,true);
// if(hit->z() > 8020. && hit->z() < 8055.) hit->setLayer(3,zone,true);
// if(hit->z() > 8520. && hit->z() < 8555.) hit->setLayer(4,zone,true);
// if(hit->z() > 8705. && hit->z() < 8740.) hit->setLayer(7,zone,true);
// if(hit->z() > 9200. && hit->z() < 9250.) hit->setLayer(8,zone,true);
// if(hit->z() > 9390. && hit->z() < 9430.) hit->setLayer(11,zone,true);
//
// if(hit->z() >7890. && hit->z() < )
//
//
// }
