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

#include "TrackStudy.h"
#include <TH2.h>
#include <TStyle.h>
#include "LinParFit.h"
#include "Track.h"
#include "PatHit.h"
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
Double_t Chi2_ParabolaXZ;
Double_t ax_par;
Double_t bx_par;
Double_t cx_par;
//Cubic XZ
Double_t Chi2_CubicXZ;
Double_t ax_cub;
Double_t bx_cub;
Double_t cx_cub;
Double_t dx_cub;
//Line Y
Double_t Chi2_LineY;
Double_t ay_line;
Double_t by_line;
Double_t ay_par;
Double_t by_par;
Double_t cy_par;
//Parabolic Y
Double_t Chi2_ParabolaY;

//Track Flags
Double_t Momentum;
Double_t Px;
Double_t Py;
Double_t Track_eta;
Double_t Pt;
Double_t Track_Phi;
Double_t Track_Pt;
//-------------------------------------------------------------

//NHIts study
//Glancing of Z position
TH1D *zMCHitClones = new TH1D("zMCHitClones","zMCHitClones;z[mm];Counts",1000,7600,9500);
TH2D *zMCHitClonevsP = new TH2D("zMCHitClones_Vs_P","zMCHit Clones vs P;z[mm];P[MeV]",200,7600,9500,200,1000,12000);
TH2D *xvsyClones = new TH2D("xVsy_Clones","x vs y Clones;x[mm];y[mm]",300,-3000,3000,300,-2500,2500);
//TH3D *MCHitCloneXYZ = new TH3D("X Y Z MCHitClone","")

void TrackStudy::Begin(TTree * /*tree*/)
{

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
  //std::cout<<"Ciao"<<std::endl;
  // HERE DEBUG AND PRINT AT VIDEO OF HIT CONTENT
  Double_t PCut =3000;
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
  // if(isElectron) return kTRUE; //No Electrons
  // if(!Eta_in25) return kTRUE; //Eta Cut
  // //if(P>5000) return kTRUE;
  // if(!isSeed) return kTRUE;
  // if(P<2000) return kTRUE; //Study done for P>2000 !electron & LongTracks
  if(MC_Ovtx_z>4000) return kTRUE; //Study done for Tracks not generated in interactions (naively)
  //std::vector<PatHit> track(100); //Create the Vector of Hit to be analysed in terms of hit content (remove some of them? (like duplicates in layers?))
  std::vector<MCHit> track_MC(100); //Vector of MCHit going into cluster
  std::vector<MCHit> Particle_GeantHit(100);
  //Create a vector of hits to be placed in the PrSeedTrack to be analysed: check if it's the case to remove clones
  Bool_t printHit = false;
  if(printHit){
    std::cout<<"i \t Xat0 \t Zat0 \t planeCode \t zone \t dxDy \t dzDy"<<"\t Fraction"<<"\t Size"<<"\t Charge"<<std::endl;}
    for (Int_t i = 0 ;  i< CheatedSeeding_NHits; i++){
      PatHit hit = PatHit();
      hit.setHit(PrHit_Xat0[i],PrHit_Zat0[i],PrHit_dxDy[i],PrHit_dzDy[i],std::sqrt(PrHit_w2[i]),PrHit_yMin[i],PrHit_yMax[i],PrHit_zone[i],PrHit_planeCode[i],PrHit_isX[i],PrHit_LHCbID[i]);
      //hit.setCluser(); Need the corresponding branches in the ntuple
      if(printHit){
        std::cout<<i<<"\t"<<hit.x(0.)<<"\t"<<hit.z(0.)<<"\t"<<hit.planeCode()<<"\t"<<hit.cluster().charge()<<"\t"<<hit.cluster().fraction()<<"\t"<<hit.cluster().size()<<"\t"<<hit.cluster().charge()<<std::endl;
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
  bool doFit=true;
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
    Chi2_ParabolaXZ+= std::pow( ((double)MCHit_Assoc_X[i]- (solution_xz_par[0]+ solution_xz_par[1]*dz+ solution_xz_par[2]*dz*dz))/0.100,2);
    Chi2_CubicXZ+= std::pow( ((double)MCHit_Assoc_X[i]- (solution_xz_cubic[0]+ solution_xz_cubic[1]*dz+ solution_xz_cubic[2]*dz*dz +solution_xz_cubic[3]*dz*dz*dz))/0.100,2);
    Chi2_LineY+=std::pow( ((double)MCHit_Assoc_Y[i]-(solution_yz_line[0]+solution_yz_line[1]*dz))/1.14,2);
    Chi2_ParabolaY+=std::pow( ((double)MCHit_Assoc_Y[i]-(solution_yz_par[0]+solution_yz_par[1]*dz+solution_yz_par[2]*dz*dz))/1.14,2);
  }

  Momentum  = P;
  t1->Fill();

  //Study the ySlope business
  //For each track let's plot the distribution of the beta
  //I need to fit the X projection Before...load the PrHits X position


  //FitXProjection(track);
  //Direct Track Study for Fit
  Int_t zone = 0;
  if(nXUp>=4)  zone = 1;
  if(nxDown>=4) zone = 0;
  PrSeedTrack track = PrSeedTrack(zone,m_zReference);
  for(Int_t i=0;i<PrHit;i++){
    PatHit hit = PatHit();
    hit.setHit(PrHit_Xat0[i],PrHit_Zat0[i],PrHit_dxDy[i],PrHit_dzDy[i],std::sqrt(PrHit_w2[i]),PrHit_yMin[i],PrHit_yMax[i],PrHit_zone[i],PrHit_planeCode[i],PrHit_isX[i],PrHit_LHCbID[i]);
    //Select only Hits with charge ?
    track.addHit(hit);
  }
  track.PrintHits();
  //Here the study for the Fitting on PrHit



  return kTRUE;
}

// //void TrackStudy::fitXProjection(){

//   double mat[6];
//   float rhs[3];
//   //std::fill(rhs,rhs+3,0);
//   std::vector<Hit> Hits = track->hits();

//   for(int loop = 0;3>loop;++loop)
//   {
//     std::fill(mat,mat+6,0.);
//     std::fill(rhs,rhs+3,0.);
//     for( int i=0; i < Hits.size() ;i++ )
//     {
//       const float w =  Hits[i].w2();//squared
//       //std::cout<<"W\t"<<w<<std::endl;
//       const float dz = Hits[i].GetZ() - m_m_zReference;
//       float deta = 0;
//       deta = dz*dz*(1-m_dRatio*dz);
//       Hit *hit = new Hit(Hits[i].GetX(),Hits[i].GetY(),Hits[i].GetZ());
//       float dist = track->distance(hit);
//       //always()<<"Loop \t"<<loop<<"\n Distance From Hit \t"<<dist<<endmsg;
//       // if(loop>0)
//       // dist = track.distance( *itH ); //try the effect
//       mat[0]+= w     ;
//       mat[1]+= w * dz;
//       mat[2]+= w * dz * dz;
//       mat[3]+= w * deta;
//       mat[4]+= w * dz * deta;
//       mat[5]+= w * deta * deta;
//       //right hand side
//       rhs[0]+= w * dist;
//       rhs[1]+= w * dist * dz;
//       rhs[2]+= w * dist * deta;
//     }
//   }
//}
Bool_t TrackStudy::fitXProjection(PrSeedTrack * track){
  float mat[6];
  float rhs[3];
  std::vector<PatHit> Hits = track->hits();
  float dRatio=0;
  for(int loop = 0;3>loop;++loop){
    std::fill(mat,mat+6,0.);
    std::fill(rhs,rhs+3,0.);
    for( int i=0;i<Hits.size();i++){
      const float w = Hits[i].w2();
      std::cout<<"Hit w2"<<Hits[i].w2()<<std::endl;
      const float dz = Hits[i].z(0.)-m_zReference;
      float deta = 0;
      //if(m_usedRatio){
      //deta = dz*dz*(1-dRatio*dz);
      //}
      float dist = track->distance(Hits[i]);

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
      //std::cout<<"Failed to decompose matrix"<<std::endl;
      return false;
    }
    decomp.Solve(rhs);
    track->updateParameters(rhs[0],rhs[1],rhs[2],0.,0.);
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
