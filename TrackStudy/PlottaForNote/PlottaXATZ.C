#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
void CenterTitles( TH2D * h2 );
void CenterTitles( TH1D * h1 );


void PlottaXATZ(){
  TFile *f = new TFile("Analysed_Track_XatZAdd.root","READ");
  if(f==0){
    std::cout<<"Error Opening File"<<std::endl;
    return;
  }
  gROOT->ProcessLine(".x ./style.C");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TTree *t  = (TTree*)f->Get("treee");
  // Double_t alpha[3];
  // alpha[0] = 120.24;
  // alpha[1] = 510.0;
  // alpha[2] = 730.6;
  Double_t PCut[3];
  PCut[0] = 5000.;
  PCut[1] = 2000.;
  PCut[2] = 1500.;
  if( t==0){
    std::cout<<"Error Loading the tree"<<std::endl;
  }
  //Case 0 Variables L0_Selections
  Double_t  Delta_Case0, tx_inf_Case0, Delta_Case1,tx_inf_Case1,Delta_Case2, tx_inf_Case2;
  Bool_t CaseAccepted[3];
  Int_t Nev = t->GetEntries();
  Int_t nSeed1[3];
  Int_t nSeed2[3];
  Double_t DelSeed1_Case0[3];
  Double_t DelSeed2_Case0[3];
  Double_t DelSeed1_Case1[3];
  Double_t DelSeed2_Case1[3];
  Double_t DelSeed1_Case2[3];
  Double_t DelSeed2_Case2[3];
  Double_t Track_P;

  TH1D * DeltaSeed[3];
  // TH1D * DeltaSeed[3];
  for(int i = 0; i<3; i++){
    DeltaSeed[i] = new TH1D(Form("DeltaSeed1_Case%i",i),Form("Case %i, p> %.0f MeV; #Delta xAtZ [mm] ; Counts / 20 #mu m", i, PCut[i]),300,-3,3);
  //   DeltaSeed2[i] = new TH1D(Form("DeltaSeed2_Case%i",i),Form("Case %i, p> %.0f MeV , Seed Plane 7; #Delta xAtZ ; Counts;", i, PCut[i]),300,0,5);
  // }
  }
  t->SetBranchAddress("Track_P",&Track_P);

  t->SetBranchAddress("Case0Accepted",&CaseAccepted[0]);
  t->SetBranchAddress("Case1Accepted",&CaseAccepted[1]);  
  t->SetBranchAddress("Case2Accepted",&CaseAccepted[2]);

  t->SetBranchAddress("nSeed1Case0Rem",&nSeed1[0]);
  t->SetBranchAddress("nSeed2Case0Rem",&nSeed2[0]);
  
  t->SetBranchAddress("nSeed1Case1Rem",&nSeed1[1]);
  t->SetBranchAddress("nSeed2Case1Rem",&nSeed2[1]);
    
  t->SetBranchAddress("nSeed1Case2Rem",&nSeed1[2]);
  t->SetBranchAddress("nSeed2Case2Rem",&nSeed2[2]);

  

  t->SetBranchAddress("DelSeed1_Case0Rem",&DelSeed1_Case0);
  t->SetBranchAddress("DelSeed2_Case0Rem",&DelSeed2_Case0);

  t->SetBranchAddress("DelSeed1_Case1Rem",&DelSeed1_Case1);
  t->SetBranchAddress("DelSeed2_Case1Rem",&DelSeed2_Case1);
  
  t->SetBranchAddress("DelSeed1_Case2Rem",&DelSeed1_Case2);
  t->SetBranchAddress("DelSeed2_Case2Rem",&DelSeed2_Case2);

  int nEv = t->GetEntries();
  for(int i = 0; i<nEv; i++){
    t->GetEntry(i);
    for(int j = 0; j<3; j++){ //loop cases
      if(CaseAccepted[j] && Track_P>PCut[j]){
        if(j==0){//Case0
          std::cout<<"Case 0 Accepted "<<std::endl;
          for( int k = 0; k<nSeed1[j]; k++){
            DeltaSeed[j]->Fill(DelSeed1_Case0[k]);
            //std::cout<<"DelSeed1Case0 = "<<DelSeed1_Case0[k]<<std::endl;
          }
          for( int k = 0; k<nSeed2[j]; k++){
            DeltaSeed[j]->Fill(DelSeed2_Case0[k]);
            //std::cout<<"DelSeed1Case0 = "<<DelSeed1_Case0[k]<<std::endl;
          }
        }
        if(j==1) { //Case1
          std::cout<<"Case 0 Accepted "<<std::endl;
          for( int k = 0; k<nSeed1[j]; k++){
            DeltaSeed[j]->Fill(DelSeed1_Case1[k]);
            //std::cout<<"DelSeed1Case0 = "<<DelSeed1_Case0[k]<<std::endl;
          }
          for( int k = 0; k<nSeed2[j]; k++){
             DeltaSeed[j]->Fill(DelSeed2_Case1[k]);
             //std::cout<<"DelSeed1Case0 = "<<DelSeed1_Case0[k]<<std::endl;
          }
        }
        if(j==2) { //Case2
          std::cout<<"Case 0 Accepted "<<std::endl;
          for( int k = 0; k<nSeed1[j]; k++){
            DeltaSeed[j]->Fill(DelSeed1_Case2[k]);
            //std::cout<<"DelSeed1Case0 = "<<DelSeed1_Case0[k]<<std::endl;
          }
          for( int k = 0; k<nSeed2[j]; k++){
            DeltaSeed[j]->Fill(DelSeed2_Case2[k]);
            //std::cout<<"DelSeed1Case0 = "<<DelSeed1_Case0[k]<<std::endl;
          }
        }
      }
    }
  }
  
  TCanvas * c1 = new TCanvas();
  CenterTitles(DeltaSeed[0]);
  DeltaSeed[0]->SetLineColor(kBlack);
  DeltaSeed[0]->SetName("Case 0 p >5 GeV");
  CenterTitles(DeltaSeed[1]);
  DeltaSeed[1]->SetLineColor(kRed);
  DeltaSeed[1]->SetName("Case 1 p >2 GeV");
  CenterTitles(DeltaSeed[2]);
  DeltaSeed[2]->SetLineColor(kBlue);
  DeltaSeed[2]->SetName("Case 2 p >1.5 GeV");

  
  DeltaSeed[2]->Draw();
  DeltaSeed[1]->Draw("same");
  DeltaSeed[0]->Draw("same");
  gPad->BuildLegend();
  //c1->Divide(3,1);
  // for(int i = 0; i<3; i++){
  //   DeltaSeed1[i]->SetLineColor(kBlue);
  //   CenterTitles(DeltaSeed[i]);
  //   DeltaSeed2[i]->SetLineColor(kRed);
  //   CenterTitles(DeltaSeed2[i]);
  //   c1->cd(i+1);
  //   gPad->SetLogy();
  //   DeltaSeed[i]->Draw();
  //   //DeltaSeed2[i]->Draw("same");
  // }
  }
void CenterTitles(TH2D * h2){
  h2->GetXaxis()->CenterTitle();
  h2->GetYaxis()->CenterTitle();
  h2->GetXaxis()->SetNdivisions(510);
}
void CenterTitles(TH1D *h1){
  h1->GetXaxis()->CenterTitle();
  h1->GetYaxis()->CenterTitle();
  h1->SetLineColor(kBlue);
  h1->SetLineWidth(1.0);
  h1->GetXaxis()->SetNdivisions(510);
  //gPad-> SetLogy();
}
