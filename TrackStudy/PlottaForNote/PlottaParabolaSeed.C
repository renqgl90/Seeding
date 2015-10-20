#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
using namespace std;
void CenterTitles( TH2D * h2 );
void CenterTitle( TH1D * h1);
void PlottaParabolaSeed(){
  
  gROOT->ProcessLine(".x ./style.C");
  gStyle->SetOptStat(0);
  TFile *f = new TFile("Analysed_Track.root","READ");
  if(f==0){
    std::cout<<"Error Opening File"<<std::endl;
    return;
  }
  TTree *t  = (TTree*)f->Get("treee");

  Double_t x0Corr[3];
  x0Corr[0] = 0.002152;
  x0Corr[1] = 0.001534;
  x0Corr[2] = 0.001834;
  Double_t PCut[3];
  PCut[0] = 5000.;
  PCut[1] = 2000.;
  PCut[2] = 1500.;

  Bool_t CaseAccepted[3];
  Double_t x0[3];
  Double_t DelSeed1[3];
  Double_t DelSeed2[3];
  t->SetBranchAddress("Case0Accepted",&CaseAccepted[0]);
  t->SetBranchAddress("Case1Accepted",&CaseAccepted[1]);
  t->SetBranchAddress("Case2Accepted",&CaseAccepted[2]);
  t->SetBranchAddress("x0_Case0",&x0[0]);
  t->SetBranchAddress("x0_Case1",&x0[1]);
  t->SetBranchAddress("x0_Case2",&x0[2]);
  t->SetBranchAddress("DelxProjectedSeed1_Case0",&DelSeed1[0]);
  t->SetBranchAddress("DelxProjectedSeed2_Case0",&DelSeed2[0]);
  t->SetBranchAddress("DelxProjectedSeed1_Case1",&DelSeed1[1]);
  t->SetBranchAddress("DelxProjectedSeed2_Case1",&DelSeed2[1]);
  t->SetBranchAddress("DelxProjectedSeed1_Case2",&DelSeed1[2]);
  t->SetBranchAddress("DelxProjectedSeed2_Case2",&DelSeed2[2]);

  Double_t Track_P;
  TH2D * Case[3];
  TH1D * CaseOneD[3];
  Case[0] = new TH2D("DeltaSeed0","Case 0 , p > 5 GeV ; x_{0} [mm] ; #Delta^{Seed}_{Corrected} [ mm ]",1000,-2000,2000,300,-5,+5);
  Case[1] = new TH2D("DeltaSeed1", "Case 1 , p > 2 GeV ; x_{0} [mm] ; #Delta^{Seed}_{Corrected} [ mm ]",1000,-6500,6500,300,-10,+10);
  Case[2] =  new TH2D("DeltaSeed2","Case 2 , p > 1.5 GeV; x_{0} [mm] ;#Delta^{Seed}_{Corrected} [ mm ]",1000,-6500,6500,300,-15,+15);
  
  CaseOneD[0] = new TH1D("Del0","Case0 , p> 5GeV; #Delta^{Seed}_{Corrected} [mm] ; Counts/200 #mu m", 50, -5,5);
  CaseOneD[1] = new TH1D("Del1","Case1 , p> 2GeV; #Delta^{Seed}_{Corrected} [mm] ; Counts/200 #mu m", 100, -10,10);
  CaseOneD[2] = new TH1D("Del2","Case2 , p> 1.5 GeV; #Delta^{Seed}_{Corrected} [mm] ; Counts/200 #mu m",150,-15,15);
  t->SetBranchAddress("Track_P",&Track_P);
  int nEv = t->GetEntries();
  for( int i = 0; i<nEv; i++){
    t->GetEntry(i);
    for(int j = 0; j<3;j++){
      if(CaseAccepted[j] && Track_P>PCut[j]){
        Double_t DeltaNew1 = -999;
        Double_t DeltaNew2 = -999;
        if(DelSeed1[0]>-99)
        {
          DeltaNew1 = DelSeed1[j]-x0Corr[j]*x0[j];
          Case[j]->Fill(x0[j],DeltaNew1);
          CaseOneD[j]->Fill(DeltaNew1);
        }
        if(DelSeed2[j]>-99){
          DeltaNew2 = DelSeed2[j]-x0Corr[j]*x0[j];
          Case[j]->Fill(x0[j],DeltaNew2);
          CaseOneD[j]->Fill(DeltaNew2);
        }
      }
    }
  }
  TCanvas *c1 = new TCanvas();
  c1->Divide(2,1);
  CenterTitles(Case[0]);
  c1->cd(1);
  Case[0]->Draw("colz");
  c1->cd(2);
  CenterTitle(CaseOneD[0]);
  CaseOneD[0]->Draw();
  
  TCanvas *c2 = new TCanvas();
  c2->Divide(2,1);
  c2->cd(1);
  CenterTitles(Case[1]);
  Case[1]->Draw("colz");
  c2->cd(2);
  CenterTitle(CaseOneD[1]);
  CaseOneD[1]->Draw();
  TCanvas *c3 = new TCanvas();
  c3->Divide(2,1);
  c3->cd(1);
  CenterTitles(Case[2]);
  Case[2]->Draw("colz");
  c3->cd(2);
  CenterTitle(CaseOneD[2]);
  CaseOneD[2]->Draw();
}
void CenterTitles(TH2D * h2){
  h2->GetXaxis()->CenterTitle();
  h2->GetYaxis()->CenterTitle();
  h2->GetXaxis()->SetNdivisions(510);
}
void CenterTitle(TH1D *h1){
  h1->GetXaxis()->CenterTitle();
  h1->GetYaxis()->CenterTitle();
  h1->SetLineColor(kBlue);
  h1->SetLineWidth(1.0);
  h1->GetXaxis()->SetNdivisions(510);
  gPad-> SetLogy();
}
