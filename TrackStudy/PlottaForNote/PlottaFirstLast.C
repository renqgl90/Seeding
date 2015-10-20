#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
void CenterTitles( TH2D * h2 );
void CenterTitles( TH2D * h1 );
void PlottaFirstLast(){
  gROOT->ProcessLine(".x ./style.C");
  gStyle->SetOptStat(0);
  // gStyle->SetTitleX(0.5);
  // gStyle->SetTitleAlign(23); 
  TFile *f = new TFile("Analysed_Track.root","READ");
  if(f==0){
    std::cout<<"Error Opening File"<<std::endl;
    return;
  }
  TTree *t  = (TTree*)f->Get("treee");

  Double_t alpha[3];
  alpha[0] = 120.24;
  alpha[1] = 510.0;
  alpha[2] = 730.6;
  Double_t PCut[3];
  PCut[0] = 5000.;
  PCut[1] = 2000.;
  PCut[2] = 1500.;
  
  if( t==0){
    std::cout<<"Error Loading the tree"<<std::endl;
  }
  //Case 0 Variables L0_Selections
  Double_t  Delta_Case0, tx_inf_Case0, Delta_Case1,tx_inf_Case1,Delta_Case2, tx_inf_Case2;
  Double_t  Track_P;
  Bool_t Case0Accepted,Case1Accepted,Case2Accepted;
  Int_t Nev = t->GetEntries();
  
  t->SetBranchAddress("Del_x_inf_Case0",&Delta_Case0);
  t->SetBranchAddress("tx_inf_Case0",&tx_inf_Case0);
  t->SetBranchAddress("Case0Accepted",&Case0Accepted);
  t->SetBranchAddress("Track_P",&Track_P);
  
  t->SetBranchAddress("Del_x_inf_Case1",&Delta_Case1);
  t->SetBranchAddress("tx_inf_Case1",&tx_inf_Case1);
  t->SetBranchAddress("Case1Accepted",&Case1Accepted);
  
  t->SetBranchAddress("Del_x_inf_Case2",&Delta_Case2);
  t->SetBranchAddress("tx_inf_Case2",&tx_inf_Case2);
  t->SetBranchAddress("Case2Accepted",&Case2Accepted); 
  
  TH2D * Case0_FirstLast_AllP = new TH2D("Case1_AllP","Case 1 : all P;t_{x}^{#infty};#Delta = x^{last}_{true} - x^{last}_{predicted} [mm]", 1000, -0.4,0.4 , 2000, -2000,2000);
  TH2D * Case0_FirstLast_PMore5 = new TH2D("Case1_Pmore5","Case 1 : P > 5 GeV/c;t_{x}^{#infty}; #Delta = x^{last}_{true} - x^{last}_{predicted} [mm]", 1000, -0.4,0.4 , 2000, -2000,2000);
  TH2D * Case0_FirstLast_PMore5Rot = new TH2D("Case1_Pmore5Rot","Case 1 : P > 5 GeV/c;t_{x}^{#infty}; #Delta - t_{x}^{#infty} * L0_alphaCorr [mm]", 1000, -0.4, 0.4, 1000, -500,500);


  

  TH2D * Case1_FirstLast_AllP = new TH2D("Case2_AllP","Case2 : all P;t_{x}^{#infty};#Delta = x^{last}_{true} - x^{last}_{predicted} [mm]", 1000, -0.4,0.4 , 2000, -2000,2000);
  TH2D * Case1_FirstLast_PMore5 = new TH2D("Case2_Pmore3","Case2 P > 2 GeV/c;t_{x}^{#infty}; #Delta = x^{last}_{true} - x^{last}_{predicted} [mm]", 1000, -0.4,0.4 , 2000, -2000,2000);
  TH2D * Case1_FirstLast_PMore5Rot = new TH2D("Case2_Pmore5Rot","Case 2 : P > 2 GeV/c;t_{x}^{#infty};#Delta -  t_{x}^{#infty} * L0_alphaCorr [mm]", 1000, -0.4,0.4, 1000, -700,700);

  
  TH2D * Case2_FirstLast_AllP = new TH2D("Case3_AllP","Case 3 : all P;t_{x}^{#infty};#Delta = x^{last}_{true} - x^{last}_{predicted} [mm]", 1000, -0.4,0.4 , 2000, -2000,2000);
  TH2D * Case2_FirstLast_PMore5 = new TH2D("Case3_Pmore3","Case 3 : P > 1.5 GeV/c;t_{x}^{#infty};#Delta = x^{last}_{true} - x^{last}_{predicted} [mm]", 1000, -0.4,0.4 , 2000, -2000,2000); 
  TH2D * Case2_FirstLast_PMore5Rot = new TH2D("Case3_Pmore5Rot","Case 3 : P > 1.5 GeV/c;t_{x}^{#infty};#Delta - t_{x}^{#infty} * L0_alphaCorr [mm]", 1000, -0.4,0.4, 1000, -1200,1200);
  
  for(int i = 0; i< Nev ;i++){
    t->GetEntry(i);
    if(Case0Accepted){
      // std::cout<<"Deta_Case 0 = "<<Delta_Case0<<"\t tx_inf_Case0  "<<tx_inf_Case0<<"\t P "<<Track_P<<"Is Accepted "<<Case0Accepted << std::endl;
      Case0_FirstLast_AllP->Fill(tx_inf_Case0,Delta_Case0);
      if(Track_P>PCut[0]){
        Case0_FirstLast_PMore5->Fill(tx_inf_Case0,Delta_Case0);
        Case0_FirstLast_PMore5Rot->Fill(tx_inf_Case0,Delta_Case0-alpha[0]*tx_inf_Case0);
      }
    }
    if(Case1Accepted){
      // std::cout<<"Deta_Case 1 = "<<Delta_Case1<<"\t tx_inf_Case1  "<<tx_inf_Case1<<"\t P "<<Track_P<<"Is Accepted "<<Case1Accepted << std::endl;
      Case1_FirstLast_AllP->Fill(tx_inf_Case1,Delta_Case1);
      if(Track_P>PCut[1]){
        Case1_FirstLast_PMore5->Fill(tx_inf_Case1,Delta_Case1);
        Case1_FirstLast_PMore5Rot->Fill(tx_inf_Case1,Delta_Case1-alpha[1]*tx_inf_Case1);
      }
    }
    if(Case2Accepted){
      // std::cout<<"Deta_Case 2 = "<<Delta_Case2<<"\t tx_inf_Case2  "<<tx_inf_Case2<<"\t P "<<Track_P<<"Is Accepted "<<Case2Accepted << std::endl;
      Case2_FirstLast_AllP->Fill(tx_inf_Case2,Delta_Case2);
      if(Track_P>PCut[2]){
        Case2_FirstLast_PMore5->Fill(tx_inf_Case2,Delta_Case2);
        Case2_FirstLast_PMore5Rot->Fill(tx_inf_Case2,Delta_Case2 - alpha[2]*tx_inf_Case2);
      }
    }
  }
  
  TCanvas * c1  = new TCanvas();
  c1->Divide(3,1);


  c1->cd(1);
  CenterTitles(Case0_FirstLast_AllP);
  Case0_FirstLast_AllP->Draw("colz");


  c1->cd(2);
  CenterTitles(Case1_FirstLast_AllP);
  Case1_FirstLast_AllP->Draw("colz");


  c1->cd(3);
  CenterTitles(Case2_FirstLast_AllP);
  Case2_FirstLast_AllP->Draw("colz");

  TCanvas * c2 = new TCanvas();
  c2->Divide(3,1);
  c2->cd(1);
  CenterTitles(Case0_FirstLast_PMore5);
  Case0_FirstLast_PMore5->Draw("colz");


  c2->cd(2);
  CenterTitles(Case1_FirstLast_PMore5);
  Case1_FirstLast_PMore5->Draw("colz");


  c2->cd(3);
  CenterTitles(Case2_FirstLast_PMore5);
  Case2_FirstLast_PMore5->Draw("colz");

  TCanvas * c3 = new TCanvas();
  c3->Divide(3,1);
  c3->cd(1);
  CenterTitles(Case0_FirstLast_PMore5Rot);
  Case0_FirstLast_PMore5Rot->Draw("colz");
  
  c3->cd(2);
  CenterTitles(Case1_FirstLast_PMore5Rot);
  Case1_FirstLast_PMore5Rot->Draw("colz");
  c3->cd(3);
  CenterTitles(Case2_FirstLast_PMore5Rot);
  Case2_FirstLast_PMore5Rot->Draw("colz");
  
}
  
void CenterTitles( TH2D * h2){
  h2->GetXaxis()->CenterTitle();
  h2->GetYaxis()->CenterTitle();
}
