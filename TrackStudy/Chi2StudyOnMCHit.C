#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "iostream"

void Chi2StudyOnMCHit()
{
  Double_t Phi;
  Double_t Pt;
  Double_t P;
  Double_t Eta;
  Double_t Chi2ParXZ;
  Double_t Chi2CubXZ;
  Double_t Chi2LineY;
  Double_t Chi2ParY;
  Double_t Chi2dRatioXZ;
  Double_t Chi2dRatioXZFix;
  using namespace std;
  //Macro to study dRatio(x(zRef) , y(zRef) );
  TFile *file = new TFile("Analysed_Track.root");
  TTree* t1 = (TTree*)file->Get("treee");
  Int_t nEntries = (Int_t)t1->GetEntries();
  cout<<"Loaded"<<endl;
  t1->SetBranchAddress("Track_P",&P);
  t1->SetBranchAddress("Track_Pt",&Pt);
  t1->SetBranchAddress("Track_Phi",&Phi);
  t1->SetBranchAddress("Track_eta",&Eta);

  t1->SetBranchAddress("Chi2_CubicXZ",&Chi2CubXZ);
  t1->SetBranchAddress("Chi2_ParabolaXZ",&Chi2ParXZ);
  t1->SetBranchAddress("Chi2_ParabolaY",&Chi2ParY);
  t1->SetBranchAddress("Chi2_LineY",&Chi2LineY);

  t1->SetBranchAddress("Chi2_XZDRATIO",&Chi2dRatioXZ);
  t1->SetBranchAddress("Chi2_XZDRATIOFixed",&Chi2dRatioXZFix);

  gEnv->GetValue("Canvas.SavePrecision", -1);
  gEnv->SetValue("Canvas.SavePrecision", 16);
  //Phi
  Int_t NBins = 100;
  Double_t minPhi = -3.2;
  Double_t maxPhi = 3.2;
  TH1D * Phi_ParXZ = new TH1D("Phi_Chi2ParXZ","Phi_ParXZ",NBins,minPhi,maxPhi);
  TH1D * Phi_ParXZ_norm = new TH1D("Phit_Chi2ParXZ_norm","Phit_ParXZ_n",NBins,minPhi,maxPhi);

  TH1D * Phi_CubXZ = new TH1D("Phi_Chi2CubXZ","Phi_CubXZ",NBins,minPhi,maxPhi);
  TH1D * Phi_CubXZ_norm = new TH1D("Phit_Chi2CubXZ_norm","Phit_CubXZ_n",NBins,minPhi,maxPhi);

  TH1D * Phi_dRatioXZ = new TH1D("Phi_Chi2dRatioXZ","Phi_dRatioXZ",NBins,minPhi,maxPhi);
  TH1D * Phi_dRatioXZ_norm = new TH1D("Phit_Chi2dRatioXZ_norm","Phit_dRatioXZ_n",NBins,minPhi,maxPhi);
  TH1D * Phi_dRatioXZFix = new TH1D("Phi_Chi2dRatioXZFix","Phi_dRatioXZFix",NBins,minPhi,maxPhi);
  TH1D * Phi_dRatioXZFix_norm = new TH1D("Phit_Chi2dRatioXZFix_norm","Phit_dRatioXZFix_n",NBins,minPhi,maxPhi);
//
  TH1D * Phi_LineY = new TH1D("Phi_Chi2LineY","Phi_LineY",NBins,minPhi,maxPhi);
  TH1D * Phi_LineY_norm = new TH1D("Phit_Chi2LineY_norm","Phit_LineY_n",NBins,minPhi,maxPhi);

  TH1D * Phi_ParY = new TH1D("Phi_Chi2ParY","Phi_ParY",NBins,minPhi,maxPhi);
  TH1D * Phi_ParY_norm = new TH1D("Phi_Chi2ParY_norm","Phit_ParY_n",NBins,minPhi,maxPhi);

  //P
  Double_t minP = 1000;
  Double_t maxP = 10000;
  TH1D * P_ParXZ = new TH1D("P_Chi2ParXZ","P_ParXZ",NBins,minP,maxP);
  TH1D * P_ParXZ_norm = new TH1D("Pt_Chi2ParXZ_norm","Pt_ParXZ_n",NBins,minP,maxP);

  TH1D * P_CubXZ = new TH1D("P_Chi2CubXZ","P_CubXZ",NBins,minP,maxP);
  TH1D * P_CubXZ_norm = new TH1D("Pt_Chi2CubXZ_norm","Pt_CubXZ_n",NBins,minP,maxP);

  TH1D * P_dRatioXZ = new TH1D("P_Chi2dRatioXZ","P_dRatioXZ",NBins,minP,maxP);
  TH1D * P_dRatioXZ_norm = new TH1D("Pt_Chi2dRatioXZ_norm","Pt_dRatioXZ_n",NBins,minP,maxP);
  TH1D * P_dRatioXZFix = new TH1D("P_Chi2dRatioXZFix","P_dRatioXZFix",NBins,minP,maxP);
  TH1D * P_dRatioXZFix_norm = new TH1D("Pt_Chi2dRatioXZFix_norm","Pt_dRatioXZFix_n",NBins,minP,maxP);


  TH1D * P_LineY = new TH1D("P_Chi2LineY","P_LineY",NBins,minP,maxP);
  TH1D * P_LineY_norm = new TH1D("Pt_Chi2LineY_norm","Pt_LineY_n",NBins,minP,maxP);

  TH1D * P_ParY = new TH1D("P_Chi2ParY","P_ParY",NBins,minP,maxP);
  TH1D * P_ParY_norm = new TH1D("Pt_Chi2ParY_norm","Pt_ParY_n",NBins,minP,maxP);

  //Pt
  Double_t minPt = 0;
  Double_t maxPt = 5000;
  TH1D * Pt_ParXZ = new TH1D("Pt_Chi2ParXZ","Pt_ParXZ",NBins,minPt,maxPt);
  TH1D * Pt_ParXZ_norm = new TH1D("Pt_Chi2ParXZ_norm","Pt_ParXZ_n",NBins,minPt,maxPt);

  TH1D * Pt_CubXZ = new TH1D("Pt_Chi2CubXZ","Pt_CubXZ",NBins,minPt,maxPt);
  TH1D * Pt_CubXZ_norm = new TH1D("Pt_Chi2CubXZ_norm","Pt_CubXZ_n",NBins,minPt,maxPt);

  TH1D * Pt_dRatioXZ = new TH1D("Pt_Chi2dRatioXZ","Pt_dRatioXZ",NBins,minPt,maxPt);
  TH1D * Pt_dRatioXZ_norm = new TH1D("Pt_Chi2dRatioXZ_norm","Pt_dRatioXZ_n",NBins,minPt,maxPt);
  TH1D * Pt_dRatioXZFix = new TH1D("Pt_Chi2dRatioXZFix","Pt_dRatioXZFix",NBins,minPt,maxPt);
  TH1D * Pt_dRatioXZFix_norm = new TH1D("Pt_Chi2dRatioXZFix_norm","Pt_dRatioXZFix_n",NBins,minPt,maxPt);

  TH1D * Pt_LineY = new TH1D("Pt_Chi2LineY","Pt_LineY",NBins,minPt,maxPt);
  TH1D * Pt_LineY_norm = new TH1D("Pt_Chi2LineY_norm","Pt_LineY_n",NBins,minPt,maxPt);

  TH1D * Pt_ParY = new TH1D("Pt_Chi2ParY","Pt_ParY",NBins,minPt,maxPt);
  TH1D * Pt_ParY_norm = new TH1D("Pt_Chi2ParY_norm","Pt_ParY_n",NBins,minPt,maxPt);



  Double_t minEta = 0.5;
  Double_t maxEta = 7;
  TH1D * Eta_ParXZ = new TH1D("Eta_Chi2ParXZ","Eta_ParXZ",NBins,minEta,maxEta);
  TH1D * Eta_ParXZ_norm = new TH1D("Eta_Chi2ParXZ_norm","Eta_ParXZ_n",NBins,minEta,maxEta);

  TH1D * Eta_CubXZ = new TH1D("Eta_Chi2CubXZ","Eta_CubXZ",NBins,minEta,maxEta);
  TH1D * Eta_CubXZ_norm = new TH1D("Eta_Chi2CubXZ_norm","Eta_CubXZ_n",NBins,minEta,maxEta);

  TH1D * Eta_dRatioXZ = new TH1D("Eta_Chi2dRatioXZ","Eta_dRatioXZ",NBins,minEta,maxEta);
  TH1D * Eta_dRatioXZ_norm = new TH1D("Eta_Chi2dRatioXZ_norm","Eta_dRatioXZ_n",NBins,minEta,maxEta);
  TH1D * Eta_dRatioXZFix = new TH1D("Eta_Chi2dRatioFixXZ","Eta_dRatioXZFix",NBins,minEta,maxEta);
  TH1D * Eta_dRatioXZFix_norm = new TH1D("Eta_Chi2dRatioXZFix_norm","Eta_dRatioXZFix_n",NBins,minEta,maxEta);



  TH1D * Eta_LineY = new TH1D("Eta_Chi2LineY","Eta_LineY",NBins,minEta,maxEta);
  TH1D * Eta_LineY_norm = new TH1D("Eta_Chi2LineY_norm","Eta_LineY_n",NBins,minEta,maxEta);

  TH1D * Eta_ParY = new TH1D("Eta_Chi2ParY","Eta_ParY",NBins,minEta,maxEta);
  TH1D * Eta_ParY_norm = new TH1D("Eta_Chi2ParY_norm","Eta_ParY_n",NBins,minEta,maxEta);

  for (Int_t i= 0; i<nEntries; i++) {
    t1->GetEntry(i);
    if(Chi2ParXZ>800) continue;
    if(Chi2LineY>800) continue;
    Phi_ParXZ->Fill(Phi,Chi2ParXZ);
    Phi_ParXZ_norm->Fill(Phi);

    Phi_CubXZ->Fill(Phi,Chi2CubXZ);
    Phi_CubXZ_norm->Fill(Phi);

    Phi_dRatioXZ->Fill(Phi,Chi2dRatioXZ);
    Phi_dRatioXZ_norm->Fill(Phi);
    Phi_dRatioXZFix->Fill(Phi,Chi2dRatioXZFix);
    Phi_dRatioXZFix_norm->Fill(Phi);

    Phi_LineY->Fill(Phi,Chi2LineY);
    Phi_LineY_norm->Fill(Phi);

    Phi_ParY->Fill(Phi,Chi2ParY);
    Phi_ParY_norm->Fill(Phi);

    P_ParXZ->Fill(P,Chi2ParXZ);
    P_ParXZ_norm->Fill(P);

    P_CubXZ->Fill(P,Chi2CubXZ);
    P_CubXZ_norm->Fill(P);

    P_dRatioXZ->Fill(P,Chi2dRatioXZ);
    P_dRatioXZ_norm->Fill(P);
    P_dRatioXZFix->Fill(P,Chi2dRatioXZFix);
    P_dRatioXZFix_norm->Fill(P);


    P_LineY->Fill(P,Chi2LineY);
    P_LineY_norm->Fill(P);

    P_ParY->Fill(P,Chi2ParY);
    P_ParY_norm->Fill(P);


    Pt_ParXZ->Fill(Pt,Chi2ParXZ);
    Pt_ParXZ_norm->Fill(Pt);

    Pt_CubXZ->Fill(Pt,Chi2CubXZ);
    Pt_CubXZ_norm->Fill(Pt);

    Pt_dRatioXZ->Fill(Pt,Chi2dRatioXZ);
    Pt_dRatioXZ_norm->Fill(Pt);
    Pt_dRatioXZFix->Fill(Pt,Chi2dRatioXZFix);
    Pt_dRatioXZFix_norm->Fill(Pt);

    Pt_LineY->Fill(Pt,Chi2LineY);
    Pt_LineY_norm->Fill(Pt);

    Pt_ParY->Fill(Pt,Chi2ParY);
    Pt_ParY_norm->Fill(Pt);

    Eta_ParXZ->Fill(Eta,Chi2ParXZ);
    Eta_ParXZ_norm->Fill(Eta);

    Eta_CubXZ->Fill(Eta,Chi2CubXZ);
    Eta_CubXZ_norm->Fill(Eta);
    Eta_dRatioXZ->Fill(Eta,Chi2dRatioXZ);
    Eta_dRatioXZ_norm->Fill(Eta);
    Eta_dRatioXZFix->Fill(Eta,Chi2dRatioXZFix);
    Eta_dRatioXZFix_norm->Fill(Eta);


    Eta_LineY->Fill(Eta,Chi2LineY);
    Eta_LineY_norm->Fill(Eta);

    Eta_ParY->Fill(Eta,Chi2ParY);
    Eta_ParY_norm->Fill(Eta);


  }
  Phi_ParY->Divide(Phi_ParY,Phi_ParY_norm);
  Eta_ParY->Divide(Eta_ParY,Eta_ParY_norm);
  P_ParY->Divide(P_ParY,P_ParY_norm);
  Pt_ParY->Divide(Pt_ParY,Pt_ParY_norm);

  Phi_LineY->Divide(Phi_LineY,Phi_LineY_norm);
  Eta_LineY->Divide(Eta_LineY,Eta_LineY_norm);
  P_LineY->Divide(P_LineY,P_LineY_norm);
  Pt_LineY->Divide(Pt_LineY,Pt_LineY_norm);

  Phi_CubXZ->Divide(Phi_CubXZ,Phi_CubXZ_norm);
  Eta_CubXZ->Divide(Eta_CubXZ,Eta_CubXZ_norm);
  P_CubXZ->Divide(P_CubXZ,P_CubXZ_norm);
  Pt_CubXZ->Divide(Pt_CubXZ,Pt_CubXZ_norm);

  Phi_dRatioXZ->Divide(Phi_dRatioXZ,Phi_dRatioXZ_norm);
  Eta_dRatioXZ->Divide(Eta_dRatioXZ,Eta_dRatioXZ_norm);
  P_dRatioXZ->Divide(P_dRatioXZ,P_dRatioXZ_norm);
  Pt_dRatioXZ->Divide(Pt_dRatioXZ,Pt_dRatioXZ_norm);

  Phi_dRatioXZFix->Divide(Phi_dRatioXZFix,Phi_dRatioXZFix_norm);
  Eta_dRatioXZFix->Divide(Eta_dRatioXZFix,Eta_dRatioXZFix_norm);
  P_dRatioXZFix->Divide(P_dRatioXZFix,P_dRatioXZFix_norm);
  Pt_dRatioXZFix->Divide(Pt_dRatioXZFix,Pt_dRatioXZFix_norm);

  Phi_ParXZ->Divide(Phi_ParXZ,Phi_ParXZ_norm);
  Eta_ParXZ->Divide(Eta_ParXZ,Eta_ParXZ_norm);
  P_ParXZ->Divide(P_ParXZ,P_ParXZ_norm);
  Pt_ParXZ->Divide(Pt_ParXZ,Pt_ParXZ_norm);


  TCanvas *c1 = new TCanvas("YFit","YFit");
  c1->Divide(2,2);
  c1->cd(1);
  Phi_LineY->SetLineColor(kBlue);
  Phi_LineY->Draw();
  Phi_ParY->SetLineColor(kRed);
  Phi_ParY->Draw("same");
  //c1->BuildLegend();
  c1->cd(2);
  Eta_LineY->SetLineColor(kBlue);
  Eta_LineY->Draw();
  Eta_ParY->SetLineColor(kRed);
  Eta_ParY->Draw("same");
  //c1->BuildLegend();
  c1->cd(3);
  Pt_LineY->SetLineColor(kBlue);
  Pt_LineY->Draw();
  Pt_ParY->SetLineColor(kRed);
  Pt_ParY->Draw("same");
  //c1->BuildLegend();
  c1->cd(4);
  P_LineY->SetLineColor(kBlue);
  P_LineY->Draw();
  P_ParY->SetLineColor(kRed);
  P_ParY->Draw("same");

  TCanvas *c2 = new TCanvas("XZFit", "XZFit");
  c2->Divide(2,2);
  c2->cd(1);
  Phi_ParXZ->SetLineColor(kBlack);
  Phi_ParXZ->Draw();
  Phi_CubXZ->SetLineColor(kBlue);
  Phi_CubXZ->Draw("same");
  Phi_dRatioXZ->SetLineColor(kRed);
  Phi_dRatioXZ->Draw("same");
  Phi_dRatioXZFix->SetLineColor(kGreen);
  Phi_dRatioXZFix->Draw("same");

  //->BuildLegend();
  c2->cd(2);
  Eta_ParXZ->SetLineColor(kBlack);
  Eta_ParXZ->Draw();
  Eta_CubXZ->SetLineColor(kBlue);
  Eta_CubXZ->Draw("same");
  Eta_dRatioXZ->SetLineColor(kRed);
  Eta_dRatioXZ->Draw("same");
  Eta_dRatioXZFix->SetLineColor(kGreen);
  Eta_dRatioXZFix->Draw("same");

  //c2->BuildLegend();
  c2->cd(3);
  Pt_ParXZ->SetLineColor(kBlack);
  Pt_ParXZ->Draw();
  Pt_CubXZ->SetLineColor(kBlue);
  Pt_CubXZ->Draw("same");
  Pt_dRatioXZ->SetLineColor(kRed);
  Pt_dRatioXZ->Draw("same");
  Pt_dRatioXZFix->SetLineColor(kGreen);
  Pt_dRatioXZFix->Draw("same");

  c2->cd(4);
  P_ParXZ->SetLineColor(kBlack);
  P_ParXZ->Draw();
  P_CubXZ->SetLineColor(kBlue);
  P_CubXZ->Draw("same");
  P_dRatioXZ->SetLineColor(kRed);
  P_dRatioXZ->Draw("same");
  P_dRatioXZFix->SetLineColor(kGreen);
  P_dRatioXZFix->Draw("same");


}
