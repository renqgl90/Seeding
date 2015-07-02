#define FakedSeeding_cxx

#include "FakedSeeding.h"
//#include <TH2.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
//#include <stdio>
#include <iostream>
//#include "params.h"
#include "TF1.h"
#include <vector>
#include "LinParFit.h"
//===============
//L0
//==============
//L0 Params
const double m_ConstC = 2.34e8;
//const double m_alphaCorrection = 2900.;
const double m_alphaCorrection = 120.64;
const double L0_Up = +265.*5.5;
const double L0_Down = -265.*5.5;
const double m_x0Corr = 0.0023;
//X0 correction (See L0 Plots)
//x0<0 : res_Corr < m_m1(x+m_x0Offset) +m_yOff1
//       res_Corr > m_m2(x+m_x0Offset) +m_yOff2

//x0>0 : res_Corr > m_m2(x-m_x0Offset) -m_yOff2
//       res_Corr < m_m1(x+m_x0Offset) -m_yOff1
//L1 Params
const double m_m1 = 0;//-0.00272;
const double m_m2 = 0;//0.001;
const double m_x0Offset = 50;
const double m_yOff1 = +3.5;
const double m_yOff2 = -3.5;
//L2 Params
//Delta From Parab 1 Extrapolation
const double maxDelta_par1 = 3.0;
const double maxDelta_par2 = 3.0;

//TH2D * L0_0_P = new TH2D("L0","X_{last}^{true}-X_{infinite}^{projected} Vs X_{infite}^{projected};X_{infinite}^{projected}[mm];X_{last}^{true}-X_{infinite}^{projected}",100,-3000,3000,100,-500,500);



TH2D * L0_0 = new TH2D("L0","X_{last}^{true}-X_{infinite}^{projected} Vs X_{infite}^{projected} (P>5GeV & OVTX < 9m);X_{infinite}^{projected}[mm];X_{last}^{true}-X_{infinite}^{projected}",100,-3000,3000,100,-4000,4000);
TH2D * L0_1 = new TH2D("L0_1","X_{last}^{true}-X_{infinite}^{projected} Vs t^{infinite}_{x};t_{x}^{infinite};X_{last}^{true}-X_{infinite}^{projected}[mm] ",100,-0.5,0.5,100,-500,500);
TH2D * L0_Corrected = new TH2D("L0_Corrected","X_{last}^{true}-X_{infinite}^{projected} Corrected Vs t^{infinite}_{x};t_{x}^{infinite};X_{last}^{true}-X_{infinite}^{projected} Corrected [mm]",100,-0.5,0.5,100,-4000,4000);
TH1D *L0_1d = new TH1D("Residual_noCorrection","X_{last}-X_{infinite}^{projected};X_{last}^{true}-X_{infinite}^{projected} [mm]",300,-600,600);
TH1D *L0_Corrected_1d = new TH1D("Residual_Correction","X_{last}^{true}-X_{infinite}^{projected} Corrected;X_{last}^{true}-X_{infinite}^{projected} [mm] ; Counts ",300,-4000,4000);
TF1 *L0_Upper_Bound = new TF1("L0_Upper_Bounds","x*[0]+[1]",-0.5,0.5);
TF1 *L0_Lower_Bound = new TF1("L0_Lower_Bounds","x*[0]+[1]",-0.5,0.5);
//L0_Lower_Bound->FixParameters(m_alphaCorrection,L0_Down);


//===============
//L1
//==============
TH1D * L1_0_1D = new TH1D("L1_0_1D","#Delta(X_{3rd}-x^{projected});#Delta(X_{3rd}-x^{projected}) [mm];Counts/0.12 [mm]",100,-6,6);
TH1D * L1_0_1D_Corr = new TH1D("L1_0_1D_Corr","#Delta(X_{3rd}-x^{projected}) Corrected;#Delta(X_{3rd}-x^{projected}) [mm];Counts/0.12 [mm]",100,-6,6);
TH2D * L1_0 = new TH2D("L1_0","#Delta(X_{3rd}-x^{projected}) Vs t_{x}^{picked};x_{0};X_{T2-2}-x^{projected}[mm]",100,-10000,10000,100,-6.,6.);
TH2D * L1_0_corr = new TH2D("L1_0_Corr","#Delta(X_{3rd}-x^{projected}) Vs t_{x}^{picked} Corrected;x_{0};X_{T2-2}-x^{projected} Corrected [mm]",100,-1000,10000,100,-6.,6.);
TH2D * L1_0_x0pos = new TH2D("L1_0_x0pos","#Delta(X_{3rd}-x^{projected}) Vs t_{x}^{picked} (x_{0} >0);t_{x}^{picked};X_{3rd}-x^{projected}[mm]",100,-0.5,0.5,100,-6.,6.);
TH2D * L1_0_x0neg = new TH2D("L1_0_x0neg","#Delta(X_{3rd}-x^{projected}) Vs t_{x}^{picked} (x_{0} <0);t_{x}^{picked};X_{3rd}-x^{projected}[mm]",100,-0.5,0.5,100,-6.,6.);
TH2D * L1_1_x0pos = new TH2D("L1_1_x0pos","#Delta(X_{3rd}-x^{projected}) Vs #Delta(t_{x}^{picked} - t_{x}^{inf}) (x_{0} >0);#Delta(t_{x}^{picked} - t_{x}^{inf});X_{3rd}-x^{projected} [mm]}",100,-0.2,0.2,100,-6.,6.);
TH2D * L1_1_x0neg = new TH2D("L1_1_x0neg","#Delta(X_{3rd}-x^{projected}) Vs #Delta(t_{x}^{picked} - t_{x}^{inf}) (x_{0} <0);#Delta(t_{x}^{picked} - t_{x}^{inf});X_{3rd}-x^{projected} [mm]",100,-0.2,0.2,100,-6.,6.);
TH2D * L1_2_x0pos = new TH2D("L1_2_x0pos","#Delta(X_{3rd}-x^{projected}) Vs x0 (x0 >0);x_{0} [mm] ;X_{3rd}-x^{projected}[mm]",100,-200,1500,100,-6.,6.);
TH2D * L1_2_x0neg = new TH2D("L1_2_x0neg","#Delta(X_{3rd}-x^{projected}) Vs x0 (x0 <0);x_{0} [mm] ;X_{3rd}-x^{projected}[mm]",100,-1500,200,100,-6.,6.);

TH2D * L1_0_x0posCorr = new TH2D("L1_0_x0posCorr","#Delta(X_{3rd}-x^{projected}) Vs t_{x}^{picked} (x_{0} >0);t_{x}^{picked};X_{3rd}-x^{projected}[mm]",100,-0.5,0.5,100,-6.,6.);
TH2D * L1_0_x0negCorr = new TH2D("L1_0_x0negCorr","#Delta(X_{3rd}-x^{projected}) Vs t_{x}^{picked} (x_{0} <0);t_{x}^{picked};X_{3rd}-x^{projected}[mm]",100,-0.5,0.5,100,-6.,6.);

TH2D * L1_1_x0posCorr = new TH2D("L1_1_x0posCorr","#Delta(X_{3rd}-x^{projected}) Vs #Delta(t_{x}^{picked} - t_{x}^{inf}) (x_{0} >0);#Delta(t_{x}^{picked} - t_{x}^{inf});X_{3rd}-x^{projected} [mm]}",100,-0.2,0.2,100,-6.,6.);
TH2D * L1_1_x0negCorr = new TH2D("L1_1_x0negCorr","#Delta(X_{3rd}-x^{projected}) Vs #Delta(t_{x}^{picked} - t_{x}^{inf}) (x_{0} <0);#Delta(t_{x}^{picked} - t_{x}^{inf});X_{3rd}-x^{projected} [mm]",100,-0.2,0.2,100,-6.,6.);
TH2D * L1_2_x0posCorr = new TH2D("L1_2_x0posCorr","#Delta(X_{3rd}-x^{projected}) Vs x0 (x0 >0);x_{0} [mm] ;X_{3rd}-x^{projected}[mm]",100,-200,1500,1000,-6.,6.);
TH2D * L1_2_x0negCorr = new TH2D("L1_2_x0negCorr","#Delta(X_{3rd}-x^{projected}) Vs x0 (x0 <0);x_{0} [mm] ;X_{3rd}-x^{projected}[mm]",100,-1500,200,1000,-4.,6.);

//x0<0
TF1 * L1_x0Upper_x0neg = new TF1("L2_x0Upper_x0pos","[0]*(x+[1])+[2]",-1500,1500);
TF1 * L1_x0Lower_x0neg = new TF1("L2_x0Lower_x0pos","[0]*(x+[1])+[2]",-1500,1500);
//x0>0
TF1 * L1_x0Upper_x0pos = new TF1("L2_x0Upper_x0pos","[0]*(x-[1])-[2]",-1500,1500);
TF1 * L1_x0Lower_x0pos = new TF1("L2_x0Lower_x0pos","[0]*(x-[1])-[2]",-1500,1500);




//==========================================
//L2
//=========================================
TH1D * L2_0_T1_par1 = new TH1D("L2_0_T1_par1","#Delta (X_{true} - X_{extrap}) in T1;#Delta(X_{true}-X_{extrapol});Counts",100,-2,2);
TH1D * L2_0_T2_par1 = new TH1D("L2_0_T2_par1","#Delta (X_{true} - X_{extrap}) in T2;#Delta(X_{true}-X_{extrapol});Counts",100,-2,2);
TH1D * L2_0_T3_par1 = new TH1D("L2_0_T3_par1","#Delta (X_{true} - X_{extrap}) in T3;#Delta(X_{true}-X_{extrapol});Counts",100,-2,2);

TH1D * L2_0_T1_par2 = new TH1D("L2_0_T1_par2","#Delta (X_{true} - X_{extrap}) in T1;#Delta(X_{true}-X_{extrapol});Counts",100,-2,2);
TH1D * L2_0_T2_par2 = new TH1D("L2_0_T2_par2","#Delta (X_{true} - X_{extrap}) in T2;#Delta(X_{true}-X_{extrapol});Counts",100,-2,2);
TH1D * L2_0_T3_par2 = new TH1D("L2_0_T3_par2","#Delta (X_{true} - X_{extrap}) in T3;#Delta(X_{true}-X_{extrapol});Counts",100,-2,2);
TH1D * Chi2_par1 = new TH1D("Chi2_par1","#Chi^{2} (parabola) on 6 Hit;#Chi^{2};Counts",100,0,30);
TH1D * Chi2_par2 = new TH1D("Chi2_par2","#Chi^{2} (parabola corr) on 6 Hit;#Chi^{2};Counts",100,0,30);
//L2_0_T1_par1
//L2_0_T2_par1
//L2_0_T3_par1
//L2_0_T1_par2
//L2_0_T2_par2
//L2_0_T3_par2
Double_t xFirst,xLast,yFirst,yLast,yBack;
Double_t XSlopeT3,XSlopeT1;
Double_t OrVTX_Z,OrVTX_X,OrVTX_y;
Int_t RegionT1,RegionT3;
Int_t XZType;
Int_t XYType;
Double_t P,ax_0,bx_0,cx_0,dx_0,ay_0,by_0,cy_0;
Int_t ID;

char nameFile[100]="data/ProcessedTrack.root";
TFile *fFile = new TFile(nameFile,"RECREATE");
TTree *t1 = new TTree("t1","myData");
void FakedSeeding::Begin(TTree * /*tree*/)
{
  m_Total = 0;
  m_L0Selected =0;
  m_L1Selected =0;
  m_L2_Par1_Selected = 0;
  m_L2_Par2_Selected = 0;
  m_dRatio = 0.000244;
  m_zReference = 8520.;
  doFit = true;
  TString option = GetOption();
  std::cout << "Begin selector ... " <<std::endl;
  t1->Branch("P",&P,"P/D");	//energia nel centro di massa
  t1->Branch("ax",&ax_0,"ax/D");
  t1->Branch("bx",&bx_0,"bx/D");
  t1->Branch("cx",&cx_0,"cx/D");
  t1->Branch("ay",&ay_0,"ay/D");
  t1->Branch("by",&by_0,"by/D");
  t1->Branch("cy",&cy_0,"cy/D");
  t1->Branch("dx",&dx_0,"dx/D");
  t1->Branch("ID",&ID,"ID/I");
  t1->Branch("xFirst",&xFirst,"xFirst/D");
  t1->Branch("xLast",&xLast,"xLast/D");
  t1->Branch("yFirst",&yFirst,"yFirst/D");
  t1->Branch("yLast",&yLast,"yLast/D");
  t1->Branch("yBack",&yBack,"yBack/D");
  t1->Branch("XZType",&XZType,"XZType/I");
  t1->Branch("XYType",&XYType,"XYType/I");
  t1->Branch("XSlopeT1",&XSlopeT1,"XSlopeT1/D");
  t1->Branch("XSlopeT3",&XSlopeT3,"XSlopeT3/D");
  t1->Branch("OVTX_Z",&OrVTX_Z,"OVTX_Z/D");
  t1->Branch("OVTX_Y",&OrVTX_Y,"OVTX_Y/D");
  t1->Branch("OVTX_X",&OrVTX_X,"OVTX_X/D");
  //t1->Branch("xBackward",&xBack,"xBack/D");
}
void FakedSeeding::SlaveBegin(TTree * /*tree*/){
  float coao;
  TString option = GetOption();
}
Bool_t FakedSeeding::Process(Long64_t entry)
{
  fChain->GetTree()->GetEntry(entry);
  //Map
  float zref = (float)m_zReference;
  PrSeedTrack *Track = new PrSeedTrack(zref);
  Track->setZref(zref);
  Track->setdRatio((float)m_dRatio);
  PrSeedTrack *Track_All = new PrSeedTrack(zref);
  Track_All->setZref(zref);
  Track_All->setdRatio((float)m_dRatio);
  double Z1  = 7855.6;
  double Z2  = 8040.4;
  double Z3  = 8537.6;
  const float dxDy = (5./360.00000)*2.*TMath::Pi();
  double Z4  = 8722.4;
  double Z5  = 9222.6;
  double Z6  = 9407.4;
  double ZU1 = 7917.2;
  double ZV1 = 7978.8;
  double ZU2 = 8599.2;
  double ZV2 = 8660.8;
  double ZU3 = 9284.8;
  double ZV3 = 9345.8;
  double X1 = X1st;
  double X2 = X2nd;
  double X3 = X3rd;
  double X4 = X4th;
  double X5 = X5th;
  double X6 = X6th;
  double XU1 = Ux1st;
  double XU2 = Ux2nd;
  double XU3 = Ux3rd;
  double XV1 = Vx1st;
  double XV2 = Vx2nd;
  double XV3 = Vx3rd;
  Hit hit1  =Hit(X1,Y1st,Z1);
  Hit hit2  =Hit(XU1,Uy1st,ZU1);
  hit2.SetDxDy(dxDy);
  Hit hit3  =Hit(XV1,Vy1st,ZV1);
  hit3.SetDxDy(-dxDy);
  Hit hit4  =Hit(X2,Y2nd,Z2);
  Hit hit5  =Hit(X3,Y3rd,Z3);
  Hit hit6  =Hit(XU2,Uy2nd,ZU2);
  hit6.SetDxDy(dxDy);
  Hit hit7  =Hit(XV2,Vy2nd,ZV2);
  hit7.SetDxDy(-dxDy);
  Hit hit8  =Hit(X4,Y4th,Z4);
  Hit hit9  =Hit(X5,Y5th,Z5);
  Hit hit10 =Hit(XU3,Uy3rd,ZU3);
  hit10.SetDxDy(-dxDy);
  Hit hit11 =Hit(XV3,Vy3rd,ZV3);
  hit11.SetDxDy(dxDy);
  Hit hit12 =Hit(X6,Y6th,Z6);


  std::vector<Hit> AllHits;
  std::random_device rd_1;
  std::mt19937 gen_1(rd_1());
  std::discrete_distribution<> d_1({4, 24, 33, 39});
  hit1.SetSize((double)d_1(gen_1));
  hit2.SetSize((double)d_1(gen_1));
  hit3.SetSize((double)d_1(gen_1));
  hit4.SetSize((double)d_1(gen_1));
  hit5.SetSize((double)d_1(gen_1));
  hit6.SetSize((double)d_1(gen_1));
  hit7.SetSize((double)d_1(gen_1));
  hit8.SetSize((double)d_1(gen_1));
  hit9.SetSize((double)d_1(gen_1));
  hit10.SetSize((double)d_1(gen_1));
  hit11.SetSize((double)d_1(gen_1));
  hit12.SetSize((double)d_1(gen_1));
  AllHits.push_back(hit1);
  AllHits.push_back(hit2);
  AllHits.push_back(hit3);
  AllHits.push_back(hit4);
  AllHits.push_back(hit5);
  AllHits.push_back(hit6);
  AllHits.push_back(hit7);
  AllHits.push_back(hit8);
  AllHits.push_back(hit9);
  AllHits.push_back(hit10);
  AllHits.push_back(hit11);
  AllHits.push_back(hit12);
  Track_All->addHits(AllHits);
  //Track_All->printHits();
  //in a for loop Size = (double)d_1(gen_1);

  //FITTA(Track_All);
  std::vector<Hit> Hitss = Track_All->hits();
  double solutions[4];
  double solutions_y[3];
  std::vector<double> dsolution;
  std::vector<double> dsolution_y;
  std::fill(solutions,solutions+4,0.);
  std::fill(solutions_y,solutions_y+3,0.);
  double res = 0.;
  double resy = 0.;
  double dz = 0;
  // if(P<5000) return kTRUE;
  // if(OVTX_Z>4000) return kTRUE;
  // if( std::abs(PID)==11 || (std::abs(PID)!=211 && std::abs(PID)!=2212 && std::abs(PID)!=321 && std::abs(PID)!=13))
  //   return kTRUE;
  // if ((Y1st>=0 && Y6th<=0)||(Y1st>=0 && Y6th<=0))
  //   return kTRUE;
  // if(P<500)
  //   return kTRUE;
  int regionT1= -1;
  int regionT3 = -1;
  if(std::fabs(Y1st)>500){
    if(std::fabs(X1st)>500)
      regionT1 = 3;
    if(std::fabs(X1st)<=500)
      regionT1 = 2;
  }
  if(std::fabs(Y1st)<500){
    if(std::fabs(X1st>500))
      regionT1 = 1;
    if(std::fabs(X1st<500))
      regionT1 = 0;
  }
  if(std::fabs(Y6th)>500){
    if(std::fabs(X6thx)>500)
      regionT3 = 3;
    if(std::fabs(X6th)<=500)
      regionT3 = 2;
  }
  if(std::fabs(Y6th)<500){
    if(std::fabs(X6th>500))
      regionT3 = 1;
    if(std::fabs(X1st<500))
      regionT3 = 0;
  }
  XYType=-1;
  XZType=-1;
  //Region  motion With Y info, without Y info
  //2 / 3  0-->0 XYType = 0              0/2 --> 0/2 XZType = 0   2*regionT1%2 + regionT3%2
  //0 / 1  0-->1 XYType = 1              0/2 --> 1/3 XZType = 1
  //       0-->2 XYType = 2              1/3 --> 0/2 XZType = 2
  //       0-->3 XYType = 3              1/3 --> 1/3 XYType = 3
  //       1-->0 XYType = 4
  //       1-->1 XYType = 5
  //       1-->2 XYType = 6
  //       1-->3 XYType = 7
  //       2-->0 XYtype = 8
  //       2-->1 XYType = 9
  //       2-->2 XYType = 10
  //       2-->3 XYType = 11
  //       3-->0 XYType = 12
  //       3-->1 XYType = 13
  //       3-->2 XYType = 14
  //       3-->3 XYType = 15
  XYType = regionT1*4+regionT3;
  XZType = 2*regionT2%2+regionT3%2;

  if(doFit){
    for(int j=0;j<3;j++){//loop
      //std::fill(dsolutions,dsolutions+4,0.);
      LinParFit<double> fit(4);
      LinParFit<double> fit_y(3);
      for (int i = 0;i<Hitss.size();i++){
        dz =(double)(Hitss[i].GetZ()-m_zReference);
        res=(double)(Hitss[i].GetX() - (solutions[0]+solutions[1]*dz+solutions[2]*dz*dz+solutions[3]*dz*dz*dz));//*dz*dz*dz+solutions[4]*dz*dz*dz*dz));
        resy = (double)(Hitss[i].GetY()-(solutions_y[0]+solutions_y[1]*dz + solutions_y[2]*dz*dz));
        //std::cout<<"Loop\t"<<j<<"\t Residual \t"<<res<<std::endl;
        fit_y.accumulate(resy,0.100,dz);
        fit.accumulate(res,0.100,dz);
      }
      //std::vector<double> solution;
      dsolution_y= fit_y.solution();
      dsolution=fit.solution();
      for (int k=0;k<4;k++){
        if(k<3){
          solutions_y[k]+=dsolution_y[k];
        }
        solutions[k]+=dsolution[k];
      }
    }
    //Fit for y(z) = [0]+[1]*(z-zRef)
    //So y(zRef) = [0]

    double ay = solutions_y[0];// - m_zReference*solutions_y[1];
    double by = solutions_y[1];
    double cy = solutions_y[2];
    //std::cout<<"Par \t Value "<<std::endl;
    double a = solutions[0];
    double b = solutions[1];
    double c = solutions[2];
    double e = solutions[3]/solutions[2];
    double Constant_C;
    ay_0=ay;
    by_0=by;
    cy_0=cy;
    ax_0=a;
    bx_0=b;
    cx_0=c;
    dx_0=e;
    xFirst=X1st;
    xLast=X6th;
    yFirst=Y1st;
    yLast=Y6th;
    t1->Fill();
    //  double f = solutions[4]/(solutions[3]*solutions[2]);
    // std::cout<<" a   \t"<<a
    //          <<"\n b \t"<<b
    //          <<"\n c \t"<<c
    //          <<"\n e \t"<<e<<std::endl;
    for (int i = 0;i<Hitss.size();i++){
      dz = (double)(Hitss[i].GetZ()-m_zReference);
      res = (double)(Hitss[i].GetX() - (a+b*dz+c*dz*dz*(1.+e*dz)));//+solutions[3]));//*dz*dz*dz+solutions[4]*dz*dz*dz*dz));
      //   std::d::cout<<"Hit Index\t Residual"<<"\n "<<i<<"\t"<<res<<std::endl;
    }
    //          <<"\n f \t"<<f
    //          <<std::endl;
    double dRatioy = solutions_y[2]/solutions_y[1];
    //std::cout<<"Fitted dRatio y " <<solutions_y[2]/solutions_y[1]<<std::endl;
    Fitted_dRatioy->Fill(dRatioy*10e6);
    Fitted_dRatioyVsC->Fill(dRatioy*10e6,c);
    Fitted_dRatioyVsP->Fill(dRatioy*10e6,P);
    Fitted_dRatioyVsXFirst->Fill(dRatioy*10e6,a);
    Constant_C =( b * m_zReference - a)/c;
    C_Constant->Fill(Constant_C);
    Fitted_dRatio->Fill((double)e);
    Fitted_dRatio_Vs_p->Fill((double)P,(double)e);
    Fitted_dRatio_VsXAtZ->Fill((double)a,(double)e);
    //Vs aY
    Fitted_dRatio_VsYatZ->Fill((double)(ay),(double)e);
    Fitted_dRatio_VsbyatZ->Fill((double)by,(double)e);
    Fitted_dRatio_VsbyOveray->Fill((double)(ay/by),(double)e);
    // std::cout<<"ay = "<<ay
    //          <<"by = "<<by
    //          <<"by/ay ="<<ay/by<<std::endl;
    C_Constant_Vs_P->Fill((double)P,(double)Constant_C);
    if(P>1000 && OVTX_Z<9000 && std::abs(PID)!=11)
    {
      //Vs y
      Fitted_dRatio_VsYatZSel->Fill((double)(ay),(double)e);
      Fitted_dRatio_VsbyatZSel->Fill((double)(by),(double)e);
      Fitted_dRatio_VsbyOveraySel->Fill((double)(ay/by),(double)e);
      Fitted_dRatio_Sel->Fill((double)e);
      Fitted_dRatio_Vs_p_Sel->Fill((double)P,(double)e);
      C_Constant_Vs_P_Sel->Fill((double)P,(double)Constant_C);
      C_Constant_Sel->Fill(Constant_C);
      C_Constant_Vs_a->Fill(Constant_C,std::abs(a));
      C_Constant_Vs_b->Fill(Constant_C,b);
      C_Constant_Vs_c->Fill(Constant_C,c);
    }
  }

  //----------------
  // Only X Layers
  //----------------
  std::vector<double> Z;
  Z.push_back((float)Z1);
  Z.push_back((float)Z2);
  Z.push_back((float)Z3);
  Z.push_back((float)Z4);
  Z.push_back((float)Z5);
  Z.push_back((float)Z6);

  //--------------
  //Remove Electrons , Keep Muon(13), Pions(211), Protons(2212),Kaons(321)
  if(P<500) return kTRUE;
  if(OVTX_Z>9000) return kTRUE;
  if( std::abs(PID)==11 || (std::abs(PID)!=211 && std::abs(PID)!=2212 && std::abs(PID)!=321 && std::abs(PID)!=13))
    return kTRUE;
  if ((Y1st>=0 && Y6th<=0)||(Y1st>=0 && Y6th<=0))
    return kTRUE;
  //Fit result plotted for this selection

  m_Total++;
  std::vector<double> X_true;
  std::vector<double> X_smear;
  X_true.push_back(X1st);
  X_true.push_back(X2nd);
  X_true.push_back(X3rd);
  X_true.push_back(X4th);
  X_true.push_back(X5th);
  X_true.push_back(X6th);
  Int_t Case=0;
  double X_T3_2 = X6th +gRandom->Gaus(0,m_sigmaSmearing); //  X in last  Layer
  double X_T1_1 = X1st +gRandom->Gaus(0,m_sigmaSmearing); // X in first Layer
  double Z_T3_2 = Z6; //  Z of last  Layer
  double Z_T1_1 = Z1; // Z of first Layer


  double X_T1_2 = X2nd +gRandom->Gaus(0,m_sigmaSmearing);
  double Z_T1_2 = Z2;
  double X_T2_1 = X3rd +gRandom->Gaus(0,m_sigmaSmearing); // X of 1st Layer in T2
  double Z_T2_1 = Z3; // Z of 1st Layer in T2
  double X_T2_2 = X4th+gRandom->Gaus(0,m_sigmaSmearing); // X of 2nd Layer in T2
  double Z_T2_2 = Z4; // Z of 2nd Layer in T2
  double X_T3_1 = X5th+gRandom->Gaus(0,m_sigmaSmearing);
  double Z_T3_1 = Z5th;
  if(Case==0){
    Xlast= X_T3_2;
    Zlast= Z_T3_2;
    Xfirst= X_T1_2;
    Zfirst= Z_T1_2;
  }
  if(Case==1){
    Xfirst=X_T1_1;
    Zfirst=Z_T1_1;
    Xlast=X_T3_1;
    Zlast=Z_T3_1;
  }
  if(Case==2){
    Xfirst =X_T1_2;
    Xlast  =X_T3_2;
    Zfirst =Z_T1_1;
    Xlast  =Z_T3_2;
  }
  if(Case==0){
    X_smear.push_back((float)Xfirst);
    X_smear.push_back((float)X_T1_2);
    X_smear.push_back((float)X_T2_1);
    X_smear.push_back((float)X_T2_2);
    X_smear.push_back((float)X_T3_1);
    X_smear.push_back((float)Xlast);
  }
  if(Case==1){
    X_smear.push_back((float)Xfirst);
    X_smear.push_back((float)X_T1_2);
    X_smear.push_back((float)X_T2_1);
    X_smear.push_back((float)X_T2_2);
    X_smear.push_back((float)Xlast);
    X_smear.push_back((float)X_T3_2);
  }
  if(Case==1){
    X_smear.push_back((float)XT1_1);
    X_smear.push_back((float)Xfirst);
    X_smear.push_back((float)X_T2_1);
    X_smear.push_back((float)X_T2_2);
    X_smear.push_back((float)X_T3_1);
    X_smear.push_back((float)Xlast);
  }
  for (int i =0;i<6;i++)
  {
    Hit_res->Fill(X_smear[i]-X_true[i]);
  }
  std::vector<Hit> Hits;

  double t_inf = Xfirst/Zfirst;
  double X_last_Projected = Xfirst/Zfirst*Zlast;
  double Distance_LastToInfProjected = Xlast-X_last_Projected;
  //LastInf
  //LastInfDistance = Xlast-X_last_Projected
  double Correction = m_alphaCorrection*(t_inf);
  double Distance_LastToInfProjected_Corrected = Distance_LastToInfProjected - Correction; // Xlast -(m_alphaCoorection*t_inf+Xfirst/Zfirst*Zlast)



  L0_Corrected_1d->Fill(Distance_LastToInfProjected_Corrected);
  L0_1d->Fill(Distance_LastToInfProjected);
  //----------
  //L0 Plots (selection for last layer dependengly on first layer picked position)
  //L0_1d->Write();
  //---------
  L0_0->Fill(X_last_Projected,Distance_LastToInfProjected);
  L0_1->Fill(t_inf,Distance_LastToInfProjected);
  L0_Corrected->Fill(t_inf,Distance_LastToInfProjected_Corrected);
  //To D0 : check for cases different from this and if one can gain from selecting only hits in central region applying another rotation (diagonal one)



  //-----------------L1 step Do the selection
  if( std::abs(Distance_LastToInfProjected_Corrected) > L0_Up) return kTRUE;
  m_L0Selected++;
  double tx_picked =  (Xlast-Xfirst)/(Zlast-Zfirst) ;
  double X_projected_3rd = Xfirst + tx_picked *(Z3-Zfirst);
  double Deltatx = tx_picked - t_inf;
  double Delta = X_T2_1 - X_projected_3rd;
  double x0 = -1*(tx_picked)*Zfirst + Xfirst;

  double CorrX0 = m_x0Corr*x0;
  double DeltaCorr = Delta - CorrX0  ;
  x0vsTxInf->Fill(x0,t_inf);
  L1_0->Fill(x0,Delta);
  L1_0_corr->Fill(x0,DeltaCorr);
  L1_0_1D ->Fill(Delta);
  L1_0_1D_Corr->Fill(DeltaCorr);
  if(x0>0)
  {
    L1_0_x0pos->Fill(tx_picked,Delta);
    L1_1_x0pos->Fill(Deltatx,Delta);
    L1_2_x0pos->Fill(x0,Delta);

    L1_0_x0posCorr->Fill(tx_picked,DeltaCorr);
    L1_1_x0posCorr->Fill(Deltatx,DeltaCorr);
    L1_2_x0posCorr->Fill(x0,DeltaCorr);
  }
  if(x0<0)
  {
    L1_0_x0neg->Fill(tx_picked,Delta);
    L1_1_x0neg->Fill(Deltatx,Delta);
    L1_2_x0neg->Fill(x0,Delta);
    L1_0_x0negCorr->Fill(tx_picked,DeltaCorr);
    L1_1_x0negCorr->Fill(Deltatx,DeltaCorr);
    L1_2_x0negCorr->Fill(x0,DeltaCorr);

  }


  //Do L1 Selection:
  //x0<0 : res_Corr < m_m1(x+m_x0Offset) +m_yOff1 //Upper
  //       res_Corr > m_m2(x+m_x0Offset) +m_yOff2 //Lower
  //x0>0 : res_Corr > m_m2(x-m_x0Offset) -m_yOff2
  //       res_Corr < m_m1(x+m_x0Offset) -m_yOff1
  if(x0<0.)
  {
    double upperBound = m_m1*(x0+m_x0Offset)+m_yOff1;
    double lowerBound = m_m2*(x0+m_x0Offset)+m_yOff2;
    if(DeltaCorr> upperBound || DeltaCorr<lowerBound) return kTRUE;
  }
  if(x0>0.)
  {
    double lowerBound = m_m1*(x0-m_x0Offset)-m_yOff1;
    double upperBound = m_m2*(x0-m_x0Offset)-m_yOff2;
    if(DeltaCorr> upperBound || DeltaCorr<lowerBound) return kTRUE;
  }
  m_L1Selected++;

  //Be carefull in real Life these parameters has to be optimised and look to occupancy (1Hit /8mm expected)

  //================================
  //Selection L2 In this step we will have a 3 Hit Combination (case0_1) 1st,3rd,6th station
  //===============================
  //We Fit the parabola with 3 hits (or we compute the parabola)
  //2 Methods :
  //Compute Parabola in 2 ways
  //Parabola simple
  //Parabola with 3rd order correction
  //Create Method SolveParabola(Hit1,Hit2,Hit3,a,b,c);
  //Create Method SolveParabola2(Hit1,Hit2,Hit3,a,b,c);
  //Try to plot the q/p vs c with SolveParabola1 and SolveParabola2
  //Look to Delta ParabolaExtrapolated Hit_X2, Hit_X4,Hit_X5
  Hit *Hit1 = new Hit(Xfirst,0,Zfirst);
  Hit *Hit2 = new Hit(X_T2_1,0,Z3);
  Hit *Hit3 = new Hit(Xlast,0,Zlast);
  double a_par1=0;
  double b_par1=0;
  double c_par1=0;
  double a_par2=0;
  double b_par2=0;
  double c_par2=0;


  SolveParabola(Hit1,Hit2,Hit3,a_par1,b_par1,c_par1);

  //std::cout<<"a \t "<<a_par1<<"\t b \t "<<b_par1<<"\t c \t "<<c_par1<<std::endl;
  SolveParabola2(Hit1,Hit2,Hit3,a_par2,b_par2,c_par2);
  //std::cout<<"a \t "<<a<<"\t b \t "<<b<<"\t c \t "<<c<<std::endl;
  //std::cout<<"Next Track"<<std::endl;

  Hit *HitT1 = new Hit(X_T1_2,0.,Z_T1_2 );
  Hit *HitT2 = new Hit(X_T2_2,0.,Z_T2_2 );
  Hit *HitT3 = new Hit(X_T3_1,0.,Z_T3_1 );


  //Compute distance from Parabola 1:
  //Extrapolate Hit
  double ExtrapolX_T1_par1 = Extrapol(HitT1->GetZ(),a_par1,b_par1,c_par1);
  double ExtrapolX_T2_par1 = Extrapol(HitT2->GetZ(),a_par1,b_par1,c_par1);
  double ExtrapolX_T3_par1 = Extrapol(HitT3->GetZ(),a_par1,b_par1,c_par1);
  double DeltaX_T1_par1 =  HitT1->GetX()-ExtrapolX_T1_par1;
  double DeltaX_T2_par1 =  HitT2->GetX()-ExtrapolX_T2_par1;
  double DeltaX_T3_par1 =  HitT3->GetX()-ExtrapolX_T3_par1;
  // //Compute distance from Parabola 2:
  double ExtrapolX_T1_par2 =  Extrapol2(HitT1->GetZ(),a_par2,b_par2,c_par2);
  double ExtrapolX_T2_par2 =  Extrapol2(HitT2->GetZ(),a_par2,b_par2,c_par2);
  double ExtrapolX_T3_par2 =  Extrapol2(HitT3->GetZ(),a_par2,b_par2,c_par2);
  double DeltaX_T1_par2 =  HitT1->GetX()-ExtrapolX_T1_par2;
  double DeltaX_T2_par2 =  HitT2->GetX()-ExtrapolX_T2_par2;
  double DeltaX_T3_par2 =  HitT3->GetX()-ExtrapolX_T3_par2;


  L2_0_T1_par1->Fill(DeltaX_T1_par1);
  L2_0_T2_par1->Fill(DeltaX_T2_par1);
  L2_0_T3_par1->Fill(DeltaX_T3_par1);
  L2_0_T1_par2->Fill(DeltaX_T1_par2);
  L2_0_T2_par2->Fill(DeltaX_T2_par2);
  L2_0_T3_par2->Fill(DeltaX_T3_par2);
  // std::cout<<DeltaX_T1_par1<<"\t"<<DeltaX_T2_par1<<"\t"<<DeltaX_T3_par1<<std::endl;
  // std::cout<<DeltaX_T1_par2<<"\t"<<DeltaX_T2_par2<<"\t"<<DeltaX_T3_par2<<std::endl;

  if( std::abs(DeltaX_T1_par1) <  maxDelta_par1 && std::abs(DeltaX_T2_par1) <  maxDelta_par1 && std::abs(DeltaX_T3_par1) <  maxDelta_par1) m_L2_Par1_Selected ++;
  if( std::abs(DeltaX_T1_par2) <  maxDelta_par2 && std::abs(DeltaX_T2_par2) <  maxDelta_par2 && std::abs(DeltaX_T3_par2) <  maxDelta_par2) m_L2_Par2_Selected ++;

  //---------------------
  //L2//Parabola 1 Selection
  //---------------------

  if( std::abs(DeltaX_T1_par1) >  maxDelta_par1 && std::abs(DeltaX_T2_par1) <  maxDelta_par1 && std::abs(DeltaX_T3_par1) >  maxDelta_par1) return kTRUE;

  //LinParFit<double> fit_quartic(5);
  LinParFit<double> fit_parabola(3);
  std::vector<double> x,errx,z;

  //Create the vector of Hits
  // double Xlast  = X6th; //  X in last  Layer
  // double Xfirst = X1st; // X in first Layer
  // double Zlast  = Z6; //  Z of last  Layer
  // double Zfirst = Z1; // Z of first Layer
  // double X_T1_2 = X2nd;
  // double Z_T1_2 = Z2;
  // double X_T2_1 = X3rd; // X of 1st Layer in T2
  // double Z_T2_1 = Z3; // Z of 1st Layer in T2
  // double X_T2_2 = X4th; // X of 2nd Layer in T2
  // double Z_T2_2 = Z4; // Z of 2nd Layer in T2
  // double X_T3_1 = X5th;
  // double Z_T3_1 = Z5;
  x.push_back(Xfirst);
  x.push_back(X_T1_2);
  x.push_back(X_T2_1);
  x.push_back(X_T2_2);
  x.push_back(X_T3_1);
  x.push_back(Xlast);
  z.push_back(Zfirst-m_zReference);
  z.push_back(Z2-m_zReference);
  z.push_back(Z3-m_zReference);
  z.push_back(Z4-m_zReference);
  z.push_back(Z5-m_zReference);
  z.push_back(Zlast-m_zReference);
  std::random_device rd;
  std::mt19937 gen(rd());
  std::discrete_distribution<> d({4, 24, 33, 39});
  // for (int i=0;i<6;i++)
  //  {
  //    Hit Hit(X_smear[i],0.,Z[i]);
  //    Hits.push_back(Hit);
  //  }
  // Track->addHits(Hits);

  // Track->printTrack();


  double ClSize = 0;
  for (int i=0;i<6;i++)
  {
    Hit Hit(X_smear[i],0.,Z[i]);

    ClSize = (double)d(gen);
    ClSize_1->Fill(ClSize+1);
    //GeneratedSize->Fill(ClSize+1.);
    errx.push_back(0.05 + 0.03*(ClSize+1.));
    Hit.SetSize((int)(ClSize+1));
    Hits.push_back(Hit);
  }
  Track->addHits(Hits);
  for (auto ity = std::begin(x), its = std::begin(errx), itx = std::begin(z);
       std::end(x) != ity; ++ity, ++its, ++itx) {
    fit_parabola.accumulate(*ity, *its, *itx);
  }
  std::vector<double> solution;
  solution = fit_parabola.solution();
  double Chi2 = 0;
  for (auto ity = std::begin(x), its = std::begin(errx), itx = std::begin(z);
       std::end(x) != ity; ++ity, ++its, ++itx)
  {
    Chi2+=std::pow(*ity - (solution[0]+solution[1]*(*itx)+solution[2]*(*itx)*(*itx)),2)/((*its)*(*its));
  }
  Chi2_par1->Fill(Chi2);


  bool Ok = fitXProjection(Track);
  if(Ok)
  {
    float C = (Track->b() * m_zReference - Track->a())/(Track->c());
    float X0 = Track->a() - Track->b()*m_zReference +Track->c()*m_ConstC;
    C_par->Fill(C);
    XBackProj->Fill(X0);
    Chi2_par2->Fill(Track->chi2());
  }
  if(Ok)
  {
    m_L2_Selection_ParCubic++;
  //std::cout<<"C\t"<<C<<std::endl;
  //Track->printTrack();

  }
  //Track->printTrack();
  //  std::cout<<"Chi2 = \t "<<Chi2<<std::endl;

  //std::cout<<fit_parabola<<std::endl;


  return kTRUE;

}

void FakedSeeding::SlaveTerminate()
{

}

void FakedSeeding::Terminate()
{
  fFile->Write();
  fFile->Close();
  //TCanvas *c1 = new TCanvas();
  //c1->Print("./Plots/DeltaFakedSeeding_P_5Gev.pdf");
  TFile *f = new TFile("Histogram_FakeSeeding_Phigh.root","recreate");
  std::cout<<"L0 Selection efficiency \t "<<(double)m_L0Selected/m_Total<<std::endl;
  std::cout<<"L1 Selection efficiency \t "<<(double)m_L1Selected/m_Total<<std::endl;
  std::cout<<"L2 Selection Par1 efficiency \t "<<(double)m_L2_Par1_Selected/m_Total<<std::endl;
  std::cout<<"L2 Selection Par2 efficiency \t "<<(double)m_L2_Par2_Selected/m_Total<<std::endl;
  C_par->Write();
  Fitted_dRatio_VsXAtZ->Write();
  XBackProjVsChi2->Write();
  track_distance->Write();
  track_chi2->Write();
  track_pullHitsVsP->Write();
  track_pullHits->Write();
  L0_1d->Write();
  if(doFit)
  {
    //ALL
    Fitted_dRatio_VsYatZ->Write();
    Fitted_dRatio_VsYatZSel->Write();
    Fitted_dRatio_VsbyOveray->Write();
    Fitted_dRatio_VsbyatZ->Write();
    Fitted_dRatio_VsbyatZSel->Write();
    Fitted_dRatio_VsbyOveraySel->Write();


    Fitted_dRatioy->Write();
    Fitted_dRatioyVsC->Write();
    Fitted_dRatioyVsP->Write();
    Fitted_dRatioyVsXFirst->Write();
    Fitted_dRatio_Vs_p->Write();
    Fitted_dRatio->Write();
    C_Constant->Write();
    C_Constant_Vs_c->Write();
    C_Constant_Vs_b->Write();
    C_Constant_Vs_a->Write();
    C_Constant_Vs_P->Write();
    //No electron P>5 Gev
    Fitted_dRatio_Vs_p_Sel->Write();
    Fitted_dRatio_Sel->Write();
    C_Constant_Sel->Write();
    C_Constant_Vs_P_Sel->Write();
    XBackProj->Write();

  }
  TCanvas *c0 = new TCanvas();
  Hit_res->Draw();
  c0->Print("./Plots/ResidualDistribution_Generated.pdf");
  Hit_res->Write();
  x0vsTxInf->Write();
  ClSize_1->Write();
  //=================
  //L1
  //=================
  std::vector<TH2D*> L0_Histos;
  L0_Histos.push_back(L0_0);
  L0_Histos.push_back(L0_1);
  L0_Histos.push_back(L0_Corrected);
  TCanvas *c1 = new TCanvas();
  int size = L0_Histos.size();
  c1->Divide(2,2);
  for(Int_t i=0; i < L0_Histos.size(); i++)
  {
    c1->cd(1+i);
    L0_Histos[i]->Draw("colz");
    L0_Histos[i]->Write();
  }
  c1->cd(2);
  L0_1->Draw("colz");
  std::cout<<"Here_3"<<std::endl;
  L0_Upper_Bound->FixParameter(0,m_alphaCorrection);
  L0_Upper_Bound->FixParameter(1,L0_Up);
  L0_Lower_Bound->FixParameter(0,m_alphaCorrection);
  L0_Lower_Bound->FixParameter(1,L0_Down);
  L0_Upper_Bound->Draw("same");
  L0_Lower_Bound->Draw("same");
  c1->cd(4);
  //L0_Corrected_1d->Sumw2(true);
  //gPad->SetLogy();
  L0_Corrected_1d->Draw();
  L0_Corrected->Write();
  //L0_Corrected_1d->Write();

  std::string outputLocation("./Plots/");
  std::string extension (".pdf");
  std::string FileNameL0 ("L0Plots");
  std::string outL0 = outputLocation+FileNameL0+extension;
  const char *cstrL0 = outL0.c_str();
  c1->Print(cstrL0);
//====================
//L1
//====================
  L1_0_1D->Write();
  L1_0_corr->Write();
  L1_0->Write();
  L1_0_1D_Corr->Write();
  std::vector<TH2D*> L1_x0pos_Histos;
  L1_x0pos_Histos.push_back(L1_0_x0pos);
  L1_x0pos_Histos.push_back(L1_1_x0pos);
  L1_x0pos_Histos.push_back(L1_2_x0pos);
  std::vector<TH2D*> L1_x0neg_Histos;
  L1_x0neg_Histos.push_back(L1_0_x0neg);
  L1_x0neg_Histos.push_back(L1_1_x0neg);
  L1_x0neg_Histos.push_back(L1_2_x0neg);
  //Correction
  std::vector<TH2D*> L1_x0pos_Corr_Histos;
  L1_x0pos_Corr_Histos.push_back(L1_0_x0posCorr);
  L1_x0pos_Corr_Histos.push_back(L1_1_x0posCorr);
  L1_x0pos_Corr_Histos.push_back(L1_2_x0posCorr);
  std::vector<TH2D*> L1_x0neg_Corr_Histos;
  L1_x0neg_Corr_Histos.push_back(L1_0_x0negCorr);
  L1_x0neg_Corr_Histos.push_back(L1_1_x0negCorr);
  L1_x0neg_Corr_Histos.push_back(L1_2_x0negCorr);
  TCanvas *c2 = new TCanvas();
  c2->Divide(3,2);
  for(Int_t i=0; i < L1_x0pos_Histos.size(); i++)
  {
    c2->cd(1+i);
    //L1_x0pos_Histos[i]->SetMarkerColor(kBlue);
    L1_x0pos_Histos[i]->Draw("colz");
    L1_x0pos_Histos[i]->Write();
  }
  for(Int_t i=0; i < L1_x0neg_Histos.size(); i++)
  {
    c2->cd(4+i);
    //L1_x0pos_Histos[i]->SetMarkerColor(kBlue);
    L1_x0neg_Histos[i]->Draw("colz");
    L1_x0neg_Histos[i]->Write();
  }

  TCanvas *c3 = new TCanvas();
  c3->Divide(3,2);
  for(Int_t i=0; i < L1_x0pos_Corr_Histos.size(); i++)
  {
    c3->cd(1+i);
    //L1_x0pos_Histos[i]->SetMarkerColor(kBlue);
    L1_x0pos_Corr_Histos[i]->Draw("colz");
    L1_x0pos_Corr_Histos[i]->Write();
  }
  for(Int_t i=0; i < L1_x0neg_Corr_Histos.size(); i++)
  {
    c3->cd(4+i);
    //L1_x0pos_Histos[i]->SetMarkerColor(kBlue);
    L1_x0neg_Corr_Histos[i]->Draw("colz");
    L1_x0neg_Corr_Histos[i]->Write();
  }
  c3->cd(6);//
  ////X0 correction
  //x0<0 : res_Corr < m_m1(x+m_x0Offset) +m_yOff1 //Upper
  //       res_Corr > m_m2(x+m_x0Offset) +m_yOff2 //Lower
  //x0>0 : res_Corr > m_m2(x-m_x0Offset) -m_yOff2
  //       res_Corr < m_m1(x-m_x0Offset) -m_yOff1
  L1_x0Upper_x0neg->FixParameter(0,m_m1);
  L1_x0Upper_x0neg->FixParameter(1,m_x0Offset);
  L1_x0Upper_x0neg->FixParameter(2,m_yOff1);
  L1_x0Lower_x0neg->FixParameter(0,m_m2);
  L1_x0Lower_x0neg->FixParameter(1,m_x0Offset);
  L1_x0Lower_x0neg->FixParameter(2,m_yOff2);
  L1_x0Upper_x0neg->Draw("same");
  L1_x0Lower_x0neg->Draw("same");
  L1_x0Upper_x0neg->Write();
  L1_x0Lower_x0neg->Write();
  //x0>0
  L1_x0Upper_x0pos->FixParameter(0,m_m1);
  L1_x0Upper_x0pos->FixParameter(1,m_x0Offset);
  L1_x0Upper_x0pos->FixParameter(2,m_yOff1);
  L1_x0Lower_x0pos->FixParameter(0,m_m2);
  L1_x0Lower_x0pos->FixParameter(1,m_x0Offset);
  L1_x0Lower_x0pos->FixParameter(2,m_yOff2);
  c3->cd(3);
  L1_x0Upper_x0pos->Draw("same");
  L1_x0Lower_x0pos->Draw("same");
  L1_x0Lower_x0pos->Write();
  L1_x0Upper_x0pos->Write();
  std::string FileNameL1 ("L1Plots");
  std::string outL1 = outputLocation+FileNameL1+extension;
  const char *cstrL1 = outL1.c_str();
  std::string FileNameL1_Corr ("L1Plots_Corr");
  std::string outL1_Corr = outputLocation+FileNameL1_Corr+extension;
  const char *cstrL1_Corr = outL1_Corr.c_str();
  c2->Print(cstrL1);
  c3->Print(cstrL1_Corr);

  //=======================
  //L2
  //======================
  std::vector<TH1D*> L2_0_Histos;
  L2_0_Histos.push_back(L2_0_T1_par1);
  L2_0_Histos.push_back(L2_0_T1_par2);
  L2_0_Histos.push_back(L2_0_T2_par1);
  L2_0_Histos.push_back(L2_0_T2_par2);
  L2_0_Histos.push_back(L2_0_T3_par1);
  L2_0_Histos.push_back(L2_0_T3_par2);
  TCanvas *c4 = new TCanvas();
  c4->Divide(3,2);
  int j=0;
  for(Int_t i=0; i < (L2_0_Histos.size())/2; i++)
  {
    c4->cd(1+i);
    //L1_x0pos_Histos[i]->SetMarkerColor(kBlue);
    gPad->SetLogy();
    L2_0_Histos[j+1]->SetName("With Par2");
    L2_0_Histos[j+1]->SetLineColor(kBlue);
    L2_0_Histos[j+1]->Draw("");
    L2_0_Histos[j+1]->Write();
    L2_0_Histos[j]->SetLineColor(kRed);
    L2_0_Histos[j]->SetName("With Par1");
    L2_0_Histos[j]->Draw("same");
    L2_0_Histos[j]->Write();

    //gPad->BuildLegend();
    j+=2;
  }
  c4->cd(4);
  gPad->SetLogy();
  Chi2_par1->SetLineColor(kRed);
  Chi2_par1->Draw();
  Chi2_par2->SetLineColor(kBlue);
  Chi2_par2->Draw("same");
  c4->cd(5);
  Hit_res->Draw();
  Chi2_par1->Write();
  Chi2_par2->Write();
  c4->Print("./Plots/L2_0_ResidualsPlots.pdf");
  f->Close();
}




void FakedSeeding::SolveParabola(Hit * hit1,Hit * hit2,Hit * hit3,double &a, double &b,double &c)
{
  const float z1 = hit1->GetZ() - m_zReference;
  const float z2 = hit2->GetZ() -  m_zReference;
  const float z3 = hit3->GetZ() -  m_zReference;
  const float x1 = hit1->GetX();
  const float x2 = hit2->GetX();
  const float x3 = hit3->GetX();
  const float det = (z1*z1)*z2 + z1*(z3*z3) + (z2*z2)*z3 - z2*(z3*z3) - z1*(z2*z2) - z3*(z1*z1);
  if( std::abs(det) < 1e-8 )
  {
    a = 0.0;
    b = 0.0;
    c = 0.0;
    return;
  }
  const float det1 = (x1)*z2 + z1*(x3) + (x2)*z3 - z2*(x3) - z1*(x2) - z3*(x1);
  const float det2 = (z1*z1)*x2 + x1*(z3*z3) + (z2*z2)*x3 - x2*(z3*z3) - x1*(z2*z2) - x3*(z1*z1);
  const float det3 = (z1*z1)*z2*x3 + z1*(z3*z3)*x2 + (z2*z2)*z3*x1 - z2*(z3*z3)*x1 - z1*(z2*z2)*x3 - z3*(z1*z1)*x2;
  a = det1/det;
  b = det2/det;
  c = det3/det;


}

void FakedSeeding::SolveParabola2(Hit * hit1,Hit * hit2,Hit * hit3,double &a, double &b,double &c)
{
  const float z1 = hit1->GetZ() - m_zReference;
  const float z2 = hit2->GetZ() - m_zReference;
  const float z3 = hit3->GetZ() - m_zReference;
  const float x1 = hit1->GetX();
  const float x2 = hit2->GetX();
  const float x3 = hit3->GetX();
  const float e = m_dRatio;
  const float corrZ1 = 1-e*z1;
  const float corrZ2 = 1-e*z2;
  const float corrZ3 = 1-e*z3;
  const float det = (z1*z1)*corrZ1*z2 + z1*(z3*z3)*corrZ3 + (z2*z2)*corrZ2*z3 - z2*(z3*z3)*corrZ3 - z1*(z2*z2)*corrZ2 - z3*(z1*z1)*corrZ1;
  if( std::abs(det) < 1e-8 )
  {
    a = 0.0;
    b = 0.0;
    c = 0.0;
    return;
  }
  const float det1 = (x1)*z2 + z1*(x3) + (x2)*z3 - z2*(x3) - z1*(x2) - z3*(x1);
  const float det2 = (z1*z1)*corrZ1*x2 + x1*(z3*z3)*corrZ3 + (z2*z2)*corrZ2*x3 - x2*(z3*z3)*corrZ3 - x1*(z2*z2)*corrZ2 - x3*(z1*z1)*corrZ1;
  const float det3 = (z1*z1)*corrZ1*z2*x3 + z1*(z3*z3)*corrZ3*x2 + (z2*z2)*corrZ2*z3*x1 - z2*(z3*z3)*corrZ3*x1 - z1*(z2*z2)*corrZ2*x3 - z3*(z1*z1)*corrZ1*x2;
  a = det1/det;
  b = det2/det;
  c = det3/det;

}

double FakedSeeding::Extrapol(double Z,double a,double b,double c)
{
  double zmod = Z-m_zReference;
  //x(z) = a(dz*dz)+b*(dz)+c
  double x_proj = a*(zmod*zmod)+b*zmod+c;
  return x_proj;
}
double FakedSeeding::Extrapol2(double Z1,double a1,double b1,double c1)
{

  double zmod = Z1-m_zReference;
  const float e = m_dRatio;
  double x_proj = a1*(zmod*zmod)*(1-e*zmod)+b1*zmod+c1;
  return x_proj;
}


bool FakedSeeding::fitXProjection(PrSeedTrack * track)
{
  float mat[6];
  float rhs[3];
  //std::fill(rhs,rhs+3,0);
  std::vector<Hit> Hits = track->hits();

  for(int loop = 0;3>loop;++loop)
  {
    std::fill(mat,mat+6,0.);
    std::fill(rhs,rhs+3,0.);
    for( int i=0; i < Hits.size() ;i++ )
    {
      const float w =  Hits[i].w2();//squared
      //std::cout<<"W\t"<<w<<std::endl;
      const float dz = Hits[i].GetZ() - m_zReference;
      float deta = 0;
      deta = dz*dz*(1-m_dRatio*dz);
      Hit *hit = new Hit(Hits[i].GetX(),Hits[i].GetY(),Hits[i].GetZ());
      float dist = track->distance(hit);
      //always()<<"Loop \t"<<loop<<"\n Distance From Hit \t"<<dist<<endmsg;
      // if(loop>0)
      // dist = track.distance( *itH ); //try the effect
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
    if(!decomp)
    {
      //std::cout<<"Failed to decompose matrix"<<std::endl;
      return false;
    }
    decomp.Solve(rhs);
    track->updateParameters(rhs[0],rhs[1],rhs[2],0.,0.);
  }
  float chi2_track = 0.;
  float maxChi2 = 0.;
  float maxDistance = 0.;
  // int i = 0; unused
  //compute Chi2 of the fitted track and the single Hit chi2
  for ( int i=0;i<Hits.size(); i++)
  {
    Hit *hit = new Hit(Hits[i].GetX(),Hits[i].GetY(),Hits[i].GetZ());
    hit->SetW2(Hits[i].w2());

    float distance = track->distance(hit);
    float chi2_onHit = track->chi2( hit ); //\frac{dist^{2}}{\sigma^{2}}
    if (chi2_onHit>maxChi2)
    {
      maxChi2 = chi2_onHit;
    }
    track_distance->Fill(distance);
    track_pullHits->Fill(chi2_onHit);
    chi2_track += track->chi2( hit );
    track_pullHitsVsP->Fill(chi2_onHit,(float)P);
    // //All Hits in 3 sigma? to check externally too
  }
  float constanteC = (track->b() * m_zReference - track->a())/(track->c());
  //backward extrapolation
  float X0 = track->a() - track->b()*m_zReference +track->c()*m_ConstC;
  //track_chi2PerDoFVs
  track->setChi2(chi2_track,3);
  XBackProjVsChi2->Fill(track->chi2(),X0);

  track_chi2->Fill(track->chi2());
  track_chi2PerDoF->Fill(track->chi2()/(3.));
  //std::cout<<"Delta Chi2 (should be 0) \t"<<chi2_track-track->chi2()<<std::endl;
  //if(std::abs(maxDistance) < 2 && std::sqrt(maxChi2) < 4) return true;
  return true;
}
// void FakedSeeding::FITTA(PrSeedTrack * track)
// {
//   std::vector<Hit> Hits = track->hits();
//   std::vector<double> solution;
//   std::vector<double> dsolution;

//   float res = 0.;
//   float dz = 0;
//   for(int j=0;j<3;j++)
//   {
//     LinParFit<double> fit(5);
//     for (int i = 0;i<Hits.size();i++)
//     {
//       dz = Hits[i].GetZ()-m_zReference;
//       res = Hits[i].GetX() - (solution[0]+solution[1]*dz+solution[2]*dz*dz+solution[3]*dz*dz*dz+solution[4]*dz*dz*dz*dz);
//       std::cout<<"Loop\t"<<j<<"\t Residual \t"<<res<<std::endl;
//       fit.accumulate(res,(double)Hits[i].Err(),(double)Hits[i].GetZ()-m_zReference);
//     }
//     //std::vector<double> solution;
//     dsolution=fit.solution();
//     for (int k=0;k<5;k++)
//     {
//       solution[k]=solution[k]+dsolution[k];
//     }
//   }

//   std::cout<<"Par \t Value "<<std::endl;
//   double a = solution[0];
//   double b = solution[1];
//   double c = solution[2];
//   double e = solution[3]/solution[2];
//   double f = solution[4]/(solution[3]*solution[2]);
//   std::cout<<" a   \t"<<a
//            <<"\n b \t"<<b
//            <<"\n c \t"<<c
//            <<"\n e \t"<<e
//            <<"\n f \t"<<f
//            <<std::endl;
//   // for (int i = 0;i<solution.size();i++)
//   // {
//   //   std::cout<<i<<"\t"<<solution[i]<<std::endl;
//   // }
// }
