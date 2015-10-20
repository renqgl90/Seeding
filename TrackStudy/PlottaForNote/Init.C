#include "TString.h"

void Init(bool rebuild=true){
  const int nLibraries=1;
  TString libraries[nLibraries]={TString("CountDoubleCandidate")/*TString("CountDoubleCandidate_DST")*/};
  gROOT->Reset();
  gROOT->SetMacroPath("./src");
  gSystem->SetBuildDir("build_dir",kTRUE);
  gSystem->AddIncludePath("-I./include");
  TString extension=rebuild?TString(".C++"):TString(".C+");
  for (int i=0; i<nLibraries; ++i){
    std::cout<<"Loading "+libraries[i]<<std::endl;
    gROOT->LoadMacro(libraries[i]+extension);
  }
  std::cout<<"Done!"<<endl;
}             
