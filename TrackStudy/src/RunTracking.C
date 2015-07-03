#include <iostream>
#include <fstream>
#include "TChain.h"
#include "TSystem.h"
#include "TString.h"
#include "TrackStudy.h"
#include "PatHit.h"
#include "MCHit.h"
//#include "TSelector.h"
#include <TROOT.h>

int main() {

  // TChain ch1("PrClustersResidual/TrackStudy");
  // ch1.Add("./data/SciFi-Tuple-Debug_1000_AllNew.root");
  // TrackStudy *perfectSeeding = new TrackStudy();
  // perfectSeeding->Init( (TTree*)  &ch1);
  // perfectSeeding->Begin( (TTree*) &ch1);
  // Int_t nEntries = ch1.GetEntries();
  // for (Int_t i = 0; i<nEntries;i++)
  //   {
  //     perfectSeeding->Process(i);
  //   }
  // perfectSeeding->Terminate();
  // std::cout<<"With Jacco Modification"<<std::endl;

  TChain ch2("PrClustersResidual/TrackStudy");
  ch2.Add("./data/SciFi-Tuple-Debug_1000_old.root");
  TrackStudy *perfectSeeding1 = new TrackStudy();
  perfectSeeding1->Init( (TTree*)  &ch2);
  perfectSeeding1->Begin( (TTree*) &ch2);
  Int_t nEntries1 = ch2.GetEntries();
  for (Int_t i = 0; i<nEntries1;i++)
    {
      perfectSeeding1->Process(i);
    }
  perfectSeeding1->Terminate();
  std::cout<<"Before Jacco Modification"<<std::endl;


  return 0;
}
