#include <iostream>
#include <fstream>
#include "TChain.h"
#include "TSystem.h"
#include "TString.h"
#include "FakedSeeding.h"

//#include "TSelector.h"
#include <TROOT.h>

int main() {

  TChain ch1("AllLayers");
  ch1.Add("./data/AllLayers.root");//SelectedSearchWindows.root");
  //TSelector *selector2 = TSelector::GetSelector("./src/FakedSeeding.C");
  FakedSeeding *sel2 = new FakedSeeding();
  sel2->Init( (TTree*)  &ch1);
  sel2->Begin( (TTree*) &ch1);
  Int_t nEntries = ch1.GetEntries();
  //Int_t nEntries = 10;
  //nEntries = ;
  for (Int_t i = 0; i<nEntries;i++)
    {
      sel2->Process(i);
    }
  sel2->Terminate();
  // ch1.Process(sel2);
  // std::cout<<"*****End Of The Story*****"<<std::endl;
  return 0;
  
  // if (doOccupancy)
  //   {
  //     TChain ch1("PrClustersResidual/ClusterMCHitAndTrackStudy");	    
  //     ch1.Add("./data/SciFi-Tuple.root");
  //     TSelector *selector = TSelector::GetSelector("./src/OccupancyStudy.C");
  //     OccupancyStudy *sel = new OccupancyStudy();
  //     sel->Init( (TTree*)  &ch1);
  //     sel->Begin( (TTree*) &ch1);
  //     Int_t nEntries = ch1.GetEntries();
  //     for (Int_t i = 0; i<nEntries;i++)
  // 	{
  // 	  sel->Process(i);
  // 	}
  //     sel->Terminate();
  //     //ch.Process(selector);
  //     std::cout<<"*****End Of The Story*****"<<std::endl;
  //     return 0;
  //   }
  
}
