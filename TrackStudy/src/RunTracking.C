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

  TChain ch1("PrClustersResidual/TrackStudy");
  ch1.Add("./data/SciFi-Tuple-Debug_1000.root");
  TrackStudy *perfectSeeding = new TrackStudy();
  perfectSeeding->Init( (TTree*)  &ch1);
  perfectSeeding->Begin( (TTree*) &ch1);
  Int_t nEntries = ch1.GetEntries();
  for (Int_t i = 0; i<nEntries;i++)
    {
      perfectSeeding->Process(i);
    }
  perfectSeeding->Terminate();
  return 0;
}
