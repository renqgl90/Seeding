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
#include <iostream>
#include "LinParFit.h"
#include "PatHit.h"
#include "Track.h"
#include "MCHit.h"
#include <iostream>
#include <cmath>
#include <cerrno>
#include <cfenv>
#include <cstring>
#include <vector>
#include <string>
#include <iostream>

void TrackStudy::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

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
  std::vector<PatHit> track(100);
  std::vector<MCHit> track_MC(100);
  for (Int_t i = 0 ;  i< CheatedSeeding_NHits; i++){
    PatHit *hit = PatHit();
    hit.setHit(PrHit_Xat0[i],PrHit_Zat0[i],PrHit_dxDy[i],PrHit_dzDy[i],std::sqrt(PrHit_w2),PrHit_yMin[i],PrHit_yMax[i],PrHit_zone[i],PrHit_planeCode[i],PrHit_isX[i]);

    /* code */
  }
  for (Int_t i = 0;  i< N_MCHit_Assoc; i++) {
  MCHit mcHit = MCHit();
  MCHit.setMCHit(MCHit_Assoc_X[i], MCHit_Assoc_Y[i],MCHit_Assoc_Z[i],MCHit_tx[i],MCHit_ty[i],MCHit_p[i],MCHit_pathlength[i],P,MC_px,MC_py,MC_pz);
  }


  // loaded MCParticle
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either TrackStudy::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing
  // keyed objects with PROOF, the object is already loaded and is available
  // via the fObject pointer.
  //
  // This function should contain the "body" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.


   return kTRUE;
}

void TrackStudy::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void TrackStudy::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}
