#ifndef PRCHEATEDSCIFITRACKING_H 
#define PRCHEATEDSCIFITRACKING_H 1

// Include files
// from Gaudi

#include "GaudiAlg/GaudiHistoAlg.h"
#include "PrKernel/PrHitManager.h"

/** @class PrCheatedSciFiTracking PrCheatedSciFiTracking.h
 *  
 *  Cheated track reconstruction in the SciFi. 
 *  It creates tracks by getting all Clusters associated to a reconstructible MCParticle. 
 *  Cuts can be set on the minimum number of total hits, hits in x layers and hits in stereo layers.
 *  Top and bottom modules are not mixed in the x layers.
 *
 * - HitManagerName: Name of the hit manager
 * - DecodeData: Decode the data 
 * - OutputName: Name of the output location
 * - NumZones: Number of zones (normally 2 x number of layers if no y-segmentation)
 * - MinXHits: Minimum number of required hits in x layers to make a track.
 * - MinStereoHits: Minimum number of required hits in stereo layers to make a track.
 * - MinTotHits: Minimum number of total hits to make a track.
 * - SpecialCase: If requiring 10(9) Hits, allow also the (4X+6UV) or (6X +4UV) combination
 * - SpecialCase2 : If requiring 9 Hits, allow also the (3X+6UV) or (6UV +3X) combination
 * - NoDoubleCounting: Check that the Hits are required to be in different Layers
 * - CheckReco : Check (nT-i_X>=1 & nT-i_UV>=1) for all the stations
 *  @author Michel De Cian
 *  @date   2015-03-23
 */

class PrCheatedSciFiTracking : public GaudiHistoAlg {
public:
  /// Standard constructor
  PrCheatedSciFiTracking( const std::string& name, ISvcLocator* pSvcLocator );

  virtual ~PrCheatedSciFiTracking( ); ///< Destructor

  virtual StatusCode initialize();    ///< Algorithm initialization
  virtual StatusCode execute   ();    ///< Algorithm execution
  virtual StatusCode finalize  ();    ///< Algorithm finalization

 

protected:

  /// make cheated tracks by getting the clusters matched to an MCParticle
  void makeLHCbTracks( LHCb::Tracks* result );
  void printTrack(std::vector<PrHit*>& track);
  void printHit(const PrHit* hit,std::string title="");
private:
  
  PrHitManager*   m_hitManager;
  
  std::string     m_hitManagerName;
  bool            m_decodeData;
  std::string     m_outputName;
  int             m_numZones;
  int             m_minXHits;
  int             m_minStereoHits;
  int             m_minTotHits;
  bool m_specialCase;
  bool m_noDoubleCounting;
  bool m_CheckReco;
  bool m_specialCase2;

};
#endif // PRCHEATEDSCIFITRACKING_H
