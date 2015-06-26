// Include files 

// from Gaudi
#include "GaudiKernel/AlgFactory.h"

// local
#include "Event/Track.h"
#include "Event/StateParameters.h"
#include "PrCheatedSciFiTracking.h"
#include "Event/FTCluster.h"
#include "Linker/LinkedTo.h"
#include "Event/MCTrackInfo.h"

//-----------------------------------------------------------------------------
// Implementation file for class : PrCheatedSciFiTracking
//
// 2015-03-23 : Michel De Cian
// 2015-06-26 : Renato Quagliani
//-----------------------------------------------------------------------------

// Declaration of the Algorithm Factory
DECLARE_ALGORITHM_FACTORY( PrCheatedSciFiTracking )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
PrCheatedSciFiTracking::PrCheatedSciFiTracking( const std::string& name,
                                        ISvcLocator* pSvcLocator)
: 
  GaudiHistoAlg ( name, pSvcLocator),
  m_hitManager(nullptr)
{
  declareProperty( "HitManagerName",      m_hitManagerName       = "PrFTHitManager"             );
  declareProperty( "DecodeData",          m_decodeData           = false                        );
  declareProperty( "OutputName",          m_outputName           = LHCb::TrackLocation::Seed    );
  declareProperty( "NumZones",            m_numZones             = 24                           );
  declareProperty( "MinXHits",            m_minXHits             = 5                            );
  declareProperty( "MinStereoHits",       m_minStereoHits        = 5                            );
  declareProperty( "MinTotHits",          m_minTotHits           = 10                           );
  declareProperty( "SpecialCase",         m_specialCase= false);
  declareProperty( "SpecialCase2", m_specialCase2=false);
  declareProperty( "NoDoubleCounting",    m_noDoubleCounting=false);
  declareProperty( "CheckReco", m_CheckReco= false);
}
//=============================================================================
// Destructor
//=============================================================================
PrCheatedSciFiTracking::~PrCheatedSciFiTracking() {}

//=============================================================================
// Initialization
//=============================================================================
StatusCode PrCheatedSciFiTracking::initialize() {
  StatusCode sc = GaudiHistoAlg::initialize(); // must be executed first
  
  if ( sc.isFailure() ) return sc;  // error printed already by GaudiAlgorithm

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Initialize" << endmsg;

  m_hitManager = tool<PrHitManager>( m_hitManagerName );
  m_hitManager->buildGeometry();
  
  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode PrCheatedSciFiTracking::execute() {
  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Execute" << endmsg;
  
  
  LHCb::Tracks* result = new LHCb::Tracks();
  put( result, m_outputName );

  if( m_decodeData ) m_hitManager->decodeData();
  
  makeLHCbTracks( result );

  

  
  return StatusCode::SUCCESS;
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode PrCheatedSciFiTracking::finalize() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Finalize" << endmsg;

  return GaudiHistoAlg::finalize();  // must be called after all other actions
}


//=========================================================================
//  Make the cheated tracks
//=========================================================================
void PrCheatedSciFiTracking::makeLHCbTracks ( LHCb::Tracks* result ) {

  LinkedTo<LHCb::MCParticle, LHCb::FTCluster> myClusterLink ( evtSvc(), msgSvc(), LHCb::FTClusterLocation::Default );
  LHCb::MCParticles* mcParts = getIfExists<LHCb::MCParticles> ( LHCb::MCParticleLocation::Default );

  MCTrackInfo trackInfo( evtSvc(), msgSvc() );

  for( LHCb::MCParticles::const_iterator iPart = mcParts->begin(); iPart != mcParts->end(); ++iPart){
    
    const LHCb::MCParticle* mcPart = *iPart;

    //const bool isLong = trackInfo.hasVeloAndT( mcPart );
    //const bool isDown = trackInfo.hasT( mcPart ) && trackInfo.hasTT( mcPart );
    const bool isSeed = trackInfo.hasT( mcPart );
    std::vector<PrHit*> track;

    //if( !isLong && !isDown ) continue;
    if( !isSeed ) continue;
    LHCb::Track* tmp = new LHCb::Track;
    tmp->setType( LHCb::Track::Ttrack );
    tmp->setHistory( LHCb::Track::PrSeeding );
    tmp->setPatRecStatus( LHCb::Track::PatRecIDs );
    
    std::vector<int> firedXLayers(m_numZones,0);
    std::vector<int> firedStereoLayers(m_numZones,0);
    
    int totHits = 0;
    
    // -- loop over all zones
    for(int iZone = 0; iZone < m_numZones; ++iZone){
      
      PrHitZone* zone = m_hitManager->zone( iZone );
      
      // -- loop over all hits in a zone
      for ( PrHits::iterator itH = zone->hits().begin(); zone->hits().end() != itH; ++itH){
        
        PrHit* hit = *itH;
      
        LHCb::MCParticle* mcPart1 = myClusterLink.first( hit->id().ftID() ); 
        bool found = false;
	
        while( mcPart1 != nullptr){
          if( mcPart1 == mcPart){
            found = true;
          }
          mcPart1 = myClusterLink.next();
        }
      
  
        if( hit->isX() && found ){
          if( firedXLayers[iZone] == 0){
            firedXLayers[iZone]++;
          }
        }
        
	if( !(hit->isX()) && found ){
          if( firedStereoLayers[iZone] == 0){
            firedStereoLayers[iZone]++;
          }
        }
       
        if( found ){
          totHits++;
	  track.push_back(hit);
          tmp->addToLhcbIDs( hit->id() );
        }
        
      }
    }
    
        
    int sumLowerX = 0;
    int sumUpperX = 0;
    int sumStereo = 0;
    int sumLowerX_noDoubleCounting=0;
    int sumUpperX_noDoubleCounting=0;
    int sumStereo_noDoubleCounting=0;
    bool T1_Ok=false;
    bool T2_Ok=false;
    bool T3_Ok=false;
    int nT1_X=0;
    int nT1_UV=0;
    int nT2_X=0;
    int nT2_UV=0;
    int nT3_X=0;
    int nT3_UV=0;
    for(int i = 0; i < m_numZones; i = i+2){
      //if you find at least one hit in a x-layer Lower increase the counter for non-doubleCouting of 1 unit
      if(firedXLayers[i]!=0){
	sumLowerX_noDoubleCounting++;
      }
      sumLowerX += firedXLayers[i];
    }
    
    for(int i = 1; i < m_numZones; i = i+2){
      //if you find at leas one hit in a x-layer Upper increase the counter for non-doubleCounting of 1 unit
      if(firedXLayers[i]!=0){
	sumUpperX_noDoubleCounting++;
      }
      sumUpperX += firedXLayers[i];
    }

    for(int i = 0; i < m_numZones; i++){
      //check both X and UV Layers and depending on which Zone you are increase the counter
      if(firedXLayers[i]!=0 || firedStereoLayers[i]!=0){
	if(i==0|| i==1 || i==6 || i==7) nT1_X++;
	if(i==8|| i==9|| i==14 || i==15) nT2_X++;
	if(i==16 || i==17 || i== 22 ||i==23) nT3_X++;
	if(i==2||i==3||i==4||i==5) nT1_UV++;
	if(i==10 ||i==11 || i==12 ||i==13) nT2_UV++;
	if(i==18 || i==19 || i==20 || i==21) nT3_UV++;
      }
      //if you find at least one hit in a Stereo-Layer increase the counter for non-doubleCounting of 1 unit
      if(firedStereoLayers[i]!=0){
	sumStereo_noDoubleCounting++;
      }
      sumStereo += firedStereoLayers[i];
    }
    
    
    debug() << "sumLowerX: " << sumLowerX 
            << " sumUpperX " << sumUpperX 
            << " sumStereo " << sumStereo 
            << " totHits "   << totHits << endmsg;
    //Is Particle reconstructible? In T1 at least one X and one UV, same for T2,T3
    T1_Ok= (nT1_X>=1 && nT1_UV>=1);
    T2_Ok= (nT2_X>=1 && nT2_UV>=1);
    T3_Ok= (nT3_X>=1 && nT3_UV>=1);
    //Total Number of Hits without counting duplicates
    int TotalHits_noDoubleCounting = nT1_X + nT2_X + nT3_X + nT1_UV + nT2_UV + nT3_UV;
    
    //Check Total Hits without counting twice the X,UV hit if >1 Cluster per layer (Renato)
    //Check with double Counting that the totalNHits is not less than the required one
    //delete the track if the non-doubleCounting is On and the TotalHits_without double counting is less than the minimal number of Hits required
    if(m_noDoubleCounting && TotalHits_noDoubleCounting < m_minTotHits){
      always()<<"Will Remove this track becaust total Layers is less than \t"<<m_minTotHits<<endmsg;
      printTrack(track);
      delete tmp;
      continue;
    }
    //Check the reconstructibily of a track: NX>=1 & NUV>=1 per layer & no DoubleCounting! (Renato)
    if( !(T1_Ok && T2_Ok && T3_Ok) && m_CheckReco){
      always()<<"Will Remove this track because do not satisfy Reconstructible Criteria"<<endmsg;                         
      printTrack(track);
      delete tmp;
      continue;
    }
    
    //Standard counter (Michel)
    if( !m_noDoubleCounting && ((sumLowerX < m_minXHits && sumUpperX < m_minXHits) || sumStereo < m_minStereoHits || totHits < m_minTotHits)){
      always()<<"Will Remove this track because of Michel Conditions"<<endmsg;
      printTrack(track);
      delete tmp;
      continue;
    }

    //No SpecialCases but check you are not doubleCounting a layer
    if( m_noDoubleCounting && !m_specialCase && !m_specialCase2 && ((sumLowerX_noDoubleCounting < m_minXHits && sumUpperX_noDoubleCounting < m_minXHits) || sumStereo < m_minStereoHits || totHits < m_minTotHits)){
      always()<<"Will Remove This Track because of Michel Condition but considering just Duplicates"<<endmsg;
      printTrack(track);
      delete tmp;
      continue;
    }
   
    int MinProduct=m_minXHits*m_minStereoHits;
    if(m_specialCase && m_minTotHits==9) MinProduct =20; //5*4
    if(m_specialCase2 && m_minTotHits==9) MinProduct =18; //6*3
    if(m_specialCase && m_minTotHits==10) MinProduct =24; //6*4
    bool keepTrack =( (sumLowerX_noDoubleCounting * sumStereo_noDoubleCounting>=MinProduct) || (sumUpperX_noDoubleCounting*sumStereo_noDoubleCounting >=MinProduct ));
    //allow for 10 Hit requirement , i.e, 10 Layers no repetition, SumLowerX>=4 & SumStereo ==6 case 
    if( m_noDoubleCounting && (m_specialCase || m_specialCase2) && !keepTrack){
      always()<<"Will Remove This Track because specialCases are not satisfied"<<endmsg;
      printTrack(track);
      delete tmp;
      continue;
    }
    
    
    // -- these are obviously useless numbers, but it should just make the checkers happy
    double qOverP = 0.01;

    LHCb::State tState;
    double z = StateParameters::ZEndT;
    tState.setLocation( LHCb::State::AtT );
    tState.setState( 100, 50, z, 0.1, 0.1, qOverP );

    //tState.setCovariance( m_geoTool->covariance( qOverP ) );
    tmp->addToStates( tState );
    
    result->insert( tmp );
    
    
  }
    
    
    

 

}

void PrCheatedSciFiTracking::printTrack(std::vector<PrHit*> &track){
  for(PrHits::iterator itH= track.begin(); track.end()!=itH;++itH){
    printHit(*itH);
  }
}
void PrCheatedSciFiTracking::printHit( const PrHit* hit,std::string title){
  info() << "  " << title << " "
	 << format( " Plane%3d zone%2d z0 %8.2f x0 %8.2f  size%2d charge%3d coord %8.3f used%2d ",
		    hit->planeCode(), hit->zone(), hit->z(), hit->x(),
		    hit->size(), hit->charge(), hit->coord(), hit->isUsed() ),
    info()<<endmsg;
}
