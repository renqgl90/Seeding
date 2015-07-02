// Include files 

// from Gaudi
#include "GaudiKernel/AlgFactory.h"
// local
#include "PrClustersResidual.h"
#include "GaudiKernel/AlgFactory.h" 
#include "Event/MCParticle.h"
#include "Event/MCHit.h"
//#include "Event/MCParticle.h"
#include "Event/MCProperty.h"
// local
#include "Event/MCProperty.h"
#include "Linker/LinkedFrom.h"
#include "Event/MCTrackInfo.h"
#include "Linker/AllLinks.h"
#include "Kernel/FTChannelID.h"
//-----------------------------------------------------------------------------
// Implementation file for class : PrClustersResidual
//
// 2015-01-29 : renato quagliani
//-----------------------------------------------------------------------------

// Declaration of the Algorithm Factory
DECLARE_ALGORITHM_FACTORY( PrClustersResidual )


//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
PrClustersResidual::PrClustersResidual( const std::string& name,
                                        ISvcLocator* pSvcLocator)
: GaudiTupleAlg ( name , pSvcLocator ),
  m_ftHitManager(nullptr),
  m_zone(24)
{
  declareProperty("MCHitsLocation",m_mcHitLocation = "/Event/MC/FT/Hits");
  declareProperty("HitManagerName",m_hitManagerName = "PrFTHitManager");
  always()<<"***In the Constructor \n"
          <<"m_zone"<<m_zone<<endmsg;
}
//=============================================================================
// Destructor
//=============================================================================
PrClustersResidual::~PrClustersResidual() {} 

//=============================================================================
// Initialization
//=============================================================================
StatusCode PrClustersResidual::initialize() {
  StatusCode sc = GaudiTupleAlg::initialize(); // must be executed first
  if ( sc.isFailure() ) return sc;  // error printed already by GaudiAlgorithm
  
  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Initialize" << endmsg;
  m_ftHitManager = tool<PrHitManager>(m_hitManagerName);
  always()<<"Initialize ...."
          <<"Loading PrFTHitManager"<<endmsg;
  
  //====from PrCounter::InitEvent
  // m_validData = false;
  // LHCb::Tracks* tracks = getIfExists<LHCb::Tracks>( m_container );
  // if ( NULL == tracks ) {
  //   if( msgLevel(MSG::DEBUG) ) debug() << "Track container '" << m_container << "' does not exist" <<endmsg;
  //   return;
  // }
  // if ( NULL == m_link ) m_link = new MyAsct( evtSvc(), m_container );
  // m_nbGhost = 0;
  // m_nbTrack = tracks->size();
  // const Table* table = m_link->direct();
  // if ( NULL == table ) { 
  //   Warning( "Problem with MC associations for " + m_container ).ignore();
  //   return; 
  // }
  
  
  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode PrClustersResidual::execute() {
  
  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Execute" << endmsg;
  
  //return StatusCode::SUCCESS;
  LinkedTo<LHCb::MCParticle, LHCb::Track> myForwardLink ( evtSvc(), msgSvc(),LHCb::TrackLocation::Forward);
  LinkedTo<LHCb::MCParticle, LHCb::Track> mySeedLink ( evtSvc(), msgSvc(),LHCb::TrackLocation::Seed);
  LinkedTo<LHCb::MCParticle, LHCb::FTCluster> myClusterLink ( evtSvc(), msgSvc(), LHCb::FTClusterLocation::Default );
  LinkedTo<LHCb::MCHit, LHCb::FTCluster> myFTCluster2MCHitLink ( evtSvc(),msgSvc(), LHCb::FTClusterLocation::Default + "2MCHits");
  MCTrackInfo trackInfo( evtSvc(), msgSvc() );
  Tuples::Tuple tuple = GaudiTupleAlg::nTuple("ClusterMCHitAndTrackStudy","Events");
  
  for (unsigned int zone = 0; zone < m_zone; ++zone)
  {
    bool isX = false;
    bool isU = false;
    bool isV = false;
    bool isT1 = false;
    bool isT2 = false;
    bool isT3 = false;
    
    if(zone==0 || zone == 2 || zone ==6 || zone ==7 || zone ==8 || zone || 9 || zone == 14 || zone == 15 || zone == 16 || zone == 17 || zone == 23)
    {
      isX = true;
      if(zone==0 || zone ==1 || zone == 6 || zone == 7)
        isT1 = true;
      if(zone==8 || zone ==9 || zone == 14 || zone == 15)
        isT2 = true;
      if(zone==16 || zone ==17 || zone == 22 || zone == 23)
        isT3 = true;
    }
    if(zone == 2 || zone ==3 ||zone ==10 || zone==11 || zone == 18|| zone == 19)
    {
      isU = true;
      if (zone >0 && zone <8)
        isT1 = true;
      if (zone >5 && zone <16)
        isT2 = true;
      if (zone >15 && zone <23)
        isT3 = true; 
    }
    
    if(zone == 4 || zone ==5 ||zone ==12 || zone==13 || zone == 20|| zone == 21)
    {
      isV = true;
      if (zone >0 && zone <8)
        isT1 = true;
      if (zone >5 && zone <16)
        isT2 = true;
      if (zone >15 && zone <23)
        isT3 = true;
    }
    
    std::string m_mcHitLocation;
    //int m_zone;
    
    HitRange range = m_ftHitManager->hits(zone);
    for (HitRange::const_iterator iHit = range.begin(); range.end() != iHit; ++iHit) {
      const LHCb::MCParticle* mcPart1 = myClusterLink.first( (*iHit)->id().ftID() );
      if( mcPart1 == nullptr ) continue; //Considering Only hits linked to firt hit in MyClusterLink
      LHCb::MCHit* mcHit = myFTCluster2MCHitLink.first( (*iHit)->id().ftID() );
      Int_t numberMCHitToCluster =0 ;
      while(mcHit != nullptr)
      {
        Gaudi::XYZPoint pMid = mcHit->midPoint();
        //here fill the tuple with watever you want
        numberMCHitToCluster++;
        
        //get Cluster associated to the Hit
        LHCb::FTLiteCluster litecluster = getLiteCluster( (*iHit)->id());
        tuple->column("numberMCHitToCluster",numberMCHitToCluster);
        tuple->column("ClusterCharge",litecluster.charge());
        tuple->column("ClusterSize",litecluster.size());
        tuple->column("ClusterFraction",litecluster.fraction());
        tuple->column("ClusterChannelID",litecluster.channelID());
        tuple->column("ClusterChannelIDSipmCell",litecluster.channelID().sipmCell());
        tuple->column("ClusterChannelIDSipmID",litecluster.channelID().sipmId());
        tuple->column("ClusterChannelIDMat",litecluster.channelID().mat());
        tuple->column("ClusterChannelIDModule",litecluster.channelID().module());
        tuple->column("ClusterChannelLayer",litecluster.channelID().layer());
        tuple->column("ClusterChannelQuarter",litecluster.channelID().quarter());
        tuple->column("isX",isX);
        tuple->column("isU",isU);
        tuple->column("isV",isV);
        tuple->column("isT1",isT1);
        tuple->column("isT2",isT2);
        tuple->column("isT3",isT3);
        // -- As the key of an FTCluster is its channelID, we can link LHCbID and MCParticle directly!
        // -- Caveat: Only take the first link. This might not be fully correct if high precision is needed.
        
        bool isLong = trackInfo.hasVeloAndT( mcPart1 );
        bool isSeed  = trackInfo.hasT( mcPart1 );
        bool accT = trackInfo.accT(mcPart1);
        bool accTT = trackInfo.accTT(mcPart1);
        
        
        tuple->column("MCParticlePID",mcPart1->particleID().pid());
        tuple->column("MCParticleIsLong",isLong);
        tuple->column("MCParticleIsSeed",isSeed);
        tuple->column("MCParticleP",mcPart1->p());
        tuple->column("MCParticlePt",mcPart1->pt());
        tuple->column("MCParticleGamma",mcPart1->gamma());
        tuple->column("MCParticleBeta",mcPart1->beta());
        tuple->column("MCParticleVirtualMass",mcPart1->virtualMass());
        tuple->column("MCParticleCharge",mcPart1->particleID().threeCharge()/3.);  
        tuple->column("MCParticlePseudoRapidity",mcPart1->pseudoRapidity());
        tuple->column("MCParticleAccT",accT);
        tuple->column("MCParticleAccTT",accTT);
        tuple->column("zone",zone);
        //MCHitInfos
        tuple->column("MCHit_X",pMid.X());
        tuple->column("MCHit_Y",pMid.Y());
        tuple->column("MCHit_Z",pMid.Z());
        tuple->column("PrHit_X",(*iHit)->x(pMid.Y()));
        tuple->column("PrHit_Z",(*iHit)->z(pMid.Y()));
        tuple->column("XResidual",(*iHit)->x(pMid.Y())-pMid.X());
        tuple->column("ZResidual",(*iHit)->z(pMid.Y())-pMid.Z());
        tuple->column("Hit_dxDy",(*iHit)->dxDy());
        tuple->column("Hit_werr",(*iHit)->werr());
        tuple->column("Hit_w",(*iHit)->w());
        tuple->column("Hit_coord",(*iHit)->coord());
        tuple->column("Hit_isUsed",(*iHit)->isUsed());
        tuple->column("Hit_yMin",(*iHit)->yMin());
        tuple->column("Hit_yMax",(*iHit)->yMax());
        tuple->column("Hit_Zone",(*iHit)->zone());
        
        
        
        
        //if( !isLong && !isSeed ) continue;//only care about tracks that have at `least a T-Station Segment
        
        //Get the list of tracks generated by the Hit i am lookign to 
        std::vector<const LHCb::Track*> tracksFwd = getTrack( (*iHit)->id(), LHCb::TrackLocation::Forward );
        std::vector<const LHCb::Track*> tracksSeed = getTrack( (*iHit)->id(), LHCb::TrackLocation::Seed );
        //when assocFwd is true we consider it as reconstructed
        //when assocFwd is false we hit belongs to 
        bool isUsedByFwd = false;
        bool assocFwd = false;
        double Chi2NDofFwd = -10;
        int NDoFFwd = -1;
        
        if( !tracksFwd.empty() ){
          isUsedByFwd = true;
          assocFwd = false;
          // -- If one of the tracks associated to this LHCbID is associated to the same MCParticle 
          // -- as the LHCbID is associated, then it is associated...
          // -- (aka: if LHCbID->MCParticle == LHCbID->Track->MCParticle, then the hit was "efficient")
          for(  const LHCb::Track* track : tracksFwd){
            //LHCb::MCParticle* mcPart2 = mySeedLink.first( track->key() ); 
            Chi2NDofFwd = track->chi2PerDoF();
            NDoFFwd = track->nDoF();
            LHCb::MCParticle* mcPart2 = myForwardLink.first( track->key() );
            if( mcPart1 == mcPart2 ){
              assocFwd = true;
              break;
            } 
          }
        }
        double Chi2NDofSeed = -10;
        int NDoFSeed = -1;
        
        bool isUsedBySeed = false;
        bool assocSeed = false;
        if( !tracksSeed.empty() ){
          isUsedBySeed = true;
          assocSeed = false;
          // -- If one of the tracks associated to this LHCbID is associated to the same MCParticle 
          // -- as the LHCbID is associated, then it is associated...
          // -- (aka: if LHCbID->MCParticle == LHCbID->Track->MCParticle, then the hit was "efficient")
          for(  const LHCb::Track* track : tracksSeed){
            Chi2NDofSeed = track->chi2PerDoF();
            NDoFSeed = track->nDoF();
            //LHCb::MCParticle* mcPart2 = mySeedLink.first( track->key() ); 
            LHCb::MCParticle* mcPart2 = mySeedLink.first( track->key() );
            if( mcPart1 == mcPart2 ){
              assocSeed = true;
              break;
            }
          }
        }
        //mcHit = myMCHitLink.next();
        //continue to fill tuples
        
        
        
        tuple->column("TrackChi2NDOFSeed",Chi2NDofSeed);
        tuple->column("TrackNDoFSeed",NDoFSeed);
        
        tuple->column("TrackChi2NDoFFwd",Chi2NDofFwd);
        tuple->column("TrackNDoFFwd",NDoFFwd);
        
        tuple->column("assocSeed",assocSeed);
        tuple->column("isUsedBySeed",isUsedBySeed);
        tuple->column("assocFwd",assocFwd);
        tuple->column("isUsedByFwd",isUsedByFwd);
        StatusCode status = tuple->write();
        mcHit = myFTCluster2MCHitLink.next();
      }
    }
  }
  return StatusCode::SUCCESS;
}

      
      //access all the infos of the MCParticle
      //we have the cluster too
      //Now we should see the MCHit from the FTCluster
      
      

      
      
      
      
      
      
//=============================================================================
//  Finalize
//=============================================================================
StatusCode PrClustersResidual::finalize() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Finalize" << endmsg;

  return GaudiTupleAlg::finalize();  // must be called after all other actions
}



//=============================================================================
//=============================================================================
//  Get the Cluster corresponding to the LHCbID
//=============================================================================
LHCb::FTLiteCluster PrClustersResidual::getLiteCluster(const LHCb::LHCbID id)
{
  //  if( !id.isFT() ) ; 
  LHCb::FTLiteCluster cluster ;
  //  typedef FastClusterContainer<LHCb::FTLiteCluster,int> FTLiteClusters;
  FTLiteClusters* clusters = getIfExists<FTLiteClusters>( LHCb::FTLiteClusterLocation::Default );
  
  if( clusters == nullptr && msgLevel( MSG::ERROR ))
  {
    error() << "Could not find FTLite clusters at: " << LHCb::FTLiteClusterLocation::Default << endmsg;   
  }
  //loop over Hits
  for( FTLiteClusters::const_iterator it = clusters->begin(); it != clusters->end(); it++)
  {
    if( (*it).channelID() == id.ftID())
    {
      cluster = (*it);
      break;
      
    }
  }
  return cluster;
}


//=============================================================================
//  Get the Track(s) corresponding to the LHCbID
//=============================================================================
std::vector<const LHCb::Track*> PrClustersResidual::getTrack(const LHCb::LHCbID id, const std::string location){
  
  LHCb::Track::LHCbIDContainer idCont;
  idCont.push_back( id );
  
  std::vector< const LHCb::Track* > idTracks;
  idTracks.clear();
  
  const LHCb::Tracks* tracks = getIfExists<LHCb::Tracks>( location );
  if( tracks == nullptr ) return idTracks;
  
  for( LHCb::Tracks::const_iterator it = tracks->begin(); it != tracks->end(); it++){

    const LHCb::Track* track = *it;

    if( track->containsLhcbIDs( idCont ) ) idTracks.push_back( track );
    
  }

  return idTracks;
}
