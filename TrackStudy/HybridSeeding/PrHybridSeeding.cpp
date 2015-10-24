// Include files
// from bost
#include <boost/assign/list_of.hpp>
#include <boost/array.hpp>
#include <boost/foreach.hpp>
// from Gaudi
#include "GaudiKernel/AlgFactory.h"
#include "Event/Track.h"
#include "Event/StateParameters.h"
#include "Math/CholeskyDecomp.h"
// local
#include "PrLineFitterY.h"
#include "PrHybridSeeding.h"
#include "PrPlaneCounter2.h"
#include "Event/MCTrackInfo.h"
#include "Event/MCProperty.h"
#include "Event/LinksByKey.h"
#include "Event/MCHit.h"
#include "Linker/LinkedTo.h"
#include "Linker/LinkedFrom.h"
#include "Linker/AllLinks.h"
#include "Event/FTCluster.h"
#include "Event/FTLiteCluster.h"
#include "Event/MCParticle.h"

#ifdef TRUTH_MATCH_Histos
#include "GaudiAlg/GaudiTupleAlg.h"
#else
#include "GaudiAlg/GaudiAlgorithm.h"
#endif

//-----------------------------------------------------------------------------
// Implementation file for class : PrHybridSeeding
//
// 2015-03-11 : renato quagliani
//-----------------------------------------------------------------------------

// Declaration of the Algorithm Factory
DECLARE_ALGORITHM_FACTORY( PrHybridSeeding )

//=============================================================================
// Standard constructor, initializes variable//=============================================================================

PrHybridSeeding::PrHybridSeeding( const std::string& name,
                                  ISvcLocator* pSvcLocator)
:
#ifdef TRUTH_MATCH_Histos
GaudiTupleAlg(name,pSvcLocator),
#else 
  GaudiAlgorithm(name,pSvcLocator),
#endif
  m_hitManager(nullptr)
  ,m_geoTool(nullptr)
  ,m_debugTool(nullptr)
  ,m_zones(24)
  ,m_timerTool(nullptr)
{
#ifdef TRUTH_MATCH_Histos
  declareProperty("EtaCut",m_etaCut = true);
  declareProperty("NoEle",m_noElectrons = true);
  // declareProperty("RecoverTrack", m_recoverTrack = false);
  // declareProperty("ChiDoFRecover", m_chi2Recover = 1.0);
  
#endif
  
  declareProperty( "RecoverTrack", m_recoverTrack = false);
  declareProperty( "ChiDoFRecover", m_chi2Recover = 1.0);
  declareProperty( "UseCorrPosition" , m_useCorrPos = false);
  declareProperty( "UseCorrSlopes" , m_useCorrSlopes = false);
  declareProperty( "InputName" ,              m_inputName= LHCb::TrackLocation::Forward);
  declareProperty( "OutputName",             m_outputName=LHCb::TrackLocation::Seed);
  declareProperty( "HitManagerName",       m_hitManagerName= "PrFTHitManager");
  declareProperty( "DecodeData",              m_decodeData= false);
  declareProperty( "XOnly",                       m_xOnly= false);
  declareProperty( "MinXPlanes",                m_minXPlanes = 4);
  declareProperty( "MaxNHits",                m_maxNHits =13);
  declareProperty( "NCases" ,                    m_nCases = 3); //Max
  declareProperty( "TimingMeasurement",   m_doTiming= false);
  declareProperty( "PrintSettings",               m_printSettings= true);
  declareProperty( "UseCubicCorrection",     m_useCubic = true);
  declareProperty( "RemoveClones" ,          m_removeClones = true); // to be optimised : track.sorthits LHCbId
  declareProperty( "RemoveClonesUpDown" , m_ClonesUpDown = false);
  declareProperty( "RemoveClonesX",         m_removeClonesX = true); // to be optimised
  declareProperty( "FlagHits",                      m_FlagHits = true);//to be impoved
  declareProperty( "RemoveFlagged",          m_removeFlagged = false);// to be improved
  {
    std::vector<unsigned int> tmp = boost::assign::list_of( 12)(11)(11);
    declareProperty( "SizeToFlag" ,                 m_SizeFlag = tmp);
  }
  //-------------------Flag Hits Settings
  {
    std::vector<double> tmp = boost::assign::list_of( 1.0)(1.0)(1.0);
    declareProperty( "Flag_MaxChi2DoF_11Hits" ,      m_MaxChi2Flag = tmp); //Chi2 Contribution of the Hit to flag
  }
  {
    std::vector<double> tmp = boost::assign::list_of( 400.)(400.)(400.);
    declareProperty( "Flag_MaxX0_11Hits"  , m_MaxX0Flag = tmp ); //Track Backward projection max value for flagging
  }
  //--------------------X-Search Parametrisation
  //1st / Last Layer search windows  
  declareProperty( "removeHp", m_removeHighP = false);
  {
    std::vector<double> tmp=boost::assign::list_of(120.64)(120.64)(120.64);
    declareProperty( "L0_AlphaCorr" , m_alphaCorrection = tmp);
  }
  {
    std::vector<double> tmp=boost::assign::list_of(280.0)(540.0)(840.0);
    declareProperty( "L0_tolHp" , m_TolFirstLast = tmp);
  }
  
  {
    std::vector<double> tmp = boost::assign::list_of(1.2)(1.2)(1.2);
    declareProperty( "TolXRemaining", m_tolRemaining = tmp);
  }
  {
    std::vector<double> X0Rotation = boost::assign::list_of(0.002134)(0.001534)(0.001834);
    declareProperty( "x0Corr",  m_x0Corr = X0Rotation);  // Rotation angle
  }
  {
    std::vector<double> X0SlopeChange = boost::assign::list_of(400.)(500.)(500.);
    declareProperty( "X0SlopeChange", m_x0SlopeChange = X0SlopeChange);
  }
  {
    std::vector<double> TolX0SignUp = boost::assign::list_of(0.75)(0.75)(0.75);
    declareProperty( "ToleranceX0Up", m_TolX0SameSign = TolX0SignUp);
  }
  {
    std::vector<double> x0Cut = boost::assign::list_of( 4000.)( 4000.)(4000.);
    declareProperty( "x0Cut", m_x0Cut = x0Cut);
  }
  {
    std::vector<double> TolAtX0Cut = boost::assign::list_of(12.0)(12.0)(12.0);
    declareProperty( "TolAtX0Cut" , m_tolAtX0Cut = TolAtX0Cut);
  }
  {
    std::vector<double> x0SlopeChange2 = boost::assign::list_of(1500.)(1500.)(1500.);
    declareProperty( "X0SlopeChangeDown" , m_x0SlopeChange2 = x0SlopeChange2);
  }
  {
    //Tolerance inferior for |x0| > m_x0SlopeChange2 when x0 = m_x0Cut
    std::vector<double> tolAtx0CutOpp = boost::assign::list_of(1.5)(2.0)(7.0);
    declareProperty( "TolAtX0CutOpp" , m_tolAtx0CutOppSig = tolAtx0CutOpp);
  }
  {
    std::vector<double> tolOpp = boost::assign::list_of(0.75)(0.75)(0.75);
    declareProperty( "ToleranceX0Down", m_tolX0Oppsig = tolOpp);
  }
  
  
  declareProperty( "maxParabolaSeedHits", m_maxParabolaSeedHits = 8); // for Case 1
  //Look up in remaining X layers and Track Fitting  X
  declareProperty( "dRatio",                    m_dRatio= -0.000262);
  declareProperty( "CConst" ,                    m_ConstC = 2.458e8); //Backward Projection
  //--------------------UV search  Parametrisation (inherit from PrSeedingXLayers ( to be modified)
  declareProperty("TolCoord",                    m_coord = 0.005);
  declareProperty("RemoveHole",              m_removeHole = true);
  //---Added
  declareProperty("yMin"   ,         m_yMin = -1.0 *Gaudi::Units::mm);
  declareProperty("yMin_TrFix",      m_yMin_TrFix = -2.0 * Gaudi::Units::mm);  
  declareProperty("yMax"   ,         m_yMax = 2700.*Gaudi::Units::mm);
  declareProperty("yMax_TrFix",      m_yMax_TrFix = +30.0 * Gaudi::Units::mm);
  //---Line Parameters
  declareProperty("UseLineY", m_useLineY = false);
  {
    std::vector<double> Chi2LowLine = boost::assign::list_of( 5.0)( 6.0)( 7.0);
    declareProperty("Chi2LowLine" ,  m_Chi2LowLine = Chi2LowLine);
  }
  {
    std::vector<double> Chi2HighLine = boost::assign::list_of( 30.0)(50.0)(80.0);
    declareProperty("Chi2HighLine", m_Chi2HighLine = Chi2HighLine);
  }
  {  
    //Unused for the moment (use the backward projection to define Hough cluster)
    std::vector<double> tmp = boost::assign::list_of(200.)(200.)(200.);
    declareProperty( "X0ChangeCoord", m_X0ChangeCoord = tmp);
  }
  {
    std::vector<double> tmp = boost::assign::list_of(0.002)(0.002)(0.0035);
    declareProperty( "TolTyOffset", m_tolTyOffset = tmp);
  }
  {
    std::vector<double> tmp = boost::assign::list_of(0.0)(0.0)(0.015);
    declareProperty( "TolTySlope", m_tolTySlope = tmp);
  }

  declareProperty("DoAsymm"       ,           m_doAsymmUV = true);
  //---
  declareProperty( "TriangleFix" ,              m_useFix               = true);
  declareProperty( "TriangleFix2ndOrder",   m_useFix2ndOrder = true);
  
  {
    std::vector<double> maxChi2X = boost::assign::list_of( 5.5)(5.5)(5.5);
    declareProperty( "maxChi2HitsX",          m_maxChi2HitsX = maxChi2X );
  }
  {
    std::vector<double> maxChi2DoFX = boost::assign::list_of(4.0)(5.0)(6.0);
    declareProperty( "maxChi2DoFX" , m_maxChi2DoFX = maxChi2DoFX);
  }
  
  {
    std::vector<double> maxChi2FullFit_11and12Hit = boost::assign::list_of(5.5)(5.5)(5.5);
    declareProperty( "maxChi2Hits_11and12Hit" , m_maxChi2HitFullFitHigh = maxChi2FullFit_11and12Hit);
  }
  {
    std::vector<double> maxChi2FullFit_less11Hit = boost::assign::list_of(2.5)(2.5)(2.5);
    declareProperty( "maxChi2Hits_less11Hit" , m_maxChi2HitFullFitLow = maxChi2FullFit_less11Hit);
  }
  {
    std::vector<double> maxY0 = boost::assign::list_of( 50.)(50.)(50.);
    declareProperty( "maxYatZeroLow" , m_maxY0Low = maxY0);
  }
  {
    std::vector<double> maxyRef = boost::assign::list_of( 500.)(500.)(500.);
    declareProperty( "maxYatzRefLow" , m_maxYZrefLow = maxyRef);
  }
  {
    std::vector<double> chi = boost::assign::list_of( 8.0)(8.0)(8.0);
    declareProperty( "maxChi2PerDoF",       m_maxChi2PerDoF = chi );
  }
  
  declareProperty( "UseCorrectionDz" ,     m_corrLeftSide = false);
  
  //-------------------Clone Removal Settings
  {
    std::vector<unsigned int> nCommonX = boost::assign::list_of( 1)( 2 )( 3);
    declareProperty( "minCommonX", m_nCommonX = nCommonX);
  }
  declareProperty( "minNCommonUV",     m_nCommonUV=7); //for clone removal in Look Up UV Hits
  declareProperty( "minNCommonUVUpDown" , m_nCommonUVTriangle = 2);
  // Parameters for debugging
}
//=============================================================================
// Destructor
//=============================================================================
PrHybridSeeding::~PrHybridSeeding() {}

//=============================================================================
// Initialization
//=============================================================================
StatusCode PrHybridSeeding::initialize() {
#ifdef TRUTH_MATCH_Histos
  StatusCode sc = GaudiTupleAlg::initialize();
#else
  StatusCode sc = GaudiAlgorithm::initialize();
#endif
  if ( sc.isFailure() ) return sc;  // error printed already by GaudiAlgorithm
  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Initialize" << endmsg;
  m_hitManager= tool<PrHitManager>( m_hitManagerName );
  m_geoTool = tool<PrGeometryTool>("PrGeometryTool");
  m_hitManager->buildGeometry();
  m_debugTool = 0;
  if ( "" != m_debugToolName ) {
    m_debugTool = tool<IPrDebugTool>( m_debugToolName );
    info()<<"Debug tool "<<m_debugToolName<<" loaded."<<endmsg;
  } else {
    m_wantedKey = -100;  //no debug
  }
  if ( m_doTiming) {
    m_timerTool = tool<ISequencerTimerTool>( "SequencerTimerTool/Timer", this );
    m_timeTotal   = m_timerTool->addTimer( "PrSeeding total" );
    m_timerTool->increaseIndent(); //+1
    m_timeFromForward = m_timerTool->addTimer( "Time from Forward" );
    m_timerTool->increaseIndent(); //+2
    m_timeXProjeUp[0] = m_timerTool->addTimer( "Case 0 : X Projections Up");
    m_timerTool->increaseIndent(); //+3
    m_timeCloneXUp[0] = m_timerTool->addTimer( "Case 0 : Clones X Up");
    m_timerTool->increaseIndent(); //+4
    m_timeStereoUp[0] = m_timerTool->addTimer( "Case 0 : AddStereo Up");
    m_timerTool->increaseIndent();//+5
    m_timeFlagUp[0] = m_timerTool->addTimer( " Case 0 : Flag Up");
    m_timerTool->decreaseIndent(); // +4
    m_timerTool->decreaseIndent(); // + 3
    m_timerTool->decreaseIndent(); //+ 2
    m_timeXProjeUp[1] = m_timerTool->addTimer( "Case 1 : X Projections Up");
    m_timerTool->increaseIndent(); //+ 3
    m_timeCloneXUp[1] = m_timerTool->addTimer( "Case 1 : Clones X Up");                                                               
    m_timerTool->increaseIndent();      //+4                  
    m_timeStereoUp[1] = m_timerTool->addTimer( "Case 1 : AddStereo Up");                                                                                                                                                    
    m_timerTool->increaseIndent();          //+5                                                
    m_timeFlagUp[1] = m_timerTool->addTimer( "Case 1 : Flag Up");
    m_timerTool->decreaseIndent();//+4
    m_timerTool->decreaseIndent();//+3                                                           
    m_timerTool->decreaseIndent();//+2
    m_timeXProjeUp[2] = m_timerTool->addTimer( "Case 2 : X Projections Up");                                                               
    m_timerTool->increaseIndent();//+3
    m_timeCloneXUp[2] = m_timerTool->addTimer( "Case 2 : Clones X Up");                                                               
    m_timerTool->increaseIndent();//+4
    m_timeStereoUp[2] = m_timerTool->addTimer( "Case 2 : AddStereo Up");                                                             
    m_timerTool->increaseIndent();//+5
    m_timeFlagUp[2] = m_timerTool->addTimer( "Case 2 : Flag Up");
    m_timerTool->decreaseIndent();//+4
    m_timerTool->decreaseIndent();//+3
    m_timerTool->decreaseIndent();//+2
    
    m_timeXProjeDo[0] = m_timerTool->addTimer( "Case 0 : X Projections Down");                                                                 
    m_timerTool->increaseIndent();//+3
    m_timeCloneXDo[0] = m_timerTool->addTimer( "Case 0 : Clones X Down");                                                                 
    m_timerTool->increaseIndent();//+4
    m_timeStereoDo[0] = m_timerTool->addTimer( "Case 0 : AddStereo Dowm");                                                               
    m_timerTool->increaseIndent();    //+5                           
    m_timeFlagDo[0] = m_timerTool->addTimer( " Case 0 : Flag Down");
    m_timerTool->decreaseIndent();        //+4                                                      
    m_timerTool->decreaseIndent();            //+3
    m_timerTool->decreaseIndent();//+2
    m_timeXProjeDo[1] = m_timerTool->addTimer( "Case 1 : X Projections Down");                                                                 
    m_timerTool->increaseIndent();    //+3
    m_timeCloneXDo[1] = m_timerTool->addTimer( "Case 1 : Clones X Down");
    m_timerTool->increaseIndent();        //+4
    m_timeStereoDo[1] = m_timerTool->addTimer( "Case 1 : AddStereo Down");   
    m_timerTool->increaseIndent(); //+5
    m_timeFlagDo[1] = m_timerTool->addTimer( "Case 1 : Flag Down");
    m_timerTool->decreaseIndent(); //+4
    m_timerTool->decreaseIndent(); //+3
    m_timerTool->decreaseIndent(); //+2
    m_timeXProjeUp[2] = m_timerTool->addTimer( "Case 2 : X Projections Up");
    m_timerTool->increaseIndent(); //+3
    m_timeCloneXUp[2] = m_timerTool->addTimer( "Case 2 : Clones X Up");
    m_timerTool->increaseIndent(); //+4
    m_timeStereoUp[2] = m_timerTool->addTimer( "Case 2 : AddStereo Up");
    m_timerTool->increaseIndent(); //+5                                                  
    m_timeFlagUp[2] = m_timerTool->addTimer( "Case 2 : Flag Up");
    m_timerTool->decreaseIndent(); //+4
    m_timerTool->decreaseIndent(); //+3
    m_timerTool->decreaseIndent(); //+2
    
    m_timerTool->decreaseIndent();//+1
    m_timeClone = m_timerTool->addTimer( "Remove Clones");
    m_timeConvert = m_timerTool->addTimer("Convert Tracks");
    
    // m_timeStereo      = m_timerTool->addTimer( "Add stereo" );
    // m_timeFinal       = m_timerTool->addTimer( "Convert tracks" );
    m_timerTool->decreaseIndent();//0
  }
  if( m_decodeData )         info() << "Will decode the FT clusters!" << endmsg;
  if( m_FlagHits)                 info()<<"Will Flag the Hits" << endmsg;
  if( m_removeFlagged)      info()<<"Will Not re-use Flagged"<<endmsg;
  if( m_inputName == "")    info()<<"Standalone Seeding"<<endmsg;
  if( !(m_inputName == "")) info()<<"Forward tracks as input"<<endmsg;
  if( m_removeClones)      info()<<"Will Remove Clones"<<endmsg;
  if( m_xOnly) info()<<"Will use Only X Layers" <<endmsg;
  if( m_useFix && !m_xOnly) info()<<"WIll Use the triangle Fix"<<endmsg;
  if( m_recoverTrack) info()<<"Will Recover the track at the end (to be implemented)"<<endmsg;
  
  if( m_printSettings){
    info() <<"==================================================="<<endmsg
           << "===============GLOBAL SETTINGS===================="<<endmsg
           <<" InputName                                       = "<<  m_inputName              << endmsg
           << " OutputName                                     = "<<  m_outputName             << endmsg
           << " HitManagerName                                 = "<<  m_hitManagerName         << endmsg
           << " DecodeData                                     = "<<  m_decodeData             << endmsg
           << " XOnly                                          = "<<  m_xOnly                  << endmsg
           << " MinXPlanes                                     = "<<  m_minXPlanes             << endmsg
           << " MaxNHits                                       = "<<  m_maxNHits               << endmsg
           << " NCases                                         = "<<  m_nCases                 << endmsg
           << " TimingMeasurement                              = "<<  m_doTiming               << endmsg
           << " dRatiio                                        = "<< m_dRatio                 << endmsg
           << " CConst BackProj                                = "<< m_ConstC<<endmsg
           << "=================ADD UV Layer Settings============"<< endmsg
           << "doAsymmUV                                       = "<< m_doAsymmUV  <<endmsg
           << "Use Triangle Fix                                = " << m_useFix                  << endmsg
           << "Use SecondOrder `Triangle Fixing                = " << m_useFix2ndOrder <<endmsg
           << "*******Fit Settings"<<endmsg
           << "==================Clone Removal and Flag Hits Settings=========="<<endmsg
           << "Remove Clones after X searh                    = " <<m_removeClonesX<<endmsg
           << "RemoveClones after add stereo UV               = " << m_removeClones<<endmsg
           << "Min Hits Common Total                         = " << m_nCommonUV <<endmsg
           << "Compare Up and Down                           = " << m_ClonesUpDown << endmsg
    
           << "Min Hits Common Total when Up and Down compared = " << m_nCommonUVTriangle<< endmsg
           << "Flag the hits                                  = " << m_FlagHits <<endmsg
           << "Remove All Flagged                             = " << m_removeFlagged<<endmsg;
    // << "Flag_MaxChi2 Hit                               = " << m_MaxChi2Flag<<endmsg
    // << "Flag_MaxX0                                     = " <<m_MaxX0Flag <<endmsg;
    
    info()<<"  Will Run N Cases = " << m_nCases<< endmsg;
    
    info()<<" ============================ Find X Projection Settings ======================="<<endmsg;
    info()<<" 1 - First-Last Layer Settings  (2 Hit Combo)  "<<endmsg;
    for(unsigned int kk =0;m_alphaCorrection.size()>kk;kk++){
      info()<<"\t Case "<<kk<<"   Rotation                  ="<<m_alphaCorrection[kk]<<endmsg;
      info()<<"\t Case "<<kk<<"   Tolerance after Rotation  ="<<m_TolFirstLast[kk]<<endmsg;
    }
    
    info()<<" 2 - Add Hit in T2 - XLayers Stations  (3 Hit Combo) "<<endmsg;
    info()<<" Allow N Parabola Seed Hits = " << m_maxParabolaSeedHits;
    for(unsigned int kk =0;m_x0Corr.size() > kk; kk++)
    {
      info()<<"\t Case "<<kk<<"       x0Rotation = "<<m_x0Corr[kk]<<endmsg;
      info()<<"\t Case "<<kk<<"    x0SlopeChange = "<<m_x0SlopeChange[kk]<<"    Tolerance at x0SlopeChange "<<  m_TolX0SameSign[kk]<<endmsg;
      info()<<"\t Case "<<kk<<"   x0SlopeChange2 = "<<m_x0SlopeChange2[kk]<<"    Tolerance at x0SlopeChange2"<< m_tolX0Oppsig[kk]<<endmsg;
      info()<<"\t Case "<<kk<<"             x0Max = "<<m_x0Cut[kk]<<"    Tolerance Up at x0Max"<< m_tolAtX0Cut[kk]<< " Tolerance Down at x0Max "<<m_tolAtx0CutOppSig[kk]<<endmsg;
    }
    info()<<" 3 - Add Hit in remaining Layers (4/5/6 Hit Combo) ---- "<<endmsg;
    for( unsigned int kk =0; m_tolRemaining.size()>kk;kk++){
      info()<<"\t Case"<<kk<<"    Tolerance = "<<m_tolRemaining[kk]<<endmsg;
    }
    info()<<" 4 - Fit XZ Projection  "<<endmsg;
    info()<<"\t minXPlanes "<< m_minXPlanes;
    info()<<"\t dRatio " << m_dRatio<<endmsg;
    for( unsigned int kk =0; m_maxChi2HitsX.size()>kk;kk++){      
      info()<<"\t Case "<<kk<<"      MaxChi2Hit X Fit    = "<<m_maxChi2HitsX[kk]<<endmsg;
      info()<<"\t Case "<<kk<<"      MaxChi2Track X Fit  = "<<m_maxChi2DoFX[kk]<<endmsg;
    }
    info()<<" Remove Clones X " << m_removeClonesX<< endmsg;
    
    info()<<" ========================= Add The UV part =================================== "<<endmsg;
    info()<<" 1 - Hit Selection pre-Hough Clustering "<<endmsg;
    info()<<"\t Y Min                 "<< m_yMin      <<endmsg
          <<"\t Y Max                 "<< m_yMax      <<endmsg
          <<"\t UseTrFix              "<<m_useFix     <<endmsg
          <<"\t Y Min TrFix           "<< m_yMin_TrFix<<endmsg
          <<"\t Y Max TrFix           "<< m_yMax_TrFix<<endmsg
          <<"\t Asymmetric UV search  "<< m_doAsymmUV  <<endmsg
          <<"\t RemoveHole            "<<m_removeHole <<endmsg
          <<"\t TrFix 2nd Order       "<<m_useFix2ndOrder<<endmsg
          <<endmsg;
    
    info()<<" 2 - Hough Cluster in UV Layers (y/z) hits "<<endmsg;
    for( unsigned int kk = 0; m_tolTyOffset.size()>kk; kk++){
      info()<<"\t Case "<<kk<<" maxCoord-minCoord < " <<m_tolTyOffset[kk]<<" + "<< m_tolTySlope[kk] <<" * minCoord"<<endmsg;
    }
    
    info()<<" 3 - Do preliminary Y line fit "<< m_useLineY;
    if( m_useLineY){
      for( unsigned int kk = 0;  m_Chi2LowLine.size()>kk; kk++){
        info()<<"\t Case"<<kk<<"max   Chi2DoFX + Chi2Line (<11) hit  "<<m_Chi2LowLine[kk]<<endmsg;
        info()<<"\t Case"<<kk<<"max   Chi2DoFX + Chi2Line (>10) hit  "<<m_Chi2HighLine[kk]<<endmsg;
      }
    }
    info()<<" 4 - Fit Full Track "<<endmsg;
    info()<<" Force track to have maxNHits ( if UseLineY, auto 12)  = " << m_maxNHits<<endmsg;
    info()<<" Use dRatioCorrection with XY position " << m_useCorrPos<<endmsg;
    info()<<" Use dRatioCorrection with Slopes      " << m_useCorrSlopes<<endmsg;
    info()<<"\t \t Hits <11 "<<endmsg;
    for( unsigned int kk = 0;  m_maxY0Low.size() >kk; kk++){
      info()<<"\t Case"<<kk<<" maxY at 0 after Fit    "<< m_maxY0Low[kk]<<endmsg;
      info()<<"\t Case"<<kk<<" maxY at zRef after Fit "<<m_maxYZrefLow[kk]<<endmsg;
      info()<<"\t Case"<<kk<<" MaxChi2Hit             "<<m_maxChi2HitFullFitLow[kk]<<endmsg;
    }
    
    info()<<"\t \t Hits >=11 "<<endmsg;
    for( unsigned int kk =0; m_maxChi2HitFull.size()>kk;kk++){
      info()<<"\t Case"<<kk<<" MaxChi2Hit             "<<m_maxChi2HitFullFitHigh[kk]<<endmsg;
    }
    info()<<" 5 - FinalChi2DoF      "<<endmsg;
    for( unsigned int kk = 0; m_maxChi2FullFit.size()>kk;kk++){
      info()<<"\t Case"<<kk<<" MaxChi2PerDoF Total   "<< m_maxChi2FullFit[kk]<<endmsg;
    }
    info()<<" ====================== Flag Hits ============================== " << endmsg;
    for( unsigned int kk = 0; m_SizeFlag.size()>kk;kk++){
      info()<<"\t Case"<<kk<<" Size Flag"<< m_SizeFlag[kk]<<endmsg;
    }
  }

  
#ifdef TRUTH_MATCH_Histos
  setHistoTopDir("FT/");
#endif
  return StatusCode::SUCCESS;
}


//=============================================================================
// Main execution
//=============================================================================
StatusCode PrHybridSeeding::execute() {
#ifdef TRUTH_MATCH_Histos
  using namespace Tuples;
#endif
  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Execute" << endmsg;
  if ( m_doTiming )
  {
    m_timerTool->start( m_timeTotal );
    m_timerTool->start( m_timeFromForward );
  }
  LHCb::Tracks* result = new LHCb::Tracks();
  put( result, m_outputName );
  if( m_decodeData ){
    debug()<<"Decoding Data"<<endmsg;
    m_hitManager->decodeData();
    debug()<<"Decoding Done"<<endmsg;
  }
  //  int multiplicity[24];
  int multiplicityTot = 0;
  //char zoneName[100];
  // UNFLAG USED TO ALL HITS
  debug()<<"UNFLAGGING ALL HITS"<<endmsg;
  for ( unsigned int zone = 0; m_hitManager->nbZones() > zone; ++zone ) {
    //multiplicity[zone]=0;
    //sprintf (zoneName, "Multiplicity_InZone_%i",zone);
    for ( PrHits::const_iterator itH = m_hitManager->hits( zone ).begin();
          m_hitManager->hits( zone ).end() != itH; ++itH ) {
      (*itH)->setUsed( false );
      //multiplicity[zone]++;
      multiplicityTot++;
    }
  }
  //========================================================
  // Remove Seed segment if we pass the  Forward as Input
  //========================================================
  if ( "" != m_inputName){
    debug()<<"Removing Seed Segment from Forward Tracking"<<endmsg;
    for(int i = 0; i < 24; i++){
      PrHitZone* zone = m_hitManager->zone(i);
      std::stable_sort( zone->hits().begin(),  zone->hits().end(), compLHCbID());
    }
    LHCb::Tracks* forward = get<LHCb::Tracks>( m_inputName );
    //Loop over forward tracks container
    for ( LHCb::Tracks::iterator itT = forward->begin(); forward->end() != itT; ++itT ) {
      //Vector of LHCbID
      std::vector<LHCb::LHCbID> ids;
      ids.reserve(20);
      // Loop over LHCbIDs of the forward track
      for ( std::vector<LHCb::LHCbID>::const_iterator itId = (*itT)->lhcbIDs().begin();
            (*itT)->lhcbIDs().end() != itId; ++itId ) {
        if ( (*itId).isFT() && (*itId).ftID().layer() < 12 ) {
          LHCb::FTChannelID ftId =(*itId).ftID();
          int zoneNb = 2 * ftId.layer() + ftId.mat(); //zones top are even (0, 2, 4, ....,22)  and zones bottom are odd
          //Load the PrHitZone
          PrHitZone* zone = m_hitManager->zone(zoneNb);

          // -- The hits are sorted according to LHCbID, we can therefore use a lower bound to speed up the search

          PrHits::iterator itH = std::lower_bound(  zone->hits().begin(),  zone->hits().begin(), *itId, lowerBoundLHCbID() );

          for ( ; zone->hits().end() != itH; ++itH ) {
            if( *itId < (*itH)->id() ) break;
            if ( (*itH)->id() == *itId ) (*itH)->setUsed( true );
          }
          ids.push_back( *itId );
        }
      }

      // Forward output loaded to TTrack container with History flag
      LHCb::Track* seed = new LHCb::Track;
      seed->setLhcbIDs( ids );
      seed->setType( LHCb::Track::Ttrack );
      seed->setHistory( LHCb::Track::PrSeeding );
      seed->setPatRecStatus( LHCb::Track::PatRecIDs );
      seed->addToStates( (*itT)->closestState( 9000. ) );
      result->insert( seed );
    }



    //=======================================
    // Sort hits according to x for each zone
    //=======================================
    for(int i = 0; i < 24; i++){
      PrHitZone* zone = m_hitManager->zone(i);
      std::stable_sort( zone->hits().begin(),  zone->hits().end(), compX());
    }
  }
  //=======================================
  // Fill zones
  //=======================================
  m_zones.clear();
  for(int i = 0; i < 24; i++){
    m_zones.push_back( m_hitManager->zone(i) );
  }
  //==========================================================
  //END FLAGGING HITS FROM FORWARD
  //==========================================================
  // Hits are ready to be processed
  //==========================================================
  m_trackCandidates.clear();
  if ( m_doTiming ) {
    m_timerTool->stop( m_timeFromForward );
  }
  //========================================================
  //------------------MAIN SEQUENCE IS HERE-----------------
  //========================================================

  // ----- Loop through lower and upper half
  //for( unsigned int part= 0; 2 > part; ++part ){
    //----- Loop For difference Cases
  for( unsigned int part = 0 ; 2>part; ++part){
    for(unsigned int icase = 0; m_nCases>icase ; ++icase){
      
      if( m_doTiming){
        if( part ==0){m_timerTool->start( m_timeXProjeUp[icase]);}
        if( part ==1){m_timerTool->start( m_timeXProjeDo[icase]);} 
      }
      
      // Find The X Projection
      findXProjections(part,icase);
      if( m_doTiming){
        if( part ==0){
          m_timerTool->stop( m_timeXProjeUp[icase]);
          m_timerTool->start( m_timeCloneXUp[icase]);
        }
        if( part ==1){
          m_timerTool->stop( m_timeXProjeDo[icase]);
          m_timerTool->start(m_timeCloneXDo[icase]);  
        }
      }
      std::sort(m_xCandidates.begin(), m_xCandidates.end(), PrSeedTrack2::GreaterBySize());
      if(m_removeClonesX) removeClonesX( m_nCommonX[icase] , part, icase, m_xOnly);

      if( m_doTiming){
        if( part == 0){
          m_timerTool->stop( m_timeCloneXUp[icase]);
          m_timerTool->start(m_timeStereoUp[icase]);
        }
        if( part == 1){
          m_timerTool->stop( m_timeCloneXDo[icase]);
          m_timerTool->start(m_timeStereoDo[icase]);
        }
      }
      
      
      //Add The stereo Part
      if(!m_xOnly){
        addStereo( part, icase ); 
      }
      if(m_doTiming){
        if(part == 0){
          m_timerTool->stop( m_timeStereoUp[icase]);
          m_timerTool->start(  m_timeFlagUp[icase]);
        }
        if( part ==1){
          m_timerTool->stop( m_timeStereoDo[icase]);
          m_timerTool->start( m_timeFlagDo[icase]);
        }
      }
      
      //Flag found Hits at the end of each single case ( exclude the latest one )
      
      if(m_FlagHits && (icase ==0 || icase ==1) && !m_xOnly){
        flagHits(icase,part);
      }
      if( m_doTiming){
        if( part ==0)
          m_timerTool->stop( m_timeFlagUp[icase]);
        if( part ==1)
          m_timerTool->stop(m_timeFlagDo[icase]);
      }
    }
    if( m_xOnly) m_xCandidates.clear();
  }
  if(m_doTiming) m_timerTool->start( m_timeClone);

  // Global Clone Removal
  if(m_removeClones) removeClones(m_nCommonUV);
  
  if( m_doTiming){ 
    m_timerTool->stop( m_timeClone); 
    m_timerTool->start( m_timeConvert);
  }
  
  // Convert tracks to LHCb objects
  makeLHCbTracks( result );
  if( m_doTiming){
    m_timerTool->stop(m_timeConvert);
  }
  
  if( msgLevel(MSG::DEBUG)) debug()<<"Making LHCb Tracks Done"<<endmsg;
  if( m_doTiming){
    m_timerTool->stop( m_timeTotal);
  }
  return StatusCode::SUCCESS;
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode PrHybridSeeding::finalize() {

  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Finalize" << endmsg;
#ifdef TRUTH_MATCH_Histos
  return GaudiTupleAlg::finalize();
#else
  return GaudiAlgorithm::finalize();
#endif
}
//==================
//Add Stereo void
//==================
void PrHybridSeeding::addStereo(unsigned int part, unsigned int iCase)
{
  // use the named IDs of zones (should be useless considering the following consideration):
  // firstZone = part + s_T1U
  // lastZone  = part + s_T3V
  unsigned int firstZone = part+2; //1st station U Layers
  unsigned int lastZone = part+22;
  if(m_useFix){
    firstZone = 2;
    lastZone = 22;
  }
  for( PrSeedTrack2s::iterator itT =m_xCandidates.begin();m_xCandidates.end() != itT; ++itT){

    // (*itT).Case() can and must be written itT->Case()  (easier to read)
    if( (*itT).Case() != iCase) continue;
    if( (*itT).zone() != part) continue;
    if( !(*itT).valid()) continue;
    PrHits myStereo;
    myStereo.reserve(60);
#ifdef TRUTH_MATCH_Histos
    bool AssocX = false;
    LHCb::MCParticle* partic = nullptr;
    double efficiency = -1;
    int nAssocHits = -10;
    AssocX = AssocTrack( (*itT), efficiency, partic, nAssocHits);
#endif
    const unsigned int stepSize = m_useFix ? 1:2;

    // It is possible to consider to store UV zones exactly like in findXProjection
    // ex:

    // std::vector<PrHitZone*> uvZones;
    // uvZones.reserve(6);
    // for (unsigned int uvZoneId : {s_T1U, s_T1V, s_T2U, s_T2V, s_T3U, s_T3V}){
    //   uvZones.push_back(m_zones[uvZoneId|part]);
    // }
    //
    // for (PrHitZone* zone : uvZones) {



    // With the trick, it could be something like that:

    // std::vector<PrHitZone*> uvZones;
    // uvZones.reserve(12);
    // for (unsigned int detector_part : { 0, 1 }) {
    //   for (unsigned int uvZoneId : {s_T1U, s_T1V, s_T2U, s_T2V, s_T3U, s_T3V}) {
    //     uvZones.push_back(m_zones[uvZoneId|detector_part]);
    //   }
    // }
    //
    // for (PrHitZone* zone : uvZones) {


    for(unsigned int kk = firstZone; lastZone >kk;kk+=stepSize){
      if(m_zones[kk]->isX()) continue; // can be avoid with the previous trick
      double yMin = -1*std::fabs(m_yMax);
      double yMax = +1*std::fabs(m_yMax);
      double dxDy = m_zones[kk]->dxDy();
      double zPlane = m_zones[kk]->z();
      double xPred = (*itT).x(zPlane);
      if(yMin >yMax){
        double tmp = yMax;
        yMax = yMin;
        yMin = tmp;
      }
      if(yMax <yMin) always()<<"Stereo Problem"<<std::endl;
      if(m_doAsymmUV){

        // It is probably better to merge the two parts
        // the inner condition can written like this: (kk % 2 != part && m_useFix)
        if(part==0){ //Low
          //yMax = 1.0; //[-2700, 1.0]
          yMax = -m_yMin;
          if(kk%2==1 && m_useFix){ //[-30, 2.0]
            // yMax = 2.0;
            // yMin = -30.0; //mm
            yMax = -m_yMin_TrFix;
            yMin = -m_yMax_TrFix;
          }
        }
        if(part==1){ //Upper tracks
          //yMin = -1.0; // [-1.0 , 2500]
          yMin = m_yMin;
          if(kk%2==0 && m_useFix){//[-2.0,30.0]
            // yMin = -2.0;
            // yMax = 30.0;
            yMin = m_yMin_TrFix;
            yMax = m_yMax_TrFix;
          }
        }
      }
      double xMin = xPred + yMin*dxDy;
      double xMax = xPred + yMax*dxDy;
      if(xMin > xMax){
        double tmp = xMax;
        xMax = xMin;
        xMin = tmp;
      }
      PrHits& hits = m_zones[kk]->hits();
      PrHits::iterator itH = std::lower_bound(hits.begin(),hits.end(),xMin, lowerBoundX());
      // PrHits::iterator itEnd = std::upper_bound(hits.begin(),hits.end(),xMax, upperBoundX());
      //PrHits::iterator itH = m_zones[kk]->hits().begin();
      for(  ; hits.end()!= itH; ++itH ){

        // It is possible to avoid (*itH)-> syntax defining the following variable:
        // PrHit* hit = *itH
        //
        // Then (*itH)-> is replaced by hit->

        if( (*itH)->isUsed() && m_removeFlagged) continue;
        if( (*itH)->x() > xMax) break;
        double y = ((*itH)->x()-xPred)/dxDy;
        if(y >yMax) continue;
        if(y <yMin) continue;
#ifdef TRUTH_MATCH_Histos
        bool skip = false;
#endif
        if(m_useFix2ndOrder){
          if(kk%2==1){
            if( part==0 && (*itH)->yMax()<0.) continue;
            if( part==0 && (*itH)->yMax()>0.){ yMax = 1.0; yMin = -2.0 -std::fabs((*itH)->yMax()); }
            
            if( part==1 && (*itH)->yMax()>0.){ yMin = -1.0; yMax = std::fabs(m_yMax); }
            if( part==1 && (*itH)->yMax()<0.){ yMin = -2.0 + std::fabs((*itH)->yMax()); yMax = std::fabs(m_yMax);}
          }
          if(kk%2==0){ //part = 0 natural ; part = 1 opposite
            if( part==0 && (*itH)->yMin()<0.){ yMax = 1.0; yMin = -m_yMax;}
            if( part==0 && (*itH)->yMin()>0.){ yMax = 2.0 - std::fabs( (*itH)->yMin());   yMin = -std::fabs(m_yMax);}
            
            if( part==1 && (*itH)->yMin()>0.) continue;
            if( part==1 && (*itH)->yMin()<0.){ yMax = 2.0 + std::fabs((*itH)->yMin()); yMin = -1.0;}
          }
        }
        if(y<yMin) continue;
        if(y>yMax) continue;
        if(m_removeHole){
          // It is better to avoid the square root computation like this:
          // double radius2 = y*y + xPred8xPred;
          // if (radius2 < 87.0*87.0) continue;

          double radius = std::sqrt( y*y + xPred*xPred);
          if(  radius < 87.0) continue;
        }
        (*itH)->setCoord(  std::fabs(((*itH)->x() - xPred ) / dxDy  / zPlane  ));
#ifdef TRUTH_MATCH_Histos
        // if(partic!=nullptr){
        //   if( AssocX ){
        //     Tuple nTupleUV = nTuple("Seeding/CollectUV/AddStereo","AddStereoHit");
        //     nTupleUV->column("isWanted",isWanted(partic));
        //     nTupleUV->column("AssHit",matchKey( (*itH), partic->key());
        //     nTupleUV->column("P",partic->p());
        //     nTupleUV->column("nX",(*itT).hits().size());
        //     nTupleUV->column("Chi2DoF_X",(*itT).chi2PerDoF());
        //     nTupleUV->column("MaxChi2_X",(*itT).MaxChi2());
        //     nTupleUV->column("X0Back",(*itT).X0());
        //     nTupleUV->column("x_uorV",(*itH)->x());
        //     nTupleUV->column("xMin",xMin);
        //     nTupleUV->column("yMin",yMin);
        //     nTupleUV->column("radius",radius);
        //     nTupleUV->column("yMax",yMax);
        //     nTupleUV->column("xMax",xMax);
        //     nTupleUV->column("part",part);
        //     nTupleUV->column("dxDy",dxDy);
        //     nTupleUV->column("xPred",xPred);
        //     nTupleUV->column("Hit_yMin",(*itH)->yMin());
        //     nTupleUV->column("Hit_yMax",(*itH)->yMax());
        //     nTupleUV->column("y", ((*itH)->x()-xPred)/dxDy);
        //     nTupleUV->column("plane",((*itH)->planeCode()));
        //     nTupleUV->column("kk",kk);
        //     nTupleUV->column("Coord",((*itH)->x()-xPred)/dxDy/zPlane);
        //     nTupleUV->write(); 
        //   }
        // }
#endif
        if(m_removeFlagged && (*itH)->isUsed()) continue;
        myStereo.push_back( (*itH) );
      }
    }
    std::sort( myStereo.begin(), myStereo.end(), PrHit::LowerByCoord() );
    
    PrPlaneCounter2 plCount;
    //Save position of this x candidate, for later use
    //My Stereo is a collection of Hits in the stereo Layers
    //int minHough  = 5;
    unsigned int minUV = 4; //minimal different UV layer Hits
    unsigned int minTot = 10; //minimal different number of Hits
    if(iCase == 0){
      if((*itT).hits().size() == 6 ) { minUV = 4; minTot = 9;}
      if((*itT).hits().size() == 5 ) { minUV = 5; minTot = 9;}
      if((*itT).hits().size() == 4 ) { minUV = 6; minTot = 9;}
      // minHough = 6;
    }
    if(iCase == 1){
      if((*itT).hits().size() == 6 ) { minUV = 4; minTot = 9;}
      if((*itT).hits().size() == 5 ) { minUV = 5; minTot = 9;}
      if((*itT).hits().size() == 4 ) { minUV = 6; minTot = 9;}
      // minHough = 6;
    }
    if(iCase == 2){
      if((*itT).hits().size() == 6 ){ minUV = 4; minTot = 9;}
      if((*itT).hits().size() == 5 ){ minUV = 4; minTot = 9;}
      if((*itT).hits().size() == 4 ){ minUV = 5; minTot = 9;} 
    }
    // unsigned int minHough  = 5;
    unsigned int firstSpace = m_trackCandidates.size();
    PrHits::iterator itBeg = myStereo.begin();  //first hit in U-V layer with small Ty

    // itEnd seems to be one element earlier than expected
    PrHits::iterator itEnd = itBeg + minUV-1; //go ahead of minUV hits
    //Long tracks >5 GeV tolTyOffset = 0.002
    //additional hough cluster window based on the backward projection of the x-z plane?
    //Long tracks > 5 GeV:
    // double signSlope = (std::fabs((*itT).X0()) >200.)? +1. :  -1. ;  //true
    // double TyOffset = m_tolTyOffset;
    // double TySlope = m_tolTySlope;
    // //double minTy = 0.;
    // // Case 0 : All P>5GeV TyOffset = 0.002 is OK
    // // Case 1 : All P TyOffset + Slope 0.0035 ????
    // // Case 2 : All P TyOffset + Slope 0.0035 ????
    // //What do we have to do here?
    // //Case 0 : useImprovedStereo 0.002 + slope
    PrLineFitterY BestLine(m_geoTool->zReference(), (*itT));
    while( itEnd < myStereo.end()) {
      double tolTy = m_tolTyOffset[iCase] + m_tolTySlope[iCase]*std::fabs( (*itBeg)->coord() );
      if( ( (*(itEnd-1))->coord() - (*itBeg)->coord() )< tolTy){//there was -1
        while( itEnd+1 < myStereo.end() &&
               ((*itEnd)->coord() - (*itBeg)->coord() < tolTy) ){ // 
          ++itEnd; //extend last hit until you don't reach out of tolerance
        }
        plCount.set( itBeg, itEnd );
        if( (minUV-1) < plCount.nbDifferentUV() && plCount.isOKUV()){ // if inside the hough cluster process you have minUV u-v layers + at least one station with 2 hits you can start building up the track
          PrSeedTrack2 temp( *itT ) ;//maybe generate it later on?;
            if(m_useLineY){ //Preselection with line on Y with 1 hit per layer. Avoid due to occupancy to fit tracks with >1 hit per layer
              //bool fitDone = false;
              bool fit = false;
              if( plCount.nbSingleUV() == plCount.nbDifferentUV() && plCount.nbSingleUV() > minUV-1){ //if the case is nUVLayer = nHits
                fit = BestLine.fit( itBeg,itEnd ); //fit for the line will set the chi2 for BestLine
                //fitDone = true;
                if(fit && LineOK( m_Chi2LowLine[iCase] ,m_Chi2HighLine[iCase],  BestLine , temp) ){ //criteria satisfied => add Hits on track
                  // curly brackets
                  for( PrHits::iterator hit = itBeg; itEnd!= hit; ++hit)
                    temp.addHit( (*hit));
                }
              }
              if( plCount.nbDifferentUV() != plCount.nbSingleUV() && plCount.nbDifferentUV() > minUV-1){ //more than 1 hit per layer
                std::vector<std::vector<PrHit*> > PlanesMultiple; //Contains all the planes having > 1 layer
                PlanesMultiple.clear();

                // it is possible to consider a loop on layers instead of copying multiple times
                // the computation on every layer

                // hitsT1U or hits1U ?
                std::vector<PrHit*> hits1; //hits in layer 1
                hits1.clear(); // it is useless to clear an array just after its declaration
                std::vector<PrHit*> hits2;
                hits2.clear();
                std::vector<PrHit*> hits5;
                hits5.clear();
                std::vector<PrHit*> hits6;
                hits6.clear();
                std::vector<PrHit*> hits9;
                hits9.clear();
                std::vector<PrHit*> hits10;
                hits10.clear();
                std::vector<PrHit*> PlanesSingle;
                PlanesSingle.clear();
                PlanesMultiple.reserve(6);
                PlanesSingle.reserve(6);
                //Improve this by doing std::vector< vector<PrHit*>> hits where hits[0] is the vector of his in hits1.
                for(PrHits::iterator itH = itBeg; itEnd!=itH;++itH){
                  if(plCount.nbInPlane( (*itH)->planeCode()) ==1 ){//if 1 hit in plane of processed hit push back the hit in the container of PlaneSingle
                    PlanesSingle.push_back((*itH));
                    continue;
                  }
                  if((*itH)->planeCode()==1 && plCount.nbInPlane((*itH)->planeCode())>1){
                    hits1.push_back((*itH));
                  }
                  if((*itH)->planeCode()==2 && plCount.nbInPlane((*itH)->planeCode())>1){
                    hits2.push_back((*itH));
                  }
                  if((*itH)->planeCode()==5 && plCount.nbInPlane((*itH)->planeCode())>1){
                    hits5.push_back((*itH));
                  }
                  if((*itH)->planeCode()==6 && plCount.nbInPlane((*itH)->planeCode())>1){
                    hits6.push_back((*itH));
                  }
                  if((*itH)->planeCode()==9 && plCount.nbInPlane((*itH)->planeCode())>1){
                    hits9.push_back((*itH));
                  }
                  if((*itH)->planeCode()==10 && plCount.nbInPlane((*itH)->planeCode())>1){
                    hits10.push_back((*itH));
                  }
                }
                if(hits1.size()!=0) PlanesMultiple.push_back(hits1);
                if(hits2.size()!=0) PlanesMultiple.push_back(hits2);
                if(hits5.size()!=0) PlanesMultiple.push_back(hits5);
                if(hits6.size()!=0) PlanesMultiple.push_back(hits6);
                if(hits9.size()!=0) PlanesMultiple.push_back(hits9);
                if(hits10.size()!=0) PlanesMultiple.push_back(hits10);
                PrSeedTrack2s Tracks; //combinatorical tracks container when >1 hit in at least 1 layer.
                //PrHit* hit1 = nullptr;

                // Is it really necessary to have 6 nested loops? (maybe it is...)
                for( PrHit* hit1: PlanesMultiple[0]){
                  if(PlanesMultiple.size()==1){

                    // tem doesn't seem to be an explicit name
                    PrSeedTrack2 tem(*itT);
                    if(PlanesSingle.size()!=0){
                      tem.addHits2(PlanesSingle);
                    }
                    tem.addHit(hit1);
                    Tracks.push_back(tem);
                    continue;
                  }
                  if(PlanesMultiple.size()>1){
                    for( PrHit* hit2 : PlanesMultiple[1]){
                      if(PlanesMultiple.size() ==2){
                        PrSeedTrack2 tem(*itT);
                        if(PlanesSingle.size()!=0){
                          tem.addHits2(PlanesSingle);
                        }
                        tem.addHit(hit1);
                        tem.addHit(hit2);
                        Tracks.push_back(tem);
                        continue;
                      }
                      if(PlanesMultiple.size()>2){
                        for(PrHit* hit3 : PlanesMultiple[2]){
                          if(PlanesMultiple.size() == 3){
                            PrSeedTrack2 tem(*itT);
                            if(PlanesSingle.size()!=0){
                              tem.addHits2(PlanesSingle);
                            }
                            tem.addHit( hit1);
                            tem.addHit( hit2);
                            tem.addHit( hit3);
                            Tracks.push_back(tem);
                            continue;
                          }
                          if(PlanesMultiple.size()>3){
                            for(PrHit* hit4 : PlanesMultiple[3]){
                              if(PlanesMultiple.size() == 4){
                                PrSeedTrack2 tem(*itT);
                                if(PlanesSingle.size()!=0){
                                  tem.addHits2(PlanesSingle);
                                }
                              tem.addHit(hit1);
                              tem.addHit(hit2);
                              tem.addHit(hit3);
                              tem.addHit(hit4);
                              Tracks.push_back(tem); 
                              continue;
                              }
                              if(PlanesMultiple.size()>4){
                                for(PrHit* hit5 : PlanesMultiple[4]){
                                  if(PlanesMultiple.size() == 5){
                                    PrSeedTrack2 tem(*itT);
                                    if(PlanesSingle.size()!=0){
                                      tem.addHits2(PlanesSingle);
                                    }
                                    tem.addHit(hit1);
                                    tem.addHit(hit2);
                                    tem.addHit(hit3);
                                    tem.addHit(hit4);
                                    tem.addHit(hit5);
                                    Tracks.push_back(tem);
                                    continue; 
                                  }
                                  if(PlanesMultiple.size()>5){
                                    for(PrHit* hit6: PlanesMultiple[5])
                                    { 
                                      if(PlanesMultiple.size()==6)
                                      {
                                        PrSeedTrack2 tem(*itT);
                                        if(PlanesSingle.size()!=0){
                                          tem.addHits2(PlanesSingle);
                                        }
                                        tem.addHit(hit1);
                                        tem.addHit(hit2);
                                        tem.addHit(hit3);
                                        tem.addHit(hit4);
                                        tem.addHit(hit5);
                                        tem.addHit(hit6);
                                        Tracks.push_back(tem);
                                        continue; // useless continue
                                      }//end if planes size 6
                                    }//end loop hit in multiple 6
                                  }//end if plane size >5
                                }//end loop hit in multiple 5
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
                if(msgLevel(MSG::DEBUG)){
                  always()<<"UVSegments candidates"<<Tracks.size()<<endmsg;
                  always()<<"XZ segment start"<<endmsg; printTrack( (*itT));
                }
              // Tracks contains the Tracks with 1 hit per layer
                bool fit = false;
                for(PrSeedTrack2s::iterator tr = Tracks.begin(); Tracks.end()!=tr; ++tr){ //for each combinatorics (1 hit in 1 layer), choose the best one based on the LineChi2DoF + XZProjectionChi2DoF
                  fit = BestLine.fit( (*tr).hits().begin(), (*tr).hits().end() );
                  if(fit)
                    (*tr).setChi2LineY( BestLine.Chi2DoF() , BestLine.nHitsLine());
                  BestLine.reset();
                }
                // stable_sort?
                std::sort( Tracks.begin(), Tracks.end() , [temp]( const PrSeedTrack2& track1, const PrSeedTrack2& track2)->bool{
                    return( (track1.chi2DoFLine() + temp.chi2PerDoF()) < (track2.chi2DoFLine() + temp.chi2PerDoF()));
                  });
                //sort the vector based on the chi2 and get the first element (best chi2)
                if(( Tracks.front().chi2PerDoF() + Tracks.front().chi2DoFLine() < m_Chi2HighLine[iCase] 
                     && Tracks.front().hits().size()>10) 
                   ||( Tracks.front().chi2PerDoF() + Tracks.front().chi2DoFLine() < m_Chi2LowLine[iCase] 
                       && Tracks.front().hits().size()<11)){
                  for( PrHits::iterator it = Tracks.front().hits().begin(); Tracks.front().hits().end() != it; ++it){
                    if( (*it)->isX()) continue;
                    temp.addHit( (*it));
                  }
                }
                Tracks.clear();
                hits1.clear();
                hits2.clear();
                hits5.clear();
                hits6.clear();
                hits9.clear();
                hits10.clear();
              }
              BestLine.reset();
              if(msgLevel(MSG::DEBUG)){
                temp.sortbyz();
                // debug()?
                always()<<"temp track with the best line added"<<endmsg; printTrack(temp);
              }
              if(msgLevel(MSG::DEBUG)){
                always()<<"Will Fit the following track with N Hits = "<<temp.hits().size()<<endmsg;
                temp.sortbyz();
                
                printTrack(temp);
              }
              //for all hit added minCoord and so on and cut on that 
            }//end Use Liney
            if(!m_useLineY){
              for(PrHits::iterator itH = itBeg; itEnd!= itH;++itH){ 
                temp.addHit((*itH));
              }
            }
            bool ok = false;
            if(temp.hits().size()>12 && m_useLineY){
              always()<<"Error on number of hits  = "<<temp.hits().size()<<endmsg; printTrack(temp);
            }
            if(temp.hits().size()>=minTot){
              ok = fitSimultaneouslyXY(temp,iCase); 
            }
            if(temp.hits().size()>m_maxNHits){
              ok = false; //Max N hits is 13 here
            }
            while ( !ok && temp.hits().size() > minTot){
              //debug()?
              if( msgLevel(MSG::DEBUG) ) always()<<"RemoveWorst and Refit UV"<<endmsg;
              ok = removeWorstAndRefit( temp   , iCase);
              //step++
              if(temp.hits().size()>m_maxNHits) ok = false;
            }
            // #ifdef TRUTH_MATCH_Histos
            //           if( temp.hits().size() > 6){
            //             int nHitsAssociated = 0;
            //             LHCb::MCParticle *mcPart = nullptr;
            //             double effic = -1.;
            //             bool wanted = false;
            //             bool Assoc = AssocTrack(temp, effic, mcPart, nHitsAssociated);
            //             Tuple nTupleAfterFit = nTuple("Seeding/AddStereo/TupleafterFit","AfterFitBeforeChi2");
            //             if(Assoc){
            //               wanted = isWanted(mcPart);
            //             }
            //             PrPlaneCounter2 counterTuple;
            //             PrLineFitterY liner(m_geoTool->zReference(),temp);
            //             liner.set( temp.hits());
            //             counterTuple.set(temp.hits().begin(), temp.hits().end());
            //             nTupleAfterFit->column("Line_minCoord",liner.minCoord());
            //             nTupleAfterFit->column("Line_maxCoor",liner.maxCoord());
            //             nTupleAfterFit->column("Full_nXSingle", counterTuple.nbSingleX());
            //             nTupleAfterFit->column("Full_nUVSingle", counterTuple.nbSingleUV());
            //             nTupleAfterFit->column("Full_nUVDiff",counterTuple.nbDifferentUV());
            //             nTupleAfterFit->column("Full_isOKUV",counterTuple.isOKUV());
            //             nTupleAfterFit->column("Full_isOKX",counterTuple.isOKX());
            //             nTupleAfterFit->column("Full_isOK",counterTuple.isOK());
            //             nTupleAfterFit->column("Full_total",(int) temp.hits().size());
            //             if(ok){
            //               setChi2(temp);
            //             }
            //             nTupleAfterFit->column("Full_Chi2DoF", temp.chi2PerDoF());
            //             nTupleAfterFit->column("MaxChi2", temp.MaxChi2());
            //             nTupleAfterFit->column("Full_Chi2",temp.chi2());
            //             nTupleAfterFit->column("Full_ay",temp.ay());
            //             nTupleAfterFit->column("Full_by",temp.by());
            //             nTupleAfterFit->column("Full_X0",temp.X0());
            //             nTupleAfterFit->column("Full_Case",iCase);
            //             nTupleAfterFit->column("Full_nIter",step);
            //             nTupleAfterFit->column("Full_OK",ok);
            //             nTupleAfterFit->column("Full_Assoc", Assoc);
            //             if(Assoc){
            //               nTupleAfterFit->column("P",mcPart->p());
            //               nTupleAfterFit->column("isWanted",isWanted(mcPart));
            //             }
            //             if(!Assoc){
            //               nTupleAfterFit->column("P",-1000.);
            //               nTupleAfterFit->column("isWanted",isWanted(mcPart));
            //             }
            //             nTupleAfterFit->column("Eff",effic);
            //             nTupleAfterFit->column("nAss",nHitsAssociated);
            //             nTupleAfterFit->column("xSlope_9000",temp.xSlope(9000.));
            //             nTupleAfterFit->column("Ass_100",nHitsAssociated == temp.hits().size());
            //             nTupleAfterFit->column("Ass_OneOff",nHitsAssociated ==( temp.hits().size()-1)) ;
            //             nTupleAfterFit->write();
            //           }        
            //           setChi2(temp);
            // #endif
            if( ok  && temp.hits().size()>=minTot){
              setChi2(temp);
              double maxChi2 = m_maxChi2PerDoF[iCase];
              //number of hits hard coded??? More than in PatSeeding?
              if( temp.hits().size() >= minTot &&
                  temp.chi2PerDoF() < maxChi2 ){
                temp.setCase(iCase); 
                if(m_removeClones) std::sort(temp.hits().begin(), temp.hits().end(), compLHCbID());
                m_trackCandidates.push_back( temp );
              }
              // PROBABLE BUG:
              // the following line is equivalent to
              // minUV = 1;
              // itBeg++;
              itBeg += minUV=1 ;//-1 because 
              //itBeg += minUV-2;// was always 4
            }
        }//nUV check
      }//not in tolerance
      ++itBeg; // Move Forward the itBeg
      itEnd = itBeg + minUV-1 ; //was always 5
    }
    //=== Remove bad candidates: Keep the best one for this input track
    // FirstSpace is the number of track candidates before the finding of the UV Hits
    if( msgLevel(MSG::DEBUG)) always()<<"track Candidates Size"<<m_trackCandidates.size()<<"firstSpace"<<firstSpace<<endmsg;
    if( m_trackCandidates.size() > firstSpace+1 ){
      // it is better to avoid unsigned int loops
      for( unsigned int kk = firstSpace; m_trackCandidates.size() > kk -1 ; ++kk ){ 
        //maybe just use where you find more UV? it's automatic because you start from one X-Z candidate
        if(!m_trackCandidates[kk].valid()) continue;
        for( unsigned int ll = kk + 1; m_trackCandidates.size() > ll; ++ll ){
          if( !m_trackCandidates[ll].valid() )continue;
          if( m_trackCandidates[ll].hits().size() < m_trackCandidates[kk].hits().size()){
            m_trackCandidates[ll].setValid( false );
          }else if ( m_trackCandidates[ll].hits().size() > m_trackCandidates[kk].hits().size()){
            m_trackCandidates[kk].setValid( false );
            ////equal size, take the one with the better chi2
          }else if( m_trackCandidates[kk].chi2() < m_trackCandidates[ll].chi2() ){
            m_trackCandidates[ll].setValid( false );
          }else{
            m_trackCandidates[kk].setValid( false );
          }
        }
      }//end loop track Candidates for removal bad tracks
    }//loop candidates removal
  }//end loop xProjections
  m_xCandidates.clear(); //At the end of each case delete it or you
}
bool PrHybridSeeding::LineOK( double minChi2Low, double minChi2High, PrLineFitterY line, PrSeedTrack2& xProje){
  const int nHits = line.nHitsLine() + xProje.hits().size();
  const double Chi2DoFLineXProj = xProje.chi2PerDoF() + line.Chi2DoF();
  
  if( nHits > 10 && Chi2DoFLineXProj < minChi2High) return true;
  if( nHits < 11 && Chi2DoFLineXProj < minChi2Low) return true;
  return false;
}


void PrHybridSeeding::removeClonesX(unsigned int maxCommon, unsigned int part, unsigned int icase, bool xOnly)
{
  for( PrSeedTrack2s::iterator itT1 = m_xCandidates.begin(); m_xCandidates.end() !=itT1; ++itT1 ){
    if((*itT1).zone() != part) continue;
    if(!(*itT1).valid()) continue;
    for( PrSeedTrack2s::iterator itT2 = itT1 + 1; m_xCandidates.end() !=itT2; ++itT2 ){
      if( !m_removeClonesX && xOnly)break;
      if( (*itT2).zone() != part) continue;
      if (!(*itT2).valid()) continue;
      if( (*itT2).Case() != (*itT1).Case() ) continue;

      // looks weird, but does work
      int Compare = (*itT1).hits().size()*(*itT2).hits().size();
      switch(Compare){
      case 36: //6 vs 6
        maxCommon = 3;
        break;
      case 30: //6 vs 5
        maxCommon = 3;
        break;
      case 24: //6 vs 4
        maxCommon = 2;
        break; 
      case 25: //5 vs 5
        maxCommon= 2;
        break;
      case 20: //5 vs 4
        maxCommon = 1;
        break;
      case 16: /// 4 vs 4
        maxCommon = 1;
        break;
      }
      unsigned int nCommon = 0;
      PrHits::iterator itH1 = (*itT1).hits().begin();
      PrHits::iterator itH2 = (*itT2).hits().begin();
      PrHits::iterator itEnd1 = (*itT1).hits().end();
      PrHits::iterator itEnd2 = (*itT2).hits().end();
      //count number of common hits between track 1 and track 2
      while( itH1 != itEnd1 && itH2 != itEnd2 ){
        if ( (*itH1)->id() == (*itH2)->id() ){
          ++nCommon;

          // trailing spaces
          ++itH1;                                                                 
          ++itH2; 
        }
        else if( (*itH1)->id() < (*itH2)->id() ){ 
          ++itH1;
        }
        else{
          ++itH2;                                                          
        }
      }
      if ( nCommon >= maxCommon ){
        if ( (*itT1).hits().size() > (*itT2).hits().size() ) {
          (*itT2).setValid( false );
        } else if ( (*itT1).hits().size() < (*itT2).hits().size() ) {
          (*itT1).setValid( false );
        } else if ( (*itT1).chi2PerDoF() < (*itT2).chi2PerDoF() ){
          (*itT2).setValid( false );
        } else {
          (*itT1).setValid( false );
        }
      }
    }
    // curly brackets
    if(xOnly && (*itT1).valid() && (*itT1).zone() == part && icase==(m_nCases-1) )
      m_trackCandidates.push_back(*itT1);
  }
}

void PrHybridSeeding::removeClones(unsigned int maxCommon){
  std::sort(m_trackCandidates.begin(), m_trackCandidates.end(),PrSeedTrack2::GreaterBySize());
  for ( PrSeedTrack2s::iterator itT1 = m_trackCandidates.begin(); m_trackCandidates.end() !=itT1; ++itT1 ){
    if( !(*itT1).valid()) continue;
    for ( PrSeedTrack2s::iterator itT2 = itT1 + 1; m_trackCandidates.end() !=itT2; ++itT2 ) {
      if ( !(*itT2).valid()) continue;      
      bool oppositeZone = (*itT2).zone() != (*itT1).zone();
      if( oppositeZone && !m_ClonesUpDown) continue;
      unsigned int nCommon = 0;
      unsigned int nCommonUV = 0;
      PrHits::iterator itH1 = (*itT1).hits().begin();
      PrHits::iterator itH2 = (*itT2).hits().begin();
      PrHits::iterator itEnd1 = (*itT1).hits().end();
      PrHits::iterator itEnd2 = (*itT2).hits().end();
      while( itH1 != itEnd1 && itH2 != itEnd2 ){
        if ( (*itH1)->id() == (*itH2)->id() ){
          ++nCommon;
          ++itH1;
          ++itH2;                                                    
        }
        else if( (*itH1)->id() < (*itH2)->id() ){
          ++itH1;
        }
        else{                                           
          ++itH2;                                         
        }               
      }
      unsigned int maxCommonUV = maxCommon;
      if( ( nCommon>=maxCommonUV ) || (oppositeZone && nCommonUV>=m_nCommonUVTriangle)){
        if((*itT1).hits().size() > (*itT2).hits().size()){
          (*itT2).setValid( false);
        }
        if( (*itT1).hits().size() < (*itT2).hits().size())
        {
          (*itT1).setValid( false);
        }
        else if((*itT1).chi2PerDoF() < (*itT2).chi2PerDoF()){
          (*itT2).setValid(false);
        }else
        {
          (*itT1).setValid(false);
        }
      }
    }
  }
}


void PrHybridSeeding::flagHits(unsigned int icase, unsigned int part)
{
  std::sort( m_trackCandidates.begin() , m_trackCandidates.end() , PrSeedTrack2::LowerBySize()); //bigger size is in front
  for(PrSeedTrack2s::iterator track = m_trackCandidates.begin(); m_trackCandidates.end()!=track ; ++track){
    if( (*track).hits().size() < m_SizeFlag[icase]) break; // Important the sorting of before
    if( (*track).Case() != icase || (*track).zone()!=part) continue;
#ifdef TRUTH_MATCH_Histos
    //PrPlaneCounter2 plCount;
    //plCount.set(track);
    // if(!track.valid()) continue;
    // PrPlaneCounter2 plCount;
    // plCount.set(track.hits().begin(), track.hits().end());
    // Tuple tupleHitFlag = nTuple("Seeding/FlagHits/Hits","Flagging");
    // tupleHitFlag->column("nHits",track.hits().size());
    // tupleHitFlag->column("nx",plCount.nbSingleX());
    // tupleHitFlag->column("nUV",plCount.nbSingleUV());
    // tupleHitFlag->column("Chi2_DoF",track.chi2PerDoF());
    // tupleHitFlag->column("X0Back",track.X0());
    // tupleHitFlag->column("ay",track.ay());
    // tupleHitFlag->column("by",track.by());
    // tupleHitFlag->column("bx",track.bx());
    // tupleHitFlag->column("cx",track.cx());
    // tupleHitFlag->column("maxCHi2",track.MaxChi2());
    // bool Assoc = false;
    // LHCb::MCParticle* partic = nullptr;
    // double efficiency = -1;                         
    // int nAssocHits = 0;
    // Assoc = AssocTrack( track, efficiency, partic, nAssocHits);
    // tupleHitFlag->column("Assoc",Assoc);
    // tupleHitFlag->column("Case",(int)icase);  
    // if(partic != nullptr)
    // {
    //   tupleHitFlag->column("P",partic->p());
    // }
    // if( partic == nullptr)
    // {
    //   tupleHitFlag->column("P", -1000.);
    // }
    // bool isGhost = false;
    // isGhost = !isWanted(partic) || ( nAssocHits/track.hits().size() < 0.79);
    // bool PureSignal = isWanted(partic) && (-nAssocHits+track.hits().size())<1;
    // tupleHitFlag->column("isGhost", isGhost);
    // tupleHitFlag->column("Signal", PureSignal);
    // tupleHitFlag->column("isWanted",isWanted(partic));
    // tupleHitFlag->column("Hit_Eff", (double)nAssocHits/track.hits().size());
    // tupleHitFlag->write();
    //Look for tracks
#endif
    if(!(*track).valid())continue;
    if(! (((*track).hits().size()==11
           && (*track).chi2PerDoF()< m_MaxChi2Flag[icase]
           && std::fabs((*track).X0()) < m_MaxX0Flag[icase]) || 
          ((*track).hits().size()==12)) ) continue;
    for(PrHits::iterator it = (*track).hits().begin();(*track).hits().end()!=it; ++it){
      (*it)->setUsed(true);
    }
  }
}

//=========================================================================
//  Convert to LHCb tracks
//=========================================================================
void PrHybridSeeding::makeLHCbTracks ( LHCb::Tracks* result ) {
  for ( PrSeedTrack2s::iterator itT = m_trackCandidates.begin();m_trackCandidates.end() != itT; ++itT ) {


#ifdef TRUTH_MATCH_Histos
    Tuple tupleFinalTrack = nTuple("Seeding/UVFit/StoreTrack","UVEvents");
    PrPlaneCounter2 plcont;
    plcont.set( (*itT).hits().begin(),(*itT).hits().end());
    if( plcont.nbSingleX() <4 || plcont.nbSingleUV()<4)
    {
      //printTrack( *itT);
    }
    tupleFinalTrack->column("nX",(int)plcont.nbSingleX());
    tupleFinalTrack->column("nUV",(int)plcont.nbSingleUV());
    PrLineFitterY lineY(m_geoTool->zReference(), (*itT));
    lineY.set( (*itT).hits().begin(), (*itT).hits().end() );
    lineY.fit();
    tupleFinalTrack->column("chi2_line",lineY.Chi2());
    tupleFinalTrack->column("nHits", (int)(*itT).hits().size());
    tupleFinalTrack->column("valid", (*itT).valid());
    tupleFinalTrack->column("chi2", (*itT).chi2());
    tupleFinalTrack->column("chi2Dof", (*itT).chi2PerDoF());
    tupleFinalTrack->column("X0",(*itT).X0());
    tupleFinalTrack->column("ay",(*itT).ay());
    tupleFinalTrack->column("by",(*itT).by());
    tupleFinalTrack->column("bx",(*itT).bx());
    tupleFinalTrack->column("slope9000",(*itT).xSlope(9000.));
    LHCb::MCParticle * partic = nullptr;
    double efficiency = -1;
    int nAssocHits = -1;
    bool Assoc = AssocTrack((*itT),efficiency,partic,nAssocHits);
    tupleFinalTrack->column("eff",efficiency);    
    if(Assoc){
      tupleFinalTrack->column("P",partic->p());
    }
    if(!Assoc){
      tupleFinalTrack->column("P",-1000.);
    }
    tupleFinalTrack->column("Ass",Assoc);
    tupleFinalTrack->column("nAssocHits",nAssocHits);
    tupleFinalTrack->column("isWanted",isWanted(partic));
    tupleFinalTrack->column("Case",(*itT).Case());
    tupleFinalTrack->write();
#endif
    
    if( m_recoverTrack && ((*itT).nx()>5 || (*itT).ny()>5)){
      if((*itT).chi2PerDoF() < m_chi2Recover ) (*itT).setValid(true);
    }
    if ( !(*itT).valid() ) continue;
    if ( msgLevel(MSG::DEBUG) ) debug()<<"Creating LHCb Track"<<endmsg;
    LHCb::Track* tmp = new LHCb::Track;
    tmp->setType( LHCb::Track::Ttrack );
    tmp->setHistory( LHCb::Track::PrSeeding );
    double qOverP = m_geoTool->qOverP( *itT );
    LHCb::State tState;
    double z = StateParameters::ZEndT;
    tState.setLocation( LHCb::State::AtT );
    tState.setState( (*itT).x( z ), (*itT).y( z ), z, (*itT).xSlope( z ), (*itT).ySlope( ), qOverP );
    //== overestimated covariance matrix, as input to the Kalman fit
    tState.setCovariance( m_geoTool->covariance( qOverP ) );
    tmp->addToStates( tState );
    //== LHCb ids.
    tmp->setPatRecStatus( LHCb::Track::PatRecIDs );
    for ( PrHits::iterator itH = (*itT).hits().begin(); (*itT).hits().end() != itH; ++itH ) {
      tmp->addToLhcbIDs( (*itH)->id() );
    }
    tmp->setChi2PerDoF( (*itT).chi2PerDoF() );
    tmp->setNDoF(       (*itT).nDoF() );
    result->insert( tmp );
  }
}






void PrHybridSeeding::solveParabola(const PrHit* hit1, const PrHit* hit2, const PrHit* hit3, double& a, double& b, double& c){
  const double x3 = hit3->x();
  
  //New Parabola Method (be adapted to new parametrisation)
  const double z1_PB = hit1->z() - hit3->z();
  const double z2_PB = hit2->z() - hit3->z();
  const double x1_PB = hit1->x() - hit3->x();
  const double x2_PB = hit2->x() - hit3->x();

  const double det_PB = z1_PB*z2_PB*(z1_PB-z2_PB);

  if( std::fabs(det_PB) < 1e-8 ){
    a = 0.0;
    b = 0.0;
    c = 0.0;
    return;

  }

  a=0; b=0;c=0;
  a = (z2_PB*x1_PB-z1_PB*x2_PB)/det_PB;
  b = (-z2_PB*z2_PB*x1_PB+z1_PB*z1_PB*x2_PB)/det_PB;
  const double z3_PB = hit3->z() - m_geoTool->zReference();
  c = x3 + a * z3_PB * z3_PB - b *z3_PB;
  b -= 2 * a * z3_PB;
}

void PrHybridSeeding::solveParabola2(const PrHit* hit1,const PrHit* hit2,const PrHit* hit3,double& a1, double& b1,double& c1){
  const double z1 = hit1->z() - m_geoTool->zReference();
  const double z2 = hit2->z() - m_geoTool->zReference();
  const double z3 = hit3->z() - m_geoTool->zReference();
  const double x1 = hit1->x();
  const double x2 = hit2->x();
  const double x3 = hit3->x();
  //const double e = m_dRatio;
  const double corrZ1 = 1.+m_dRatio*z1;
  const double corrZ2 = 1.+m_dRatio*z2;
  const double corrZ3 = 1.+m_dRatio*z3;

  // As you never use corrZ* alone,
  // You can help the compiler by defining the following product:
  // double zzcorrZ* = z* * z* * corrZ*;
  const double det = (z1*z1)*corrZ1*z2 + z1*(z3*z3)*corrZ3 + (z2*z2)*corrZ2*z3 - z2*(z3*z3)*corrZ3 - z1*(z2*z2)*corrZ2 - z3*(z1*z1)*corrZ1;
  if( std::fabs(det) < 1e-8 )
  {
    a1 = 0.0;
    b1 = 0.0;
    c1 = 0.0;
    return;
  }
  const double det1 = (x1)*z2 + z1*(x3) + (x2)*z3 - z2*(x3) - z1*(x2) - z3*(x1);
  const double det2 = (z1*z1)*corrZ1*x2 + x1*(z3*z3)*corrZ3 + (z2*z2)*corrZ2*x3 - x2*(z3*z3)*corrZ3 - x1*(z2*z2)*corrZ2 - x3*(z1*z1)*corrZ1;
  const double det3 = (z1*z1)*corrZ1*z2*x3 + z1*(z3*z3)*corrZ3*x2 + (z2*z2)*corrZ2*z3*x1 - z2*(z3*z3)*corrZ3*x1 - z1*(z2*z2)*corrZ2*x3 - z3*(z1*z1)*corrZ1*x2;
  a1 = det1/det;
  b1 = det2/det;
  c1 = det3/det;
}
//=========================================================================
//  Fit the track, return OK if fit sucecssfull
//=========================================================================
bool PrHybridSeeding::fitSimultaneouslyXY( PrSeedTrack2& track , unsigned int iCase){
  double mat[15];
  double rhs[5];
  unsigned int nHitsX = 0;
  unsigned int nHitsStereo = 0;
  const double zRef = m_geoTool->zReference();
  for ( int loop = 0; 3 > loop ; ++loop ){
    if(loop ==1 && m_useCubic && (m_useCorrSlopes || m_useCorrPos)){
      // it could be interesting to keep in memory the squared radius for further computation
      double RadiusPosition = std::sqrt( (track.ax()*track.ax()*std::fabs(track.ax())/2000.) +
                                         (track.y(zRef)*track.y(zRef)*std::fabs(track.y(zRef))/1000.));
      double RadiusSlopes = std::sqrt( track.bx()*track.bx()*std::fabs(track.bx())/0.3 +
                                       track.by()*track.by()*std::fabs(track.by())/0.1);

      // Magic Numbers !
      double dRatioPos = -1.*(2.622e-4 +1.943e-8*RadiusPosition + 1.08e-11*RadiusPosition*RadiusPosition);
      double dRatioSlopes = -1.*( 2.6098e-4+ 6.31e-5*RadiusSlopes  -0.000156778*RadiusSlopes*RadiusSlopes + 0.000134126*RadiusSlopes*RadiusSlopes*RadiusSlopes);
      if(m_useCorrPos){
        track.setdRatio(dRatioPos);
      }
      if(m_useCorrSlopes){
        track.setdRatio(dRatioSlopes);
      }
    }
    std::fill(mat,mat+15,0.);
    std::fill(rhs,rhs+5,0.);
    // std::fill(mat, mat + 20, 0.);
    for( PrHits::iterator itH = track.hits().begin(); track.hits().end() != itH; ++itH ){
      if( (*itH)->isX()){
        nHitsX++;
      }else{
        nHitsStereo++;
      }  
      const double w = (*itH)->w();
      const double dxdy = (*itH)->dxDy();
      // const double dz = ((*itH)->z()-zRef)*0.001;
      const double yOnTrack = track.yOnTrack( (*itH) ) ;
      const double   dz = 0.001*((*itH)->z( yOnTrack ) - zRef);
      //it need the z at
      const double dRatio = track.dRatio();
      const double deta = dz*dz*(1. + dz*dRatio);
      const double wdz = w * dz;
      const double weta = w * deta;
      const double wdxdy = w * dxdy;
      const double wdxdydz = wdxdy * dz;
      const double dist = track.distance( *itH );
      //Fill Matrix
      mat[0] += w;
      mat[1] += wdz; 
      mat[2] += wdz * dz;
      mat[3] += weta;
      mat[4] += weta * dz; 
      mat[5] += weta * deta;
      mat[6] -= wdxdy;
      mat[7] -= wdxdydz;   
      mat[8] -= wdxdy * deta;
      mat[9] += wdxdy * dxdy;
      mat[10] -= wdxdydz;
      mat[11] -= wdxdydz * dz;
      mat[12] -= wdxdydz * deta;  
      mat[13] += wdxdydz * dxdy;
      mat[14] += wdxdydz * dz * dxdy;

      // fill right hand side
      rhs[0] += w * dist;
      rhs[1] += wdz * dist;
      rhs[2] += weta * dist;
      rhs[3] -= wdxdy * dist;
      rhs[4] -= wdxdydz * dist;
    }//Loop over Hits to fill the matrix
    // decompose matrix, protect against numerical trouble
    // track.setnXnY( nHitsX, nHitsStereo );
    if(nHitsX < 4 || nHitsStereo < 4) return false;
    ROOT::Math::CholeskyDecomp<double, 5> decomp(mat);
    if (!decomp) return false;
    decomp.Solve(rhs);
    rhs[1]*=1.e-3;
    rhs[2]*=1.e-6;
    rhs[4]*=1.e-3;    
    rhs[3]-=rhs[4]*zRef;
    if( loop >0 && (std::fabs(rhs[0]) > 1e4 || std::fabs(rhs[1]) > 5. ||
                    std::fabs(rhs[2]) > 1e-3 || std::fabs(rhs[3]) > 1e4 || std::fabs(rhs[4]) > 1.)) return false;
    track.updateParameters(rhs[0],rhs[1],rhs[2],rhs[3],rhs[4]);
  }
  track.setnXnY( nHitsX, nHitsStereo);
  //double chi2_track = 0.;
  double maxChi2 =0.;
  // double maxDistance = 0.;
  // double absdistanceSum = 0.;
  // double distanceSum = 0;
  // unsigned int NHits = (unsigned int)track.hits().size();
  //PrHit *worst = nullptr;
  for( PrHits::iterator itH = track.hits().begin(); track.hits().end() != itH; ++itH ){
    double chi2_onHit = track.chi2( *itH);
    //chi2_track+=chi2_onHit;
    if ( chi2_onHit > maxChi2 ){
      maxChi2 = chi2_onHit;
    }
  }//Set Max Chi2DoF
  double X0 = track.ax() - track.bx()*m_geoTool->zReference()+track.cx()*m_ConstC;
  track.setX0(X0);
  track.setMaxChi2(maxChi2);

  
#ifdef TRUTH_MATCH_Histos
  PrLineFitterY bestLine(m_geoTool->zReference(), track);   
  bestLine.set(track.hits());                                                        
  bool fit = false;                                                      
  // if(temp.hits().size()>minTot){
  //fit = bestLine.fit();                          
  // }
  Tuple nTupleLineY = nTuple("Seeding/AddStereo/TupleLineY","LineY");                                         
  LHCb::MCParticle* part = nullptr;
  double eff = -1;                                                                                     
  int nass =-10;                
  bool assoc = AssocTrack(track,eff,part,nass);
  bool isGhost = true;                                                                                 
  int nassLine = 0;                                                                             
  double minY = -1e5;                                                                               
  double maxY = 1e5;                                                                       
  double minCoord = 1e5;                                                                     
  double maxCoord = -1e5;                                                                                           
  int nUVLine = 0;
  for( PrHits::iterator hit = track.hits().begin(); track.hits().end()!=hit; ++hit){
    if( (*hit)->isX()) continue;
    nUVLine++;
    if(matchKey( (*hit) , part->key())) nassLine++;
    double y  =  ( (*hit)->x() - track.x( (*hit)->z())) / (*hit)->dxDy();
    if(y>minY){
      minY = y;
    }                                                                               
    if(y<maxY){                                                                   
      maxY = y;
    }                                                                                 
    if( (*hit)->coord() < minCoord){
      minCoord = (*hit)->coord();               
    }                                                                                 
    if( (*hit)->coord() >maxCoord){
      maxCoord = (*hit)->coord();
    }
  }
  Tuple nTupleFitMethodY = nTuple("Seeding/AddStereo/TupleFit","TupleFitMethod");
  nTupleFitMethodY->column("iCase",iCase);
  nTupleFitMethodY->column("maxChi2",maxChi2);
  nTupleFitMethodY->column("Chi2",chi2_track);
  nTupleFitMethodY->column("nHits",track.hits().size());
  nTupleFitMethodY->column("nUV",nUVLine);
  nTupleFitMethodY->column("nAssLine",nassLine);
  nTupleFitMethodY->column("nX",track.hits().size()-nUVLine);
  nTupleFitMethodY->column("EffUV",(double)nassLine/nUVLine);
  nTupleFitMethodY->column("X0",X0);
  nTupleFitMethodY->column("part",track.zone());
  nTupleFitMethodY->column("minY",minY);
  nTupleFitMethodY->column("maxY",maxY);
  nTupleFitMethodY->column("minCoord",minCoord);
  nTupleFitMethodY->column("maxCoord",maxCoord);
  nTupleFitMethodY->column("nUV",nUVLine);
  nTupleFitMethodY->column("nX",track.hits().size()-nUVLine);
  nTupleFitMethodY->column("nIter",refit);
  nTupleFitMethodY->column("Chi2DoF",chi2_track/(track.hits().size()-5.));
  nTupleFitMethodY->column("ay",track.ay());
  nTupleFitMethodY->column("y0",track.y(0.));
  nTupleFitMethodY->column("by",track.by());
  if(eff>0 && eff >0.7){
    nTupleFitMethodY->column("P",part->p());
  }else{
    nTupleFitMethodY->column("P",-1000.);
  }
  nTupleFitMethodY->column("isWanted",isWanted(part));
  nTupleFitMethodY->column("ax",track.ax());
  nTupleFitMethodY->column("ax",track.bx());
  nTupleFitMethodY->column("cx",track.cx());
  nTupleFitMethodY->column("xSlope8000",track.xSlope(8000.));
  nTupleFitMethodY->column("xSlope10000",track.xSlope(10000.));
  nTupleFitMethodY->column("yZref",track.y(m_geoTool->zReference()));
  nTupleFitMethodY->write();
#endif
  
  if( ( track.hits().size()>10) && maxChi2<m_maxChi2HitFullFitHigh[iCase]) return true;
  if(std::fabs(track.y(0.))< m_maxY0Low[iCase] &&
     maxChi2<m_maxChi2HitFullFitLow[iCase] && 
     track.hits().size()<11 &&
     std::fabs(track.y(m_geoTool->zReference())) < m_maxYZrefLow[iCase])
    return true;
  return false;
}
//=======================================
//Fit Only X Projection
//=======================================

bool PrHybridSeeding::fitXProjection(PrSeedTrack2& track, unsigned int iCase ){
  if (msgLevel(MSG::DEBUG)) debug()<<"Fitting"<<endmsg;
  if(track.hits().size()<m_minXPlanes) return false;  
  double mat[6];
  double rhs[3];
  //track.setdRatio(m_dRatio);
  for(int loop = 0;3>loop;++loop)
  {
    std::fill(mat,mat+6,0.);
    std::fill(rhs,rhs+3,0.);
    for( PrHits::iterator itH = track.hits().begin(); track.hits().end() != itH; ++itH ){
      const double dRatio = track.dRatio();
      const double w = (*itH)->w();//squared
      const double dz= 0.001*((*itH)->z() - m_geoTool->zReference());
      const double deta = dz*dz*(1. + dRatio*dz);
      const double dist = track.distance( *itH );
      if (msgLevel(MSG::DEBUG)) debug()<<"Loop \t"<<loop<<"\n Distance From Hit \t"<<dist<<endmsg;
      mat[0]+= w;
      mat[1]+= w * dz;   mat[2]+= w * dz * dz;
      mat[3]+= w * deta; mat[4]+= w * dz * deta;  mat[5]+= w * deta * deta;
      rhs[0]+= w * dist;
      rhs[1]+= w * dist * dz;
      rhs[2]+= w * dist * deta;
    }
    ROOT::Math::CholeskyDecomp<double,3> decomp(mat);
    if(!decomp){
      return false;
    }
    //Solve linear system
    decomp.Solve(rhs);
    rhs[1]*=1.e-3;
    rhs[2]*=1.e-6;
    if (msgLevel(MSG::DEBUG)) debug()<<"Loop \t"<<loop<<"\n a = \t"<<rhs[0]<<"\n b = \t"<<rhs[1]<<"\n c = \t"<<rhs[2]<<endmsg;
    // protect against unreasonable track parameter corrections
    // (check that out)
    if(std::abs(rhs[0]) > 1.e4 || std::abs(rhs[1]) > 5. ||
       std::abs(rhs[2]) > 1.e-3 ) return false;
    //Small corrections
    track.updateParameters(rhs[0],rhs[1],rhs[2],0.,0.);
    if(loop==0 && m_useCubic){
      track.setdRatio(m_dRatio);
    }
    //Put back later faster maybe
    if(loop >0 && std::abs(rhs[0]) < 5e-5 && std::abs(rhs[1]) < 5e-8 &&
       std::abs(rhs[2]) < 5e-11){
      break;
    }
  }
  //Compute some values on the track
  double chi2_track = 0.;
  double maxChi2 = 0.;
  for ( PrHits::iterator itH = track.hits().begin(); track.hits().end() != itH; ++itH )  //Loop over all hits in PrSeedTrack2
  {
    double chi2_onHit = track.chi2( *itH );
    chi2_track += track.chi2( *itH );
    if ( chi2_onHit > maxChi2 ){
      maxChi2 = chi2_onHit;
    }
  }
  double X0 = track.ax() - track.bx()*m_geoTool->zReference()+track.cx()*m_ConstC;
  track.setX0(X0);
  track.setMaxChi2(maxChi2);
  
#ifdef TRUTH_MATCH_Histos
  PrPlaneCounter2 counter;
  counter.set(track.hits().begin(),track.hits().end());
  if(track.hits().size()!=counter.nbSingleX()) printTrack(track);
  //if(maxChi2<20.){
  //DOTUPLE
  LHCb::MCParticle * partic = nullptr;
  double efficiency = -1;
  int nAssocHits = -1;
  bool Assoc = AssocTrack(track,efficiency,partic,nAssocHits);
  char tupleTitle[100];
  //  if(iCase==1 || iCase==2){
  Tuple tupleTrackXZFit = nTuple("Seeding/FitXZProjection/FitMethod","FitMethod");
  tupleTrackXZFit->column("XT1", track.x( m_zones[s_T1X1]->z()));
  tupleTrackXZFit->column("XT3", track.x( m_zones[s_T3X2]->z()));
  
  tupleTrackXZFit->column("Case",(int)iCase);
  //tupleTrackXZFit->column("nIter",(int)Refit);
  tupleTrackXZFit->column("Part", track.zone());
  tupleTrackXZFit->column("Assoc",Assoc);
  tupleTrackXZFit->column("MaxChi2Hit",maxChi2);
  tupleTrackXZFit->column("Chi2PerDoF",chi2_track/(track.hits().size()-3));
  tupleTrackXZFit->column("nbSingleX",(int)counter.nbSingleX());
  tupleTrackXZFit->column("nXHits",(int)track.hits().size());
  tupleTrackXZFit->column("ax",track.ax());
  tupleTrackXZFit->column("bx",track.bx());
  tupleTrackXZFit->column("cx",track.cx());
  tupleTrackXZFit->column("X0Back",X0);
  track.sortbyz();
  double slope = (track.hits().back()->x() - track.hits().front()->x())/(track.hits().back()->z()-track.hits().front()->z());
  double x0 = track.hits().front()->x()-track.hits().front()->z()*slope;
  tupleTrackXZFit->column("X0_Mid",x0);
  if(partic!=nullptr){
    tupleTrackXZFit->column("P",partic->p());
    tupleTrackXZFit->column("isWanted",isWanted(partic));
    tupleTrackXZFit->column("isXOk",counter.isOKX());
    tupleTrackXZFit->column("nAssocX",nAssocHits);
  }
  if(partic==nullptr){
    tupleTrackXZFit->column("P",-999.);
    tupleTrackXZFit->column("isWanted",false);
    tupleTrackXZFit->column("isXOk",false);
    tupleTrackXZFit->column("nAssocX",nAssocHits);
  }
  tupleTrackXZFit->write();
  
#endif
  //track.setRefitX(Refit);
  // double m_maxChi2Hit = 8.;
  // double m_maxX0Track = 8000.;
  // double m_maxChi2Offset = 0.75;
  // double m_maxX0Offset = 80.;
  // double m_Cut = 200.;
  // if( std::fabs(x0) > 300. && iCase >0) return false;


  // You can just:
  // return (maxChi2 < m_maxChi2HitsX[iCase]);

  if(maxChi2 < m_maxChi2HitsX[iCase]) return true;
  // if( track.hits().size() == 6  && maxChi2 < 8. ) return true;
  // if( maxChi2 < 1.5 && track.hits().size() < 6. ) return true;
  // if(  ( maxChi2-m_maxChi2Offset)*(std::fabs(X0)-m_maxX0Offset)<m_Cut
  //      &&(  maxChi2<m_maxChi2Hit || maxChi2<m_maxChi2Offset) ) return true;
  // if( std::fabs(X0) < 200)
  // //if( std::fabs(X0) > m_maxX0Track) return false;
  // if( iCase==0 &&  ( ( (maxChi2-m_maxChi2Offset )*(std::fabs(X0)-m_maxX0Offset) < m_Cut && maxChi2 < m_maxChi2Hit) || maxChi2<m_maxChi2Offset) ) return true;
  // if(maxChi2>m_maxChi2Hit) return false;
  // if( track.hits().size() == 4 && std::fabs(X0) < 400) return false;
  // if( track.hits().size() == 4 ){ m_maxChi2Hit = 4.;}
  // if( track.hits().size() == 5 ){ m_maxChi2Hit = 2.5;}
  // if( (maxChi2-m_maxChi2Offset )*(std::fabs(X0)-m_maxX0Offset) < m_Cut && maxChi2 < m_maxChi2Hit) return true;
  //if(track.hits().size()==4) return false;
  // if( track.hits().size()==6 && maxChi2<5.5) return true;
  // if( track.hits().size()==5 && maxChi2<5.5) return true;
  // if( track.hits().size()==4 && maxChi2<5.5) return true;
  return false;
}



bool PrHybridSeeding::removeWorstAndRefit(PrSeedTrack2& track, unsigned int iCase)
{
  bool removeX = false;
  bool removeUV = false;
  if( track.nx() == 4) removeUV = true;
  if( track.ny() == 4) removeX = true;
  if(msgLevel(MSG::DEBUG)) always()<<"Removing Worst UV layer"<<endmsg;
  PrHits::iterator worst= track.hits().begin();
  // PrHits::iterator worstUV = track.hits().begin();
  // PrHits::iterator worstDupl = track.hits().begin();
  double maxChi2 = -1.;
  
  for( PrHits::iterator itH = track.hits().begin(); track.hits().end() != itH; ++itH ){
    if( removeUV && (*itH)->isX() ) continue;
    if( removeX && !(*itH)->isX() ) continue;
    if(track.chi2(*itH) > maxChi2 ){
      maxChi2 = track.chi2( *itH );
      worst = itH;
    }
  }
  track.hits().erase(worst);
  return fitSimultaneouslyXY(track, iCase);
  return false;
}

//=========================================================================
//  Remove the worst hit and refit.
//=========================================================================
bool PrHybridSeeding::removeWorstAndRefitX ( PrSeedTrack2& track , unsigned int iCase)
{
  if(track.hits().size()<=m_minXPlanes) return false;
  if (msgLevel(MSG::DEBUG)) debug()<<"Removing Worst and Refitting"<<endmsg;
  //Find maxChi2 contribution of an Hit
  double maxChi2 = 0.;
  PrHits::iterator worst = track.hits().begin();    
  for ( PrHits::iterator itH = track.hits().begin(); track.hits().end() != itH; ++itH ){
    double chi2 = track.chi2( *itH );
    if( chi2 > maxChi2 ){
      maxChi2 = chi2;
      worst = itH;
    }
  }
  track.hits().erase(worst);
  bool OK = fitXProjection(track, iCase);
  return OK;
}

void PrHybridSeeding::setChi2 ( PrSeedTrack2& track )
{
  float chi2 = 0.;
  int   nDoF = -3;
  // Fitted a parabola
  bool hasStereo = false;
  //for ( PrHits::iterator itH = track.hits().begin(); track.hits().end() != itH; ++itH )
  for (PrHit* hit : track.hits()){
    const double d = track.distance( hit );
    if ( hit->dxDy() != 0 ) hasStereo = true;
    //hasStereo = hasStereo || (hit->dxDy() != 0);
    const double w = hit->w();
    chi2 += w * d * d;
    nDoF += 1;
  }
  if (hasStereo) {
    nDoF -= 2;
  }
  track.setChi2( chi2, nDoF );
}



//=========================================================================
//  Set the chi2 of the track
//=========================================================================
void PrHybridSeeding::setChi2X ( PrSeedTrack2& track ) {
  double chi2 = 0.;
  int   nDoF = -3;  // Fitted a parabola
  //bool hasStereo = false;
  //for ( PrHits::iterator itH = track.hits().begin(); track.hits().end() != itH; ++itH ) {
  for (PrHit* hit : track.hits()) {
    if(hit->dxDy() !=0 && msgLevel(MSG::DEBUG)) debug()<<"You were picking up Stereo Layers!!!"<<endmsg;
    const double d = track.distance( hit );
    const double w = hit->w();
    chi2 += w * d * d;
    nDoF += 1;
  }
  if (msgLevel(MSG::DEBUG)) debug()<<"Chi2 Set for track = \t"<<chi2<<endmsg;
  track.setChi2( chi2, nDoF );
}


//===============================================
//  Print the whole track
//=========================================================================
void PrHybridSeeding::printTrack ( PrSeedTrack2& track ) {
  for ( PrHits::iterator itH = track.hits().begin(); track.hits().end() != itH; ++itH ) {
    info() << format( "dist %7.3f dy %7.2f chi2 %7.2f ", track.distance( *itH ), track.deltaY( *itH ), track.chi2( *itH ) );
    printHit( *itH );
  }
}


//=========================================================================
//  Print the information of the selected hit
//=========================================================================
void PrHybridSeeding::printHit ( const PrHit* hit, std::string title ) {
  info() << "  " << title << " "
         << format( " Plane%3d zone%2d z0 %8.2f x0 %8.2f  size%2d charge%3d coord %8.3f used%2d ",
                    hit->planeCode(), hit->zone(), hit->z(), hit->x(),
                    hit->size(), hit->charge(), hit->coord(), hit->isUsed() );
  //if ( m_debugTool ) m_debugTool->printKey( info(), hit->id() );
  //if ( matchKey( hit ) ) info() << " ***";
  info() << endmsg;
}


void PrHybridSeeding::findXProjections(unsigned int part, unsigned int iCase)
{
  //m_xCandidates.clear();
  PrHits parabolaSeedHits;
  std::vector<PrHits> xHitsLists; //vector of list of Hits
  PrHits xHits;
  //just do the 1st one here //1st layer and last one
  unsigned int firstZoneId = -1;
  unsigned int lastZoneId = -1;
  if(0 == iCase){
    firstZoneId = s_T1X1 | part;
    lastZoneId  = s_T3X2 | part;
  }else if ( 1 == iCase ){
    firstZoneId = s_T1X2 | part;
    lastZoneId  = s_T3X1 | part;
  }else if ( 2 == iCase ){
    firstZoneId = s_T1X1 | part;
    lastZoneId  = s_T3X1 | part;
  }else if ( 3 == iCase ){
    firstZoneId = s_T1X2 | part;
    lastZoneId  = s_T3X2 | part;
  }
  if (msgLevel(MSG::DEBUG)) debug()<<"\t Loading Case Hit in first and last Zone"<<endmsg;
  PrHitZone* fZone = m_zones[firstZoneId];
  PrHitZone* lZone = m_zones[lastZoneId];
  if (msgLevel(MSG::DEBUG)) debug()<<"Hits in last and first Zone Loaded"<<endmsg;
  //load hits in first zone and last one
  PrHits& fHits = fZone->hits();
  PrHits& lHits = lZone->hits();
  double zFirst     = fZone->z(0.); //First  Zone Z value
  double zLast      = lZone->z(0.); //First  Zone Z value
  double DeltaZ = zLast-zFirst;     //Delta First-Last
  std::vector<PrHitZone*> xZones;
  xZones.reserve(4);
  for (unsigned int xZoneId : {s_T1X1, s_T1X2, s_T2X1, s_T2X2, s_T3X1, s_T3X2}){
    xZoneId |= part;
    if (xZoneId != firstZoneId && xZoneId != lastZoneId){
      //should i have to remove the layers already used?
      xZones.push_back(m_zones[xZoneId]);
    }
  }
  //PrHits::iterator itLBeg = lHits.begin();
  if (msgLevel(MSG::DEBUG)) debug()<<"Hits in the InBetween Zones Loaded"<<endmsg;
  if (msgLevel(MSG::DEBUG)) debug()<<"Will Loop over Hits in first Zone"<<endmsg;
  //for(PrHits::iterator itF = fHits.begin(); fHits.end() != itF; ++itF){
  for (PrHit* fHit : fHits){
    if ( 0 != iCase && fHit->isUsed() && m_removeFlagged ) continue;
    //define search window as a function of the x in the first layer
    double tx_inf =  fHit->x()/zFirst;
    double xProjeInf = tx_inf*zLast;
    double tolHp = m_TolFirstLast[iCase];
    //From Branch ImproveIt we change the alphaCOrrection : case 1 : 1988.63; case 2 : 2354.0
    //For all cases except case = 0
    double maxXl = xProjeInf + tx_inf*m_alphaCorrection[iCase]  +tolHp;
    double minXl = xProjeInf + tx_inf*m_alphaCorrection[iCase]  -tolHp;;
    if(maxXl < minXl){
      double temp = maxXl;
      minXl = maxXl;
      maxXl = temp;
    }
    if (msgLevel(MSG::DEBUG)) debug()<<"iCase "<<iCase<<"\t X last Min \t "<<minXl<<"\n\t\t\t Max X last \t "<<maxXl<<endmsg;
    if (msgLevel(MSG::DEBUG)) debug()<<"Will Loop over Last Layer"<<endmsg;
    
#ifdef TRUTH_MATCH_Histos
    char min[100];
    char max[100];
    sprintf(min,"L0Selection/HighMomentum2/Case%i/minXl_vs_tinf",iCase);
    sprintf(max,"L0Selection/HighMomentum2/Case%i/maxXl_vs_tinf",iCase);
    plot2D(minXl,tx_inf,min,"minXl vs t_{x}^{inf};minXl[mm];t_{x}^{inf}",-3300.,3300.,-0.5,0.5,200,200);
    plot2D(maxXl,tx_inf,max,"maxXl vs t_{x}^{inf};maxXl[mm];t_{x}^{inf}",-3300.,3300.,-0.5,0.5,200,200);
#endif
    
    
    
    PrHits::iterator itL = std::lower_bound(lHits.begin(), lHits.end(), minXl, lowerBoundX());
    PrHits::iterator itLEnd = std::upper_bound(lHits.begin(), lHits.end(), maxXl, upperBoundX());

    // should be a little faster with itL instead of lHits.begin() (itLEnd is obvisously after itL)
    // PrHits::iterator itLEnd = std::upper_bound(itL, lHits.end(), maxXl, upperBoundX());

    PrHit* lHit;
    for (; itLEnd != itL; ++itL) {
      lHit = *itL;
      if(nullptr == lHit){ // should never happen
        if (msgLevel(MSG::DEBUG)) debug()<<"Not Picking Hits in Last layer in the good search Window"<<endmsg;
        break;
      }
      if(lHit->x() > maxXl) break;
      if (lHit->isUsed() && m_removeFlagged) continue;
      double tx_pickedcombination = (lHit->x()-fHit->x())/DeltaZ;
      // double delta = lHit->x() - xProjeInf;
      parabolaSeedHits.clear();
      parabolaSeedHits.reserve(8);
      double x0 = fHit->x()-tx_pickedcombination*zFirst;
      double CorrX0 = m_x0Corr[iCase]*x0;
      double x0new = x0*(1.+m_x0Corr[iCase]);
      // double hp4alpha = m_x0Corr[iCase];
      if( msgLevel(MSG::DEBUG) ) debug() <<" x0 " << x0 << "CorrX0" << CorrX0 << "x0new" << x0new << "slope"<< m_x0Corr[iCase]<< endmsg;
      
#ifdef TRUTH_MATCH_Histos
      char min2[100];
      sprintf(min2,"L0Selection/HighMomentum2/Case%i/last_vs_tinf",iCase+10);
      plot2D(lHit->x() , tx_inf,min2,"last x vs t_{x}^{inf};x_{last}; t_{x}^{inf}",-3300.,3300.,-0.5,0.5,200,200);
      //Plot of 2 Hit associator for this case
      double xProjeCorrected = xProjeInf + tx_inf*m_alphaCorrection[iCase];
      bool assoc2Hit = false;
      LinkedTo<LHCb::MCParticle> ftLink(evtSvc(),msgSvc(),LHCb::FTClusterLocation::Default);
      LHCb::FTChannelID idFirst = fHit->id().ftID();
      LHCb::MCParticle* partic = ftLink.first(idFirst);
      while(0!=partic){
        assoc2Hit = matchKey(lHit,partic->key());
        if(assoc2Hit) break;
        partic = ftLink.next();
      }
      char truthSelected[100];
      char titleSelected[100];
      sprintf(titleSelected,"#Delta(X_{Hit}^{Last} - x_{inf}^{proj});#Delta(X_{Hit}^{Last} - x_{inf}^{proj}) [mm];Counts");
      char truthSelectedCorrected[100];
      char titleSelectedCorrected[100];
      sprintf(titleSelectedCorrected,"#Delta(X_{Hit}^{Last} - x_{inf}^{pro}}) Corr;#Delta(X_{Hit}^{Last} - x_{inf}^{proj}) [mm] Corr;Counts");
      
      char SelectedVsTxInf[100];
      char titleSelectedVsTxInf[100];
      sprintf(titleSelectedVsTxInf,"t_{x} infinte projection vs #Delta(x_{Last} - x_{inf}^{proje}) True;t_{x}^{inf};#Delta(x_{last}-x_{inf}) [mm] ");
      char SelectedVsTxInfCorrected[100];
      char titleSelectedVsTxInfCorrected[100];
      sprintf(titleSelectedVsTxInfCorrected,"t_{x} infinte projection vs #Delta(x_{Last} - x_{inf}^{proje}) True;t_{x}^{inf};#Delta(x_{last}-x_{inf}) [mm]");
      char x0Vstx[100];
      char titlex0Vstx[100];
      sprintf(titlex0Vstx,"x_{0} Vs t_{x}^{inf} ;x_{0};t_{x}^{inf}");
      if(assoc2Hit)
      {
        //All this plots are for 2 hit associated combination
        sprintf(truthSelected,"L0Selection/HighMomentum2/Case%i/DeltaHitProjection/AllAssociated",iCase+10);
        sprintf(truthSelectedCorrected,"L0Selection/HighMomentum2/Case%i/DeltaHitProjection/AllAssociated_withCorr",iCase+10);
        sprintf(SelectedVsTxInf,"L0Selection/HighMomentum2/Case%i/DeltaHitProjectionVsTxInf/AllAssociated",iCase+10);
        sprintf(SelectedVsTxInfCorrected,"L0Selection/HighMomentum2/Case%i/DeltaHitProjectionVsTxInf/AllAssociated_withCorr",iCase+10);
        sprintf(x0Vstx,"L0Selection/HighMomentum2/Case%i/X0VsTx/AllAssociated",iCase+10);
        plot2D(x0,tx_inf,x0Vstx,titlex0Vstx , -3000. , 3000., -0.6,0.6,300,200);
        plot(lHit->x()-xProjeInf,truthSelected,    titleSelected,-500.,500.,200);
        plot(lHit->x()-xProjeCorrected,truthSelectedCorrected,    titleSelectedCorrected,-500.,500.200);
        plot2D(tx_inf,lHit->x()-xProjeInf,SelectedVsTxInf,   titleSelectedVsTxInf,-0.6,0.6,-500.,500.,200,300);
        plot2D(tx_inf,lHit->x()-xProjeCorrected,SelectedVsTxInfCorrected,  titleSelectedVsTxInfCorrected,-0.6,0.6,-500.,500.,200,300);
        
        //This plots are for associated 2 hit combination from not physics particle interest
        if( !isWanted(partic)){
          sprintf(truthSelected,"L0Selection/HighMomentum2/Case%i/DeltaHitProjection/NotWanted",iCase+10);
            sprintf(truthSelectedCorrected,"L0Selection/HighMomentum2/Case%i/DeltaHitProjection/NotWanted_withCorr",iCase+10);
            sprintf(SelectedVsTxInf,"L0Selection/HighMomentum2/Case%i/DeltaHitProjectionVsTxInf/NotWanted",iCase+10);
            sprintf(SelectedVsTxInfCorrected,"L0Selection/HighMomentum2/Case%i/DeltaHitProjectionVsTxInf/NotWanted_withCorr",iCase+10);
            sprintf(x0Vstx,"L0Selection/HighMomentum2/Case%i/X0VsTx/NotWanted_withCorr",iCase+10);
            plot2D(x0,tx_inf,x0Vstx,titlex0Vstx , -3000. , 3000., -0.6,0.6,300,200);
            plot(lHit->x()-xProjeInf,truthSelected,    titleSelected,-500.,500.,200);
            plot(lHit->x()-xProjeCorrected,truthSelectedCorrected,    titleSelectedCorrected,-500.,500.200);
            plot2D(tx_inf,lHit->x()-xProjeInf,SelectedVsTxInf,   titleSelectedVsTxInf,-0.6,0.6,-500.,500.,200,300);
            plot2D(tx_inf,lHit->x()-xProjeCorrected,SelectedVsTxInfCorrected,  titleSelectedVsTxInfCorrected,-0.6,0.6,-500.,500.,200,300);
          }
          //This plots are for associated 2 hit combination from physics particle interest divided into momentum region
          if( isWanted(partic)){
            if(partic->p()>5000.){
              sprintf(truthSelected,"L0Selection/HighMomentum2/More5000/Case%i/DeltaHitProjection/Wanted_more5GeV",iCase+10);
              sprintf(truthSelectedCorrected,"L0Selection/HighMomentum2/More5000/Case%i/DeltaHitProjection/Wanted_more5GeV_withCorr",iCase+10);
              sprintf(SelectedVsTxInf,"L0Selection/HighMomentum2/More5000/Case%i/DeltaHitProjectionVsTxInf/Wanted_more5GeV",iCase+10);
              sprintf(SelectedVsTxInfCorrected,"L0Selection/HighMomentum2/More5000/Case%i/DeltaHitProjectionVsTxInf/Wanted_more5GeV_withCorr",iCase+10);
              plot(lHit->x()-xProjeInf,truthSelected,    titleSelected,-500.,500.,200);
              plot(lHit->x()-xProjeCorrected,truthSelectedCorrected,    titleSelectedCorrected,-500.,500.200);
              plot2D(tx_inf,lHit->x()-xProjeInf,SelectedVsTxInf,   titleSelectedVsTxInf,-0.6,0.6,-500.,500.,200,300);
              plot2D(tx_inf,lHit->x()-xProjeCorrected, SelectedVsTxInfCorrected, titleSelectedVsTxInfCorrected,-0.6,0.6,-500.,500.,200,300);
              sprintf(x0Vstx,"L0Selection/HighMomentum2/More5000/Case%i/X0VsTx/Wanted_more5GeV",iCase+10);
              plot2D(x0,tx_inf,x0Vstx,titlex0Vstx , -3000. , 3000., -0.6,0.6,300,200);
            }
            if(partic->p()>3000.){
              sprintf(truthSelected,"L0Selection/HighMomentum2/More3000/Case%i/DeltaHitProjection/Wanted_more3GeV",iCase+10);
              sprintf(truthSelectedCorrected,"L0Selection/HighMomentum2/More3000/Case%i/DeltaHitProjection/Wanted_more3GeV_withCorr",iCase+10);
              sprintf(SelectedVsTxInf,"L0Selection/HighMomentum2/More3000/Case%i/DeltaHitProjectionVsTxInf/Wanted_more3GeV",iCase+10);
              sprintf(SelectedVsTxInfCorrected,"L0Selection/HighMomentum2/More3000/Case%i/DeltaHitProjectionVsTxInf/Wanted_more3GeV_withCorr",iCase+10);
              plot(lHit->x()-xProjeInf,truthSelected,    titleSelected,-500.,500.,200);
              plot(lHit->x()-xProjeCorrected,truthSelectedCorrected,    titleSelectedCorrected,-500.,500.200);
              plot2D(tx_inf,lHit->x()-xProjeInf,SelectedVsTxInf,   titleSelectedVsTxInf,-0.6,0.6,-500.,500.,200,300);
              plot2D(tx_inf,lHit->x()-xProjeCorrected,SelectedVsTxInfCorrected,  titleSelectedVsTxInfCorrected,-0.6,0.6,-500.,500.,200,300);
              sprintf(x0Vstx,"L0Selection/HighMomentum2/More3000/Case%i/X0VsTx/Wanted_more3GeV",iCase+10);
              plot2D(x0,tx_inf,x0Vstx,titlex0Vstx , -3000. , 3000., -0.6,0.6,300,200);
            }
            if(partic->p()>500.){
              sprintf(truthSelected,"L0Selection/HighMomentum2/More500/Case%i/DeltaHitProjection/Wanted_more500MeV",iCase+10);
              sprintf(truthSelectedCorrected,"L0Selection/HighMomentum2/More500/Case%i/DeltaHitProjection/Wanted_withCorr_more500MeV",iCase+10);
              sprintf(SelectedVsTxInf,"L0Selection/HighMomentum2/More500/Case%i/DeltaHitProjectionVsTxInf/Wanted_more500MeV",iCase+10);
              sprintf(SelectedVsTxInfCorrected,"L0Selection/HighMomentum2/More500/Case%i/DeltaHitProjectionVsTxInf/Wanted_more500MeV_withCorr",iCase+10);
              plot(lHit->x()-xProjeInf,truthSelected,    titleSelected,-500.,500.,200);
              plot(lHit->x()-xProjeCorrected,truthSelectedCorrected,    titleSelectedCorrected,-500.,500.,200);
              plot2D(tx_inf,lHit->x()-xProjeInf,SelectedVsTxInf,   titleSelectedVsTxInf,-0.6,0.6,-500.,500.,200,300);
              plot2D(tx_inf,lHit->x()-xProjeCorrected,SelectedVsTxInfCorrected,  titleSelectedVsTxInfCorrected,-0.6,0.6,-500.,500.,200,300);
         
              sprintf(x0Vstx,"L0Selection/HighMomentum2/More500/Case%i/X0VsTx/Wanted_more500Mev",iCase+10);
              plot2D(x0,tx_inf,x0Vstx,titlex0Vstx , -3000. , 3000., -0.6,0.6,300,200);
            }
          }
        }
        //This pot are for not associated 2 hit combination
        if(!assoc2Hit){
          sprintf(truthSelected,"L0Selection/HighMomentum2/NotAss/Case%i/DeltaHitProjection/NotAssociated",iCase+10);
          sprintf(truthSelectedCorrected,"L0Selection/HighMomentum2/NotAss/Case%i/DeltaHitProjection/NotAssociated_withCorr",iCase+10);
          sprintf(SelectedVsTxInf,"L0Selection/HighMomentum2/NotAss/Case%i/DeltaHitProjectionVsTxInf/NotAssociated",iCase+10);
          sprintf(SelectedVsTxInfCorrected,"L0Selection/HighMomentum2/NotAss/Case%i/DeltaHitProjectionVsTxInf/NotAssociated_withCorr",iCase+10);
          plot(lHit->x()-xProjeInf,truthSelected,titleSelected,-500.,500.,200);
          plot(lHit->x()-xProjeCorrected,truthSelectedCorrected,    titleSelectedCorrected,-500.,500.,200);
          plot2D(tx_inf,lHit->x()-xProjeInf,SelectedVsTxInf,   titleSelectedVsTxInf,-0.6,0.6,-500.,500.,200,300);
          plot2D(tx_inf,lHit->x()-xProjeCorrected, SelectedVsTxInfCorrected, titleSelectedVsTxInfCorrected,-0.6,0.6,-500.,500.,200,300);
          sprintf(x0Vstx,"L0Selection/HighMomentum2/NotAss/Case%i/X0VsTx/NotAssociated",iCase+10);
          plot2D(x0,tx_inf,x0Vstx,titlex0Vstx , -3000. , 3000., -0.6,0.6,300,200);
        }
#endif
        if (msgLevel(MSG::DEBUG)) debug()<<"Will loop over Parabola Seed Hits: n Layers"<<m_zones.size()<<endmsg; 
        for(PrHitZone* xZone : {m_zones[s_T2X1 | part], m_zones[s_T2X2 | part]}) {
          double xProjected = x0 + xZone->z()*tx_pickedcombination;
          double xProjectedCorrected = xProjected+CorrX0;
          double xMax =0.;
          double xMin =0.;
          double max = 0.;
          double min = 0.;
          double slope = (m_tolAtX0Cut[iCase]-m_TolX0SameSign[iCase])/(m_x0Cut[iCase]-m_x0SlopeChange[iCase]);
          double slopeopp = (m_tolAtx0CutOppSig[iCase] -m_tolX0Oppsig[iCase])/(m_x0Cut[iCase]-m_x0SlopeChange2[iCase]);
          
          if(x0>0.){
            max = m_TolX0SameSign[iCase];
            //max= 1.2; //tight ?
            //min = x0>400.? -m_hp4_slope*x0 : -1.0; // <0
            min = x0 > m_x0SlopeChange[iCase]?  -slope*( x0 - m_x0SlopeChange[iCase]) - m_TolX0SameSign[iCase] : -m_TolX0SameSign[iCase];
            xMin = xProjectedCorrected + min;
            max = x0>m_x0SlopeChange2[iCase]? slopeopp*( x0 - m_x0SlopeChange2[iCase]) + m_tolX0Oppsig[iCase] : +m_tolX0Oppsig[iCase];
            xMax = xProjectedCorrected + max;
          }
          if(x0 < 0.){
            //max = x0<-400. ? -m_hp4_slope*: 1.0; // >0
            max = x0 <-m_x0SlopeChange[iCase]? -slope*( x0 + m_x0SlopeChange[iCase]) + m_TolX0SameSign[iCase]: m_TolX0SameSign[iCase];
            //min = -1.2;
            min = x0 < - m_x0SlopeChange2[iCase]? slopeopp*( x0 + m_x0SlopeChange2[iCase]) - m_tolX0Oppsig[iCase]: -m_tolX0Oppsig[iCase] ;
            xMin = xProjectedCorrected + min;
            xMax = xProjectedCorrected + max;
          }
          
          if(xMin > xMax) always()<<"Error xMin xMax"<<endmsg;
#ifdef TRUTH_MATCH_Histos
          // bool is7 =xZone->planeCode()==7;
          char titlemin[100];
          char titlemax[100];
          sprintf(titlemin,"L1Selection/HighMomentum2/Case%i/ParSeedLayer/Min",iCase);
          sprintf(titlemax,"L1Selection/HighMomentum2/Case%i/ParSeedLayer/Max",iCase);
          plot2D(x0,xMax,titlemax,"x0 vs xMax [mm]",-10000.,10000.,-3000.,3000.,400,300);
          plot2D(x0,xMin,titlemin,"x0 vs xMin [mm]",-10000.,10000.,-3000.,3000.,400,300);
#endif
          
          if( xMax<xMin && msgLevel(MSG::DEBUG)) debug()<<"\t\t\t\t\t Wrong xMax/xMin"<<endmsg;
          if( msgLevel(MSG::DEBUG)) debug()<<"Lower bound the zones"<<endmsg;
          PrHits::iterator itH = std::lower_bound(xZone->hits().begin(), xZone->hits().end(), xMin, lowerBoundX());
          if( itH == xZone->hits().end()) continue; //next Zone
          PrHit* mHit = nullptr;;
          if(msgLevel(MSG::DEBUG) )debug()<<"Will Loop over xZones Hits"<<endmsg;
          for(; xZone->hits().end() != itH; ++itH){
            mHit = *itH;
            if( mHit == nullptr) break; // should never happen
            if( mHit->x() > xMax ) break;
            // we can try to avoid this test
            if( mHit->isUsed() && m_removeFlagged) continue; //Not re use Hits in the middle
            if( msgLevel(MSG::DEBUG)) debug()<<"Filling Parabola Seed Hits"<<endmsg;
            parabolaSeedHits.push_back(mHit);
#ifdef TRUTH_MATCH_Histos
            bool HitsAssociated = false;
            char LinearVsX0[100];
            char DeltaLinear[100]; //1D
            char DeltaLinearCorrected[100]; //1D
            char DeltaLinearVsX0[100];
            char DeltaLinearCorrectedVsX0[100];
       
            char titleLinearVsX0[100];
            sprintf(titleLinearVsX0 ,"x_{parSeed} vs x0,x_{parSeed} [mm];x0[mm]"); 
            char titleDeltaLinear[100];
            sprintf(titleDeltaLinear,"#Delta(x_{parSeed} - x_{projected} [mm]);#Delta(x_{parSeed}-x_{projected}) [mm]; Counts");
            char titleDeltaLinearCorrected[100];
            sprintf(titleDeltaLinearCorrected ,"#Delta(x_{parSeed} - x_{proj}^{corr} [mm]);#Delta(x_{parSeed}-x_{projected})[mm];Counts");
            char titleDeltaLinearvsX0[100];
            sprintf(titleDeltaLinearvsX0 , "x0 vs #Delta;x0 [mm];#Delta(x_{parSeed}-x_{proj}) [mm]");
            char titleDeltaLinearvsX0Corrected[100];
            sprintf(titleDeltaLinearvsX0Corrected,"x0 vs #Delta corrected; x0 [mm];#Delta(x_{parSeed}-x_{proj}^{corr}) [mm]");
       
            //Do the association
            LinkedTo<LHCb::MCParticle> ft1Link(evtSvc(),msgSvc(),LHCb::FTClusterLocation::Default);
            LHCb::FTChannelID idFirst = fHit->id().ftID();
            //Get MCParticle associated to first layer hit
            LHCb::MCParticle* partic = ft1Link.first(idFirst);
            while(0!=partic){
              int key = partic->key();
              HitsAssociated = ( matchKey( mHit,key) && matchKey(lHit,key));
              if(HitsAssociated) break;
              partic = ft1Link.next();
            }
            if(HitsAssociated){
              sprintf(LinearVsX0,"L1Selection/HighMomentum2/Case%i/Seed_Vs_X0_Assoc",iCase);
              sprintf(DeltaLinear,"L1Selection/HighMomentum2/Case%i/DeltaSeed_Assoc",iCase);
              sprintf(DeltaLinearCorrected,"L1Selection/HighMomentum2/Case%i/X0VsDeltaSeed_Assoc",iCase);
              sprintf(DeltaLinearVsX0,"L1Selection/HighMomentum2/Case%i/Corr/DeltaSeedCorr_Assoc",iCase);
              sprintf(DeltaLinearCorrectedVsX0,"L1Selection/HighMomentum2/Case%i/Corr/X0VsDeltaSeedCorr_Assoc",iCase);
              plot2D(mHit->x(),x0,LinearVsX0,titleLinearVsX0,-3000.,3000.,-10000,10000,300,300);
              plot(mHit->x()-xProjected,DeltaLinear,titleDeltaLinear,-10.,10.,200);
              plot(mHit->x()-xProjectedCorrected,DeltaLinearCorrected,titleDeltaLinearCorrected,-10.,10.,200);
              plot2D(x0,mHit->x()-xProjected,DeltaLinearVsX0,titleDeltaLinearvsX0,-10000.,10000.,-10.,10,300,200);
              plot2D(x0,mHit->x()-xProjectedCorrected,DeltaLinearCorrectedVsX0,titleDeltaLinearvsX0Corrected,-10000.,10000.,-10.,10,300,200);
              if(!isWanted(partic)){
                sprintf(LinearVsX0,"L1Selection/HighMomentum2/Case%i/Seed_Vs_X0_Assoc_notWanted",iCase);
                sprintf(DeltaLinear,"L1Selection/HighMomentum2/Case%i/DeltaSeed_Assoc_notWanted",iCase);
                sprintf(DeltaLinearCorrected,"L1Selection/HighMomentum2/Case%i/X0VsDeltaSeed_Assoc_notWanted",iCase);
                sprintf(DeltaLinearVsX0,"L1Selection/HighMomentum2/Case%i/Corr/DeltaSeedCorr_Assoc_notWanted",iCase);
                sprintf(DeltaLinearCorrectedVsX0,"L1Selection/HighMomentum2/Case%i/Corr/X0VsDeltaSeedCorr_Assoc_notWanted",iCase);
                plot2D(mHit->x(),x0,LinearVsX0,titleLinearVsX0,-3000.,3000.,-10000,10000,300,300);
                plot(mHit->x()-xProjected,DeltaLinear,titleDeltaLinear,-10.,10.,200);
                plot(mHit->x()-xProjectedCorrected,DeltaLinearCorrected,titleDeltaLinearCorrected,-10.,10.,200);
                plot2D(x0,mHit->x()-xProjected,DeltaLinearVsX0,titleDeltaLinearvsX0,-10000.,10000.,-10.,10,300,200);
                plot2D(x0,mHit->x()-xProjectedCorrected,DeltaLinearCorrectedVsX0,titleDeltaLinearvsX0Corrected,-10000.,10000.,-10.,10,300,200);
              }
              if(isWanted(partic)){
                sprintf(LinearVsX0,"L1Selection/HighMomentum2/Case%i/Seed_Vs_X0_Assoc_Wanted",iCase);
                sprintf(DeltaLinear,"L1Selection/HighMomentum2/Case%i/DeltaSeed_Assoc_Wanted",iCase);
                sprintf(DeltaLinearCorrected,"L1Selection/HighMomentum2/Case%i/X0VsDeltaSeed_Assoc_Wanted",iCase);
                sprintf(DeltaLinearVsX0,"L1Selection/HighMomentum2/Case%i/Corr/DeltaSeedCorr_Assoc_Wanted",iCase);
                sprintf(DeltaLinearCorrectedVsX0,"L1Selection/HighMomentum2/Case%i/Corr/X0VsDeltaSeedCorr_Assoc_Wanted",iCase);
                plot2D(mHit->x(),x0,LinearVsX0,titleLinearVsX0,-3000.,3000.,-10000,10000,300,300);
                plot(mHit->x()-xProjected,DeltaLinear,titleDeltaLinear,-10.,10.,200);
                plot(mHit->x()-xProjectedCorrected,DeltaLinearCorrected,titleDeltaLinearCorrected,-10.,10.,200);
                plot2D(x0,mHit->x()-xProjected,DeltaLinearVsX0,titleDeltaLinearvsX0,-10000.,10000.,-10.,10,300,200);
                plot2D(x0,mHit->x()-xProjectedCorrected,DeltaLinearCorrectedVsX0,titleDeltaLinearvsX0Corrected,-10000.,10000.,-10.,10,300,200);
                if(partic->p()>5000.){
                  sprintf(LinearVsX0,"L1Selection/HighMomentum2/Case%i/Seed_Vs_X0_Assoc_Wanted_more5",iCase);
                  sprintf(DeltaLinear,"L1Selection/HighMomentum2/Case%i/DeltaSeed_Assoc_Wanted_more5",iCase);
                  sprintf(DeltaLinearCorrected,"L1Selection/HighMomentum2/Case%i/X0VsDeltaSeed_Assoc_Wanted_more5",iCase);
                  sprintf(DeltaLinearVsX0,"L1Selection/HighMomentum2/Case%i/Corr/DeltaSeedCorr_Assoc_Wanted_more5",iCase);
                  sprintf(DeltaLinearCorrectedVsX0,"L1Selection/HighMomentum2/Case%i/Corr/X0VsDeltaSeedCorr_Assoc_Wanted_more5",iCase);
                  plot2D(mHit->x(),x0,LinearVsX0,titleLinearVsX0,-3000.,3000.,-10000,10000,300,300);
                  plot(mHit->x()-xProjected,DeltaLinear,titleDeltaLinear,-10.,10.,200);
                  plot(mHit->x()-xProjectedCorrected,DeltaLinearCorrected,titleDeltaLinearCorrected,-10.,10.,200);
                  plot2D(x0,mHit->x()-xProjected,DeltaLinearVsX0,titleDeltaLinearvsX0,-10000.,10000.,-50.,50,300,200);
                  plot2D(x0,mHit->x()-xProjectedCorrected,DeltaLinearCorrectedVsX0,titleDeltaLinearvsX0Corrected,-10000.,10000.,-10.,10,300,200);
                }
                if(partic->p()>3000.){
                  sprintf(LinearVsX0,"L1Selection/HighMomentum2/Case%i/Seed_Vs_X0_Assoc_Wanted_more3",iCase);
                  sprintf(DeltaLinear,"L1Selection/HighMomentum2/Case%i/DeltaSeed_Assoc_Wanted_more3",iCase);
                  sprintf(DeltaLinearCorrected,"L1Selection/HighMomentum2/Case%i/X0VsDeltaSeed_Assoc_Wanted_more3",iCase);
                  sprintf(DeltaLinearVsX0,"L1Selection/HighMomentum2/Case%i/Corr/DeltaSeedCorr_Assoc_Wanted_more3",iCase);
                  sprintf(DeltaLinearCorrectedVsX0,"L1Selection/HighMomentum2/Case%i/Corr/X0VsDeltaSeedCorr_Assoc_Wanted_more3",iCase);
                  plot2D(mHit->x(),x0,LinearVsX0,titleLinearVsX0,-3000.,3000.,-10000,10000,300,300);
             
                  plot(mHit->x()-xProjected,DeltaLinear,titleDeltaLinear,-20.,20.,200);
             
                  plot(mHit->x()-xProjectedCorrected,DeltaLinearCorrected,titleDeltaLinearCorrected,-10.,10.,200);
             
                  plot2D(x0,mHit->x()-xProjected,DeltaLinearVsX0,titleDeltaLinearvsX0,-10000.,10000.,-50.,50,300,200);
             
                  plot2D(x0,mHit->x()-xProjectedCorrected,DeltaLinearCorrectedVsX0,titleDeltaLinearvsX0Corrected,-10000.,10000.,-10.,10.,300,200);
             
                }
                if(partic->p()>500.){
                  sprintf(LinearVsX0,"L1Selection/HighMomentum2/Case%i/Seed_Vs_X0_Assoc_Wanted_more500",iCase);
                  sprintf(DeltaLinear,"L1Selection/HighMomentum2/Case%i/DeltaSeed_Assoc_Wanted_more500",iCase);
                  sprintf(DeltaLinearCorrected,"L1Selection/HighMomentum2/Case%i/X0VsDeltaSeed_Assoc_Wanted_more500",iCase);
                  sprintf(DeltaLinearVsX0,"L1Selection/HighMomentum2/Case%i/Corr/DeltaSeedCorr_Assoc_Wanted_more500",iCase);
                  sprintf(DeltaLinearCorrectedVsX0,"L1Selection/HighMomentum2/Case%i/Corr/0VsDeltaSeedCorr_Assoc_Wanted_more500",iCase);
                  plot2D(mHit->x(),x0,LinearVsX0,titleLinearVsX0,-3000.,3000.,-10000,10000,300,300);
             
                  plot(mHit->x()-xProjected,DeltaLinear,titleDeltaLinear,-10.,10.,200);
             
                  plot(mHit->x()-xProjectedCorrected,DeltaLinearCorrected,titleDeltaLinearCorrected,-10.,10.,200);
             
                  plot2D(x0,mHit->x()-xProjected,DeltaLinearVsX0,titleDeltaLinearvsX0,-10000.,10000.,-50.,50,300,200);
             
                  plot2D(x0,mHit->x()-xProjectedCorrected,DeltaLinearCorrectedVsX0,titleDeltaLinearvsX0Corrected,-10000.,10000.,-10.,10,300,200);
             
                }
              }
            }
            if(!HitsAssociated){
              sprintf(LinearVsX0,"L1Selection/HighMomentum2/Case%i/Seed_Vs_X0_NotAssoc",iCase);
              sprintf(DeltaLinear,"L1Selection/HighMomentum2/Case%i/DeltaSeed_NotAssoc",iCase);
              sprintf(DeltaLinearCorrected,"L1Selection/HighMomentum2/Case%i/X0VsDeltaSeed_NotAssoc",iCase);
              sprintf(DeltaLinearVsX0,"L1Selection/HighMomentum2/Case%i/Corr/DeltaSeedCorr_NotAssoc",iCase);
              sprintf(DeltaLinearCorrectedVsX0,"L1Selection/HighMomentum2/Case%i/Corr/X0VsDeltaSeedCorr_NotAssoc",iCase);
              plot2D(mHit->x(),x0,LinearVsX0,titleLinearVsX0,-3000.,3000.,-10000,10000,300,300);
              plot(mHit->x()-xProjected,DeltaLinear,titleDeltaLinear,-10.,10.,200);
              plot(mHit->x()-xProjectedCorrected,DeltaLinearCorrected,titleDeltaLinearCorrected,-10.,10.,200);
              plot2D(x0,mHit->x()-xProjected,DeltaLinearVsX0,titleDeltaLinearvsX0,-10000.,10000.,-10.,10,300,200);
              plot2D(x0,mHit->x()-xProjectedCorrected,DeltaLinearCorrectedVsX0,titleDeltaLinearvsX0Corrected,-10000.,10000.,-10.,10,300,200);
            }
#endif
          }
          // if (parabolaSeedHits.size() > 0 && msgLevel(MSG::DEBUG)) debug()<<"ParabolaSeedHits Size \t ="<<parabolaSeedHits.size()<<endmsg;
          //Look for another Hit in last layer
          //end loop to pick up Hits in the 2 inner Layers (was only)
        }
        if(parabolaSeedHits.size()==0) continue; //go next last layer hit
        //if we don't fine any parabola Seed Hits in the middle 2 Layers then search for another XLast Hit
        // sort the parabola seed hits wrt to distance to the linear projection
        // merged parabolaSeedHits T2-1 & T2-2
        //=======================================================
        // We have 1 Hit in 1st 1 Hit in last and a
        // vector of Hits for in-between
        //=======================================================
        //std::vector<PrHits> xHitsLists; //vector of list of Hits
        xHitsLists.clear();
   
        //=======================================================
        //Sort the ParabolaSeedHits for in-between layers in increasing distance from the Projected Corrected position only when we have more than 1 ParabolaSeedHit
        //=======================================================
        
        if(parabolaSeedHits.size()>1){
          //Principle of the Lambda funtion, Hits sorted wrt distance from linear Projection 1st-3rd layer
          std::sort( parabolaSeedHits.begin(),parabolaSeedHits.end(),
                            [x0new,tx_pickedcombination](const PrHit* lhs, const PrHit* rhs)
                            ->bool{
                              // double lhsx0 =0 ;
                              // double rhsx0 =0 ;
                              // lhsx0 = x0new;
                              // rhsx0 = x0new;
                              
                              // Wrapping the condition in a functor with an explicit name could be a good idea (just for clarity)
                              // At least an explanation of what this means
                              return std::fabs( lhs->x() - ( x0new + lhs->z()*tx_pickedcombination)) < std::fabs( rhs->x() - (x0new + rhs->z()*tx_pickedcombination));} );
        }
        if (msgLevel(MSG::DEBUG)) debug()<<"The Lambda Function Sorting end"<<endmsg;
        unsigned int maxParabolaSeedHits = m_maxParabolaSeedHits;
        if(parabolaSeedHits.size()<m_maxParabolaSeedHits){
          maxParabolaSeedHits = parabolaSeedHits.size();
        }
        // if( parabolaSeedHits.size()>m_maxParabolaSeedHits){
        //   maxParabolaSeedHits = parabolaSeedHits.size();
        // }

        // int loop (unsigned int loop must be avoided)
        for (unsigned int i = 0; i<maxParabolaSeedHits;++i) //build a parabola for each 3 hit combination
        {
          //if (maxParabolaSeedHits==0) break; //maybe a repetition
          double a = 0;
          double b = 0;
          double c = 0;
          //PrHits xHits;
          xHits.clear();
          if (m_useCubic){
            solveParabola2(fHit,parabolaSeedHits[i],lHit,a,b,c); //Extrapolation with dRatio
          }else{
            solveParabola(fHit,parabolaSeedHits[i],lHit,a,b,c); //Extrapolation without dRatio
          }
#ifdef TRUTH_MATCH_Histos
          char a_par[100];
          char b_par[100];
          char c_par[100];
          char title_apar[100];
          sprintf(title_apar,"Solve Parabola a;a;Counts");
          char title_bpar[100];
          sprintf(title_bpar,"Solve Parabola b;b;Counts");
          char title_cpar[100];
          sprintf(title_cpar,"Solve Parabola c;c;Counts");
     
          //do the 3 Hit association again
          bool HitsAssociated = false;
          LinkedTo<LHCb::MCParticle> ft1Link(evtSvc(),msgSvc(),LHCb::FTClusterLocation::Default);
          LHCb::FTChannelID idFirst = fHit->id().ftID();
          //Get MCParticle associated to first layer hit
          LHCb::MCParticle* partic = ft1Link.first(idFirst);
          while(0!=partic){
            int key = partic->key();
            HitsAssociated = ( matchKey( parabolaSeedHits[i],key) && matchKey(lHit,key));
       
            if(HitsAssociated) break;
            partic = ft1Link.next();
          }
          if(HitsAssociated){
            sprintf(a_par,"L2Selection/HighMomentum2/Case%i/Apar_assoc",iCase+10);
            plot(a,a_par,title_apar,-10e-6,10e-6,300);
            sprintf(b_par,"L2Selection/HighMomentum2/Case%i/Bpar_assoc",iCase+10);
            plot(b,b_par,title_bpar,-1.5,1.5,300);
            sprintf(c_par,"L2Selection/HighMomentum2/Case%i/Cpar_assoc",iCase+10);
            plot(c,c_par,title_cpar,-3000,3000,300);
            if(isWanted(partic)){
              sprintf(a_par,"L2Selection/HighMomentum2/Case%i/Apar_assoc_wanted",iCase+10);
              plot(a,a_par,title_apar,-10e-6,10e-6,300);
              sprintf(b_par,"L2Selection/HighMomentum2/Case%i/Bpar_assoc_wanted",iCase+10);
              plot(b,b_par,title_bpar,-1.5,1.5,300);
              sprintf(c_par,"L2Selection/HighMomentum2/Case%i/Cpar_assoc_wanted",iCase+10);
              plot(c,c_par,title_cpar,-3000,3000,300);
              if(partic->p()>5000){
                sprintf(a_par,"L2Selection/HighMomentum2/Case%i/Apar_assoc_wanted_more5",iCase+10);
                plot(a,a_par,title_apar,-10e-6,10e-6,300);
                sprintf(b_par,"L2Selection/HighMomentum2/Case%i/Bpar_assoc_wanted_more5",iCase+10);
                plot(b,b_par,title_bpar,-1.5,1.5,300);
                sprintf(c_par,"L2Selection/HighMomentum2/Case%i/Cpar_assoc_wanted_more5",iCase+10);
                plot(c,c_par,title_cpar,-3000,3000,300);
              }
              if(partic->p()>3000.){
                sprintf(a_par,"L2Selection/HighMomentum2/Case%i/Apar_assoc_wanted_more3",iCase+10);
                plot(a,a_par,title_apar,-10e-6,10e-6,300);
                sprintf(b_par,"L2Selection/HighMomentum2/Case%i/Bpar_assoc_wanted_more3",iCase+10);
                plot(b,b_par,title_bpar,-1.5,1.5,300);
                sprintf(c_par,"L2Selection/HighMomentum2/Case%i/Cpar_assoc_wanted_more3",iCase+10);
                plot(c,c_par,title_cpar,-3000,3000,300);
              }
              if(partic->p()>500.){
                sprintf(a_par,"L2Selection/HighMomentum2/Case%i/Apar_assoc_wanted_more500",iCase+10);
                plot(a,a_par,title_apar,-10e-6,10e-6,300);
                sprintf(b_par,"L2Selection/HighMomentum2/Case%i/Bpar_assoc_wanted_more500",iCase+10);
                plot(b,b_par,title_bpar,-1.5,1.5,300);
                sprintf(c_par,"L2Selection/HighMomentum2/Case%i/Cpar_assoc_wanted_more500",iCase+10);
                plot(c,c_par,title_cpar,-3000,3000,300);
              }
            }
            if(!isWanted(partic)){
              sprintf(a_par,"L2Selection/HighMomentum2/Case%i/Apar_assoc_Notwanted",iCase+10);
              plot(a,a_par,title_apar,-10e-6,10e-6,300);
              sprintf(b_par,"L2Selection/HighMomentum2/Case%i/Bpar_assoc_Notwanted",iCase+10);
              plot(b,b_par,title_bpar,-1.5,1.5,300);
              sprintf(c_par,"L2Selection/HighMomentum2/Case%i/Cpar_assoc_Notwanted",iCase+10);
              plot(c,c_par,title_cpar,-3000,3000,300);
            }
          }
          if(!HitsAssociated){
            sprintf(a_par,"L2Selection/HighMomentum2/Case%i/Apar_Notassoc",iCase+10);
            plot(a,a_par,title_apar,-10e-6,10e-6,300);
            sprintf(b_par,"L2Selection/HighMomentum2/Case%i/Bpar_Notassoc",iCase+10);
            plot(b,b_par,title_bpar,-1.5,1.5,300);
            sprintf(c_par,"L2Selection/HighMomentum2/Case%i/Cpar_Notassoc",iCase+10);
            plot(c,c_par,title_cpar,-3000,3000,300);
          }
#endif
          if (msgLevel(MSG::DEBUG)) debug()<<"Parabola Par"
                                           <<"\n a \t"<<a
                                           <<"\n b \t"<<b
                                           <<"\n c \t"<<c<<endmsg;
          //===================================================
          // Look in all the other layers except the
          // 1st/last/zone except the parabolaSeedHit
          //===================================================
          //Loop on all the xZones
          //for ( std::vector<PrHitZone*>::iterator itZ = xZones.begin(); xZones.end() != itZ; ++itZ )
          for (PrHitZone* xZone : xZones){
            if (msgLevel(MSG::DEBUG)) debug()<<"Selecting ParSeedHits"<<endmsg;
            if((int)xZone->planeCode() == (int)parabolaSeedHits[i]->planeCode()) continue;
            double dz   = xZone->z() - m_geoTool->zReference();
            double xAtZ = a * dz * dz + b * dz + c; //Parabolic computation
            if (m_useCubic){
              xAtZ= a * dz * dz * (1. + m_dRatio* dz) + b * dz + c; //with Cubic Correction
            }
            // double xMaxAtZ = xAtZ + 2.*fabs(tx_pickedcombination) + m_tolRemaining[iCase];
            // double xMinAtZ = xAtZ - 2.*fabs(tx_pickedcombination) - m_tolRemaining[iCase];
            double xMaxAtZ = xAtZ + m_tolRemaining[iCase];
            // // //std::fabs(tx_pickedcombination)+0.5;
            double xMinAtZ = xAtZ - m_tolRemaining[iCase];
            //std::fabs(tx_pickedcombination)-0.5;
            PrHit* bestProj = nullptr;
            double  bestDist = 10.0; //2.0 mm at the moment (Larger)? tighter? (see offline Seeding)
            if (xMinAtZ > xMaxAtZ){ // should never happen
              if (msgLevel(MSG::DEBUG)) debug()<<"Bad Settings!!!!!!"<<endmsg;
            }
            PrHits::iterator itH = std::lower_bound(xZone->hits().begin() ,xZone->hits().end(),xMinAtZ,lowerBoundX());
            PrHit* hit;
            for (; xZone->hits().end() != itH; ++itH){
              hit = *itH;
              if (hit->isUsed() && m_removeFlagged) continue; //try to remove it ? allow for the remaining layers to pick up flagged hits?
              if (hit->x() < xMinAtZ ) continue;
              if (hit->x() > xMaxAtZ ) break;
              //Find Hit with Best distance <2.0mm
              if(std::fabs(hit->x() - xAtZ)  <  bestDist){
                bestDist = std::fabs(hit->x() - xAtZ);
                if (msgLevel(MSG::DEBUG)) debug()<<"I found an Hit from projection"<<endmsg;
                bestProj = hit;
              }
            }
            if(bestProj != nullptr){
              xHits.push_back(bestProj);
            }
          }//end loop xZones
          //in xHits are not present the first layer and last + parabola seed hits
          if (msgLevel(MSG::DEBUG)) debug()<<"End Loop in between zones to pick up Projection of parabola"<<endmsg;
          // xHits.push_back( parabolaSeedHits[i]);
          // Add to the xHits Vector the remaining 3 Hits not considered
          if(xHits.size() < 2 ) continue; //next parabolaSeedHits ; you must find 2 hits
          xHits.push_back(fHit);
          xHits.push_back(parabolaSeedHits[i]);
          xHits.push_back(lHit);
          // xHits.push_back( parabolaSeedHits[i]);
          // xHits.push_back( fHit);
          // xHits.push_back( lHit);
          if(xHits.size()>6){
            always()<<"Smething goes wrong!!!! in the creation of the xHits list"<<endmsg;
            always()<<"xHits is bigger than 6 : ERROR"<<endmsg;
          }
          //end parabola Seed Hits loop in other Layers
          //Still in the L0 loop (selection last layer)
          //at this step we have 1 Hit in 1st Layer
          //at this step we have 1 Hit in last Layer
          //at this step we have 1 Hit in Middle Layer
          //at this step we have Hit in remaining X layers at the
          //best distance to extrapolated parabola All of them are
          //inside xHits i want to have at least min_HitXSize
          //UNDER STUDY CASE 0 Only Keep Tracks found with 6 Hits (reduce ghost rate if add UV too
          //1st Case keep only 6 Hits on found track
          //2nd Case keep tracks with 4/5/6 hits
          //if( xHits.size() <=  m_minXPlanes)  continue; //Require at least m_minXPlanes ! no for first Fit
          //std::stable_sort(xHits.begin(), xHits.end(), compX());
          //is it really needed?

          // stable_sort?
          std::sort(xHits.begin(), xHits.end(), compLHCbID());
          bool isEqual = false;
          // Remove xHits in the xHitsLists which are basically the same
          for(PrHits& hits : xHitsLists){
            if(msgLevel(MSG::DEBUG)) debug()<<"looping on xHitsLists"<<endmsg;
            if(hits == xHits){
              isEqual = true;
              break;
            }
          }
          if(!isEqual){
            if (msgLevel(MSG::DEBUG)) debug()<<"Pushing Back xHits List"<<endmsg;
            xHitsLists.push_back( xHits);
          }
        }//End loop parabolaSeedHits
        if (msgLevel(MSG::DEBUG)) debug()<<"End Loop For pick up Parabola Hits and build the xHitsLists"<<endmsg;
        //End loop Parabola Seed Hits
        //-------- Remove Duplicates from search in parabolaSeedHits
        if (msgLevel(MSG::DEBUG)) debug()<<"xHitsLists size before removing duplicates: "<<xHitsLists.size()<<endmsg;
        if(xHitsLists.size() == 0){
          continue;
        }

        // Why removing duplicates only if we have found more than 2 HitsList?
        // the 2 HitsList found could be the same, or couldn't they?
        if(xHitsLists.size() > 2){
          //---Remove Duplicates in the HitsList
          std::stable_sort( xHitsLists.begin(), xHitsLists.end());
          xHitsLists.erase( std::unique(xHitsLists.begin(), xHitsLists.end()), xHitsLists.end());
        }
        if (msgLevel(MSG::DEBUG)) debug()<<"xHitsLists size after removing duplicates: "<<xHitsLists.size()<<endmsg;
        //Now let's fit the track
        for(PrHits& xHits : xHitsLists){
          if(msgLevel(MSG::DEBUG)) debug()<<"Fit Track"<<endmsg;
          //Create the track
          PrSeedTrack2 temp_track( part , m_geoTool->zReference() , xHits); //Create the track
          //Setters for it: usefull later to parametrise
          //I load in the track these info which are then plotted
          // if(m_useCubic){
          //   temp_track.setdRatio(m_dRatio);
          // }
          //-----------------------------------------------------
          //----------------O-_The Fit_-O------------------
          //-----------------------------------------------------
          int nIter = 0;
          //bool doRefit = true;
          //temp_track.setRefitX(nIter);
          if (msgLevel(MSG::DEBUG) ){ debug()<<"Attempting to Fit the following Track"<<endmsg; printTrack(temp_track);}
          bool OK = false;
          if(temp_track.hits().size()>m_minXPlanes){ //no 4 hits at first fit
            OK = fitXProjection(temp_track,iCase);
          }
          while(!OK && temp_track.hits().size()>m_minXPlanes){
            if(temp_track.hits().size() <=m_minXPlanes){
              OK = false;
              break;
            }
            nIter++;
            if( nIter==1 && temp_track.hits().size() == 5){ OK = false; break;}
            if( temp_track.hits().size() > m_minXPlanes){
              OK = removeWorstAndRefitX(temp_track,iCase);
            }
          }
          if( OK ){
            setChi2X(temp_track);
          }
#ifdef TRUTH_MATCH_Histos
          Tuple nTupleXZ = nTuple("Seeding/FitXZProjection/HighMomentum2/BeforeStoringXCand","BeforeStoringXCandHighP2");
          PrPlaneCounter2 counter;
          counter.set(temp_track.hits().begin(),temp_track.hits().end());
          LHCb::MCParticle* partic = nullptr;
          double efficiency = -1;
          int nAssocHits = -10;
          bool assoc = AssocTrack(temp_track,efficiency,partic,nAssocHits);
          nTupleXZ->column("assoc",assoc);
          nTupleXZ->column("Case",(int)iCase+10);
          nTupleXZ->column("OK",OK);
          nTupleXZ->column("nIter",temp_track.RefitX());
          nTupleXZ->column("MaxChi2Hit",temp_track.MaxChi2());
          nTupleXZ->column("Chi2PerDoF",temp_track.chi2()/(temp_track.hits().size()-3));
          nTupleXZ->column("nbSingleX",(int)counter.nbSingleX());
          nTupleXZ->column("nXHits",(int)temp_track.hits().size());
          nTupleXZ->column("ax",temp_track.ax());
          nTupleXZ->column("bx",temp_track.bx());
          nTupleXZ->column("cx",temp_track.cx());
          nTupleXZ->column("X0Back",temp_track.X0());
          if(partic!=nullptr){
            nTupleXZ->column("P",partic->p());
            nTupleXZ->column("isWanted",isWanted(partic));
            nTupleXZ->column("isXOk",counter.isOKX());
            nTupleXZ->column("nAssocX",nAssocHits);
          }
          if(partic==nullptr){
            nTupleXZ->column("P",-999.);
            nTupleXZ->column("isWanted",false);
            nTupleXZ->column("isXOk",counter.isOKX());
            nTupleXZ->column("nAssocX",nAssocHits);
          }
          nTupleXZ->write();
#endif
          if( OK &&
              temp_track.hits().size() >= m_minXPlanes
              && (( temp_track.chi2PerDoF() < m_maxChi2DoFX[iCase]))){
            temp_track.setCase( iCase );
            m_xCandidates.push_back(temp_track); //The X Candidate is created
          }
        }//end Loop xHist:xHitsLists
    }//end loop Last Zone given a firsZone selected
  }//end loop first zone
}





#ifdef TRUTH_MATCH_Histos
bool PrHybridSeeding::matchKey( PrHit* hit, int key){
  LinkedTo<LHCb::MCParticle> fLink( evtSvc() , msgSvc(), LHCb::FTClusterLocation::Default );
  LHCb::FTChannelID id0 = hit->id().ftID();
  LHCb::MCParticle* partic = fLink.first(id0);
  while(0!=partic){
    if(key == partic->key() ) return true;
    partic=fLink.next();
  }
  return false;
}
bool PrHybridSeeding::isWanted(LHCb::MCParticle* mcPart){
  if(mcPart ==nullptr) return false;
  MCTrackInfo trackInfo( evtSvc(),msgSvc());
  if( !trackInfo.hasT(mcPart)) return false;
  if( m_etaCut && (mcPart->pseudoRapidity()>5. || mcPart->pseudoRapidity()<2.)) return false;
  if( m_noElectrons && (std::fabs(mcPart->particleID().pid() )==11)) return false;
  if( mcPart->originVertex()->isPrimary() || (mcPart->originVertex()->isDecay())) return true;

  return false;
}



bool PrHybridSeeding::AssocTrack(PrSeedTrack2 track, double& efficiency,LHCb::MCParticle*& particle, int& nHits){
  std::map<int,int> Counter; //map key and number of associated
  std::map<LHCb::MCParticle*,int> CounterParticle;
  Counter.clear();
  for(PrHits::iterator hit = track.hits().begin();hit!=track.hits().end();++hit){
    LHCb::FTChannelID id = (*hit)->id().ftID();
    LinkedTo<LHCb::MCParticle> ftLink1(evtSvc(),msgSvc(),LHCb::FTClusterLocation::Default);
    LHCb::MCParticle* part = ftLink1.first(id);
    while(0!=part){
      int key = part->key();
      Counter[key]++;
      CounterParticle[part]++;
      part=ftLink1.next();
    }
  }
  efficiency = -1;
  nHits = -1;
  if(!CounterParticle.empty()){
    std::vector<std::pair<LHCb::MCParticle*,int> > vpar(std::begin(CounterParticle),std::end(CounterParticle));
    std::sort(std::begin(vpar),std::end(vpar),[](const std::pair<LHCb::MCParticle*,int>& a, const std::pair<LHCb::MCParticle*,int>& b){return a.second<b.second;});
    std::pair<LHCb::MCParticle*,int> parmax = vpar.back();
    particle = parmax.first;
    efficiency = (double)parmax.second/(double)track.hits().size();
    nHits = parmax.second;  
  }
  Counter.clear();
  if( efficiency>0.9 )return true;
  if( efficiency<0.9 )return false;
  return false;
   
}
#endif
