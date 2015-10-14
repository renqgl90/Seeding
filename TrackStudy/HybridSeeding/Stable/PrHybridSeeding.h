#ifndef PRHYBRIDSEEDING_H 
#define PRHYBRIDSEEDING_H 1



// Include files 
// from Gaudi
#include "PrHybridSeedTrack.h"

//uncomment this line if you want to do truth matching and do plots
//#include "GaudiAlg/GaudiAlgorithm.h"

#include "GaudiAlg/ISequencerTimerTool.h"
#include "PrKernel/IPrDebugTool.h"
#include "PrKernel/PrHitManager.h"
#include "PrGeometryTool.h"
#include "TfKernel/RecoFuncs.h"
#include "Event/MCParticle.h"

/** @class PrHybridSeeding PrHybridSeeding.h
 *  
 *
 *  @author Renato Quagliani
 *  @date   2015-03-11
 */
//#define TRUTH_MATCH_Histos
#ifdef TRUTH_MATCH_Histos
#include "GaudiAlg/GaudiTupleAlg.h"
#else
#include "GaudiAlg/GaudiAlgorithm.h"
#endif

#ifdef TRUTH_MATCH_Histos
class PrHybridSeeding : public GaudiTupleAlg{
#else
  class PrHybridSeeding : public GaudiAlgorithm{
#endif
  public: 
    /// Standard constructor
    PrHybridSeeding( const std::string& name, ISvcLocator* pSvcLocator );
    virtual ~PrHybridSeeding( ); ///< Destructor
    virtual StatusCode initialize();    ///< Algorithm initialization
    virtual StatusCode execute   ();    ///< Algorithm execution
    virtual StatusCode finalize  ();    ///< Algorithm finalization
    
  protected:
    bool matchKey( const PrHit* hit ) {
      if ( m_debugTool ) return m_debugTool->matchKey( hit->id(), m_wantedKey );
      return false;
    };
    bool matchKey( const PrHybridSeedTrack& track ) {
      if ( !m_debugTool ) return false;
    for ( PrHits::const_iterator itH = track.hits().begin(); track.hits().end() != itH; ++itH ) {
      if ( m_debugTool->matchKey( (*itH)->id(), m_wantedKey ) ) return true;
    }
    return false;
    };
    
    void setKey( int key ) {
      m_wantedKey = key;
      return;
    };
    /** @brief Collect Hits in X layers for High momentum tracks
     * @param part (if 1, y<0 ; if 0 , y>0)
     */
    void findXProjections( unsigned int part, unsigned int iCase);
    /** @brief Collect hits in the stereo-layers.
     *   @param part lower (1) or upper (0) half
     */
    void addStereo( unsigned int part , unsigned int iCase);
    

    
    bool LineOK( double chi2low, double chi2high, PrLineFitterY line, PrHybridSeedTrack& temp);
    
    
    /** @brief Remove Clones of produced tracks
     */
    void removeClonesX(unsigned int maxCommon, unsigned int part, unsigned int iCase , bool xOnly);
    void removeClones(unsigned int maxCommon);
    
    /** @brief Set the Chi2 value and Chi2DoF for the track after the X+U-V search
     */
    void setChi2(PrHybridSeedTrack& track);
   
    /** @brief Flag Hits under conditions at the end of the iCase loop 
     */
    void flagHits(unsigned int icase, unsigned int part);
    /** @bried Flag Hits under conditions (similar to removal of clones ) sort of hybrid
     */
    void flagHits2();
    /** @brief Fit the track combining the XZ and YZ projections
     *  @param track The track to fit
     *  @param Refit Iteration in the Refitting after removal worst hit
     *  @return bool Success of the XY+XZ Fit
     **/
    bool fitSimultaneouslyXY( PrHybridSeedTrack& track ,unsigned int iCase);
    
    /** @brief Fit the track combining the only in the XZ plane
     *  @param track The track to fit
     *  @param Refit Iteration in the Refitting after removal worst hit
     *  @return bool Success of the XZ Fit
     **/
    bool fitXProjection( PrHybridSeedTrack & track ,unsigned int iCase );

    /** @brief Fit Y line given a list of UV hits
     *  @param hits List of UV hits from UV layers
     *  @param plCount PlaneCoutner object ot handle the removal of UV hits
     *  @return vector of hits and planeCounter updated
    **/



    /** @brief Remove the hit which gives the largest contribution to the chi2 and refit XZ
     *  @param track The track to fit
     *  @return bool Success of the fit
     */
    bool removeWorstAndRefitX( PrHybridSeedTrack& track , unsigned int iCase );
    
     /** @brief Remove the hit which gives the largest contribution to the chi2 and refit XZ + YZ
      *  @param track The track to fit
     *  @return bool Success of the fit
     */
    bool removeWorstAndRefit( PrHybridSeedTrack& track , unsigned int iCase );
    
     /** @brief Set the chi2 of the track
      *  @param track The track to set the chi2 of 
      */
    void setChi2X( PrHybridSeedTrack& track );
    /** @brief Transform the tracks from the internal representation into LHCb::Tracks
     *  @param tracks The tracks to transform
     */
    void makeLHCbTracks( LHCb::Tracks* result );
    
    /** @brief Print some information of the hit in question
     *  @param hit The hit whose information should be printed
     *  @param title Some additional information to be printed
     */
    void printHit( const PrHit* hit, std::string title="" );
    
    /** @brief Print some information of the track in question
     *  @param hit The track whose information should be printed
     */
    void printTrack( PrHybridSeedTrack& track );
    
    /** @brief Internal method to construct parabolic parametrisation out of three hits, using Cramer's rule.
   *  @param hit1 First hit
   *  @param hit2 Second hit
   *  @param hit3 Third hit
   *  @param a quadratic coefficient
   *  @param b linear coefficient
   *  @param c offset
   */
    void solveParabola(const PrHit* hit1, const PrHit* hit2, const PrHit* hit3, double& a, double& b, double& c);
    /** @brief Internal method to construct parabolic parametrisation + cubic correction included out of three hits, using Cramer's rule.
   *  @param hit1 First hit
   *  @param hit2 Second hit
   *  @param hit3 Third hit
   *  @param a quadratic coefficient
   *  @param b linear coefficient
   *  @param c offset
   */
    void solveParabola2(const PrHit* hit1, const PrHit* hit2, const PrHit* hit3, double& a1, double& b1, double& c1);
    
    /// Classe to find lower bound of x of PrHits
    class lowerBoundX {
    public:
      bool operator() (const PrHit* lhs, const double testval ) const { return lhs->x() < testval; }
    };
    /// Classe to find upper bound of x of PrHits
    class upperBoundX {
    public:
      bool operator() (const double testval, const PrHit* rhs) const { return testval < rhs->x(); }
    };
    /// Class to compare x positions of PrHits
    class compX {
    public:
      bool operator() (const PrHit* lhs, const PrHit* rhs ) const { return lhs->x() < rhs->x(); }
    };
    ///Class to find lower bound of LHCbIDs of PrHits
    class lowerBoundLHCbID {
    public:
      bool operator() (const PrHit* lhs, const LHCb::LHCbID id ) const { return lhs->id() < id; }
    };
    ///Class to compare LHCbIDs of PrHits
    class compLHCbID {
    public:
      bool operator() (const PrHit* lhs, const PrHit* rhs ) const { return lhs->id() < rhs->id(); }
    };
    
  private:
    //-------------Names for input container(if forward as input), output container name, HitManager name
    bool m_removeHighP;
    bool m_recoverTrack;
    double m_chi2Recover ;
    bool m_SlopeCorr;
    
    std::string        m_inputName;
    std::string        m_outputName;
    std::string        m_hitManagerName;
    //-------------Global configuration of the algorithm
    bool               m_decodeData;
    bool               m_xOnly;
    unsigned int   m_minXPlanes;
    unsigned int   m_maxNHits;
    unsigned int   m_nCases;
    bool               m_doTiming;
    bool               m_printSettings;
    bool               m_useCubic;
    bool               m_removeClones;
    bool               m_removeClonesX;
    bool               m_FlagHits;
    bool               m_removeFlagged;
    //------------X-search parametrisation
    //1st / Last Layer search windows
    
    std::vector<double> m_alphaCorrection;
    std::vector<double> m_TolFirstLast;
    

    //Add of the third hit in middle layers (p and Pt dependent, i.e., case dependent)
    std::vector<double> m_x0Corr;
    std::vector<double> m_x0SlopeChange;
    std::vector<double> m_TolX0SameSign;
    std::vector<double> m_x0Cut;
    std::vector<double> m_tolAtX0Cut;
    
    std::vector<double> m_tolX0Oppsig;
    std::vector<double> m_x0SlopeChange2;
    std::vector<double> m_x0CutOppSig;
    std::vector<double> m_tolAtx0CutOppSig;
    
    //Add of remaining Hits in remaining X Layers
    std::vector<double> m_tolRemaining;
    //-----------
    //Add of third hit in T2
    unsigned int m_maxParabolaSeedHits;
    //Look up in remaining X layers
    //Look up in remaining X layers and Track Fitting parameters
    double           m_dRatio;
    double           m_ConstC;
    //--------_Fit X parametrisation
    std::vector<double>            m_maxChi2HitsX;
    std::vector<double>            m_maxChi2DoFX;
    
    //--------_Full Fit parametrisation
    std::vector<double>            m_maxChi2FullFit;
    std::vector<double>            m_maxChi2HitFullFitHigh;
    std::vector<double>            m_maxChi2HitFullFitLow;
    std::vector<double>            m_maxY0Low;
    std::vector<double>            m_maxYZrefLow;
    std::vector<double>            m_maxChi2HitLow;
    
    std::vector<double>            m_maxChi2HitFull ;
    
    //Add the parameters for the Chi2 Cut
    unsigned int     m_nSigmaOffX;
    unsigned int     m_X0Max;
    
    
    //-----------------UV search parametrisation
    //Added
    double          m_yMin;
    double          m_yMin_TrFix;
    double          m_yMax;
    double          m_yMax_TrFix;
    double          m_doAsymmUV;
    //
    double          m_coord;
    //Triangle Fix
    bool            m_useFix; 
    bool            m_removeHole;
    bool            m_useFix2ndOrder;

    //LineY 
    bool m_useLineY;
    std::vector<double> m_Chi2LowLine;
    std::vector<double> m_Chi2HighLine;
    
    
    std::vector<double> m_X0ChangeCoord;
    
    std::vector<double> m_tolTyOffset;
    std::vector<double> m_tolTySlope;
    

    //X+Y fit configure
    std::vector<double>          m_maxChi2PerDoF;
    bool            m_corrLeftSide;
    //Clone removal setting
    std::vector<unsigned int>              m_nCommonX;
    unsigned int              m_nCommonUV;
    bool                      m_ClonesUpDown;
    unsigned int              m_nCommonUVTriangle;
    unsigned int              m_nUsed;
    //Flag Hits Settings
    
    std::vector<double>         m_MaxChi2Flag;
    std::vector<double>         m_MaxX0Flag;
    std::vector<unsigned int>   m_SizeFlag;
    
    //dRatio correction to use (temporary)
    bool m_useCorrPos;
    bool m_useCorrSlopes;
    bool m_useImproved;
    
    
    std::map<LHCb::LHCbID , std::vector< LHCb::MCParticle*>> m_map;
    

    //--------------------Global things

    PrHitManager*   m_hitManager;
    PrGeometryTool* m_geoTool;
    
    //== Debugging controls
    std::string     m_debugToolName;
    int             m_wantedKey;
    IPrDebugTool*   m_debugTool;
    
    
    //-------------------Containers
    std::vector<PrHybridSeedTrack>       m_trackCandidates;
    std::vector<PrHybridSeedTrack>       m_xCandidates;
    std::vector<PrHitZone*>        m_zones;
    ISequencerTimerTool* m_timerTool;
    int            m_timeTotal;
    int            m_timeFromForward;
    
    int            m_timeXProjeUp[3];
    int            m_timeStereoUp[3];
    int            m_timeCloneXUp[3];
    int            m_timeFlagUp[3];
    
    int            m_timeXProjeDo[3];
    int            m_timeStereoDo[3];
    int            m_timeCloneXDo[3];
    int            m_timeFlagDo[3];
    
    
    int            m_timeClone;
    int            m_timeConvert;
    
    //-------------------Algorithms Constants

    // Stations Zones numbering
    static const unsigned int s_T1X1 = 0;
    static const unsigned int s_T1U  = 2;
    static const unsigned int s_T1V  = 4;
    static const unsigned int s_T1X2 = 6;
    static const unsigned int s_T2X1 = 8;
    static const unsigned int s_T2U  = 10;
    static const unsigned int s_T2V  = 12;
    static const unsigned int s_T2X2 = 14;
    static const unsigned int s_T3X1 = 16;
    static const unsigned int s_T3U  = 18;
    static const unsigned int s_T3V  = 20;
    static const unsigned int s_T3X2 = 22;
    static const unsigned int s_down = 0;
    static const unsigned int s_up   = 1;

#ifdef TRUTH_MATCH_Histos
  public:
    bool matchKey( PrHit* hit,int key);
    bool isWanted( LHCb::MCParticle* mcPart);
    bool AssocTrack(PrHybridSeedTrack track,double& efficiency,LHCb::MCParticle*& particle, int& nHits);
  private:
    bool m_etaCut;
    bool m_noElectrons;
#endif


  };
  
#endif // PRHYBRIDSEEDING_H
