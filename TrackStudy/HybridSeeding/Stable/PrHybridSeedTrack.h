#ifndef PRHYBRIDSEEDTRACK_H 
#define PRHYBRIDSEEDTRACK_H 1

// Include files
#include "Event/Track.h"
#include "PrKernel/PrHit.h"

/** @class PrHybridSeedTrack PrHybridSeedTrack.h
 *  This is the working class inside the T station pattern
 *
 *  @author Renato Quagliani
 *  @date   2015-05-13
 */
//Comment it if you don't want to do truth matching
class PrHybridSeedTrack {
public: 
  /// Constructor with the z reference
  PrHybridSeedTrack( unsigned int zone, double zRef ) 
    : m_zone( zone ), m_zRef( zRef )
  {
    init();
  };
  
  PrHybridSeedTrack( unsigned int zone, double zRef, PrHits& hits )
    : m_zone( zone ), m_zRef( zRef ), m_hits( hits ) 
  {
    init();
  }; 
  
  virtual ~PrHybridSeedTrack( ) {}; ///< Destructor
  //Handling the hits on the track
  PrHits& hits()                   { return m_hits; }
  const PrHits& hits()       const{ return m_hits;  }
  
  //Add hit
  void addHit( PrHit* hit )        { m_hits.push_back( hit );}
  void addHits( PrHits& hits ) {
    m_hits.reserve( m_hits.size() + hits.size() );
    m_hits.insert( m_hits.end(), hits.begin(), hits.end() );
  }
  void addHits2( std::vector<PrHit*>& hits){
    m_hits.reserve( m_hits.size() + hits.size());
    m_hits.insert( m_hits.end(), hits.begin() , hits.end());
  }
  
  //Zone (0,1) Up Down
  unsigned int zone() const { return m_zone; }
  //====================================
  //SETTERS & Getters  Genearl One
  //====================================
  //***********Track Parameters to compute distances
  void setParameters( double ax, double bx, double cx, double ay, double by ) {
    m_ax = ax;
    m_bx = bx;
    m_cx = cx;
    m_ay = ay;
    m_by = by;
  }
  double ax()                      const { return m_ax;}
  double bx()                      const { return m_bx;}
  double cx()                      const { return m_cx;}
  double ay()                      const { return m_ay;}
  double by()                      const { return m_by;}
  
  //**********Update the parameters of a track iteratively in the Fit
  void updateParameters( double dax, double dbx, double dcx,
                         double day=0., double dby= 0.  ) {
    m_ax += dax;
    m_bx += dbx;
    m_cx += dcx;
    m_ay += day;
    m_by += dby;
  }
  void setYParam( double ay=0., double by=0.){
    m_ay = ay;
    m_by = by;
  }
  // void setnDiff( unsigne int i){
  //   m_nDiff = i;
  // }
  
  void setnXnY( unsigned int nx, unsigned int ny){
    m_nx = nx;
    m_ny = ny;
  }
  unsigned int nx() const{return m_nx;}
  
  unsigned int ny() const{return m_ny;}
  
  
  //validity
  void setValid( bool v )               { m_valid = v; }
  bool   valid()                 const { return m_valid; }
  //dRatio
  void setdRatio( double dRatio )       { m_dRatio = dRatio;}
  double dRatio()                 const{return m_dRatio;}
  //Chi2 & DoF
  void setChi2LineY( double chi2DoF, int nUV){
    m_chi2Line = chi2DoF;
    m_nUVLine = nUV;
  }
  double chi2DoFLine() const
  {
    if( m_nUVLine>1)
      return m_chi2Line;
    return 1e30;
  }
  int nUV() const{
    return m_nUVLine;
  }
  
  
  
  
  void setChi2(double chi2, int nDoF) { m_chi2 = chi2;m_nDoF = nDoF; }
  double chi2()                   const { return m_chi2; }
  double chi2PerDoF()         const { return m_chi2 / m_nDoF;}
  int nDoF()                        const { return m_nDoF; }
  //MeanDy
  double meanDy()             const { return m_meanDy; }
  void setMeanDy( double meanDy )        { m_meanDy = meanDy;}
  //MaxChi2OnHit
  void setMaxChi2(double val1)                { m_maxChi2 = val1;}
  double MaxChi2()              const { return m_maxChi2;}
 
  //BackProjection
  void setX0(double X0)                           { m_X0=X0;}
  double X0()                      const { return m_X0;}

  //RefitX status
  int RefitX()                        const { return m_RefitX;}
  void setRefitX(int Refit)                         { m_RefitX = Refit;}
  
  //dXCoord
  void setDXCoord( double dxCoord )        { m_dXCoord = dxCoord;}  
  double dXCoord()             const { return m_dXCoord; }
  //y(at a given z)
  double y( double z )          const { return (m_ay + m_by*z);}
  //Slope by
  double ySlope( )               const { return m_by; }
  
  //X at a given Z
  double x( double z )const{
    const double dz = z-m_zRef;
    return m_ax + dz * ( m_bx + dz* m_cx * (1.+m_dRatio*dz) );
  }
  
  //Slope X-Z plane track at a given z
  double xSlope( double z )const{ 
    const double dz = z-m_zRef;
    return m_bx + 2.*dz*m_cx + 3.*dz*dz*m_cx*m_dRatio; //is it need?
  }
  void setCase(unsigned int Case){
    m_case=Case ;
  }
  unsigned int Case() const{
    return m_case;
  }
  
  
  //y positon on Track
  double yOnTrack( PrHit* hit ) const{ return hit->yOnTrack( m_ay , m_by ); }
  
  //distance from Hit of the track
  double distance( PrHit* hit ) const{ 
    double yTra = yOnTrack(hit);
    return hit->x( yTra ) - x( hit->z(yTra) );
  }
  
  //Chi2 contribution from a single Hit
  double chi2( PrHit* hit )           const{ 
    const double d = distance( hit );
    double w = hit->w();
    if(m_slopeCorr){
      double yTra = yOnTrack(hit);
      const double Corr = std::cos(xSlope( (double) hit->z(yTra) ));
      w = hit->w()/(Corr*Corr);
    }
    return d * d * w;
  }
  
  void setSlopeCorr( bool corr){
    m_slopeCorr = corr;
  }
  bool slopeCorr() const{
    return m_slopeCorr;
  }
  
  //DeltaY
  double deltaY( PrHit* hit )        const{
    if ( hit->isX() ) return 0.;
    return distance( hit ) / hit->dxDy();
  }
  void sortbyz(){
    std::sort(m_hits.begin(),m_hits.end(),[](const PrHit* hit1, const PrHit* hit2)->bool{return hit1->z()<hit2->z();});
  }
  struct GreaterBySize {//there was no & in rhs
    bool operator() (const PrHybridSeedTrack& lhs, const PrHybridSeedTrack& rhs ) const { 
      if( lhs.hits().size() != rhs.hits().size() )
        return lhs.hits().size() < rhs.hits().size(); //< by default
      return lhs.chi2PerDoF() < rhs.chi2PerDoF();
    }
  };
  struct LowerBySize{
    bool operator() ( const PrHybridSeedTrack& lhs, const PrHybridSeedTrack& rhs) const{
      if(lhs.hits().size() == rhs.hits().size()) return lhs.chi2PerDoF() > rhs.chi2PerDoF();
      return lhs.hits().size() > rhs.hits().size();
    }
  };
  
  
    
  // struct LowerBySize{
  //   bool operator() ( const PrHybridSeedTrack& lhs, const PrHybridSeedTrack rhs) const { return lhs.hits().size() < rhs.hits().size(); }
  // }


  int size() const{
    return m_hits.size();
  };
  void setxPredT1Mid( double xPred)
  {
    m_xPredT1 = xPred;
  }
  double xPredT1()  const
  {
    return m_xPredT1;
  }
  
  void setxPredT2Mid( double xPred)
  {
    m_xPredT2 = xPred;
  }
  double xPredT2()  const
  {
    return m_xPredT2;
  } 
  void setxPredT3Mid( double xPred)
  {
    m_xPredT3 = xPred;
  }
  
  double xPredT3()  const
  {
    return m_xPredT3;
  }
  
protected:
  
private:

  void init(){
    m_slopeCorr = false;
    m_valid = true;
    m_hits.reserve( 32 );
    m_case = -1;
    m_nx = 0;
    m_ny = 0;
    m_nUVLine = 0;
    m_chi2Line = 1.e30;
    //track params
    m_ax = 0.;
    m_bx = 0.;
    m_cx = 0.;
    m_by = 0.;
    m_ay = 0.;
    //tracl chi2 settings
    m_chi2 = 0.;
    m_nDoF = -1;
    m_dXCoord = 0.;
    m_meanDy  = 0.;
    m_dRatio = 0.;
    m_X0 = -99999.;
    m_maxChi2 = -99999.;
    m_xPredT1 = 0.;
    m_xPredT2 = 0.;
    m_xPredT3 = 0.;
  }
  //Global Private variables
private:
  bool m_slopeCorr;
  unsigned int m_nx;
  unsigned int m_ny;
  int m_nUVLine;
  double m_chi2Line;
  
  unsigned int m_zone;
  double       m_zRef;
  PrHits       m_hits;
  bool         m_valid;
  int          m_nDoF;
  double       m_chi2;
  double       m_ax;
  double       m_bx;
  double       m_cx;
  double       m_ay;
  double       m_by;
  double       m_dRatio;
  double       m_meanDy;
  double       m_maxChi2;
  double       m_X0;
  int          m_RefitX;
  double       m_dXCoord;
  double       m_xPredT1;
  double       m_xPredT2;
  double       m_xPredT3;
  unsigned int m_case;
};

typedef std::vector<PrHybridSeedTrack> PrHybridSeedTracks;

#endif // PRSEEDTRACK2_H
