#ifndef PRLINEFITTERY_H 
#define PRLINEFITTERY_H 1

// Include files
#include "Math/CholeskyDecomp.h"
#include "PrKernel/PrHit.h"
#include "PrHybridSeedTrack.h"
#include "PrKernel/PrHitManager.h"
/** @class PrLineFitterY PrLineFitterY.h
 *  
 *
 *  @author renato quagliani
 *  @date   2015-08-01
*/
class PrLineFitterY {
public: 
  /// Standard constructor
  PrLineFitterY(double zRef , PrHybridSeedTrack track):
    m_zRef(zRef),
    m_nHitsUV(0),
    m_ax(track.ax()),
    m_bx(track.bx()),
    m_cx(track.cx()),
    m_dRatio(track.dRatio()),
    m_SlopeCorr(track.slopeCorr())
  {
    m_done = false;
    m_minCoord = 1.;
    m_maxCoord = -1.;
    m_ay = 999.;
    m_by = 999.;
    m_chi2line = 1.e14;
  };
  double zref() const
  {
    return m_zRef;
  }
  double minCoord() const{ return m_minCoord;}
  double maxCoord() const{ return m_maxCoord;}
  void setdone( bool done){
    m_done = done;
  }
  
  void reset(){
    m_ay=999.;
    m_by=999.;
    m_chi2line = 1.e14;
    m_minCoord= 1.;
    m_maxCoord=-1.;
    m_nHitsUV = 0;
  }
  // void addHit( PrHit* hit){
  //   if( !hit->isX())
  //   {
  //     if( std::fabs( hit->coord())> m_maxCoord) m_maxCoord = std::fabs(hit->coord());
  //     if( std::fabs( hit->coord())<m_minCoord) m_minCoord = std::fabs(hit->coord());
  //     //m_hits.push_back( hit);
  //     const double w = hit->w();
  //     const double dz = hit->z() - m_zRef;
  //     m_mat[0]+=w;
  //     m_mat[1]+=w*dz;
  //     m_mat[2]+=w*dz*dz;
  //   }
  // }
  // void clearHits(){
  //   m_hits.clear();
  // }
  
  // void setHits( PrHits::iterator itBeg, PrHits::iterator itEnd){
  //   //   std::fill(m_mat,m_mat+3,0.);
  //   //   for(PrHits::iterator hit = itBeg;itEnd!= hit; ++hit){
  //   if( (*hit)->isX() ) continue;
  //   //     m_nHitsUV++;
  //   m_hits.push_back( (*hit));
    //     if( std::fabs((*hit)->coord()) > m_maxCoord) m_maxCoord = std::fabs((*hit)->coord());
    //     if( std::fabs((*hit)->coord()) < m_minCoord) m_minCoord = std::fabs((*hit)->coord());
    //     const double w =  (*hit)->w();
    //     const double dz = (*hit)->z() - m_zRef;
    //     //just fill the matrix on the left hand side of the fit
    //     m_mat[0]+=w;
    //     m_mat[1]+=w*dz;
    //     m_mat[2]+=w*dz*dz;
    //   }
    // }
  //}
  
  double xSlope( const double z) const{
    const double dz = z-m_zRef;
    return m_bx + 2*m_cx*dz + 3*m_cx*m_dRatio*dz*dz;
  }
  
  double x(const double z) const{
    const double dz = z - m_zRef;
    return ( m_ax + m_bx*dz + m_cx*dz*dz*(1.+m_dRatio*dz));
  }
  double distance( PrHit *hit) const{
    const double z = hit->z();
    return ( hit->x(0.) - x( z ) )/hit->dxDy() - y( z );
  }
  
  double chi2( PrHit* hit) const{
    double erry = hit->w()*std::pow(hit->dxDy(),2);
    const double dist = distance( hit);
    if(m_SlopeCorr){
      double z = hit->z();
      double xSlopeVal = xSlope( z ) ;
      double cos = std::cos( xSlopeVal );
      erry = erry/(cos*cos);
    }
    return dist*dist*erry;
  }
  double yOnTrack( PrHit *hit) const{ return hit->yOnTrack( m_ay - m_zRef*m_by , m_by);}
  double y( double z) const{ return  m_ay + m_by*(z-m_zRef);}
  bool fit(PrHits::const_iterator itBeg, PrHits::const_iterator itEnd){
    m_ay = 0.;
    m_by = 0.;
    for( int i = 0; i<3; i++){
      std::fill(m_rhs,m_rhs+2,0.);
      std::fill(m_mat,m_mat+3,0.);
      for(PrHits::const_iterator hit = itBeg; itEnd!= hit; ++hit){
        if( (*hit)->isX()) continue;
        const double dz = (*hit)->z()-m_zRef;
        const double dist = distance( (*hit));
        double w = (*hit)->w();
        if(m_SlopeCorr){
          double cos = std::cos( xSlope( (double)(*hit)->z()));
          w = w/(cos*cos);
        }
        m_mat[0]+=w;
        m_mat[1]+=w*dz;
        m_mat[2]+=w*dz*dz;
        
        m_rhs[0]+=w*dist;
        m_rhs[1]+=w*dist*dz;
      }
      ROOT::Math::CholeskyDecomp<double,2> decomp(m_mat);
      if(!decomp){ //if(msgLevel(MSG::DEBUG)) std::cout<<"UNvable to Fit"<<std::endl;
        m_chi2line = 1e14;
        return false;
      }
      decomp.Solve(m_rhs);
      m_ay+=m_rhs[0];
      m_by+=m_rhs[1];
    }
    m_chi2line = 0;
    for( PrHits::const_iterator hit = itBeg; itEnd!=hit ;++hit){
      if( (*hit )->isX()) continue;
      m_nHitsUV++;
      m_chi2line+= chi2(*hit);
    }
    return true;
  }
  unsigned int nHitsLine(){
    return m_nHitsUV;
  }
  double ay() const{    return m_ay;}
  double ay0() const{ return m_ay-m_by*m_zRef;}
  
  double by() const{    return m_by;}

  void setChi2Line( double chi2){
    m_chi2line = chi2;
  }
  double Chi2DoF() const{
    if(m_nHitsUV >0)
      return m_chi2line/((double)(m_nHitsUV-2.));
    return 1.e30;
  }
  virtual ~PrLineFitterY( ) {} ///< Destructor
  protected:
private:
  double m_zRef;
  unsigned int m_nHitsUV;
  double m_minCoord;
  double m_maxCoord;
  double m_mat[3];
  double m_rhs[2];
  double m_ay;
  double m_by;
  double m_chi2line;
  // PrHits m_hits;
  double m_ax;
  double m_bx;
  double m_cx;
  double m_dRatio;
  bool m_done;
  bool m_SlopeCorr;
  //PrHybridSeedTrack m_track;
};
#endif // PRLINEFITTERY_H
