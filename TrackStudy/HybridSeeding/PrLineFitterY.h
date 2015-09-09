#ifndef PRLINEFITTERY_H 
#define PRLINEFITTERY_H 1

// Include files
#include "Math/CholeskyDecomp.h"
#include "PrKernel/PrHit.h"
#include "PrSeedTrack2.h"
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
  PrLineFitterY(double zRef , PrSeedTrack2 track):
    m_zRef(zRef),
    //m_track(track),
    m_ax(track.ax()),
    m_bx(track.bx()),
    m_cx(track.cx()),
    m_dRatio(track.dRatio())
  {
    std::fill(m_mat,m_mat+3,0.);
    std::fill(m_rhs,m_rhs+2,0.);
    m_done = false;
    m_minCoord = 1.;
    m_maxCoord = -1.;
    m_hits.reserve(6);
    m_ay = 0.;
    m_by = 0.;
    m_chi2line = 1.e10;
  };
  double zref() const
  {
    return m_zRef;
  }
  double minCoord() const{ return m_minCoord;}
  double maxCoord() const{ return m_maxCoord;}
  // const PrLineFitterY& operator=(const PrLineFitterY& other){
  //   m_dRatio = other.m_dRatio;
  //   m_ay = other.m_ay;
  //   m_by = other.m_by;
  //   m_ax = other.m_ax;
  //   m_bx = other.m_bx;
  //   m_cx = other.m_cx;
  //   m_zRef = other.m_zRef;
  //   m_minCoord = other.m_minCoord;
  //   m_maxCoord = other.m_maxCoord;
  //   m_hits = other.m_hits;
  //   m_done = other.m_done;
  //   // m_mat = other.m_mat;
  //   // m_rhs = other.m_rhs;
  //   return *this;
  // }
  void setdone( bool done){
    m_done = done;
  }
  
  void reset(){
    m_hits.clear();
    m_ay=0.;
    m_by=0.;
    m_chi2line = 1.e10;
    std::fill(m_mat,m_mat+3,0.);
    std::fill(m_rhs,m_rhs+2,0.);
    m_minCoord= 1.;
    m_maxCoord=-1.;
  }
  void addHit( PrHit* hit){
    if( !hit->isX())
    {
      if( std::fabs( hit->coord())> m_maxCoord) m_maxCoord = std::fabs(hit->coord());
      if( std::fabs( hit->coord())<m_minCoord) m_minCoord = std::fabs(hit->coord());
      m_hits.push_back( hit);
      const double w = hit->w();
      const double dz = hit->z() - m_zRef;
      m_mat[0]+=w;
      m_mat[1]+=w*dz;
      m_mat[2]+=w*dz*dz;
    }
  }
  void set( std::vector<PrHit*> hits){
    std::fill(m_mat,m_mat+3,0.);
    for(PrHits::const_iterator hit = hits.begin();hits.end()!= hit; ++hit){
      if( (*hit)->isX() ) continue;
      m_hits.push_back( (*hit));
      if( std::fabs((*hit)->coord()) > m_maxCoord) m_maxCoord = std::fabs((*hit)->coord());
      if( std::fabs((*hit)->coord()) < m_minCoord) m_minCoord = std::fabs((*hit)->coord());
      const double w =  (*hit)->w();
      const double dz = (*hit)->z() - m_zRef;
      //just fill the matrix on the left hand side of the fit
      m_mat[0]+=w;
      m_mat[1]+=w*dz;
      m_mat[2]+=w*dz*dz;
    }
  }
  double x(const double z) const{
    const double dz = z - m_zRef;
    return ( m_ax + m_bx*dz + m_cx*dz*dz*(1.+m_dRatio*dz));
  }
  double distance( PrHit *hit) const{
    // double yTra = yOnTrack( hit);
    const double z = hit->z();
    return ( hit->x() - x( z ))/hit->dxDy() - y( z );
  }
  
  double chi2( PrHit* hit) const{
    const double erry = hit->w()*std::pow(hit->dxDy(),2);
    const double dist = distance( hit);
    return dist*dist*erry;
  }
  double yOnTrack( PrHit *hit) const{ return hit->yOnTrack( m_ay - m_zRef*m_by , m_by);}
  double y( double z) const{ return  m_ay + m_by*(z-m_zRef) ; }
  PrHits& hits(){return m_hits;}
  bool fit(){
    for( int i = 0; i<3; i++){
      std::fill(m_rhs,m_rhs+2,0.);
      // PrHits::iterator itBeg = m_track.hits().begin();
      // PrHits::iterator itEnd = m_track.hits().end();
      for(PrHits::const_iterator hit = m_hits.begin(); m_hits.end()!= hit; ++hit){
        if( (*hit)->isX()) continue;
        const double dz = (*hit)->z()-m_zRef;
        const double dist = distance( (*hit));
        const double w = (*hit)->w();
        m_rhs[0]+=w*dist;
        m_rhs[1]+=w*dist*dz;
      }
      ROOT::Math::CholeskyDecomp<double,2> decomp(m_mat);
      if(!decomp){ //if(msgLevel(MSG::DEBUG)) std::cout<<"UNvable to Fit"<<std::endl;
        return false;
      }
      decomp.Solve(m_rhs);
      m_ay+=m_rhs[0];
      m_by+=m_rhs[1];
    }
    //m_ay-=m_by*m_zRef;
    return true;
  }
  void removeHit( PrHit *hit){
    PrHits::iterator hittoremove = m_hits.begin();
    bool found = false;
    for(PrHits::iterator itH = m_hits.begin(); m_hits.end()!=itH; ++itH){
      if( (*itH)->id() == hit->id()){
        found = true;
        hittoremove = itH;
      }
    }
    if(found){  
      m_hits.erase(hittoremove);
      const double dz = hit->z()-m_zRef;
      const double w  = hit->w(); 
      m_mat[0]-= w;
      m_mat[1]-= w*dz;
      m_mat[2]-= w*dz*dz;
    }
  }
  void SortByZ(){
    std::sort(m_hits.begin(),m_hits.end(),[](const PrHit* hit1, const PrHit* hit2)->bool{return hit1->z()<hit2->z(); });     
  }
  double ay() const{    return m_ay;}
  double ay0() const{ return m_ay-m_by*m_zRef;}
  
  double by() const{    return m_by;}
  double Chi2() const{
    //setdone(true);
    //PrHit * hit= nullptr ;
    double chi2_tot=0.;
    for(PrHits::const_iterator itH = m_hits.begin(); m_hits.end()!=itH; ++itH){ 
      if( (*itH)->isX()) continue;
      chi2_tot+= chi2((*itH));
    }
    return chi2_tot;
  }
  void setChi2Line( double chi2){
    m_chi2line = chi2;
  }
  double Chi2DoF() const{
    return m_chi2line/((double)(m_hits.size()-2.));
  }
  virtual ~PrLineFitterY( ) {} ///< Destructor
  protected:
private:
  double m_minCoord;
  double m_maxCoord;
  double m_zRef;
  double m_mat[3];
  double m_rhs[2];
  double m_ay;
  double m_by;
  double m_chi2line;
  PrHits m_hits;
  double m_ax;
  double m_bx;
  double m_cx;
  double m_dRatio;
  bool m_done;
  
  //PrSeedTrack2 m_track;
};
#endif // PRLINEFITTERY_H
