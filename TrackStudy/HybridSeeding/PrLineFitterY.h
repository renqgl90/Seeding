#ifndef PRLINEFITTERY_H
#define PRLINEFITTERY_H 1
// Include files

/** @class PrLineFitterY.h
 * Small class to perform a fit for a straight line in y, returning the Chi2/ Chi2 PerDoF/ ay, by
 * @author Renato Quagliani
 * @date 2015-07-30
 */

class PrLineFitterY{
public:
  //Constructor
  PrLineFitterY( double zref){
    m_s0 = 0.;
    m_sz = 0.;
    m_sz2 = 0.;
    m_zref = zref;

  }

  void addHit(){}
  void removeHit(){}
  void Chi2(){}

  void setAndSolve( PrHits::iterator itBeg, PrHits::iterator itEnd , PrSeedTrack2 track){
    PrHit *hit = nullptr;
    for( int i = 0; i<3; i++){
    for( PrHits::const_iterator hit = itBeg; itEnd != hit; ++hit){
      double y =(hit->x()-track.x(hit->z()))/hit->dxDy();
      double dz = ( hit->z() - m_zref);
      double dist = y - (m_ay + dz* m_by);
      double w = hit->w(); //need als the dxDy
      m_mat[0]+=err;
      m_mat[1]+=err*dz;
      m_mat[2]+=err*dz*dz;
      m_rhs[0]+=err*dist;
      m_rhs[1]+=err*dist*dz;
    }
  }
  bool Fit(){
     for()
 }
  bool solve(){
     ROOT::Math::CholeskyDecomp<double,2> decomp(m_mat);
     if(!decomp){
        return false;
     }
     decomp.Solve(m_rhs);
     m_ay

 }
  m_zref = zref;
  private:
     double m_chi2;
     double m_chi2perDof;
     double m_ay;
     double m_by;
     std::vector<PrHit> m_hits;

  /*| 1 ,   dz | [ay] = \  dy \
  /*| dz, dz*dz| [by] = \ dydz\
}
