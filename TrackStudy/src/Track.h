#ifndef TRACK_H
#define TRACK_H

#include <iostream>
#include "./PatHit.h"
/** @class PrSeedTrack PrSeedTrack.h
*  This is the working class inside the T station pattern reco
*  Adapeted to cubic Fit and also to work locally with the Hit.h
*  @author Renato Quagliani
*  @date   2015-03-20
*/

class PrSeedTrack {
public:
  PrSeedTrack(Int_t zone, Float_t zRef):
  m_zone(zone),m_zRef(zRef){
    init();
  };
  PrSeedTrack(Int_t zone, Float_t zRef ,std::vector<PatHit>& hits):
  m_zone(zone),
  m_zRef(zRef),
  m_hits(hits){
    init();
  };
  virtual ~PrSeedTrack() {};

  std::vector<PatHit>& hits() {return m_hits;}
  const std::vector<PatHit>& hits() const {return m_hits;}
  void addHit(PatHit hit) {m_hits.push_back(hit);}
  void addHits(std::vector<PatHit>& hits){
    m_hits.reserve(m_hits.size()+hits.size());
    m_hits.insert(m_hits.end(),hits.begin(),hits.end());
  }
  Int_t zone() const{return m_zone;}
  void setParameters(Float_t ax, Float_t bx, Float_t cx, Float_t ay, Float_t by){
    m_ax=ax;
    m_bx=bx;
    m_cx=cx;
    m_ay=ay;
    m_by=by;
  }
  void updateParameters(Float_t dax, Float_t dbx, Float_t dcx, Float_t day=0., Float_t dby=0.){
    m_ax+=dax;
    m_bx+=dbx;
    m_cx+=dcx;
    m_ay+=day;
    m_by+=dby;
  }
  Float_t x(Float_t z) const{
    Float_t dz = z-m_zRef;
    return m_ax + dz*(m_bx + dz*m_cx*(1+m_dRatio*dz));
  }
  Float_t xSlope( Float_t z) const{
    Float_t dz = z-m_zRef;
    return m_bx + 2. * dz * m_cx + 3.*dz*dz*m_cx*m_dRatio;
  }
  Float_t y(Float_t z) const{
    return m_ay+(z-m_zRef)*m_by;
  }
private:
  float  m_zRef;
  std::vector<PatHit> m_hits;
  bool  m_valid;
  float m_ax;
  float m_bx;
  float m_cx;
  float m_ay;
  float m_by;
  float m_chi2;
  int   m_nDoF;
  // //   float m_dXCoord;
  float m_meanDy;
  float m_dXCoord;
  float m_dRatio;
  Int_t m_zone;
private:
  void init(){
    m_hits.reserve(32);
    m_ax=0.;
    m_ax = 0.;
    m_bx = 0.;
    m_cx = 0.;
    m_by = 0.;
    m_ay = 0.;
    m_chi2 = 0.;
    m_nDoF = -1;
    m_dXCoord = 0.;
    m_meanDy  = 0.;
  }


};
#endif // #ifdef FakedSeeding_cxx
