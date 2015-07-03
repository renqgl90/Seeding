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
  m_zone(zone),
  m_zRef(zRef){
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
  void setParameters(Float_t ax, Float_t bx, Float_t cx, Float_t ay, Float_t by , Float_t dx=0.){
    m_ax=ax;
    m_bx=bx;
    m_cx=cx;
    m_ay=ay;
    m_by=by;
    m_dx=dx;
  }
  void updateParameters(Float_t dax, Float_t dbx, Float_t dcx, Float_t day=0., Float_t dby=0., Float_t ddx=0.){
    m_ax+=dax;
    m_bx+=dbx;
    m_cx+=dcx;
    m_ay+=day;
    m_by+=dby;
    m_dx+=ddx;
  }

  Float_t xSlope( Float_t z) const{
    Float_t dz = z-m_zRef;
    return m_bx + 2. * dz * m_cx + 3.*dz*dz*m_dx;
  }
  Float_t y(Float_t z) const{
    return m_ay+(z-m_zRef)*m_by;
  }

  Float_t x(Float_t z) const{
    Float_t dz = z-m_zRef;
    return m_ax+m_bx*dz+m_cx*dz*dz + m_dx*dz*dz;
  }
  Float_t distance(PatHit hit)const{
    Float_t y0 = m_ay - m_zRef*m_by;
    Float_t dyDz = m_by;
    Float_t yTra = (y0 + dyDz*hit.z(0.))/(1.- dyDz*hit.dzDy());
    return hit.x(yTra)-x(hit.z(yTra));
  }
  void PrintHits()
  {
    std::cout.precision(4);
    std::cout<<"Nb. Hits on the track = "<<m_hits.size();
    std::cout<<"i \t Xat0 \t Zat0 \t planeCode \t zone \t dxDy \t dzDy"<<"\t Fraction"<<"\t Size"<<"\t Charge"<<std::endl;
    for(Int_t i =0; i<m_hits.size(); i++){
      m_hits[i].PrintHit();
      std::cout<<i<<"\t"<<m_hits[i].x(0.)<<"\t"<<m_hits[i].z(0.)<<"\t"<<m_hits[i].planeCode()<<"\t"<<m_hits[i].cluster().charge()<<"\t"<<m_hits[i].cluster().fraction()<<"\t"<<m_hits[i].cluster().size()<<"\t"<<m_hits[i].cluster().charge()<<std::endl;
    }
  }
private:
  float  m_zRef;
  std::vector<PatHit> m_hits;
  bool  m_valid;
  float m_ax;
  float m_bx;
  float m_cx;
  float m_dx; //not dRatio (just Cubic Term)
  //float m_dRatio;
  float m_ay;
  float m_by;
  float m_chi2;
  int   m_nDoF;
  // //   float m_dXCoord;
  float m_meanDy;
  float m_dXCoord;
  float m_dRatio;
  Int_t m_zone;
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
