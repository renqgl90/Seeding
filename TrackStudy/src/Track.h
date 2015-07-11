#ifndef TRACK_H
#define TRACK_H 1

#include <iostream>
#include <iomanip>
#include "./PatHit.h"
#include "./FTCluster.h"
/** @class PrSeedTrack PrSeedTrack.h
*  This is the working class inside the T station pattern reco
*  Adapeted to cubic Fit and also to work locally with the Hit.h
*  @author Renato Quagliani
*  @date   2015-03-20
*/
using std::cout;
using std::endl;
using std::setw;
using std::internal;
class PrSeedTrack{
//public:
  //PrSeedTrack(zone,zRef,std::vector<PatHit>)
private:
  Float_t m_zRef;
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


public:
  PrSeedTrack(const Int_t zone, const Float_t zRef){
    m_zone = zone;
    m_zRef = zRef;
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
  PrSeedTrack(const Int_t zone, const Float_t zRef ,const std::vector<PatHit> hhits){
    m_zone = zone;
    m_zRef = zRef;
    m_hits.reserve(hhits.size());
    m_hits=hhits;
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
  virtual ~PrSeedTrack() {};

  const std::vector<PatHit> hits() const {return m_hits;}
  //const std::vector<PatHit>& hits() const {return m_hits;}
  void addHit(PatHit hit){m_hits.push_back(hit);}
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
  Float_t ax() const{return m_ax;}
  Float_t bx() const{return m_bx;}
  Float_t cx() const{return m_cx;}
  Float_t xSlope( Float_t z) const{
    Float_t dz = z-m_zRef;
    return m_bx + 2. * dz * m_cx + 3.*dz*dz*m_dx;
  }
  Float_t y(Float_t z) const{
    return m_ay+(z-m_zRef)*m_by;
  }

  Float_t x(Float_t z) const{
    Float_t dz = z-m_zRef;
    return m_ax+m_bx*dz+m_cx*dz*dz;//+ m_dx*dz*dz;
  }
  Float_t distance(PatHit hit)const{
    Float_t y0 = m_ay - m_zRef*m_by;
    Float_t dyDz = m_by;
    Float_t yTra = (y0 + dyDz*hit.z(0.))/(1.- dyDz*hit.dzDy());
    return hit.x(yTra)-x(hit.z(yTra));
  }
  Float_t Chi2Hit(PatHit hit) const{
    Float_t d = distance(hit);
    return d*d*hit.w2();
  }
  Float_t Chi2() const{
    Float_t Chi2=0;
    for(Int_t i = 0;i<m_hits.size(); i++){
      Chi2+= Chi2Hit(m_hits[i]);
    }
    //Chi2;
    return Chi2;
  }




  void PrintTrack(){
    std::cout.precision(6);
    std::cout<<"i"<<setw(20)<<"Hit X"<<setw(20)<<"Track X"<<setw(20)<<"Chi2 Contrib"<<std::endl;
    for(Int_t i =0; i<m_hits.size(); i++){
      Float_t distance2oe = m_hits[i].w2()*std::pow((x(m_hits[i].z(0.))-m_hits[i].x(0.)),2);
      std::cout<<i<<setw(20)<<m_hits[i].x(0.)<<setw(20)<<x(m_hits[i].z(0.))<<setw(20)<<Chi2Hit(m_hits[i])<<std::endl;
    }
  }
  void PrintHits()
  {
    std::cout.precision(3);
    std::cout<<" Nb. Hits on the track =   "<<m_hits.size()<<endl;
    std::cout<<"i"<<setw(20)<<"1/w"<<setw(20)<<"Xat0"<<setw(20)<<"Zat0"<<setw(20)<<"planeCode"<<setw(20)<<"zone"<<setw(20)<<"dxDy"<<setw(20)<<"dzDy"<<setw(20)<<"Fraction"<<setw(20)<<"Size"<<setw(20)<<"Charge \n"<<std::scientific;
    for(Int_t i =0; i<m_hits.size(); i++){
      //m_hits[i].PrintHit();
      std::cout<<i<<setw(20)<<sqrt(1./m_hits[i].w2())<<setw(20)<<m_hits[i].x(0.)<<setw(20)<<m_hits[i].z(0.)<<setw(20)<<m_hits[i].planeCode()<<setw(20)<<m_hits[i].zone()<<setw(20)<<m_hits[i].dxDy()<<setw(20)<<m_hits[i].dzDy()<<setw(20)<<m_hits[i].cluster().charge()<<setw(20)<<m_hits[i].cluster().fraction()<<setw(20)<< m_hits[i].cluster().size() <<setw(20)<<m_hits[i].cluster().charge()<<std::endl;
    }
  }
protected:
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
