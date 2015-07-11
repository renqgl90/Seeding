#ifndef PatHit_H
#define PatHit_H 1
#include <TROOT.h>
#include "./FTCluster.h"
class PatHit{

public:
  PatHit( ):
  m_xatyEq0(0),
  m_zatyEq0(0),
  m_coord(0),
  m_w(0),
  m_w2(0),
  m_yMin(0),
  m_yMax(0),
  m_zone(-1),
  m_planeCode(-1),
  m_id(-1.),
  m_dxDy(-999.),
  m_dzDy(-999.),
  m_isX(false),
  m_Cluster(FTCluster())
  {};
  virtual ~PatHit(){};

  void setHit(const Float_t xat0 , const Float_t zat0, const Float_t dxDy, const Float_t dzDy, const Float_t w ,const Float_t yMin, const Float_t yMax, const int zone, const int planeCode, Bool_t isX , Float_t id)
  {
    m_xatyEq0=xat0;
    m_zatyEq0=zat0;
    m_dxDy=dxDy;
    m_dzDy=dzDy;
    m_w=w;
    m_w2=w*w;
    m_yMin = yMin;
    m_yMax = yMax;
    m_planeCode = planeCode;
    m_zone = zone;
    m_id = id;
    m_isX = isX;
    m_coord = xat0;
  }
  Float_t yOnTrack(float y0, float dyDz) const{
    return (y0+dyDz*m_zatyEq0)/(1.-dyDz*m_dzDy) ;
  }
  Float_t coord() const{ return m_coord;}
  void setCoord(Float_t coord){
    m_coord = coord;
    }
  void PrintHit(){
    //std::cout<<i<<"\t"<<m_xa0t <<"\t"<<z(0.)<<"\t"<<m_hits[i].planeCode()<<"\t"<<m_hits[i].cluster().charge()<<"\t"<<m_hits[i].cluster().fraction()<<"\t"<<m_hits[i].cluster().size()<<"\t"<<m_hits[i].cluster().charge()<<std::endl;
  }
  void SetCluster(FTCluster cluster) {m_Cluster=cluster;}
  Float_t x(Float_t y=0) const{return m_xatyEq0+y*m_dxDy;}
  Float_t z(Float_t y=0) const{return m_zatyEq0+m_dzDy*y;}
  Int_t planeCode() const{return m_planeCode;}
  Int_t zone() const{return m_zone;}
  Float_t dxDy() const{return m_dxDy;}
  Float_t dzDy() const{return m_dzDy;}
  FTCluster cluster() const{return m_Cluster;}
  Float_t w2() const{ return m_w2;}
  Float_t w() const{ return m_w;}
private:

  Float_t m_w;
  Float_t m_w2;
  Float_t m_dxDy;
  Float_t m_dzDy;
  Float_t m_xatyEq0;
  Float_t m_zatyEq0;
  Float_t m_coord;
  Float_t m_yMin;
  Float_t m_yMax;
  Int_t m_zone;
  Float_t m_id;
  Int_t m_planeCode;
  Bool_t m_isX;
  FTCluster m_Cluster;



};
#endif
