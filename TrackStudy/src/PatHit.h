#ifndef PatHit_H
#define PatHit_H 1
#include <TROOT.h>

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
  m_id(-1),
  m_dxDy(-999.),
  m_dzDy(-999.),
  m_isX(false)
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

  Float_t x(Float_t y) const{return m_xatyEq0+y*m_dxDy;}
  Float_t z(Float_t y) const{return m_zatyEq0+m_dzDy*y;}
private:
  //Float_t m_xatyEq0;
  //Float_t m_zatyEq0;
  //Float_t m_z;
  //Int_t   m_size;
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


};
#endif
