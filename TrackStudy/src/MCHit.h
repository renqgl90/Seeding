#ifndef MCHit_H
#define MCHit_H
#include <TROOT.h>

class MCHit{

public:
  MCHit( ) :
  m_X(0),
  m_Y(0),
  m_Z(0),
  m_ty(0),
  m_tx(0),
  m_p(0),
  m_pathLength(0),
  m_pParticle(0),
  m_pxParticle(0),
  m_pyParticle(0),
  m_pzParticle(0),
  m_planeCode(-1),
  m_zone(-1),
  m_isX(false)
  {};
  virtual ~MCHit(){};

  void setMCHit(const Float_t X,Float_t Y,Float_t Z,Float_t tx,Float_t ty,Float_t patz, Float_t pathlenght,Float_t MCpart_p,Float_t MCpart_px,Float_t MCpart_py,Float_t MCpart_pz){
    m_X=X;
    m_Y=Y;
    m_Z=Z;
    m_tx=tx;
    m_ty=ty;
    m_p=patz;
    m_pathLength=pathlenght;
    m_pParticle=MCpart_p;
    m_pxParticle=MCpart_px;
    m_pyParticle=MCpart_py;
    m_pzParticle=MCpart_pz;
  }
  void setLayer(Int_t planeCode, Int_t Zone, Bool_t isX){ m_planeCode=planeCode; m_zone=Zone;m_isX=isX;}
  Float_t x() const{return m_X;};
  Float_t y() const{return m_Y;};
  Float_t z() const{return m_Z;};
  Float_t tx() const{return m_tx;};
  Float_t ty() const{return m_ty;};
  Float_t p() const{return m_p;};
  Float_t pxPart() const{return m_pxParticle;};
  Float_t pyPart() const{return m_pyParticle;};
  Float_t pzPart() const{return m_pzParticle;};
  Float_t pTPart() const{return std::sqrt(m_pxParticle*m_pxParticle+m_pyParticle*m_pyParticle);};
  Float_t pPart() const{return m_pParticle;};
  Int_t Layer() const{return m_planeCode;}
  Int_t Zone() const{return m_zone;}
private:
  Float_t m_X;
  Float_t m_Y;
  Float_t m_Z;
  Float_t m_ty;
  Float_t m_tx;
  Float_t m_p;
  Float_t m_pathLength;
  Float_t m_MCParticle_p;
  Float_t m_pParticle;
  Float_t m_pxParticle;
  Float_t m_pyParticle;
  Float_t m_pzParticle;
  Int_t m_zone;
  Int_t m_planeCode;
  Bool_t m_isX;
};
#endif
