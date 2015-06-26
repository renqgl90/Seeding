#include <TROOT.h>
#ifdef MCHit
#define MCHit_1
class MCHit{
private:
  Float_t m_X;
  Float_t m_Y;
  Float_t m_Z;
  Float_t m_ty;
  Float_t m_tx;
  Float_t m_p;
  Float_t m_pathLength;
  Float_t m_MCParticle_p;
public:
  MCHit():
  m_X(0);
  m_Y(0);
  m_Z(0);
  m_ty(0);
  m_tx(0);
  m_p(0);
  m_pathLength(0);
  m_pParticle(0);
  m_pxParticle(0);
  m_pyParticle(0);
  m_pzParticle(0);
  void setMCHit(const Float_t X,Float_t Y,Float_t Z,Float_t tx,Float_t ty,Float_t p,
    Float_t pathlenght,Float_t MCpart_p,Float_t MCPart_px,MCpart_py,MCpart_pz ){
    m_X=X;
    m_Y=Y;
    m_Z=Z;
    m_tx=tx;
    m_ty=ty;
    m_p=p;
    m_pathLength=pathlenght;
    m_pParticle=p;
    m_pxParticle=px;
    m_pyParticle=py;
    m_pzParticle=pz;
  }
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
}
#endif
