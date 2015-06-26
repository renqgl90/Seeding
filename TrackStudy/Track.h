//#ifndef TRACK_H
//#define TRACK_H 1

#include <iostream>
#include "Hit.h"
/** @class PrSeedTrack PrSeedTrack.h
 *  This is the working class inside the T station pattern reco
 *  Adapeted to cubic Fit and also to work locally with the Hit.h
 *  @author Renato Quagliani
 *  @date   2015-03-20
 */

class PrSeedTrack {
private:
    
  float  m_zRef;
  std::vector<Hit> m_hits;
  bool  m_valid;
  float m_ax;
  float m_bx;
  float m_cx;
  float m_ay;
  float m_by;
  float m_chi2;
  int   m_nDoF;
  float m_dXCoord;
  float m_meanDy;
  float m_dRatio;
public:
  PrSeedTrack(const double zRef ) 
  {
    m_zRef = zRef;
    m_valid = true;
    m_hits.reserve( 32 );
    m_ax = 0.;
    m_bx = 0.;
    m_cx = 0.;
    m_by = 0.;
    m_ay = 0.;
    m_chi2 = 0.;
    m_nDoF = -1;
    m_dXCoord = 0.;
    m_meanDy  = 0.;
    m_dRatio = 0;
  }; 
  /// Constructor with list of hits
  PrSeedTrack( float zRef, std::vector<Hit>& hits )
  {
    m_zRef = zRef;
    m_valid = true;
    m_hits.reserve( 32 );
    m_hits=hits;
    m_ax = 0.;
    m_bx = 0.;
    m_cx = 0.;
    m_by = 0.;
    m_ay = 0.;
    m_chi2 = 0.;
    m_nDoF = -1;
    m_dXCoord = 0.;
    m_meanDy  = 0.;
  }; 

  virtual ~PrSeedTrack( ) {}; ///< Destructor
  
  /// Handling of hits: acceess, insertion
  const std::vector<Hit> hits()       const { return m_hits; }
  //void addHit( Hit* hit )        { m_hits.push_back( hit );}
  void setdRatio( float dRatio){m_dRatio = dRatio;}
  void addHits( std::vector<Hit>& hits ) {
    m_hits.reserve( m_hits.size() + hits.size() );
    m_hits.insert( m_hits.end(), hits.begin(), hits.end() );
  }
  /// Parameters of the trajectory in the T stations
  void setParameters( float ax, float bx, float cx, float ay, float by ) {
    m_ax = ax;
    m_bx = bx;
    m_cx = cx;
    m_ay = ay;
    m_by = by;
  }
  void setZref(float zRef)
  {m_zRef = zRef;}
  void printHits()
    {
      std::cout<<"Hits On Track"
               <<"\n X \t Y \t Z \t Size"<<std::endl;
      for(int i =0 ; i<m_hits.size();i++)
      {
        std::cout<<m_hits[i].GetX()<<"\t"<<m_hits[i].GetY()<<"\t"<<m_hits[i].GetZ()<<"\t"<<m_hits[i].GetSize()<<std::endl;
      }
    }
  void printTrack()
  {
    std::cout<<"Hits On Track"
               <<"\n X \t \t Y \t Z"<<std::endl;
      for(int i =0 ; i<m_hits.size();i++)
      {
        std::cout<<m_hits[i].GetX()<<"\t"<<m_hits[i].GetY()<<"\t"<<m_hits[i].GetZ()<<"\t"<<std::endl;
      }
      std::cout<<"Track Params"
               <<"\n dax \t\t dbx \t\t dcx \t\t day \t\t dby \t\t zRef \t\t dRatio \n"
               <<m_ax<<"\t"<<m_bx<<"\t"<<m_cx<<"\t"<<m_ay<<"\t"<<m_by<<"\t"<<m_zRef<<"\t"<<m_dRatio<<std::endl;
  }
  void updateParameters( float dax, float dbx, float dcx,
                         float day=0., float dby= 0.  ) {
    m_ax += dax;
    m_bx += dbx;
    m_cx += dcx;
    m_ay += day;
    m_by += dby;
  }

  float x( float z )         const { float dz = z-m_zRef; return m_ax + dz * ( m_bx + dz*(m_cx*(1-m_dRatio*dz))); }
  float xSlope( float z )    const { float dz = z-m_zRef; return m_bx + 2. * dz * m_cx; }
  float y( float z )         const { return m_ay + (z-m_zRef) *  m_by; }
  float ySlope( )            const { return m_by; }  
  float distance( Hit* hit ) const { 
    return hit->GetX() - x(hit->GetZ());
  }
  float chi2( Hit* hit )     const { float d = distance( hit ); return d * d * hit->w2(); }

  float deltaY( Hit* hit )   const {
    if ( hit->isX() ) return 0.;
    return distance( hit ) / hit->dxDy();
  }

  bool valid()                  const { return m_valid; }
  void setValid( bool v )             { m_valid = v; }

  void setChi2( float chi2, int nDoF ) { m_chi2 = chi2; m_nDoF = nDoF; }
  float chi2()                  const { return m_chi2; }
  float chi2PerDoF()            const { return m_chi2 / m_nDoF; }
  int   nDoF()                  const { return m_nDoF; }

  void setDXCoord( float dxCoord )    { m_dXCoord = dxCoord; }
  float dXCoord()               const { return m_dXCoord; }

  void setMeanDy( float meanDy )      { m_meanDy = meanDy; }
  float meanDy()                const { return m_meanDy; }
  float a() const {return m_ax;}
  float b() const {return m_bx;}
  float c() const {return m_cx;}
  struct GreaterBySize {
    bool operator() (const PrSeedTrack& lhs, const PrSeedTrack rhs ) const { return lhs.hits().size() > rhs.hits().size(); }
  };


protected:

private:
  void init() {
    m_valid = true;
    m_hits.reserve( 32 );

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
//#endif // #ifdef FakedSeeding_cxx
