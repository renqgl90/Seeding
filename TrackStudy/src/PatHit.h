#include <TROOT.h>
#ifdef PatHit
#define PatHit_1
class PatHit {
private:
  Float_t m_x;
  Float_t m_y;
  Float_t m_z;
  int   m_size;
  Float_t m_Err;
  Float_t m_w;
  Float_t m_w2;
  Float_t m_dxDy;
  Float_t m_xatyEq0;
  Float_t m_zatyEq0;

 public:
    PatHit():
    m_charge(0);
    m_xatyEq0(0);
    m_zatyEq0(0);
    m_coord(0);
    m_w(0);
    m_w2(0);
    m_yMin(0);
    m_yMax(0);
    m_zone(-1);
    m_planeCone(-1);
    m_id(-1);
    m_dxDy(-999.);
    m_dzDy(-999.);
    m_isX(false);
    void setHit(const Float_t xat0 , const Float_t zat0, const Float_t dxDy, const Float_t dzDy, const Float_t w , const Float_t yMin, const Float_t yMax, const int zone, const int planeCode, bool isX){
      m_xatyEq0=xat0;
      m_zatyEq0=zat0;
      m_dxDy=dxDy;
      m_dzDy=dzDy;
      m_Err=error;
      m_w=w;
      m_w2=w2;
      m_yMin = yMin;
      m_yMax = yMax;
      m_planeCode = planeCode;
      m_zone = zone;
      m_id = id;
      m_isX = isX;
      m_coord = xat0;
    }
    Float_t yOnTrack(float y0, float dyDz){
    return (y0+dyDz*m_zatyEq0)/(1.-dyDz*m_dzDy);
    }
    Float_t coord(){return m_coord;}
    Float_t setCoord(Float_t coord){m_coord = coord;}
// Hit(float x, float z , float )//Creator for X layers
// {
//   m_x    =x;
//   m_y    =y;
//   m_z    =z;
//   m_size =-1;
//   m_Err  = 0.;
//   m_w    = 0.;
//   m_w2   = 0.;
//   m_dxDy = 0.;
//   m_xatyEq0 = x;
//   m_zatyEq0 = z;
// }
// Hit(float x,float y, float z ,float dxDy)//Creator for u-v layers
// {
//   m_x    =x; //true
//   m_y    =y; //true
//   m_z    =z; //true
//   m_size =-1;
//   m_Err  = 0.;
//   m_w    = 0.;
//   m_w2   = 0.;
//   m_dxDy = dxDy;
//   m_xatyEq0 = x-dxDy*y;
//   m_zatyEq0 =
// }
//
// Hit() {}
// /* friend  bool operator < (Hit &Hit1,Hit &Hit2);// Minoranza */
// /* friend  bool operator > (Hit &Hit1,Hit &Hit2);// Maggioranza */
// /* friend  bool compZ (Hit &Hit1,Hit &Hit2); // Compara Z */
// /* Getters */
// float GetX() {return m_x;}// Return X
// float GetY() {return m_y;}// Return Y
// float GetZ() {return m_z;}// Return Z
// float GetXat0() {return m_xatyEq0;}
// float w() {return m_w;}
// float w2() {return m_w2;}
// float Err() {return m_Err;}
// float dxDy(){return m_dxDy;}
// float x(float y) {return m_xatyEq0+m_dxDy*y;}
// float z(float y) {return m_zatyEq0+m_dzDy*y}
//
// //float x(float y){mx+m_dxDy*y}
// int   GetSize() {return m_size;}
// bool isInside() {if (m_x>-4000 && m_x<4000//Is inside T stations
// 	       && m_y>-3000 && m_y<3000
// 	       && m_z>7000 && m_z<11000)
//     return true;
//   else return false;
// }
//
// bool isX1(){//1st X Layer
//   if (m_z>7840 && m_z<7870) return true;
// else return false;
// }
// bool isX2(){//2nd X Layer
//   if (m_z>8020 && m_z<8055) return true;
//   else return false;
// }
// bool isX3(){//3rd X Layer
//   if ( m_z>8520 && m_z<8555) return true;
//   else return false;
// }
// bool isX4(){ //4th XLayer
//   if ( m_z>8705 && m_z<8740) return true;
//   else return false;
// }
// bool isX5(){//5th Layer
//   if ( m_z>9200 && m_z<9250) return true;
//   else return false;
// }
// bool isX6(){//6th Layers
//   if  ( m_z>9390 && m_z<9430) return true;
//   else return false;
// }
  // bool isX() {
  //   if ((m_z>7840 && m_z<7870) //layer0
	// || ( m_z>8020 && m_z<8055) //layer 3
	// || ( m_z>8520 && m_z<8555) //layer 4
	// || ( m_z>8705 && m_z<8740) //layer 7
	// || ( m_z>9200 && m_z<9250) //layer 8
	// || ( m_z>9390 && m_z<9430))
  //     return true;
  //   else return false;
  // }//layer 11
  // {
  // bool isT1()//1st / 4 Stations
  //   if (m_z>7840 && m_z<8520)
  //     return true;
  //   else return false;
  // }
  // bool isT2()//2nd / 4 Stations
  // {
  //   if (m_z>8520 && m_z<9200)
  //     return true;
  //   else return false;
  // }
  // bool isT3()//3rd / 4 Stations
  // {
  //   if (m_z>9200 && m_z<9430)
  //     return true;
  //   else return false;
  // }
  //
  //
  //
  // /* Setters*/
  // void  SetSize(int size) {
  //   m_size = size;
  //   m_Err = 0.05 + 0.03*(size);
  //   m_w=1./(m_Err);
  //   m_w2=1./(m_Err*m_Err);
  // }
  // void SetXYZ(float x, float y,float z) {m_x = x; m_y=y; m_z=z;}//Set 3D infos
  // void SetXYZErr(float x, float y,float z,float err)
  // {
  //   m_x = x; m_y=y; m_z=z;
  //   m_Err=err;
  //   m_w=1./(err);
  //   m_w2=m_w*m_w;
  // }
  // void SetW2(float w){m_w2 = w;
  //   m_w = std::sqrt(w);
  //   m_Err = 1./(m_w);
  // }
  // void SetDxDy(float dxDy) {m_dxDy = dxDy; m_xatyEq0 = m_x - dxDy*m_y;}
  // struct LowerByCoord(const Hit &lhs, const Hit &rhs) const {return lhs.coord()<rhs.coord;}
  //struct LowerByZ(const Hit &lhs, const Hit &rhs) const{return lhcb}
};

  //bool isInside(Hit &Hit1)
  //{

/* bool compZ( Hit &Hit1,  Hit &Hit2)// Compare Z  */
/* { */
/*   return Hit1.GetZ() < Hit2.GetZ(); */
/* } */

/* bool operator< ( Hit &Hit1, Hit &Hit2)// Less than in Z  */
/* { */
/*   return Hit1.GetZ() <Hit2.GetZ(); */
/* } */
/* bool operator> ( Hit &Hit1, Hit &Hit2)// Greater in Z */
/* { */
/*   return Hit1.GetZ() >Hit2.GetZ(); */
/* } */

//bool sortByZ(const Hit &Hit1, const Hit &Hit2) { return Hit1.Z < Hit2.Z; }
//bool sortByX(const Hit &Hit1, const Hit &Hit2) { return Hit1.X < Hit2.X; }
//bool sortByY(const Hit &Hit1, const Hit &Hit2) { return Hit1.Y < Hit2.Y; }

//typedef vector<Hit> Hits;
#endif
