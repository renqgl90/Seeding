#include <TROOT.h>

class Hit {
private:
  float m_x;
  float m_y;
  float m_z;
  int   m_size;
  float m_Err;
  float m_w;
  float m_w2;
  float m_dxDy;
  float m_xatyEq0;
 public:
  Hit(float x,float y, float z)//Creator for X layers
  { 
    m_x    =x; 
    m_y    =y; 
    m_z    =z; 
    m_size =-1;
    m_Err  = 0.;
    m_w    = 0.;
    m_w2   = 0.;
    m_dxDy = 0.;
    m_xatyEq0 = x;
  }
  Hit(float x,float y, float z ,float dxDy)//Creator for u-v layers
  { 
    m_x    =x; //true
    m_y    =y; //true
    m_z    =z; //true
    m_size =-1;
    m_Err  = 0.;
    m_w    = 0.;
    m_w2   = 0.;
    m_dxDy = dxDy;
    m_xatyEq0 = x-dxDy*y;
  }
  
  Hit() {}
  /* friend  bool operator < (Hit &Hit1,Hit &Hit2);// Minoranza */
  /* friend  bool operator > (Hit &Hit1,Hit &Hit2);// Maggioranza */
  /* friend  bool compZ (Hit &Hit1,Hit &Hit2); // Compara Z */
  /* Getters */
  float GetX() {return m_x;}// Return X
  float GetY() {return m_y;}// Return Y
  float GetZ() {return m_z;}// Return Z
  float GetXat0() {return m_xatyEq0;}
  float w() {return m_w;}
  float w2() {return m_w2;}
  float Err() {return m_Err;}
  float dxDy(){return m_dxDy;}
  float x(float y) {return m_xatyEq0+m_dxDy*y;}
  //float x(float y){mx+m_dxDy*y}
  int   GetSize() {return m_size;}
  bool isInside() {if (m_x>-4000 && m_x<4000//Is inside T stations
		       && m_y>-3000 && m_y<3000
		       && m_z>7000 && m_z<11000)
      return true;
    else return false;
  }
  bool isX1(){//1st X Layer
    if (m_z>7840 && m_z<7870) return true;
	else return false;
  }
  bool isX2(){//2nd X Layer
    if (m_z>8020 && m_z<8055) return true;
    else return false;
  }
  bool isX3(){//3rd X Layer
    if ( m_z>8520 && m_z<8555) return true;
    else return false;
  }
  bool isX4(){ //4th XLayer
    if ( m_z>8705 && m_z<8740) return true;
    else return false;
  }
  bool isX5(){//5th Layer
    if ( m_z>9200 && m_z<9250) return true;
    else return false;
  }
  bool isX6(){//6th Layers
    if  ( m_z>9390 && m_z<9430) return true;
    else return false;
  }
  bool isX() {
    if ((m_z>7840 && m_z<7870) //layer0
	|| ( m_z>8020 && m_z<8055) //layer 3
	|| ( m_z>8520 && m_z<8555) //layer 4
	|| ( m_z>8705 && m_z<8740) //layer 7
	|| ( m_z>9200 && m_z<9250) //layer 8
	|| ( m_z>9390 && m_z<9430))
      return true;
    else return false;
  }//layer 11
  bool isT1()//1st / 4 Stations
  {
    if (m_z>7840 && m_z<8520)
      return true;
    else return false;
  }
  bool isT2()//2nd / 4 Stations
  {
    if (m_z>8520 && m_z<9200)
      return true;
    else return false;
  }
  bool isT3()//3rd / 4 Stations
  {
    if (m_z>9200 && m_z<9430)
      return true;
    else return false;
  }


  
  /* Setters*/
  void  SetSize(int size) {
    m_size = size;
    m_Err = 0.05 + 0.03*(size);
    m_w=1./(m_Err);
    m_w2=1./(m_Err*m_Err);
  }
  void SetXYZ(float x, float y,float z) {m_x = x; m_y=y; m_z=z;}//Set 3D infos
  void SetXYZErr(float x, float y,float z,float err)
  {
    m_x = x; m_y=y; m_z=z;
    m_Err=err;
    m_w=1./(err);
    m_w2=m_w*m_w;
  }
  void SetW2(float w){m_w2 = w;
    m_w = std::sqrt(w);
    m_Err = 1./(m_w);
  }
  void SetDxDy(float dxDy) {m_dxDy = dxDy; m_xatyEq0 = m_x - dxDy*m_y;}
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
//#endif
