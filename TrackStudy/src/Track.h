#ifndef TRACK_H
#define TRACK_H 1

#include <iostream>
#include <iomanip>
#include "PatHit.h"
#include "FTCluster.h"
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
   float m_chi2X;
   int m_nDoFX;
   float m_chi2Full;
   int m_nDoFFUll;
   float m_ax1 =0.;
   float m_bx1 = 0.;
   float m_cx1 = 0.;
   float m_dx1 = 0.;
   float m_chi21 =0.;
   float m_nDoF1 =0.;

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

      m_ax1 =0.;
      m_bx1 = 0.;
      m_cx1 = 0.;
      m_dx1 = 0.;
      m_chi21 =0.;
      m_nDoF1 =0.;
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
      m_chi2X=0.;
      m_nDoFX=0.;
      m_chi2Full=0.;
      m_nDoFFUll=0.;
      m_dx = 0.;

   }
   virtual ~PrSeedTrack() {};

   const  std::vector<PatHit> hits() const {return m_hits;}
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
   void setParameters2(Float_t ax1, Float_t bx1, Float_t cx1, Float_t ay, Float_t by , Float_t dx1=0.){
      m_ax1=ax1;
      m_bx1=bx1;
      m_cx1=cx1;
      m_ay=ay;
      m_by=by;
      m_dx1=dx1;
   }
   PatHit hitInPlane(int val) const{
      for(int i=0; i<m_hits.size();i++){
         if(m_hits[i].planeCode()==val){
            PatHit hit = m_hits[i];
            return hit;
         }
      }
      PatHit hit1 = PatHit();
      return hit1;
   }
   void substituteHitinSameLayer(PatHit hit){
      int planeCode = hit.planeCode();
      int index = -1;
      for(int i = 0; i<m_hits.size();i++){
         if(m_hits[i].planeCode() == planeCode){ index = i; break;}
      }
      if(index !=-1){
         m_hits.erase(m_hits.begin()+index);
         addHit(hit);
         sortbyZ();
      }
   }

   void removeAllHitInPlane(int val){
      std::vector<int> indexes;
      for(int i = 0; i<m_hits.size(); i++)
      {
         if( m_hits[i].planeCode() == val) indexes.push_back(i);
      }
      for(int j = 0;j<indexes.size();j++){
         m_hits.erase(m_hits.begin() + indexes[j] -j );
      }
   }
   void updateParameters(Float_t dax, Float_t dbx, Float_t dcx, Float_t day=0., Float_t dby=0., Float_t ddx=0.){
      m_ax+=dax;
      m_bx+=dbx;
      m_cx+=dcx;
      m_ay+=day;
      m_by+=dby;
      m_dx+=ddx;
   }
   void updateParameters2(Float_t dax, Float_t dbx, Float_t dcx, Float_t day=0., Float_t dby=0., Float_t ddx=0.){
      m_ax1+=dax;
      m_bx1+=dbx;
      m_cx1+=dcx;
      m_ay+=day;
      m_by+=dby;
      m_dx1+=ddx;
   }
   Float_t ax() const{return m_ax;}
   Float_t bx() const{return m_bx;}
   Float_t cx() const{return m_cx;}
   Float_t xSlope( Float_t z) const{
      Float_t dz = z-m_zRef;
      return m_bx + 2. * dz * m_cx + 3.*dz*dz*m_dx;
   }
   Float_t y(Float_t z) const{
      return m_ay+m_by*z;
   }
   void setdRatio( Float_t val ){
      m_dx = val;
   }
   Float_t dRatio()const {return m_dx;}
   Float_t x(Float_t z) const{
     Float_t dz = z-m_zRef;
     return m_ax+m_bx*dz+m_cx*dz*dz*(1. + m_dx*dz);//+ m_dx*dz*dz;
   }
   Float_t ay()const{  return m_ay;}
   Float_t by()const{  return m_by;}

   Float_t x_new(Float_t z)const{
      Float_t dz = z-m_zRef;
      return m_ax1+m_bx1*dz+m_cx1*dz*dz + m_dx1*dz*dz*dz;
   }
   Float_t distance2(PatHit hit)const{
      Float_t y0 = m_ay;
      Float_t dyDz = m_by;
      Float_t yTra = (y0 + dyDz *hit.z(0.))/(1.-dyDz*hit.dzDy());
      return hit.x(yTra) - x_new(hit.z(yTra));

   }
   Float_t distance(PatHit hit)const{
      Float_t y0 = m_ay;
      Float_t dyDz = m_by;
      Float_t yTra = (y0 + dyDz*hit.z(0.))/(1. -dyDz*hit.dzDy());
      return hit.x(yTra)-x(hit.z(yTra));
   }
   Float_t Chi2Hit(PatHit hit) const{
      Float_t d = distance(hit);
      return d*d*hit.w2();
   }
   Float_t Chi2_modHit(PatHit hit) const{
      Float_t d2 = distance2(hit);
      return d2*d2*hit.w2();
   }
   Float_t Chi2_mod()const{
      Float_t Chi2 =0;
      for(Int_t i =0; i<m_hits.size();i++){
         Chi2+=Chi2_modHit(m_hits[i]);
      }
      return Chi2;
   }
   Float_t ax_cub() const{
      return m_ax1;
   }
   Float_t bx_cub() const{
      return m_bx1;
   }
   Float_t cx_cub() const{
      return m_cx1;
   }
   Float_t dx_cub() const{
      return m_dx1;
   }
   void setChi2X(Float_t Chi2){ m_chi2X = Chi2; }
   void setNDOFX(Int_t val){m_nDoFX = val;}
   Int_t ndofX() const{return m_nDoFX;}
   void setChi2Full(Float_t Chi2){ m_chi2Full = Chi2; }
   void setNDOFFull(Int_t val){m_nDoFFUll = val;}
   Int_t ndofFull() const{return m_nDoFFUll;}
   

   Float_t Chi2() const{
      Float_t Chi2=0.;
      for(Int_t i = 0;i<m_hits.size(); i++){
        Chi2+= Chi2Hit(m_hits[i]);
      }
      return Chi2;
   }
   
   /* Float_t Chi2Cub() const{ */
   /*   Float_t */
   /* } */
   Float_t chi2X() const{ return m_chi2X;}
   Float_t chi2XnDoF() const{ return m_chi2X/(Float_t)m_nDoFX;}
   Float_t cHi2FullnDoF() const{ return m_chi2Full/(Float_t)m_nDoFFUll;}



   void sortbyZ(){

      std::sort(m_hits.begin(), m_hits.end(),[](PatHit val1, PatHit val2)->bool{return std::abs(val1.z())<std::abs(val2.z());});

   }
   void PrintTrack(){
      std::cout.precision(6);
      std::cout<<"--------Track Parameters   XZ   ----"<<std::endl;
      std::cout<<"ax + bx *(z-zRef) + cx*(z-zRef)^{2}(1+dRatio(z-zRef))"<<std::endl;
      std::cout<<"ax"<<setw(15)<<ax()<<std::endl;
      std::cout<<"bx"<<setw(15)<<bx()<<std::endl;
      std::cout<<"cx"<<setw(15)<<cx()<<std::endl;
      std::cout<<"dRatio"<<setw(15)<<dRatio()<<std::endl;
      std::cout<<"--------Track Parameters   YZ   ----"<<std::endl;
      std::cout<<"ay + by *(z-zRef)"<<std::endl;
      std::cout<<"ay"<<setw(15)<<ay()<<std::endl;
      std::cout<<"by"<<setw(15)<<by()<<std::endl;
      std::cout<<"y(7000)"<<setw(15)<<y(7000)<<std::endl;
      std::cout<<"y(10000)"<<setw(15)<<y(10000)<<std::endl;
      std::cout<<"i"<<setw(15)<<"Hit X"<<setw(15)<<"Track X"<<setw(15)<<"Chi2 Contrib"<<setw(15)<<"y At Z"<<setw(15)<<"coord"<<std::endl;
      for(Int_t i =0; i<m_hits.size(); i++){
         Float_t distance2oe = m_hits[i].w2()*std::pow((x(m_hits[i].z(0.))-m_hits[i].x(0.)),2);
         std::cout<<i<<setw(15)<<m_hits[i].x(0.)<<setw(15)<<x(m_hits[i].z(0.))<<setw(15)<<Chi2Hit(m_hits[i])<<setw(15)<<y(m_hits[i].z(0.))<<setw(15)<<m_hits[i].coord()<<std::endl;
      }
   }
   bool hasHitInPlane(Int_t val){
      bool returnst = false;
      for( int i =0; i<m_hits.size();i++){
         if(m_hits[i].planeCode()==val){
            returnst = true;
            break;
         }
      }
      return returnst;
   }
   int nHitsInPlane(int val) const{
     int n = 0;
     for(int i=0; i<m_hits.size();i++){
       if(m_hits[i].planeCode()==val) n++;
     }
     return n;
   }
   std::vector<PatHit> hitsinPlane( int val) const{
     std::vector<PatHit> hits;
     for(int i=0;i<m_hits.size();i++){
       if(m_hits[i].planeCode()==val){
         hits.push_back(m_hits[i]);
       }
     }
     return hits;
   }
   int worst(){
      Float_t maxChi2 = -1;
      int worsti=-10;
      //PatHit worst = PatHit();
      for(int i=0; i<m_hits.size(); i++){
         Float_t Chi2onHit = Chi2Hit(m_hits[i]);
         if(Chi2onHit>maxChi2){
            maxChi2 =Chi2onHit;
            worsti = i;
         }
      }
      return worsti;
   }
   int worstPlane(){
     Float_t maxChi2 = -1;
     int worsti = -10;
     for(int i=0; i<m_hits.size(); i++){
         Float_t Chi2onHit = Chi2Hit(m_hits[i]);
         if(Chi2onHit>maxChi2){
            maxChi2 =Chi2onHit;
            worsti = m_hits[i].planeCode();
         }
      }
     return worsti;
   }
   Float_t maxChi2Hit(){
      Float_t maxChi2 = -1;
      PatHit worst = PatHit();
      for(int i=0; i<m_hits.size(); i++){
         Float_t Chi2onHit = Chi2Hit(m_hits[i]);
         if(Chi2onHit>maxChi2){
            maxChi2 =Chi2onHit;
            worst = m_hits[i];
         }
      }
      return maxChi2;
   }
   void PrintHits()
   {
     int j = worst();
     std::cout.precision(8);
     std::cout<<" Nb. Hits on the track =   "<<m_hits.size()<<endl;
     std::cout<<"i"<<setw(15)<<"1/w"<<setw(15)<<"Xat0"<<setw(15)<<"Zat0"<<setw(15)<<"planeCode"<<setw(15)<<"zone"<<setw(15)<<"isX"<<setw(15)<<"dxDy"<<setw(15)<<"dzDy"<<setw(15)<<"Fraction"<<setw(15)<<"Size"<<setw(15)<<"Charge"<<setw(15)<<"Chi2OnHit"<<setw(15)<<"Worst"<<std::endl;
      for(Int_t i =0; i<m_hits.size(); i++){
        //m_hits[i].PrintHit();
        if(i==j){
          std::cout<<i<<setw(15)<<sqrt(1./m_hits[i].w2())<<setw(15)<<m_hits[i].x(0.)<<setw(15)<<m_hits[i].z(0.)<<setw(15)<<m_hits[i].planeCode()<<setw(15)<<m_hits[i].zone()<<setw(15)<<m_hits[i].isX()<<setw(15)<<m_hits[i].dxDy()<<setw(15)<<m_hits[i].dzDy()<<setw(15)<<m_hits[i].cluster().fraction()<<setw(15)<< m_hits[i].cluster().size()<<setw(15)<<m_hits[i].cluster().charge()<<setw(15)<<Chi2Hit(m_hits[i])<<setw(15)<<"Worst"<<std::endl;
        }else{
          std::cout<<i<<setw(15)<<sqrt(1./m_hits[i].w2())<<setw(15)<<m_hits[i].x(0.)<<setw(15)<<m_hits[i].z(0.)<<setw(15)<<m_hits[i].planeCode()<<setw(15)<<m_hits[i].zone()<<setw(15)<<m_hits[i].isX()<<setw(15)<<m_hits[i].dxDy()<<setw(15)<<m_hits[i].dzDy()<<setw(15)<<m_hits[i].cluster().fraction()<<setw(15)<< m_hits[i].cluster().size()<<setw(15)<<m_hits[i].cluster().charge()<<setw(15)<<Chi2Hit(m_hits[i])<<std::endl;
        }
      }
   }
   //friend bool operator== (PrSeedTrack cP1, PrSeedTrack cP2);

protected:
private:
   void init(){
      m_hits.reserve(32);
      m_ax=0.;
      m_ax = 0.;
      m_bx = 0.;
      m_cx = 0.;
      m_dx = -0.000262;
      m_by = 0.;
      m_ay = 0.;
      m_chi2 = 0.;
      m_nDoF = -1;
      m_dXCoord = 0.;
      m_meanDy  = 0.;

      m_chi2X=0.;
      m_nDoFX=0.;
      m_chi2Full=0.;
      m_nDoFFUll=0.;


      m_ax1 =0.;
      m_bx1 = 0.;
      m_cx1 = 0.;
      m_dx1 = 0.;
      m_chi21 =0.;
      m_nDoF1 =0.;

   }
};
typedef std::vector<PrSeedTrack> PrSeedTracks;
//bool operator==(PrSeedTrack a, PrSeedTrack b)
//{
//   if(a.hits().size() != b.hits().size()) return false;
   // a.sortbyZ();
   // b.sortbyZ();
   // bool flag = true;
   // for(int i =0;i<a.hits().size();i++){
   //    if( (a.hits()[i].x() != b.hits()[i].x()) && a.hits()[i].z()!=a.hits()[i].z()) {  flag = false; break;}
   // }
   // return flag;
//}
#endif // #ifdef FakedSeeding_cxx
