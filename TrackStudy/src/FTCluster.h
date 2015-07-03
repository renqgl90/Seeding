#ifndef FTCluster_H
#define FTCluster_H 1

class FTCluster{
  public:
    FTCluster( ):
      m_Fraction(-1),
      m_Charge(-1),
      m_Size(-1),
      m_SipmID(-1),
      m_SipmCell(-1),
      m_Module(-1),
      m_Layer(-1),
      m_Mat(-1),
      m_Quarter(-1)
  {};
  virtual ~FTCluster(){};
  //THE SETTER

  void setCluster(Float_t Fraction, Float_t ssize, Float_t Charge,Float_t SipmID, Float_t SipmCell, Float_t Module, Float_t Layer,Float_t Mat, Float_t Quarter ){
      m_Fraction = Fraction;
      m_Size = ssize;
      m_Charge = Charge;
      m_SipmID = SipmID;
      m_SipmCell = SipmCell;
      m_Module = Module;
      m_Layer = Layer;
      m_Mat = Mat;
      m_Quarter = Quarter;
    }
  //THE GETTERS
  Float_t fraction() const {return m_Fraction;}
  Float_t charge() const {return m_Charge;}
  Float_t size() const {return m_Size;}
  Float_t sipmID() const {return m_SipmID;}
  Float_t spipmCell() const {return m_SipmCell;}
  Float_t module() const {return m_Module;}
  Float_t quarter() const {return m_Quarter;}
  Float_t layer() const {return m_Layer;}

private:
  Float_t m_Fraction;
  Float_t m_Size;
  Float_t m_Mat;
  Float_t m_Charge;
  Float_t m_SipmID;
  Float_t m_SipmCell;
  Float_t m_Module;
  Float_t m_Quarter;
  Float_t m_Layer;
};
#endif
