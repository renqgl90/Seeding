#ifdef FTCluster_H
#define FTCluster_H 1

class FTCluster{
  public:
    FTCluster():
      m_Fraction(-1);
      m_Charge(-1);
      m_Size(-1);
      m_SipmID(-1):
      m_SipmCell(-1):
      m_Module(-1);
      m_Layer(-1);
      m_Mat(-1);
      m_Quarter(-1);
  {};
  virtual ~FTCluster(){};
  //THE SETTER
  void setCluster(Float_t Fraction, Float_t ssize, Float_t Charge,Float_t SipmID, Float_t SipmCell, Float_t Module, Float_t Layer,Float_t Mat, Float_t Quarter ){
      m_Fraction = Fraction;
      m_Size = ssize;
      m_Charge = Charge;
      m_SipmID = SipmID;
      m_SipmCell = SipmCell;
      m_Module = m_Module;
      m_Layer = Layer;
      m_Mat = Mat;
      m_Quarter = Quarter;
    }
  //THE GETTERS
  void fraction() const {return m_Fraction;}
  void charge() const {return m_Charge;}
  void size() const {return m_Size;}
  void sipmID() const {return m_SipmID;}
  void spipmCell() const {return m_SipmCell;}
  void module() const {return m_Module;}
  void quarter() const {return m_Quarter;}
  void layer() const {return m_Layer;}

private:
  Float_t m_Fraction;
  Float_t m_Charge;
  Float_t m_SipmID;
  Float_t m_SipmCell;
  Float_t m_Module;
  Float_t m_Quarter;
  Float_t m_Layer;
};
#endif
