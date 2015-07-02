void dumpMap(TH2D * h){
  const Int_t nx = h->GetNbinsX();
  const Int_t ny = h->GetNbinsY();

  for( Int_t ix=1;ix<=nx;ix++){
    for( Int_t iy=1;iy<=ny;iy++){
      //printf("")
      std::cout<<"Bin\t ix"<<ix<<"\t iy"<<iy<<"\t Value"<<(double)h->GetBinContent(ix,iy)<<std::endl;

    }
  }
}
