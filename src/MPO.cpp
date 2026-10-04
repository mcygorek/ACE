#include "MPO.hpp"

namespace ACE{

template <typename T> 
  void MPO_ScalarType<T>::sweep_single_d2_d4(const TruncatedSVD &trunc, int l, int verbosity){
    if(a.size()<2)return;
    if(l<0||l>=a.size()-1){
      std::cerr<<"MPO::sweep_single_d2_d4: l<0||l>=a.size()-1!"<<std::endl;
      exit(1);
    }

    Eigen::MatrixXcd A = a[l].get_Matrix_d1d2d3_d4();
    TruncatedSVD_Return ret = trunc.compress(A);
//    std::cout<<"["<<l<<"]: singular values: "<<ret.sigma.transpose()<<std::endl;
    double keep=ret.sigma(0);
    if(trunc.keep!=0.)keep=trunc.keep;
#ifdef PRINT_KEEP
std::cout<<"keep: "<<keep<<std::endl;
#endif
    a[l].set_from_Matrix_d1d2d3_d4( ret.U*keep, a[l].dim_d2, a[l].dim_d3);
    Eigen::MatrixXcd P=((ret.sigma/keep).asDiagonal()*ret.Vdagger).transpose();
    a[l+1].multiply_d2(P);
    if(verbosity>0){std::cout<<"forward sweep of site "<<l<<": "<<A.rows()<<","<<A.cols()<<" -> "<<P.cols()<<std::endl;}
  }

template <typename T> 
  void MPO_ScalarType<T>::sweep_d2_d4(const TruncatedSVD &trunc, int verbosity){
    if(a.size()<2)return;
    for(int l=0; l<(int)a.size()-1; l++){
      sweep_single_d2_d4(trunc, l, verbosity);
    }
  }

template <typename T> 
  void MPO_ScalarType<T>::sweep_single_d4_d2(const TruncatedSVD &trunc, int l, int verbosity){
    if(a.size()<2)return;
    if(l<1 || l>=a.size()){
      std::cerr<<"MPO::sweep_single_d4_d2: l<1 || l>=a.size()!"<<std::endl;
      exit(1);
    }

    Eigen::MatrixXcd A = a[l].get_Matrix_d2_d1d3d4();
    TruncatedSVD_Return ret = trunc.compress(A);
//    std::cout<<"["<<l<<"]: singular values: "<<ret.sigma.transpose()<<std::endl;
    double keep=ret.sigma(0);
    if(trunc.keep!=0.)keep=trunc.keep;
#ifdef PRINT_KEEP
std::cout<<"keep: "<<keep<<std::endl;
#endif
    a[l].set_from_Matrix_d2_d1d3d4( ret.Vdagger*keep, a[l].dim_d3, a[l].dim_d4);
    Eigen::MatrixXcd P=(ret.U*(ret.sigma/keep).asDiagonal());
    a[l-1].multiply_d4(P);
    if(verbosity>0){std::cout<<"backward sweep of site "<<l<<": "<<A.rows()<<","<<A.cols()<<" -> "<<P.cols()<<std::endl;}
  }
  
template <typename T> 
  void MPO_ScalarType<T>::sweep_d4_d2(const TruncatedSVD &trunc, int verbosity){
    if(a.size()<2)return;
    for(int l=(int)a.size()-1; l>0; l--){
      sweep_single_d4_d2(trunc, l, verbosity);
    } 
  }

template <typename T>
MPO_ScalarType<T> MPO_ScalarType<T>::Identity_d1_d3(const std::vector<int> & dims){
  MPO_ScalarType<T> mpo;
  mpo.a.resize(dims.size());
  for(int i=0; i<(int)dims.size(); i++){
    mpo.a[i] = FourLeg_ScalarType<T>(dims[i], 1, dims[i], 1);
    mpo.a[i].set_zero();
    for(int j=0; j<dims[i]; j++){
      mpo.a[i](j, 0, j, 0)=1.;
    }
  }
  return mpo;
}
template <typename T> 
  void MPO_ScalarType<T>::print_dims(std::ostream &os)const{
    os<<"MPO.size()="<<size()<<": ";
    for(int l=0; l<a.size(); l++){
      a[l].print_dims(os); os<<"; ";
    }
    os<<std::endl;
  }

//template class MPO_ScalarType<double>;
template class MPO_ScalarType<std::complex<double>>;
}//namespace
