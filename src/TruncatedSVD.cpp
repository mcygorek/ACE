#include "TruncatedSVD.hpp"
#include "QRPinv_struct.hpp"
#include <Eigen/SVD>
#include "DummyException.hpp"

namespace ACE{

template <typename T> bool TruncatedSVD_T<T>::do_compress()const{
  return (threshold>0 || maxk>0);
}
template <typename T> int TruncatedSVD_T<T>::get_truncated_dim(const Eigen::VectorXd &svals)const{
  
  if(svals.size()<1){
    std::cerr<<"TruncatedSVD::get_truncated_dim: zero singular values found!"<<std::endl;
    throw DummyException();
  }
//std::cout<<"svals: "<<svals.transpose()<<std::endl;
  int mk_dim=svals.size(); 
  if(maxk>0 && maxk<mk_dim)mk_dim=maxk;

  int dim=mk_dim;
  for(int i=1; i<dim; i++){
    if(svals(i)<threshold*svals(0)){ 
      dim=i; 
      break;
    }
  }

  if(mink>0){if(dim<mink)dim=mink; if(dim>svals.size())dim=svals.size();}
  return dim;
}
template <typename T> SelectIndices TruncatedSVD_T<T>::get_select_indices(
                       const Eigen::VectorXd & s1, const Eigen::VectorXd & s2)
                                                                        const{
  if(s1.size()<1 || s2.size()<1){
    std::cerr<<"TruncatedSVD::get_selected_indices: s1.size()<1 || s2.size()<1!"<<std::endl;
    throw DummyException();
  }

  double ss_max=s1(0)*s2(0);
  double thr=threshold; 
  int mk=maxk; 
  
  if(mk<=0 && mink<=0){
    SelectIndices k_list;
    for(int k1=0; k1<(int)s1.size(); k1++){
      for(int k2=0; k2<(int)s2.size(); k2++){
        if( s1(k1)*s2(k2) >= ss_max * thr ){
          k_list.push_back(k1,k2);
        }else{
          break;
        }
      }
    }
    return k_list;

  }else{ //maxk>0 or mink>0: have to sort products of SVs first
    std::vector<std::pair<double, std::pair<int, int> > > vec;
    for(int k1=0; k1<(int)s1.size(); k1++){
      for(int k2=0; k2<(int)s2.size(); k2++){
        if( s1(k1) * s2(k2) >= ss_max * thr ){
          vec.push_back(std::make_pair( s1(k1)*s2(k2), std::make_pair(k1,k2)));
        }else{
          break;
        }
      }
    }

    std::sort(vec.begin(), vec.end(), 
      [](const std::pair<double, std::pair<int, int> > & e1, 
         const std::pair<double, std::pair<int, int> > & e2) -> bool { 
          return e1.first>e2.first;} );

    SelectIndices k_list;
    int n_elem=vec.size(); 
    if(n_elem>mk)n_elem=mk; 
    if(mink>0){if(n_elem<mink)n_elem=mink; if(n_elem>vec.size())n_elem=vec.size();}
    for(int s=0; s<n_elem; s++){
      k_list.push_back(vec[s].second);
    }

    return k_list;
  }
}

template <typename T>
TruncatedSVD_Return_T<T> TruncatedSVD_T<T>::compress(
            const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & A,
            bool calculate_residual) const{

  Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > 
                          svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);

  int tdim=get_truncated_dim(svd.singularValues());

  TruncatedSVD_Return_T<T> ret;
  ret.U=svd.matrixU().block(0, 0, A.rows(), tdim);
  ret.sigma=svd.singularValues().head(tdim);
  ret.Vdagger=svd.matrixV().block(0, 0, A.cols(), tdim).adjoint();
  
  if(calculate_residual){
    int totSV=svd.singularValues().rows();
    int last=totSV-tdim;
    ret.Residual=svd.matrixU().block(0, tdim, A.rows(), last)
                 * svd.singularValues().tail(last)
                 * svd.matrixV().block(0, tdim, A.cols(), last).adjoint();
  }

  return ret;
}

template <typename T> PassOn_T<T> TruncatedSVD_T<T>::compress_forward(
             Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & A, 
             Eigen::VectorXd & weights) const{

  PassOn_T<T> pass_on;
  if(false){return pass_on;//TODO: use_QR here
  }else{
    TruncatedSVD_Return_T<T> ret=compress(A);
    A=ret.U;
    weights=ret.sigma;

    int dim=weights.rows();
    pass_on.P=ret.Vdagger;
    for(int i=0; i<dim; i++){
      pass_on.P.row(i)*=weights(i);
    }
    pass_on.Pinv=ret.Vdagger.adjoint();
    for(int i=0; i<dim; i++){
      double w_inv=1./weights(i);
      if(Tikhonov>0.){
        w_inv=weights(i)/(weights(i)*weights(i)+Tikhonov*weights(0)*Tikhonov*weights(0));
      }
      pass_on.Pinv.col(i)*=w_inv;
    }
    return pass_on;  
  }  
}

template <typename T> PassOn_T<T> TruncatedSVD_T<T>::compress_backward(
             Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & A, 
             Eigen::VectorXd & weights) const{

  PassOn_T<T> pass_on;
  if(use_QR){
    QRPinv_struct_T<T> QRPinv(A.transpose());
    A=QRPinv.Q.transpose();
    weights=QRPinv.weights;
    pass_on.P=QRPinv.RPinv.transpose();
    pass_on.Pinv=QRPinv.RPinv_inv.transpose();
  }else{
    TruncatedSVD_Return_T<T> ret=compress(A);
    A=ret.Vdagger;
    weights=ret.sigma;

    int dim=weights.rows();
    pass_on.P=ret.U;
    for(int i=0; i<dim; i++){
      pass_on.P.col(i)*=weights(i);
    }
    pass_on.Pinv=ret.U.adjoint();
    for(int i=0; i<dim; i++){
      double w_inv=1./weights(i);
      if(Tikhonov>0.){
        w_inv=weights(i)/(weights(i)*weights(i)+Tikhonov*weights(0)*Tikhonov*weights(0));
      }
      pass_on.Pinv.row(i)*=w_inv;
    }
  }  
  return pass_on;  
}

template <typename T> void TruncatedSVD_T<T>::scale(double threshold_ratio, double maxk_ratio){
  threshold*=threshold_ratio;
  maxk=round(maxk*maxk_ratio);
}

template <typename T> void TruncatedSVD_T<T>::print_info(std::ostream &ofs)const{
  if(use_QR){
    ofs<<"use_QR=true";
  }else{
    ofs<<"threshold="<<threshold<<" maxk="<<maxk<<" mink="<<mink; //<<" keep="<<keep;
  }
  if(Tikhonov>0.){
    ofs<<" Tikhonov="<<Tikhonov;
  }
}

template <typename T> void TruncatedSVD_T<T>::setup(Parameters &param){
  threshold = param.get_as_double("threshold",0.);
  maxk = param.get_as_int("compress_maxk",0.);
  mink = param.get_as_int("compress_mink",0.);
  keep = param.get_as_double("compress_keep",0.);
  use_QR = false; //param.get_as_bool("use_QR");
  Tikhonov = param.get_as_double("compress_Tikhonov",0.);
}

template class TruncatedSVD_Return_T<std::complex<double> >;
template class TruncatedSVD_Return_T<double>;
template class TruncatedSVD_T<std::complex<double> >;
template class TruncatedSVD_T<double>;

}//namespace
