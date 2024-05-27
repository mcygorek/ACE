#include "RankCompressor_ApproxSVD.hpp"
#include "ApproxSVD.hpp"
#include "Parameters.hpp"
#include <fstream>

namespace ACE{

  // check if parameters are set up so that any compression can happen
template <typename T>
  bool RankCompressor_ApproxSVD_ScalarType<T>::has_effect()const{
    if(threshold>0. || sum_threshold>0. || maxk>0)return true;
    else return false;
  }

  //choose how many singular values should be kept:  
template <typename T>
  int RankCompressor_ApproxSVD_ScalarType<T>::get_new_dim(const Eigen::VectorXd &sv){

    int newdim=sv.size();
    if(sv.size()<1){
      std::cerr<<"Compress (SVD): svd.singularValues().size()<1!"<<std::endl;
      exit(1);
    }

    // truncate based on given maximal inner dimension
    if(maxk>0 && maxk<newdim)newdim=maxk;

    // truncate based on (relative) magnitude of singular values:
    double this_threshold=threshold;
    if(Nrange>0){
      if(count<=0){
        this_threshold=threshold;
      }else if(count<=Nrange){
        this_threshold=threshold_to;
      }else{
        double x=(double)count/(double)Nrange;
        this_threshold=exp(log(threshold)*(1.-x)+log(threshold_to)*x);
std::cout<<"this threshold: "<<this_threshold<<std::endl;
      }
    }
    
    double cutoff=sv(0)*this_threshold;
    for(int i=1; i<newdim; i++){
      if(sv(i)<cutoff){newdim=i; break;}
    }

    // truncate based on (relative) sum of neglected SVs:
    if(sum_threshold>0){ 
      if(sum_threshold>1){
        std::cerr<<"Compress (SVD): sum_threshold>1!"<<std::endl;
        exit(1);
      }
      double sum=0.;
      for(int i=sv.size()-1; i>=0; i--){
        sum+=sv(i);
      }
      double cutoff=sum*sum_threshold;
      double newsum=0.;
      for(int i=newdim-1; i>=0; i--){
        newsum+=sv(i);
        if(newsum>=cutoff){
          newdim=i+1;
          break;
        }
      }
    }
    return newdim;
  }
   
template <typename T>
  void RankCompressor_ApproxSVD_ScalarType<T>::compress(
                         Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &A, 
                         Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &L, 
                         Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &R, 
                         bool low_to_high){

    ApproxSVD<T> svd;
    if( !forceQR && precondition_repeat && 
       A.rows()==svd_bck.L.rows() && A.cols()==svd_bck.R.cols() ) { 
      //&&  A.rows() == A.cols() ){
//      std::cout<<"Fitting dimensions (at count: "<<count<<")"<<std::endl;

      std::cout<<"precondition_repeat (at count: "<<count<<")"<<std::endl;

//      svd.initialize(Lr*Lr.adjoint()*A*Rr.adjoint()*Rr, forceQR);

      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Rr, Lr;
      Rr=svd_bck.reconstructR();
      Lr=svd_bck.reconstructL();

      svd.initialize(Lr.adjoint()*A*Rr.adjoint(), forceQR);
      svd.L=Lr*svd.L;
      if(A.rows()>A.cols())svd.QRremainder=Lr*svd.QRremainder;
      svd.R=svd.R*Rr;
      if(A.rows()<A.cols())svd.QRremainder=svd.QRremainder*Rr;

/*
      svd.initialize(A, forceQR);
      svd.M=svd_bck.L.adjoint()*svd.M*svd_bck.R.adjoint();   
      svd.L=svd_bck.L*svd.L;
      svd.R=svd.R*svd_bck.R;
*/
/*
      svd.M=svd_bck.L.adjoint()*A*svd_bck.R.adjoint();
      svd.L=svd_bck.L;
      svd.R=svd_bck.R;
*/

//      svd.initialize(Lr*Lr.adjoint()*A*Rr.adjoint()*Rr, forceQR);
    }else{
      svd.initialize(A, forceQR);
    }
    svd.calculate(eps2);


    Eigen::VectorXd singularValues=svd.singularValues();
  
    int newdim=get_new_dim(singularValues);


#ifdef PRINT_SVD_DIMS
    std::cout<<"SVD dims: "<<A.rows()<<", "<<A.cols()<<" -> "<<newdim<<std::endl;
#endif

#ifdef PRINT_SVD
    std::cout<<"SVDs("<<newdim<<"):";
    for(int i=0; i<singularValues.size(); i++){
      std::cout<<" "<<singularValues(i);
    }std::cout<<std::endl;
#elif defined(PRINT_SVD_FIRST)
   std::cout<<"SVDs("<<newdim<<"):";
   std::cout<<" first: "<<singularValues(0)<<std::endl;
#endif


    if(DUMP_SVD>=0 && count==DUMP_SVD){
//Note: Total number of compressions: (2*nmax-3)*N -> middle of chain for last step, backward direction: count ~ (2*nmax-3)*N-nmax/2-1  or (2*nmax-3)*(N-1/4)
      std::cout<<"Dumping SVD at "<<count<<"-th compression to file 'DUMP_SVD.dat'!"<<std::endl;
      std::ofstream ofs("DUMP_SVD.dat");
      for(int i=0; i<singularValues.size(); i++){
        ofs<<singularValues(i)<<std::endl;
      }
    }

    if(low_to_high){
        L = svd.L.block(0,0,A.rows(),newdim);
        R = svd.M.block(0,0,newdim,newdim)*svd.R.block(0,0,newdim,A.cols());
            
    }else{
        L = svd.L.block(0,0,A.rows(),newdim)*svd.M.block(0,0,newdim,newdim); 
        R = svd.R.block(0,0,newdim,A.cols());
    }

    svd_bck.swap(svd);
    count++;
  }

template <typename T>
  void RankCompressor_ApproxSVD_ScalarType<T>::setup(Parameters &param){
    threshold=param.get_as_double("threshold", 0);
    eps2=param.get_as_double("approxSVD", threshold > 0 ? threshold : 1e-6);
    sum_threshold=param.get_as_double("sum_threshold", 0);
    maxk=param.get_as_size_t("compress_maxk", 0);
    DUMP_SVD=param.get_as_int("compress_dump",-1);
    count=param.get_as_size_t("compress_dump_initial_count",0);
    precondition_repeat=param.get_as_bool("precondition_repeat",true);
    forceQR=param.get_as_bool("forceQR",false);

    if(param.is_specified("threshold_range")){
      threshold=param.get_as_double("threshold_range",0.,0,0);
      threshold_to=param.get_as_double("threshold_range",0.,0,1);
      Nrange=param.get_as_int("threshold_range",0.,0,2);
    }else{
      Nrange=0;
    }
  }

template class RankCompressor_ApproxSVD_ScalarType<std::complex<double> >;
template class RankCompressor_ApproxSVD_ScalarType<double>;

}//namespace
