#include "RankCompressor_SVD.hpp"
#include "Parameters.hpp"
#include <fstream>

namespace ACE{

  //choose how many singular values should be kept:  
template <typename T>
int RankCompressor_SVD_ScalarType<T>::get_new_dim(const Eigen::VectorXd &sv){

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
   
  template <typename T> //
  void RankCompressor_SVD_ScalarType<T>::compress_template(
                         Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &A, 
                         Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &L, 
                         Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &R, 
                         bool low_to_high){

    Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > 
       svd( A, Eigen::ComputeThinU | Eigen::ComputeThinV);

    int newdim=get_new_dim(svd.singularValues());


#ifdef PRINT_SVD_DIMS
    std::cout<<"SVD dims: "<<A.rows()<<", "<<A.cols()<<" -> "<<newdim<<std::endl;
#endif

#ifdef PRINT_SVD
    std::cout<<"SVDs("<<newdim<<"):";
    for(int i=0; i<svd.singularValues().size(); i++){
      std::cout<<" "<<svd.singularValues()(i);
    }std::cout<<std::endl;
#elif defined(PRINT_SVD_FIRST)
   std::cout<<"SVDs("<<newdim<<"):";
   std::cout<<" first: "<<svd.singularValues()(0)<<std::endl;
#endif


    if(DUMP_SVD>=0 && count==DUMP_SVD){
//Note: Total number of compressions: (2*nmax-3)*N -> middle of chain for last step, backward direction: count ~ (2*nmax-3)*N-nmax/2-1  or (2*nmax-3)*(N-1/4)
      std::cout<<"Dumping SVD at "<<count<<"-th compression to file 'DUMP_SVD.dat'!"<<std::endl;
      std::ofstream ofs("DUMP_SVD.dat");
      for(int i=0; i<svd.singularValues().size(); i++){
        ofs<<svd.singularValues()(i)<<std::endl;
      }
    }

#ifdef PRINT_SVD_MATRIX_DIFF
std::cout<<"A-U*s*V^+:"<<std::endl;
std::cout<<A - svd.matrixU().block(0,0,A.rows(), newdim)*
               svd.singularValues().head(newdim).asDiagonal()*
               svd.matrixV().block(0,0,A.cols(),newdim).adjoint()<<std::endl;
#endif




#ifdef PRINT_SVD_ORTHO
{
std::ofstream ofs("PRINT_SVD_ORTHO", std::ios_base::app);
std::time_t print_svd_ortho_time=std::time(NULL);
double diffU=max_diff_from_ortho(svd.matrixU());
double diffV=max_diff_from_ortho(svd.matrixV());
ofs<<diffU<<" "<<diffV<<" "<<std::asctime(std::localtime(&print_svd_ortho_time))<<std::endl;
/*
if(diffU>1e-12){
std::cout<<"SVD produces non-orthogonal U: "<<diffU<<std::endl;
print_diff_from_ortho(svd.matrixU());
}
if(diffV>1e-12){
std::cout<<"SVD produces non-orthogonal V: "<<diffV<<std::endl;
print_diff_from_ortho(svd.matrixV());
}
*/
if(diffU>1e-12||diffV>1e-12){
std::cout<<"newdim="<<newdim<<std::endl;
std::cout<<"A-U*s*V^+:"<<std::endl;
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> diffmat =  
                           A - svd.matrixU().block(0,0,A.rows(), newdim)*
                           svd.singularValues().head(newdim).asDiagonal()*
                           svd.matrixV().block(0,0,A.cols(),newdim).adjoint();
{std::ofstream ofs2("PROBLEMATIC_MATRIX.dat"); ofs2<<A;}
{std::ofstream ofs2("PROBLEMATIC_SVs.dat"); ofs2<<svd.singularValues();}
double maxelem=0;
int maxi=0, maxj=0;
for(int i=0; i<diffmat.rows(); i++){
  for(int j=0; j<diffmat.cols(); j++){
    if(abs(diffmat(i,j))>maxelem){
      maxelem=abs(diffmat(i,j));
      maxi=i; maxj=j;
    }
  }
}
std::cout<<maxi<<" "<<maxj<<" "<<diffmat(maxi,maxj)<<std::endl;
exit(1);
}
}
#endif

    if(reorthogonalize>0){
      if(low_to_high){
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> U =
                                  svd.matrixU().block(0,0,A.rows(),newdim);

        for(int run=0; run<reorthogonalize; run++){
          for(int i=0; i<U.cols(); i++){
            for(int j=0; j<i; j++){
              U.col(i)-=(U.col(j).dot(U.col(i)))*U.col(j);
            }
            U.col(i).normalize();
          }
        }
        L = U;
        R = U.adjoint() * A;
      }else{
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> V =
                                svd.matrixV().block(0,0,A.cols(),newdim);

        for(int run=0; run<reorthogonalize; run++){
          for(int i=0; i<V.cols(); i++){
            for(int j=0; j<i; j++){
              V.col(i)-=(V.col(j).dot(V.col(i)))*V.col(j);
            }
            V.col(i).normalize();
          }
        }
        R = V.adjoint();
        L = A * V;
      }
    }else{
      if(low_to_high){
        L = svd.matrixU().block(0,0,A.rows(), newdim);
        R = svd.singularValues().head(newdim).asDiagonal() 
            * svd.matrixV().block(0,0,A.cols(),newdim).adjoint();
      }else{
        L = svd.matrixU().block(0,0,A.rows(), newdim)
            *svd.singularValues().head(newdim).asDiagonal() ;
        R = svd.matrixV().block(0,0,A.cols(),newdim).adjoint();
      }
    }
    count++;
  }

  template<typename T>
  void RankCompressor_SVD_ScalarType<T>::compress(
                        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &A, 
                        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &L, 
                        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &R, 
                        bool low_to_high){

//    if(use_BDCSVD){
//      compress_template<Eigen::BDCSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > >(A, L, R, low_to_high);
//    }else{
//      compress_template<Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > >(A, L, R, low_to_high);
//    }
      compress_template(A, L, R, low_to_high);
  }
 

  template<typename T>
  void RankCompressor_SVD_ScalarType<T>::setup(Parameters &param){
    threshold=param.get_as_double("threshold", 0);
    sum_threshold=param.get_as_double("sum_threshold", 0);
    maxk=param.get_as_size_t("compress_maxk", 0);
    reorthogonalize=param.get_as_int("reorthogonalize",0);
    DUMP_SVD=param.get_as_int("compress_dump",-1);
    count=param.get_as_size_t("compress_dump_initial_count",0);
//    use_BDCSVD=param.get_as_bool("use_BDCSVD",false);

    if(param.is_specified("threshold_range")){
      threshold=param.get_as_double("threshold_range",0.,0,0);
      threshold_to=param.get_as_double("threshold_range",0.,0,1);
      Nrange=param.get_as_int("threshold_range",0.,0,2);
    }else{
      Nrange=0;
    }
  }

template class RankCompressor_SVD_ScalarType<std::complex<double> >;
template class RankCompressor_SVD_ScalarType<double>;


}//namespace
