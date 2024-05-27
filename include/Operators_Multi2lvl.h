#ifndef OPERATORS_MULTI2LVL_DEFINED_H
#define OPERATORS_MULTI2LVL_DEFINED_H

#include <Eigen/Core>
#include <vector>

namespace ACE{
template <int D> class Operators_Multi2lvl{
public:
  int fulldim(){
    int f=1;
    for(int i=0; i<D; i++)f*=D;
    return f;
  }
  Eigen::MatrixXcd set_by_occupations(const std::vector<double> &n){
    if(n.size()!=D){
      std::cerr<<"Operators_Multi2lvl::set_by_occupations: n.size()!=D"<<std::endl; 
      exit(1);
    }
    int fdim=fulldim();
    Eigen::MatrixXcd ret=Eigen::MatrixXcd::Zero(fdim, fdim);
    for(int i=0; i<fdim; i++){
      int icpy=i;
      double res=1;
      for(int s=0; s<n.size(); s++){
        if(icpy%2 == 0){ 
          res*=(1.-n[n.size()-1-s]);
        }else{
          res*=n[n.size()-1-s];
        }
        icpy/=2;
      }
      ret(i,i)=res;
    }
 
    return ret;
  }
/*
  Eigen::MatrixXcd from_single_op(int n, const Eigen::MatrixXcd &mat){
    if(n<0||n>=D){
      std::cerr<<"Operators_Multi2lvl::from_single_op: n<0||n>=D"<<std::endl;
      exit(1);
    }
    if(mat.rows()!=mat.cols() || mat.rows()!=2){
      std::cerr<<"Operators_Multi2lvl::from_single_op: mat.rows()!=mat.cols() || mat.rows()!=2"<<std::endl;
      exit(1);
    }    

    int fdim=fulldim();
    Eigen::MatrixXcd ret=Eigen::MatrixXcd::Zero(fdim, fdim);

    int offset=1.;
    for(int s=0; s<n; s++){


    for(int i=0; i<fdim; i++){
      double res=1;
      
      ret(i,i)=res;
    }


    return ret;
  }
*/
};

}//namespace
#endif
