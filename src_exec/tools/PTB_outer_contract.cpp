#include "Parameters.hpp"
#include "Reader.hpp"
#include "ProcessTensorBuffer.hpp"
#include "DummyException.hpp"
#include <iomanip>

using namespace ACE;

int main(int args, char ** argv){
 try{
  Parameters param(args, argv);

  std::string initial_PT=param.get_as_string_check("initial_PT");
  std::string add_PT=param.get_as_string_check("add_PT");


  ProcessTensorBuffer PTB0(initial_PT);
  PTB0.read_only=true;
  std::cout<<"PTB0:"<<std::endl;
  PTB0.print_info(); std::cout<<std::endl;
  if(PTB0.get_n_tot()<1){
    std::cout<<"Process Tensor empty!"<<std::endl;
    return 0;
  }

  ProcessTensorBuffer PTB1(add_PT);
  PTB1.read_only=true;
  std::cout<<"PTB1:"<<std::endl;
  PTB1.print_info(); std::cout<<std::endl;
  if(PTB1.get_n_tot()<1){
    std::cout<<"Process Tensor empty!"<<std::endl;
    return 0;
  }

  int n_max=param.get_as_size_t("n_max",PTB0.get_n_tot());
  if(PTB1.get_n_tot()<n_max){
    std::cerr<<"Mismatch in lengths!"<<std::endl;
    throw DummyException();
  }

  int N=PTB0.get(0, ForwardPreload).get_N();
  int NL=N*N; 
  double base=param.get_as_double("base",NL);

  Eigen::MatrixXcd value(1,1); value(0,0)=1.;

  for(int n=0; n<n_max; n++){
    const ProcessTensorElement & e0=PTB0.get(n, ForwardPreload);
    const ProcessTensorElement & e1=PTB1.get(n, ForwardPreload);
    if(N!=e1.get_N()){
      std::cerr<<"Mismatch in system dimensions: "<<e0.get_N()<<" vs. "<<e1.get_N()<<"!"<<std::endl;
      throw DummyException();
    }

    Eigen::MatrixXcd tmp=Eigen::MatrixXcd::Zero(e0.M.dim_d2*NL*NL,e1.M.dim_d1);
    for(int i=0; i<NL; i++){
      for(int j=0; j<NL; j++){
        int i_ind0=e0.accessor.dict.beta[i*NL+j]; if(i_ind0<0)continue;
        for(int d1=0; d1<e0.M.dim_d1; d1++){
          for(int d2=0; d2<e0.M.dim_d2; d2++){
            for(int d=0; d<e1.M.dim_d1; d++){
   tmp((d2*NL+i)*NL+j,d)+=std::conj(e0.M(i_ind0, d1, d2))*value(d1,d)/base;
            }
          }
        }
      }
    }
    value=Eigen::MatrixXcd::Zero(e0.M.dim_d2,e1.M.dim_d2);
    for(int i=0; i<NL; i++){
      for(int j=0; j<NL; j++){
        int i_ind1=e1.accessor.dict.beta[i*NL+j]; if(i_ind1<0)continue;
        for(int d1=0; d1<e1.M.dim_d1; d1++){
          for(int d2=0; d2<e1.M.dim_d2; d2++){
            for(int d=0; d<e0.M.dim_d2; d++){
    value(d, d2)+=e1.M(i_ind1, d1, d2)*tmp((d*NL+i)*NL+j,d1);
            }
          }
        }
      }
    }
/*  //old version: works but takes much longer
    Eigen::MatrixXcd new_value=Eigen::MatrixXcd::Zero(e0.M.dim_d2,e1.M.dim_d2);
    for(int i=0; i<NL; i++){
      for(int j=0; j<NL; j++){
        int i_ind0=e0.accessor.dict.beta[i*NL+j]; if(i_ind0<0)continue;
        int i_ind1=e1.accessor.dict.beta[i*NL+j]; if(i_ind1<0)continue;
        for(int d01=0; d01<e0.M.dim_d1; d01++){
          for(int d02=0; d02<e0.M.dim_d2; d02++){
            for(int d11=0; d11<e1.M.dim_d1; d11++){
              for(int d12=0; d12<e1.M.dim_d2; d12++){
                 new_value(d02,d12)+=std::conj(e0.M(i_ind0, d01, d02))*e1.M(i_ind1, d11, d12)*value(d01,d11)/base;
              }
            }
          }
        }
      }
    }
    value=new_value;
*/
  }

  if(n_max<PTB0.get_n_tot()){
    Eigen::MatrixXcd new_value=Eigen::MatrixXcd::Zero(1,value.cols());
    Eigen::VectorXcd closure=PTB0.get(n_max-1).closure;
    for(int c=0; c<value.cols(); c++){
      for(int r=0; r<value.rows(); r++){
        new_value(0, c)=value(r,c)*closure(r);
      }
    }
    value=new_value;
  }
  if(n_max<PTB1.get_n_tot()){
    Eigen::MatrixXcd new_value=Eigen::MatrixXcd::Zero(value.rows(),1);
    Eigen::VectorXcd closure=PTB1.get(n_max-1).closure;
    for(int r=0; r<value.rows(); r++){
      for(int c=0; c<value.cols(); c++){
        new_value(r, 0)=value(r,c)*closure(c);
      }
    }
    value=new_value;
  } 
  std::cout<<base<<"^"<<n_max<<"="<<pow(base,n_max)<<" times"<<std::endl;
  std::cout<<std::setprecision(16);
  std::cout<<value(0,0).real()<<" "<<value(0,0).imag()<<std::endl;
 }catch (DummyException &e){
  return 1;
 }
  return 0;
}
