#ifndef OUTPUT_OPS_DEFINED_H
#define OUTPUT_OPS_DEFINED_H

#include <Eigen/Core>
#include <Eigen/Eigenvalues> 
#include "Operators.h"
#include "TimedepMatrix.h"
#include "HilbertSpaceRotation.h"

class Output_Ops{
public:
  std::vector<Eigen::MatrixXcd> ops;
  std::vector<int> proj;

  bool use_IP;
  Eigen::MatrixXcd H_IP;
  
  size_t size()const{return ops.size();}
  int get_dim()const{if(ops.size()>0)return ops[0].rows(); else return 0;}
  Eigen::MatrixXcd & operator[](size_t i){return ops[i];}
  const Eigen::MatrixXcd & operator[](size_t i)const {return ops[i];}
  std::vector<Eigen::MatrixXcd>::iterator begin(){return ops.begin();}
  std::vector<Eigen::MatrixXcd>::iterator insert(
                     std::vector<Eigen::MatrixXcd>::iterator position, 
                                          const Eigen::MatrixXcd& val){
    return ops.insert(position, val);
  }

  void add(const Eigen::MatrixXcd &op){
    if(op.rows()!=op.cols()){
      std::cerr<<"Output_Ops::add: op.rows()!=op.cols()!"<<std::endl;
      exit(1);
    }
    if(ops.size()>0 && ops[0].rows()!=op.rows()){
      std::cerr<<"Output_Ops::add: ops[0].rows()!=op.rows()!"<<std::endl;
      exit(1);
    }
    ops.push_back(op);
  }
  void add_proj(int i){
    if(ops.size()>0 && i>=ops[0].rows()){
      std::cerr<<"Output_Ops::add_proj: i>=ops[0].rows()!"<<std::endl;
      exit(1);
    }
    proj.push_back(i);
  }

  void rotate(const HilbertSpaceRotation &hs_rot){
    for(size_t i=0; i<ops.size();i++){
      ops[i]=hs_rot.apply(ops[i]);
    }
  }
 
  Eigen::MatrixXcd trafoIP(const Eigen::MatrixXcd &A, double t) const{
    if(use_IP){
      if(A.rows()!=H_IP.cols()){
        std::cerr<<"Output_Ops::trafoIP: A.rows()!=H_IP.cols()!"<<std::endl;
        exit(1);
      }
      Eigen::MatrixXcd B=std::complex<double>(0,t/Constants::hbar_in_meV_ps)*H_IP;
      Eigen::MatrixXcd C=B.exp();
      return C*A*(C.adjoint());
    }else{
      return A;
    }
  }

  void setup(Parameters &param){
    Operators2x2 op;
    if(param.is_specified("add_Output")){
      for(int i=0; i<param.get_nr_rows("add_Output"); i++){
        std::string str=param.get_as_single_string("add_Output", i);
std::cout<<"add_Output: "<<str<<std::endl;
        add(ReadExpression(str)); 
      }
    }else if(param.is_specified("MTLS")){
      int MTLS=param.get_as_size_t("MTLS",0);
      for(int i=0; i<MTLS; i++){
        add( ExpandSingleOp(loop_pow(2, i), op.ketbra(1,1), 
                                        loop_pow(2, MTLS-1-i)) );
      }
    }else{
      add( op.ketbra(1,1) );
      add( op.ketbra(0,0) );
      add( op.ketbra(0,1) );
    }
   
    if(param.is_specified("apply_InteractionPicture_Output")){ 
      use_IP=true;
      H_IP=param.get_as_operator("apply_InteractionPicture_Output");
    }else{
      use_IP=false;
    }
  }
  Output_Ops(Parameters &param){
    setup(param);
  }
  Output_Ops(){
    use_IP=false;
  }
};

#endif
