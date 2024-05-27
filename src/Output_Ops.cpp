#include "PCH.hpp"
#include "Output_Ops.hpp"
#include <Eigen/Core>
#include <Eigen/Eigenvalues> 
#include "Operators.hpp"
#include "TimedepMatrix.hpp"
#include "HilbertSpaceRotation.hpp"
#include "Reader.hpp"
#include "Parameters.hpp"
#include "Constants.hpp"

namespace ACE{

  void Output_Ops::add(const Eigen::MatrixXcd &op){
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
  void Output_Ops::add_proj(int i){
    if(ops.size()>0 && i>=ops[0].rows()){
      std::cerr<<"Output_Ops::add_proj: i>=ops[0].rows()!"<<std::endl;
      exit(1);
    }
    proj.push_back(i);
  }

  void Output_Ops::rotate(const HilbertSpaceRotation &hs_rot){
    for(size_t i=0; i<ops.size();i++){
      ops[i]=hs_rot.apply(ops[i]);
    }
  }
 
  Eigen::MatrixXcd Output_Ops::trafoIP(const Eigen::MatrixXcd &A, double t) const{
    if(use_IP){
      if(A.rows()!=H_IP.cols()){
        std::cerr<<"Output_Ops::trafoIP: A.rows()!=H_IP.cols()!"<<std::endl;
        exit(1);
      }
      Eigen::MatrixXcd B=std::complex<double>(0,t/hbar_in_meV_ps)*H_IP;
      Eigen::MatrixXcd C=B.exp();
      return C*A*(C.adjoint());
    }else{
      return A;
    }
  }

  void Output_Ops::setup(Parameters &param, bool set_default_output){
    if(param.is_specified("add_Output")){
      for(int i=0; i<param.get_nr_rows("add_Output"); i++){
        std::string str=param.get_as_single_string("add_Output", i);
//std::cout<<"add_Output: "<<str<<std::endl;
        add(ReadExpression(str)); 
      }
    }else{
      if(set_default_output){
        Operators op(2);
        add( op.ketbra(1,1) );
        add( op.ketbra(0,0) );
        add( op.ketbra(0,1) );
      }
    }
   
    if(param.is_specified("apply_InteractionPicture_Output")){ 
      use_IP=true;
      H_IP=param.get_as_operator("apply_InteractionPicture_Output");
    }else{
      use_IP=false;
    }
  }

}//namespace
