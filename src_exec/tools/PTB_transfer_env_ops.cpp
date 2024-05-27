#include "ProcessTensorBuffer.hpp"
#include "Parameters.hpp"
#include "DummyException.hpp"
#include <iostream>
#include <fstream>
#include <Eigen/Dense>

using namespace ACE;

int main(int args, char** argv){
  Parameters param(args, argv);

  std::string read_PT=param.get_as_string_check("read_PT");
  std::string multi_PT=param.get_as_string_check("multi_PT");
  int fix_phase=param.get_as_size_t("fix_phase",1);


  ProcessTensorBuffer PT(read_PT);
  ProcessTensorBuffer PT2(multi_PT, true);

  
  if(PT.get_n_tot()!=PT2.get_n_tot()){
    std::cerr<<"PT.get_n_tot()!=PT2.get_n_tot() ("<<PT.get_n_tot()<<" vs. "<<PT2.get_n_tot()<<")!"<<std::endl;
    throw DummyException();
  }
  if(PT.get_n_tot()<1){
    std::cerr<<"PT.get_n_tot()<1!"<<std::endl;
    throw DummyException();
  }

  Eigen::MatrixXcd P=Eigen::MatrixXcd::Identity(PT.get(0).M.dim_d1, PT.get(0).M.dim_d1);
  for(int n=0; n<PT.get_n_tot(); n++){

    ProcessTensorElement &e=PT.get(n);
    const ProcessTensorElement &e2=PT2.get(n);

    if(!e.is_forwardNF()){
      std::cerr<<"First PT not in forward normal form. Please run 'PTB_sweep_forward' first!"<<std::endl;
      throw DummyException();
    }
    if(!e2.is_forwardNF()){
      std::cerr<<"Second PT not in forward normal form. Please run 'PTB_sweep_forward' first!"<<std::endl;
      throw DummyException();
    }
 
    if(e.closure.rows()!=e.M.dim_d2){
      std::cerr<<"e.closure.rows()!=e.M.dim_d2!"<<std::endl;
      throw DummyException();
    }
    if(e2.closure.rows()!=e2.M.dim_d2){
      std::cerr<<"e2.closure.rows()!=e2.M.dim_d2!"<<std::endl;
      throw DummyException();
    }
    if(P.rows()!=P.cols() || P.rows()!=e.M.dim_d1){
      std::cerr<<"P.rows()!=e.M.dim_d1!"<<std::endl;
      throw DummyException();
    }
    e.M.inner_multiply_left(P);   

    int min_dim2=e.M.dim_d2; if(e2.M.dim_d2<e.M.dim_d2)min_dim2=e2.M.dim_d2;
    int max_dim2=e.M.dim_d2; if(e2.M.dim_d2>e.M.dim_d2)max_dim2=e2.M.dim_d2;
    Eigen::VectorXcd phases=Eigen::VectorXcd::Ones(e.M.dim_d2);
    for(int r=0; r<min_dim2; r++){
      std::complex<double> max_val_e2=0;
      int max_val_e2_d1=0, max_val_e2_i=0;
      for(int d1=0; d1<e2.M.dim_d1 && d1<e.M.dim_d1; d1++){
        for(int i=0; i<e2.M.dim_i && i<e.M.dim_i; i++){
          if(abs(e2.M(i,d1,r))>abs(max_val_e2)){
            max_val_e2=e2.M(i,d1,r);
            max_val_e2_d1=d1;
            max_val_e2_i=i;
          }
        } 
      }
/*
      std::complex<double> max_val_e=0;
      int max_val_e_d1=0, max_val_e_i=0;
      for(int d1=0; d1<e.M.dim_d1; d1++){
        for(int i=0; i<e.M.dim_i; i++){
          if(abs(e.M(i,d1,r))>abs(max_val_e)){
            max_val_e=e.M(i,d1,r);
            max_val_e_d1=d1;
            max_val_e_i=i;
          }
        } 
      }
      phases(r)=exp(std::complex<double>(0, 
                      std::arg(max_val_e2)-std::arg(max_val_e) ));
*/
//      phases(r)=exp(std::complex<double>(0,  \
        std::arg(max_val_e2)-std::arg(e.M(max_val_e2_i,max_val_e2_d1,r)) ));
      phases(r)=exp(std::complex<double>(0,  \
        std::arg(e.M(max_val_e2_i,max_val_e2_d1,r)) - std::arg(max_val_e2) ));
    }
    
    P=phases.asDiagonal();
    if(n<PT.get_n_tot()-1){
      e.M.inner_multiply_right(P.adjoint());
    }
    e.closure=P*e.closure;

    int o_max=e2.env_ops.ops.size();
    e.env_ops.ops.resize(o_max);
    for(int o=0; o<o_max; o++){
      e.env_ops.ops[o]=Eigen::VectorXcd::Zero(e.M.dim_d2);
      for(int r=0; r<e.env_ops.ops[o].rows() && r<e2.env_ops.ops[o].rows() && r<phases.rows(); r++){
//        e.env_ops.ops[o](r)=phases(r)*e2.env_ops.ops[o](r);
        e.env_ops.ops[o](r)=e2.env_ops.ops[o](r);
      }
    }
  }
  //PT.calculate_closures();
/*
  std::cout<<"To forward normal form and ";
  std::cout<<"extract n="<<n_extract<<"/"<<IF.a.size()<<std::endl;
  std::cout<<"using fix_phase="<<fix_phase<<std::endl;

  std::vector<Eigen::VectorXd> forwardNF=make_forwardNF(IF, final_sweep_threshold, print_timesteps, fix_phase);

  std::vector<Eigen::VectorXd> forwardNF2=make_forwardNF(IF2, final_sweep_threshold, print_timesteps, fix_phase);


  for(int r=0; r<forwardNF[n_extract].rows(); r++){
    std::cout<<forwardNF[n_extract](r)<<" ";
  }
  std::cout<<std::endl;


  for(int r=0; r<forwardNF2[n_extract].rows(); r++){
    std::cout<<forwardNF2[n_extract](r)<<" ";
  }
  std::cout<<std::endl;

  if(outfile!=""){
    std::ofstream ofs(outfile.c_str());
    for(int r=0; r<forwardNF[n_extract].rows() || r<forwardNF2[n_extract].rows(); r++){
      if(r<forwardNF[n_extract].rows()){
        ofs<<forwardNF[n_extract](r)<<" ";
      }else{
        ofs<<"nan ";
      }
      if(r<forwardNF2[n_extract].rows()){
        ofs<<forwardNF2[n_extract](r)<<" ";
      }else{
        ofs<<"nan ";
      }
      ofs<<std::endl;
    }
  }

  if(write_PT!=""){
    if(IF.env_ops.size()<1){
      std::cerr<<"IF.env_ops.size()<1!"<<std::endl;
      exit(1);
    }
    IF2.env_ops=std::vector<std::vector<Eigen::VectorXcd> >(IF2.a.size());
    for(int n=0; n<IF.a.size(); n++){
      IF2.env_ops[n]=std::vector<Eigen::VectorXcd>(IF.env_ops[n].size());
      for(int o=0; o<IF.env_ops[n].size(); o++){
        IF2.env_ops[n][o]=Eigen::VectorXcd::Zero(IF2.a[n].dim_d2);
        for(int r=0; r<IF.env_ops[n][o].rows() && r<IF2.env_ops[n][o].rows(); r++){
          IF2.env_ops[n][o](r)=IF.env_ops[n][o](r);
        }
      }
    }
    IF2.write_binary(write_PT);
  }
*/

  return 0;
}


