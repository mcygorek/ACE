#include "OutputPrinter.hpp"
#include "Output_Ops.hpp"
#include "Which_Env_Ops.hpp"
#include "Parameters.hpp"
#include "LiouvilleTools.hpp"
#include <iomanip>
#include "DummyException.hpp"

namespace ACE{

void OutputPrinter::clear(){
  print_timestep=false;
  ofs.reset(nullptr);
  output_Op=Output_Ops();
  full_densmat=false;
  which_env_ops=Which_Env_Ops_List();
  ofs_eigenstates.reset(nullptr);
  {std::vector<size_t> tmp; eigenstate_components.swap(tmp);}
  {std::vector<Eigen::VectorXcd> tmp; rho_t.swap(tmp);}
  {std::vector<double> tmp; rho_times.swap(tmp);}
  start_extract=0;
  do_extract=false;
}
void OutputPrinter::clear_results(){
  {std::vector<Eigen::VectorXcd> tmp(0); rho_t.swap(tmp);}
  {std::vector<double> tmp(0); rho_times.swap(tmp);}
}
void OutputPrinter::set_stream(const std::string &outfile, int precision){
  if(ofs){
    if(ofs->is_open())ofs->close();
    ofs.reset(nullptr);
  }
  if(outfile!="" && outfile!="/dev/null"){
    ofs.reset(new std::ofstream(outfile.c_str()));
    if(precision>0)*ofs<<std::setprecision(precision);
  }
}
void OutputPrinter::setup(Parameters & param, int setdim){

  full_densmat=false;
  //setup output_Op and check dimensions of input
  output_Op.setup(param);
  if(setdim>=0 && output_Op.size()>0 && output_Op.get_dim() != setdim){
    std::cerr<<"Mismatching dimensions: output_Op.get_dim()="<<output_Op.get_dim()<<" instead of "<<setdim<<"!"<<std::endl;
    throw DummyException();
  }
//  FTO.setup(param);
  which_env_ops.setup(param);
  if(setdim>=0 && which_env_ops.size()>0 && which_env_ops.get_dim() != setdim){
    std::cerr<<"Mismatching dimensions: which_env_ops.get_dim()="<<which_env_ops.get_dim()<<" instead of "<<setdim<<"!"<<std::endl;
    throw DummyException();
  }


/*
  if(output_Op.size()<1 && which_env_ops.size()<1){
    std::cerr<<"No output operator specified! Try 'add_Output OPERATOR'!"<<std::endl; 
    throw DummyException();
  }
*/
  if(output_Op.size()<1){
    full_densmat=true;
  } 

  //only if that checked out, open file:
  std::string outfile = param.get_as_string("outfile", "ACE.out");

  set_stream(outfile, param.get_as_int("set_precision",-1));

  //filename to put occupations of instantaneous eigenstates:
  std::string print_eigenstate_occupations = param.get_as_string("print_eigenstate_occupations", "");
  if(print_eigenstate_occupations!="" && print_eigenstate_occupations!="/dev/null"){
    ofs_eigenstates.reset(new std::ofstream(print_eigenstate_occupations.c_str()));
    int set_precision=param.get_as_int("set_precision",-1);
    if(set_precision>0)*ofs_eigenstates<<std::setprecision(set_precision);
  }else if(ofs_eigenstates){
    if(ofs_eigenstates->is_open())ofs_eigenstates->close();
    ofs_eigenstates.reset(nullptr);
  }

  eigenstate_components=param.get_all_size_t("print_eigenstate_component");

//  std::cout<<"eigenstate_components.size()="<<eigenstate_components.size()<<" :"; for(auto o : eigenstate_components){std::cout<<" "<<o;}std::cout<<std::endl;


//Initialize density matrix storage:
  {std::vector<Eigen::VectorXcd> tmp; rho_t.swap(tmp);}
  {std::vector<double> tmp; rho_times.swap(tmp);}
  do_extract=false;
  start_extract=0;
}

void OutputPrinter::print(int n, double t, const Eigen::VectorXcd & rho_reduced,
                       const std::vector<std::complex<double> > & env_reduced){ 


  Eigen::VectorXcd exp_vals;
  if(full_densmat){
    exp_vals=rho_reduced;
  }else{
    exp_vals=Eigen::VectorXcd::Zero(output_Op.size());

    Eigen::MatrixXcd rho=L_Vector_to_H_Matrix(rho_reduced);
    Eigen::MatrixXcd rho2=output_Op.trafoIP(rho,-t);

    for(size_t o=0; o<output_Op.size(); o++){
      for(int i=0; i<rho.rows(); i++){
        for(int j=0; j<rho.cols(); j++){
          exp_vals(o)+=output_Op[o](j,i)*rho2(i,j);
        }
      }
    }
  }

  if(do_extract && n>=start_extract){
    rho_t.push_back(exp_vals);
    rho_times.push_back(t);
  }

  if(!ofs)return;
  if(!ofs->is_open()){
    std::cerr<<"OutputPrinter::print not set up!"<<std::endl;
    throw DummyException();
  }

  *ofs<<t;
  for(int o=0; o<exp_vals.rows(); o++){
    *ofs<<" "<<exp_vals(o).real();
    *ofs<<" "<<exp_vals(o).imag();
  }

  for(size_t i=0; i<env_reduced.size(); i++){
    *ofs<<" "<<env_reduced[i].real();
    *ofs<<" "<<env_reduced[i].imag();
  }

  *ofs<<std::endl;
}

void OutputPrinter::print_eigenstate_occupations(double t, const Eigen::MatrixXcd &H, const Eigen::MatrixXcd &rho){
  if(ofs_eigenstates){
    int N=H.rows();
    if(rho.rows() != N){
      std::cerr<<"OutputPrinter::print_eigenstate_occupations: rho.rows() != N!"<<std::endl;
      throw DummyException();
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(H);

    *ofs_eigenstates<<t;
    for(int i=0; i<N; i++){
      Eigen::VectorXcd v=solver.eigenvectors().col(i);
      std::complex<double> c=v.adjoint()*rho*v;
      *ofs_eigenstates<<" "<<c.real();
    }
    for(int i=0; i<N; i++){
      double e=solver.eigenvalues()(i);
      *ofs_eigenstates<<" "<<e;
    }
    for(int o=0; o<eigenstate_components.size(); ++o){
      if(eigenstate_components[o]>=N){
        std::cerr<<"OutputPrinter::print_eigenstate_occupations: eigenstate_components["<<o<<"]="<<eigenstate_components[o]<<" >= N="<<N<<std::endl;
        throw DummyException();
      }
      for(int i=0; i<N; i++){
        std::complex<double> c=solver.eigenvectors()(eigenstate_components[o],i);
        *ofs_eigenstates<<" "<<c.real()<<" "<<c.imag();
      }
    }
    *ofs_eigenstates<<std::endl;
  }
}

void OutputPrinter::finish(){
}

std::pair<Eigen::VectorXd,Eigen::MatrixXcd> OutputPrinter::extract()const{
  std::pair<Eigen::VectorXd,Eigen::MatrixXcd> res;

  if(rho_t.size()<1){ return res; }
  if(rho_times.size()!=rho_t.size()){
    std::cerr<<"Error: rho_times.size()!=rho_t.size() in OutputPrinter::extract()!"<<std::endl;
    throw DummyException();
  }

  res.first=Eigen::VectorXd::Zero(rho_times.size());
  for(int i=0; i<res.first.rows();i++){res.first(i)=rho_times[i];}

  res.second=Eigen::MatrixXcd::Zero(rho_t.size(), rho_t[0].rows());
  for(size_t i=0; i<rho_t.size(); i++){
    if(rho_t[i].rows()!=res.second.cols()){
      std::cerr<<"Error: Inconsistent number of observables in OutputPrinter::extract()!"<<std::endl;
      throw DummyException();
    }
    for(int r=0; r<rho_t[i].rows(); r++){
      res.second(i, r)=rho_t[i](r);
    }
  }
  return res;
}

}//namespace
