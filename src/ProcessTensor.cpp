#include "ProcessTensor.hpp"
#include "ProcessTensorElement.hpp"
#include "PassOn.hpp"
#include "LiouvilleTools.hpp"
#include "BinaryReader.hpp"
#include <iostream>
#include "DummyException.hpp"

namespace ACE{
 
/*
const ProcessTensorElement & ProcessTensor::operator[](size_t i)const{
  if(i>=elements.size()){
    std::cerr<<"Accessing ProcessTensor out of range: "<<i<<" vs. "<<elements.size()<<std::endl;
    exit(1);
  }
  return elements[i];
}
ProcessTensorElement & ProcessTensor::operator[](size_t i){
  if(i>=elements.size()){
    std::cerr<<"Accessing ProcessTensor out of range: "<<i<<" vs. "<<elements.size()<<std::endl;
    exit(1);
  }
  return elements[i];
}
*/
const ProcessTensorElement & ProcessTensor::operator[](size_t i)const{ 
  return elements[i]; 
}
ProcessTensorElement & ProcessTensor::operator[](size_t i){ 
  return elements[i]; 
}


int ProcessTensor::get_max_dim()const{
  int maxdim=0;
  for(const ProcessTensorElement & e : elements){
    if(e.M.dim_d2>maxdim)maxdim=e.M.dim_d2;
  }
  return maxdim;
}
 
void ProcessTensor::sweep_forward(const TruncatedSVD &trunc){
  if(size()<1)return;
  PassOn pass_on(elements[0].M.dim_d1);
  for(int n=0; n<(int)size(); n++){
    elements[n].sweep_forward(trunc, pass_on, (n==(int)size()-1));
  }
}

void ProcessTensor::sweep_backward(const TruncatedSVD &trunc){
  if(size()<1)return;
  PassOn pass_on(elements.back().M.dim_d2);
  for(int n=(int)size()-1; n>=0; n--){
    elements[n].sweep_backward(trunc, pass_on, (n==0));
  }
}

void ProcessTensor::join_halfdt(ProcessTensor &other){
  if(size()<1){
    std::cerr<<"ProcessTensor::join_halfdt: size()<1!"<<std::endl;
    exit(1);
  }
  if(other.size() != 2*size()){ // TODO: for '<': think about reduction of inner dim
    std::cerr<<"ProcessTensor::join_halfdt: other.PT.size()<2*PT.size()!"<<std::endl;
    exit(1);
  }
  for(int n=0; n<(int)size(); n++){
    elements[n].join_symmetric(other[2*n], other[2*n+1]);
  }
}

void ProcessTensor::join_and_sweep_halfdt(ProcessTensor &other,
                                          const TruncatedSVD &trunc){
  if(size()<1){
    std::cerr<<"ProcessTensor::join_halfdt: size()<1!"<<std::endl;
    exit(1);
  }
  if(other.size() != 2*size()){ // TODO: for '<': think about reduction of inner dim
    std::cerr<<"ProcessTensor::join_halfdt: other.PT.size()<2*PT.size()!"<<std::endl;
    exit(1);
  }

  PassOn pass_on(elements[0].M.dim_d1*other[0].M.dim_d1);
  for(int n=0; n<(int)size(); n++){
    elements[n].join_symmetric(other[2*n], other[2*n+1]);
    elements[n].sweep_forward(trunc, pass_on, (n==(int)size()-1));
  }
  
  sweep_backward(trunc);
}

void ProcessTensor::check_compatible(const ProcessTensor &other, int factor){
  if(size()<1){
    std::cerr<<"ProcessTensor::join_select: size()<1!"<<std::endl;
    exit(1);
  }
  if(other.size() != factor*size()){ // TODO: for '<': think about reduction of inner dim
    std::cerr<<"ProcessTensor::join_select: other.PT.size()<"<<factor<<"*PT.size()!"<<std::endl;
    exit(1);
  }
}

/*
void ProcessTensor::join_and_sweep_select(ProcessTensor &other, const TruncatedSVD &trunc, int verbosity){
  check_compatible(other);

  SelectIndices k_list;
  for(int k1=0; k1<elements[0].M.dim_d1; k1++){
    for(int k2=0; k2<other[0].M.dim_d1; k2++){
      k_list.push_back(k1,k2);
    }
  }

  PassOn pass_on1(elements[0].M.dim_d1);
  PassOn pass_on2(other[0].M.dim_d1);

  for(int n=0; n<(int)size(); n++){
    elements[n].sweep_forward(trunc, pass_on1, (n==(int)size()-1));
       other[n].sweep_forward(trunc, pass_on2, (n==(int)size()-1));
    elements[n].join_forwardNF(n, other[n], trunc, k_list);
  }
  sweep_backward(trunc);
}
*/
void ProcessTensor::join_select_and_sweep_backward(ProcessTensor &other, const TruncatedSVD &trunc){

  check_compatible(other);
  if(size()<1)return;

  SelectIndices k_list_right=elements.back().get_forwardNF_selected_indices(other.back(), trunc);
  SelectIndices k_list_left;

  PassOn pass_on(k_list_right.size());
  for(int n=(int)size()-1; n>=0; n--){
    if(n==0){
      k_list_left.set_full(elements[0].M.dim_d1, other[0].M.dim_d1);
    }else{
      k_list_left=elements[n-1].get_forwardNF_selected_indices(other[n-1], trunc);
    }
    elements[n].join_selected(n, other[n], k_list_left, k_list_right);
    elements[n].sweep_backward(trunc, pass_on, (n==0));

    k_list_right=k_list_left;
  }
}

void ProcessTensor::set_from_ModePropagator(ModePropagator &mprop, const TimeGrid &tgrid, double dict_zero){
  int N=mprop.get_N_system();
//  int NL=N*N;

  resize(tgrid.n_tot);
  if(size()<1)return;

  for(int n=0; n<tgrid.n_tot; n++){
    elements[n].set_from_ModePropagator(mprop, tgrid.get_t(n), tgrid.get_dt(n), dict_zero);
  }

  Eigen::VectorXcd bath_init = H_Matrix_to_L_Vector(mprop.get_bath_init());
  elements[0].M.inner_multiply_left(bath_init.transpose());
  
  Eigen::VectorXcd Tr = H_Matrix_to_L_Vector(Eigen::MatrixXcd::Identity(N,N));
  elements.back().closure=Eigen::VectorXcd::Ones(1);
  elements.back().env_ops.set_ill_defined();
  elements.back().M.inner_multiply_right(Tr);

  calculate_closures();
}

void ProcessTensor::set_trivial(int Nsys, const TimeGrid &tgrid){
  resize(tgrid.n_tot);
  for(int n=0; n<tgrid.n_tot; n++){
    elements[n].set_trivial(Nsys);
  }
}
void ProcessTensor::calculate_closures(){
  if(size()<1)return;
  elements.back().calculate_closure(NULL);
  for(int n=elements.size()-2; n>=0; n--){
    elements[n].calculate_closure(&elements[n+1]);
  }
}

void ProcessTensor::print_dims(std::ostream &ofs)const{
  for(const auto &e : elements){
    e.print_dims(ofs);
    ofs<<std::endl;
  }
}

void ProcessTensor::add_modes(COMBINE_MODE combine_mode,
          ModePropagatorGenerator &mpg, const TimeGrid &tgrid, 
          const TruncatedSVD & trunc, int intermediate_sweep_n, 
          double dict_zero, int verbosity){

  if(verbosity>0){
    std::cout<<"Calculating PT for Generator '"<<mpg.name()<<"'"<<std::endl;
  }

  for(int k=mpg.first(); k<mpg.get_N_modes(); k=mpg.next(k)){
    if(verbosity>0){
      std::cout<<"Mode "<<k<<"/"<<mpg.get_N_modes()<<std::endl;
    }

    ModePropagatorPtr mpp=mpg.get_ModePropagator(k);

    if(verbosity>1)std::cout<<"Maxdim: before combination: "<<get_max_dim()<<std::endl;


    if(combine_mode==mode_select){
      ProcessTensor PT2(*mpp.get(), tgrid, dict_zero);
      if(verbosity>1)std::cout<<"Maxdim: new contribution: "<<PT2.get_max_dim()<<std::endl;
      sweep_forward(trunc);
      if(verbosity>1)std::cout<<"Maxdim: after forward sweep (this): "<<get_max_dim()<<std::endl;
      PT2.sweep_forward(trunc);
      if(verbosity>1)std::cout<<"Maxdim: after forward sweep (next): "<<PT2.get_max_dim()<<std::endl;

      join_select_and_sweep_backward(PT2, trunc);

    }else{
      TimeGrid tgrid2=tgrid.construct_half_dt();
      ProcessTensor PT2(*mpp.get(), tgrid2, dict_zero);
      if(verbosity>1)std::cout<<"Maxdim: new contribution: "<<PT2.get_max_dim()<<std::endl;
      join_and_sweep_halfdt(PT2, trunc);
    }
    if(verbosity>1)std::cout<<"Maxdim: after sweep: "<<get_max_dim()<<std::endl;

    for(int sweep=0; sweep<intermediate_sweep_n; sweep++){
      if(verbosity>1)std::cout<<"intermediate sweep "<<sweep<<"/"<<intermediate_sweep_n<<std::endl;
      sweep_forward(trunc);
      if(verbosity>1)std::cout<<"Maxdim: after forward sweep: "<<get_max_dim()<<std::endl;
      sweep_backward(trunc);
      if(verbosity>1)std::cout<<"Maxdim: after backward sweep: "<<get_max_dim()<<std::endl;
    }
  }
}

std::pair<int, bool> ProcessTensor::read_header(std::ifstream &is, const std::string & context){
  if(!is.good()){
    if(context!=""){
      std::cerr<<"Cannot open PT file '"<<context<<"'!"<<std::endl;
    }else{
      std::cerr<<"ProcessTensor::read_header: !is.good()"<<std::endl;
    }
    throw DummyException();
  }
  std::string magic=binary_read_fixedSizeString(is, 4);
  bool reverse=false;
  if(magic=="PT__"){
  }else if(magic=="PTr_"){
    reverse=true;
  }else{
    std::cerr<<"Cannot read PT file: signature '"<<magic<<"' instead of 'PT__'!"<<std::endl;
    throw DummyException();
  }

  int version=binary_read_int(is);
  if(version!=3){
    std::cerr<<"Cannot read PT file: version "<<version<<" instead of 3!"<<std::endl; 
    throw DummyException();
  }
  int sz=binary_read_int(is);
  std::pair<int, bool> ret(sz, reverse);
  return ret;
}

void ProcessTensor::write_header(std::ofstream &os, int sz, bool reverse){
  if(reverse){
    binary_write_fixedSizeString(os, 4, "PTr_");
  }else{
    binary_write_fixedSizeString(os, 4, "PT__");
  }
  int version=3;
  binary_write_int(os,version);

  binary_write_int(os,sz);
}

int ProcessTensor::read_size(const std::string &fname){
  std::ifstream ifs(fname.c_str());
  std::pair<int, bool> sz_r = read_header(ifs, fname);
  return sz_r.first;
}

void ProcessTensor::read_binary(std::ifstream &is, const std::string &context){
  std::pair<int, bool> header = read_header(is, context);
  int sz=header.first; 
  bool reverse=header.second;
//std::cout<<"Reading PT of size: "<<sz<<std::endl;
  elements.resize(sz);
  if(reverse){
    for(int n=(int)elements.size()-1; n>=0; n--){
      elements[n].read_binary(is);
    }
  }else{
    for(int n=0; n<(int)elements.size(); n++){
      elements[n].read_binary(is);
    }
  }
}
void ProcessTensor::read_binary(const std::string &fname){
  std::ifstream ifs(fname.c_str());
  read_binary(ifs, fname);
}

void ProcessTensor::write_binary(std::ofstream &os, bool reverse)const{
  write_header(os, elements.size(), reverse);
  if(reverse){
    for(int n=(int)elements.size()-1; n>=0; n--){
      elements[n].write_binary(os);
    }
  }else{
    for(int n=0; n<(int)elements.size(); n++){
      elements[n].write_binary(os);
    }
  }
}
void ProcessTensor::write_binary(const std::string &fname, bool reverse)const{
  std::ofstream ofs(fname.c_str());
  write_binary(ofs, reverse);
}

}//namespace
