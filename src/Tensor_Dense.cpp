#include "Tensor_Dense.hpp"
#include "Tensor.hpp"

namespace ACE{

  /* backwards:
  int Tensor_Dense::get_block_index(const Tensor_Index &ind) const{
    if(ind.get_rank()!=dim.size()){
      std::cerr<<"Tensor_Dense: operator(): ind.get_rank()!=dim.size()!"<<std::endl;
      exit(1);
    } 
    if(dim.size()<1)return 0.;
    int b=ind[dim.size()-1];
    for(int n=(int)dim.size()-2; n>=0; n--){
      b*=dim[n];
      b+=ind[n];
    }
    return b;
  } */

  int Tensor_Dense::get_block_index(const Tensor_Index &ind) const{
    if(ind.get_rank()!=(int)dim.size()){
      std::cerr<<"Tensor_Dense: operator(): ind.get_rank()!=dim.size()!"<<std::endl;
      exit(1);
    } 
    if(dim.size()<1)return 0.;
    int b=ind[0];
    for(size_t n=1; n<dim.size(); n++){
      b*=dim[n];
      b+=ind[n];
    }
    return b;
  }

  std::complex<double> Tensor_Dense::operator() (const Tensor_Index &ind)const{
    int b=get_block_index(ind);
    return mem[b];
  }
  //In contrast to, e.g., MPS versions: write access to elements
  std::complex<double> & Tensor_Dense::operator[](const Tensor_Index &ind){
    int b=get_block_index(ind);
    return mem[b];
  }

  void Tensor_Dense::allocate(){ mem=new std::complex<double>[get_total_size()]; }
  void Tensor_Dense::deallocate(){ if(mem!=NULL)delete[] mem; mem=NULL;}
  void Tensor_Dense::resize(const std::vector<int> &list){
    deallocate();
    dim=list;
    allocate();
  }
  void Tensor_Dense::resize(int rank_, int dim_){
    deallocate();
    dim.clear();
    dim.resize(rank_,dim_);
    allocate();
  }
  void Tensor_Dense::fill(std::complex<double> c){
    for(int i=0; i<get_total_size(); i++)mem[i]=c;
  }
  void Tensor_Dense::copy(const Tensor_Dense &other){
    resize(other.dim);
    for(int i=0; i<get_total_size(); i++)mem[i]=other.mem[i];
  }

  void Tensor_Dense::swap(Tensor_Dense &other){
    if(other.dim.size()!=dim.size()){
      std::cerr<<"Tensor_Dense::swap: other.dim.size()!=dim.size()!"<<std::endl;
      exit(1);
    }
    for(size_t i=0; i<dim.size(); i++){
      if(other.dim[i]!=dim[i]){
        std::cerr<<"Tensor_Dense::swap: other.dim[i]!=dim[i] for i="<<i<<std::endl;
        exit(1);
      }
    }
    std::complex<double> *swp=mem;
    mem=other.mem; 
    other.mem=swp;
  }


}//namespace
