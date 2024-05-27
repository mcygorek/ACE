#include "EDM_OneBodyTerm_FreePropagator.hpp"
#include "ReaderBasics.hpp"
#include "DummyException.hpp"

namespace ACE{

void EDM_OneBodyTerm_FreePropagator::update(double t, double dt){
  if(!fprop)return;
  fprop->update(t, dt);
}

EDM_State EDM_OneBodyTerm_FreePropagator::apply_filtered(const EDM_State &in, int r, const EDM_Filter &filter_){
  if(!fprop)return in;
  
  EDM_Filter filter=filter_;
  filter.new_round();

  EDM_State out; 
  const Eigen::MatrixXcd &M=fprop->M;

  //correct column dimension
  in.Ldim.check_in_range(r);
  if(M.cols()!=in.Ldim[r]){
    std::cerr<<"EDM_OneBodyTerm_FreePropagator::apply_filtered: ";
    std::cerr<<"M.cols()!=in.Ldim[r]!"<<std::endl;
    throw DummyException();
  }

  out.copy_empty(in);
  out.Ldim[r]=M.rows();

  //As usual for proper coeffs
  for(const auto &elem : in.coeffs){
    EDM_Index I=elem.first;
    int j=I[r];
    for(int i=0; i<M.rows(); i++){
      if(!filter.sparse_prop_passes(M(i,j))){continue;}
      I[r]=i;
      out.add_filtered(I, M(i,j)*elem.second, filter);
    }
  }

  //Only diagonals for candidates:
  for(const auto &elem : in.candidates){
    int j=elem.first[r];
    if(!filter.sparse_prop_passes(M(j,j))){continue;}
    out.add_filtered(elem.first, M(j,j)*elem.second, filter);
  }

  out.reduce(filter);
  return out;
}

void EDM_OneBodyTerm_FreePropagator::setup(Parameters &param, int site){
  Parameters param2; 
  param2.add_from_prefix("S"+int_to_string(site), param);
  fprop.reset(new FreePropagator(param2));
//  std::cout<<" TEST: OBT_FP::setup: dim="<<fprop->get_dim()<<std::endl;
  if(fprop->get_dim()<1){ 
    fprop.reset();
  }
}

}//namespace
