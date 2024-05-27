#include "EDM_TwoBodyTerm_PT.hpp"
#include "ProcessTensorForwardList.hpp"
#include "ReaderBasics.hpp"
#include "DummyException.hpp"
#include <tuple>

namespace ACE{

void EDM_TwoBodyTerm_PT::update(int n, const TimeGrid &tgrid){
  if(!PT)return;
  PT->load(n);
}

EDM_State EDM_TwoBodyTerm_PT::apply_filtered(const EDM_State &in, const std::pair<int,int> &r, const EDM_Filter &filter_){
  if(!PT)return in;

  EDM_Filter filter=filter_;
  filter.new_round();

  in.Ldim.check_in_range(r.first);
  in.Ldim.check_in_range(r.second);

  const ProcessTensorElement *e=PT->current();
  const IF_OD_Dictionary &dict=e->accessor.dict;
  int NL=dict.get_NL();

  if(NL!=in.Ldim[r.first]){
    std::cerr<<"EDM_TwoBodyTerm_PT::apply_filtered: ";
    std::cerr<<"NL="<<NL<<"!=in.Ldim["<<r.first<<"]="<<in.Ldim[r.first]<<"!"<<std::endl;
    throw DummyException();
  }
  if(e->M.dim_d1!=in.Ldim[r.second]){
    std::cerr<<"EDM_TwoBodyTerm_PT::apply_filtered: ";
    std::cerr<<"e->M.dim_d1!=in.Ldim[r.second]!"<<std::endl;
    throw DummyException();
  }
  EDM_State out; out.copy_empty(in); 
  update_closure(out, r.second);

  for(const auto &elem : in.coeffs){
    EDM_Index I=elem.first;
    int j=I[r.first];
    int d1=I[r.second];
    for(int i=0; i<NL; i++){
      int i_ind=dict.beta[i*NL+j]; if(i_ind<0)continue;
      for(int d2=0; d2<e->M.dim_d2; d2++){
        std::complex<double> c=e->M(i_ind, d1, d2);
        if(!filter.sparse_prop_passes(c)){continue;}
        I[r.first]=i;
        I[r.second]=d2;
        out.add_filtered(I, c*elem.second, filter);
      }
    }
  }
  //Unclear how to keep only "diagonals" for environment indices as they change
  //Only account for diagonal terms in system:
  for(const auto &elem : in.candidates){
    EDM_Index I=elem.first;
    int j=I[r.first];
    int d1=I[r.second];
//    for(int i=0; i<NL; i++)
    int i=j;
    {
      int i_ind=dict.beta[i*NL+j]; if(i_ind<0)continue;
      for(int d2=0; d2<e->M.dim_d2; d2++){
        std::complex<double> c=e->M(i_ind, d1, d2);
        if(!filter.sparse_prop_passes(c)){continue;}
        I[r.first]=i;
        I[r.second]=d2;
        out.add_filtered(I, c*elem.second, filter);
      }
    }
  }
  out.reduce(filter);
  return out;
}

void EDM_TwoBodyTerm_PT::update_closure(EDM_State & state, int r){
  if(!PT)return;

  state.Ldim.check_in_range(r);
  if(r>=state.closures.size()){
    std::cerr<<"EDM_TwoBodyTerm_PT::update_closure: r>=closures.size()!"<<std::endl;  
    throw DummyException();
  }
  const ProcessTensorElement *e=PT->current();
  state.closures[r]=e->closure;
  state.Ldim[r]=state.closures[r].rows();
  
  if(e->is_forwardNF()){
    if(state.relevance.size()<=r){
      state.relevance.resize(state.Ldim.size());
    }
    state.relevance[r]=e->forwardNF;
  }
}

void EDM_TwoBodyTerm_PT::setup(Parameters &param, const std::pair<int,int> &site){
  std::string key="S"+int_to_string(site.first)+"S"+int_to_string(site.second)+"_read_PT";
  PT.reset();
  std::string fname=param.get_as_string(key);
  if(fname!=""){
    std::cout<<"EDM_TwoBodyTerm_PT::setup reading PT file '"<<fname<<"'"<<std::endl;
    PT=ProcessTensorForwardList::PTptr_from_file(fname,true);
  }
}

}//namespace
