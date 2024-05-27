#include "EDM_TwoBodyTerm_FreePropagator.hpp"
#include "ReaderBasics.hpp"
#include "DummyException.hpp"
#include <tuple>

namespace ACE{

void EDM_TwoBodyTerm_FreePropagator::update(int n, const TimeGrid &tgrid){
  if(!fprop)return;
  fprop->update(tgrid.get_t(n), tgrid.get_dt(n));
}

EDM_State EDM_TwoBodyTerm_FreePropagator::apply_filtered(const EDM_State &in, const std::pair<int,int> &r, const EDM_Filter &filter_){
  if(!fprop)return in;

  EDM_Filter filter=filter_;
  filter.new_round();

  EDM_State out; 
  const Eigen::MatrixXcd &M=fprop->M;

  //correct column dimension
  in.Ldim.check_in_range(r.first);
  in.Ldim.check_in_range(r.second);
  int D0=in.Ldim[r.first];
  int D1=in.Ldim[r.second];

  if(M.cols()!=D0*D1){
    std::cerr<<"EDM_TwoBodyTerm_FreePropagator::apply_filtered: ";
    std::cerr<<"M.cols()="<<M.cols()<<"!=D0*D1="<<D0*D1<<"!"<<std::endl;
    throw DummyException();
  }
  if(M.rows()!=M.cols()){
    std::cerr<<"EDM_TwoBodyTerm_FreePropagator::apply_filtered: ";
    std::cerr<<"M.cols()!=M.rows()!"<<std::endl;
    throw DummyException();
  }  
  int d0=sqrt(D0); 
  if(d0*d0!=D0){
    std::cerr<<"EDM_TwoBodyTerm_FreePropagator::apply_filtered: ";
    std::cerr<<"d0*d0!=D0!"<<std::endl;
    throw DummyException();
  }
  int d1=sqrt(D1);
  if(d1*d1!=D1){
    std::cerr<<"EDM_TwoBodyTerm_FreePropagator::apply_filtered: ";
    std::cerr<<"d1*d1!=D1!"<<std::endl;
    throw DummyException();
  }
  out.copy_empty(in);

  //Now: Flip indices (nu0 x nu1 x mu0 x mu1 ) -> (nu0 x mu0 x nu1 x mu1 )
  //     on both sides of matrix. 
  //     While we're at it, we create a list accounting for sparsity which
  //     maps [r0][r1] -> nonzero sets [r0'][r1'][val]
 
  typedef std::tuple<int, int, std::complex<double> > tmp_tuple;
  std::vector<std::vector<std::vector< tmp_tuple > > > ML 
                    (D0, std::vector<std::vector<tmp_tuple> >(D1)); 

  for(int nu0_=0; nu0_<d0; nu0_++){
    for(int nu1_=0; nu1_<d1; nu1_++){
      int old_nu_=nu0_*d1+nu1_;
      for(int mu0_=0; mu0_<d0; mu0_++){
        int new_j0=nu0_*d0+mu0_;
        for(int mu1_=0; mu1_<d1; mu1_++){
          int old_mu_=mu0_*d1+mu1_;
          int new_j1=nu1_*d1+mu1_;
          for(int nu0=0; nu0<d0; nu0++){
            for(int nu1=0; nu1<d1; nu1++){
              int old_nu=nu0*d1+nu1;
              for(int mu0=0; mu0<d0; mu0++){
                int new_i0=nu0*d0+mu0;
                for(int mu1=0; mu1<d1; mu1++){
                  int old_mu=mu0*d1+mu1;
                  int new_i1=nu1*d1+mu1;
  std::complex<double> c=M(old_nu*d0*d1+old_mu, old_nu_*d0*d1+old_mu_);
  if(!filter.sparse_prop_passes(c)){continue;}
  ML[new_j0][new_j1].push_back(std::make_tuple(new_i0, new_i1, c));
                }
              }
            }
          }
        }
      }
    }
  }
    

  //As usual for proper coeffs
  for(const auto &elem : in.coeffs){
    EDM_Index I=elem.first;
    int j0=I[r.first];
    int j1=I[r.second];
    for(tmp_tuple & t : ML[j0][j1]){
      I[r.first]=std::get<0>(t);
      I[r.second]=std::get<1>(t);
      out.add_filtered(I, std::get<2>(t)*elem.second, filter);
    }
  }

  //Only diagonals for candidates:
  for(const auto &elem : in.candidates){
    out.add_filtered(elem.first, elem.second, filter);
  }

  out.reduce(filter);
  return out;
}

void EDM_TwoBodyTerm_FreePropagator::setup(Parameters &param, const std::pair<int,int> &site){
  Parameters param2; 
  param2.add_from_prefix("S"+int_to_string(site.first)+"S"+int_to_string(site.second), param);
  fprop.reset(new FreePropagator(param2));
//  std::cout<<" TEST: TBT_FP::setup: dim="<<fprop->get_dim()<<std::endl;
  if(fprop->get_dim()<1){ 
    fprop.reset();
  }
}

}//namespace
