#include "EDM_State.hpp"
#include "DummyException.hpp"

namespace ACE {

void EDM_State::add_to(const std::pair<EDM_Index, std::complex<double> > & P){
    auto it=coeffs.find(P.first);
    if(it==coeffs.end()){
      coeffs.insert(P);
    }else{
      it->second+=P.second;
    }
}
double EDM_State::get_relevance(const EDM_Index &I)const{
  double rel=1.;
  for(int r=0; r<I.list.size() && r<relevance.size(); r++){
    if(relevance[r].rows()>I[r]){ rel*=relevance[r](I[r]); }
  }
  return rel;
}

void EDM_State::reduce(const EDM_Filter &filter){
  EDM_State tmp;
  tmp.copy_empty(*this);
  
  for(const auto &elem : candidates){
    add_to(elem);
  }
  double largest=0;
  for(const auto &elem: coeffs){
    double val=get_relevance(elem.first)*std::abs(elem.second);
    if(val>largest){ largest=val; }
  }
  if(largest<=0.){ return; } //nothing to do
  double inv_largest=1./largest;


  // if filtering by target number of coefficients is sought:
  if(filter.max_coeffs>0){
    //literally fixing number of coefficients requires sorting. We just determine a threshold which approximately achieves max_coeffs
    std::vector<double> thrs(filter.mc_Nsteps);
    for(int i=0; i<filter.mc_Nsteps; i++){
      thrs[i]=filter.mc_max*exp((double)i/(filter.mc_Nsteps-1.) *
                          (log(filter.mc_min)-log(filter.mc_max)));
    }
//    std::vector<double> thrs={1e-2, 3e-2, 6e-2, 1e-3, 3e-3, 6e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10};
    std::vector<int> nelem(thrs.size(),0);
    //make list
    for(const auto &elem : coeffs){ 
      for(int i=0; i<thrs.size(); i++){
        if(get_relevance(elem.first)*std::abs(elem.second)*inv_largest>thrs[i]){ nelem[i]++; }
      }
    }
//    for(int i=0; i<Nsteps; i++){std::cout<<" "<<thrs[i];}std::cout<<std::endl;
//    for(int i=0; i<Nsteps; i++){std::cout<<" "<<nelem[i];}std::cout<<std::endl;
    //determine threshold
    double thresh=filter.mc_min;
    for(int i=0; i<thrs.size()-1; i++){
      if(nelem[i]>filter.max_coeffs){ 
        thresh=thrs[i];
        break;
      }
    }
//std::cout<<"TEST: max_coeffs: thresh="<<thresh<<std::endl;
    for(const auto &elem : coeffs){
      if(get_relevance(elem.first)*std::abs(elem.second)*inv_largest>thresh){
        tmp.coeffs.insert(elem);
      }
    }

  //don't use max_coeffs:
  }else{          
    for(const auto &elem : coeffs){
      if(filter.passes(get_relevance(elem.first)*elem.second*inv_largest)){
        tmp.coeffs.insert(elem);
      }else if(filter.candidate_passes(get_relevance(elem.first)*elem.second*inv_largest)){
        tmp.candidates.insert(elem);
      }
    }
  }
  coeffs.swap(tmp.coeffs);
  candidates.swap(tmp.candidates);
}


void EDM_State::print(std::ostream &os)const{
  for(const auto & i : coeffs){
    i.first.print(os);
    os<<": "<<i.second.real()<<" "<<i.second.imag()<<std::endl;
  }
}
void EDM_State::print_closures(std::ostream &os)const{
  for(const auto & q : closures){
    os<<q.transpose()<<std::endl;
  }
}

void EDM_State::write(const std::string &fname)const{
  std::ofstream ofs(fname.c_str());
  ofs<<"EDM_State "<<Ldim.size()<<" "<<coeffs.size()<<" "<<candidates.size()<<" "<<closures.size()<<std::endl;
  Ldim.print(ofs);

  std::cerr<<"Function EDM_State::write NOT IMPLEMENTED YET!"<<std::endl;
  throw DummyException();
}

void EDM_State::set_single(const Eigen::VectorXcd & rho, double thr){
  clear();

  EDM_Index I(1);
  Ldim=I; Ldim[0]=rho.rows();

  for(int i=0; i<rho.rows(); i++){
    if(std::abs(rho(i))>thr){
      I[0]=i;
      coeffs.insert(std::make_pair(I, rho(i)));
    }
  }
}

void EDM_State::set_from_product(const std::vector<Eigen::VectorXcd> & rhos, double thr){
  clear();
  if(rhos.size()<1){ 
    std::cerr<<"EDM_State::set_from_product: rhos.size()<1!"<<std::endl;
    throw DummyException();
  }
  Ldim.resize(rhos.size());
  for(int i=0; i<rhos.size(); i++){
    Ldim[i]=rhos[i].rows();
  }

  EDM_Index I(Ldim.size());  
  while(I[0]<Ldim[0]){
    std::complex<double> val=1.;
    bool brk=false;
    for(int i=0; i<Ldim.size(); i++){
      val*=rhos[i](I[i]);
      if(abs(val)<=thr){
        brk=true;
        I.increase_at(i, Ldim);
        break;
      }
    }
    if(!brk){
      coeffs.insert(std::make_pair(I, val));
      I.increase(Ldim);
    }
  }

}

void EDM_State::set_from_product(const Eigen::VectorXcd & rho0, const Eigen::VectorXcd & rho1, double thr){
  std::vector<Eigen::VectorXcd> v; 
  v.push_back(rho0);
  v.push_back(rho1);
  set_from_product(v,thr);
}


EDM_State EDM_State::trace_over(int r)const{
 try{
  Ldim.check_in_range(r);
  check_closures_consistent();
  EDM_State res;
  res.Ldim=Ldim;
  res.Ldim.erase(r);
  Eigen::VectorXcd q=closures[r];
  res.closures=closures;
  res.closures.erase(res.closures.begin()+r);
  for(const auto &elem : coeffs){
    EDM_Index I=elem.first;
    int i=elem.first[r];
    I.erase(r);
    std::complex<double> val=q(i)*elem.second;
    res.add_to(I, val);
  }
  for(const auto &elem : candidates){
    EDM_Index I=elem.first;
    int i=elem.first[r];
    I.erase(r);
    std::complex<double> val=q(i)*elem.second;
    res.add_to(I, val);
  }
  return res;

 }catch(DummyException &e){
   std::cerr<<"called by EDM_State::trace_over("<<r<<")"<<std::endl; 
   throw e;
 }
}


Eigen::VectorXcd EDM_State::get_reduced(int r)const{
  if(r>=Ldim.size()){
    std::cerr<<"EDM_State::get_reduced(r="<<r<<"): r>=Ldim.size()="<<Ldim.size()<<"!"<<std::endl;
    throw DummyException();
  }
  EDM_State tmp=*this;
  for(int i=0; i<r; i++){
    EDM_State tmp2=tmp.trace_over(0);
    tmp.swap(tmp2);
  } 
  for(int i=r+1; i<Ldim.size(); i++){
    EDM_State tmp2=tmp.trace_over(1);
    tmp.swap(tmp2);
  } 
  Eigen::VectorXcd res=Eigen::VectorXcd::Zero(Ldim[r]);
  for(const auto &elem : tmp.coeffs){
    tmp.Ldim.check_in_range(elem.first); 
    res(elem.first[0])+=elem.second;
  }
  return res;
}

void EDM_State::check_coeffs_consistent()const{
  for(const auto & I : coeffs){
    Ldim.check_in_range(I.first);
  }
}
void EDM_State::check_closures_consistent()const{
  if(Ldim.size()!=closures.size()){
    std::cerr<<"EDM_State::check_closures_consistent: Ldim.size()!=closures.size() ("<<Ldim.size()<<" vs. "<<closures.size()<<")!"<<std::endl;
    throw DummyException();
  }
  for(int r=0; r<closures.size(); r++){
    if(closures[r].rows()!=Ldim[r]){
      std::cerr<<"EDM_State::check_closures_consistent: closures["<<r<<"].rows!=Ldim["<<r<<"]("<<closures[r].rows()<<" vs. "<<Ldim[r]<<")!"<<std::endl;
    throw DummyException();
    }
  }
}
}//namespace
