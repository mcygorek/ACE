#ifndef COMBINE_PROPAGATOR_SYM_DEFINED_H
#define COMBINE_PROPAGATOR_SYM_DEFINED_H

#include "Propagator.hpp"
#include "Smart_Ptr.h"
#include "DimensionExtender.h"

namespace ACE{
/*
class ExtendPropagator: public Propagator{
public:
  Smart_Ptr<Propagator> P;
  DimensionExtender de;

  virtual int get_dim()const{ return de.dim_out(); }
  virtual bool is_time_independent()const{ return P->is_time_independent(); }
  virtual void update(double t, double dt){
    P->update(t, dt);
    M=de.get_Liouville(P->M);
  }

  ExtendPropagator(){}
  ExtendPropagator(Smart_Ptr<Propagator> &P_, const DimensionExtender &de_)
   : P(P_), de(de_) {
    de.check_consistency(); 
    if(de.dim_in()!=P->get_dim()){
      std::cerr<<"Error: ExtendPropagator: de.dim_in()!=P->get_dim()!"<<std::endl;
      exit(1);
    }
 
    if(M.rows()==de.dim_in()){   
      M=de.get_Liouville(P->M);
    }
  }
  virtual ~ExtendPropagator(){}
};
*/

// Combine two Propagators using the symmetric Trotter-type decomposition
class CombinePropagatorSym: public Propagator{
public:

  std::vector<Smart_Ptr<Propagator> > P;
  std::vector<DimensionExtender> de;

  int dim;

  void check_consistency()const{
    if(P.size()<1){
      std::cerr<<"CombinePropagatorSym: not set up!"<<std::endl;
      exit(1);
    }
    if(P.size()!=de.size()){
      std::cerr<<"CombinePropagatorSym: P.size()!=de.size()!"<<std::endl;
      exit(1);
    }
    
    for(size_t i=0; i<P.size(); i++){
      if(P[i]->get_dim()!=de[i].dim_in()){
        std::cerr<<"CombinePropagatorSym: P["<<i<<"].get_dim()!=de["<<i<<"].dim_in()!"<<std::endl;
        exit(1);
      }
      if(de[i].dim_out()!=dim){ 
        std::cerr<<"CombinePropagatorSym:de["<<i<<"].dim_out()!=dim!"<<std::endl;
        exit(1);
      }
    }
  }

  virtual int get_dim()const{
    return dim;
  }

  virtual bool is_time_independent() const{
    for(size_t i=0; i<P.size(); i++){
      if(!P[i]->is_time_independent())return false;
    }
    return true;
  }

  virtual void update(double t, double dt){
    check_consistency();
    P[0]->update(t,dt);
    M=de[0].get_Liouville(P[0]->M);

    for(size_t i=1; i<P.size(); i++){
      P[i]->update(t,dt/2.);
      Eigen::MatrixXcd M2 = M * de[i].get_Liouville(P[i]->M);
      P[i]->update(t+dt/2.,dt/2.);
      M = de[i].get_Liouville(P[i]->M) * M2;
    }
  }

  void add(Smart_Ptr<Propagator> P_, const DimensionExtender &de_){
    if(P.size()<1)dim=de_.dim_out();
    P.push_back(P_);
    de.push_back(de_);
    check_consistency();
  }
 
  CombinePropagatorSym(){
  }
  CombinePropagatorSym(Smart_Ptr<Propagator> P_, const DimensionExtender &de_){
    add(P_,de_);
  }
  CombinePropagatorSym(Smart_Ptr<Propagator> P1, const DimensionExtender &de1,
                       Smart_Ptr<Propagator> P2, const DimensionExtender &de2){
    add(P1,de1);
    add(P2,de2);
  }

  virtual ~CombinePropagatorSym(){};
};

}//namespace
#endif
