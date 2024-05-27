#ifndef COMBINE_FREE_AND_MODE_PROPAGATOR_DEFINED_H
#define COMBINE_FREE_AND_MODE_PROPAGATOR_DEFINED_H

#include "FreePropagator.hpp"
#include "ModePropagator.hpp"
#include <memory>
#include "DimensionExtender.h"

namespace ACE{

class CombineFreeAndModePropagator: public Propagator{
public:

  std::shared_ptr<FreePropagator> fprop;
  DimensionExtender f_de;
  
  std::vector<std::shared_ptr<ModePropagator> > P;
  std::vector<DimensionExtender> de;

  int dim;

  void check_consistency()const{
    if(fprop.use_count()<1){
      std::cerr<<"CombineFreeAndModePropagator: not set up!"<<std::endl;
      exit(1);
    }
    if(fprop->get_dim()!=f_de.dim_in()){
      std::cerr<<"CombineFreeAndModePropagatorSym: fprop->get_dim()!=f_de.dim_in()!"<<std::endl;
      exit(1);
    }
    if(f_de.dim_out()!=dim){ 
      std::cerr<<"CombineFreeAndModePropagator: f_de.dim_out()!=dim!"<<std::endl;
    }
   
    if(P.size()!=de.size()){
      std::cerr<<"CombineFreeAndModePropagator: P.size()!=de.size()!"<<std::endl;
      exit(1);
    }
    
    for(size_t i=0; i<P.size(); i++){
      if(P[i]->get_dim()!=de[i].dim_in()){
        std::cerr<<"CombineFreeAndModePropagator: P["<<i<<"].get_dim()!=de["<<i<<"].dim_in()!"<<std::endl;
        exit(1);
      }
      if(de[i].dim_out()!=dim){ 
        std::cerr<<"CombineFreeAndModePropagator:de["<<i<<"].dim_out()!=dim!"<<std::endl;
        exit(1);
      }
    }
  }

  virtual int get_dim()const{
    return dim;
  }

  virtual bool is_time_independent() const{
    if(!fprop->is_time_independent())return false;
    for(size_t i=0; i<P.size(); i++){
      if(!P[i]->is_time_independent())return false;
    }
    return true;
  }

  virtual void update(double t, double dt){
    check_consistency();
    
    fprop->update(t,dt);
    M=f_de.get_Liouville(fprop->M);

    for(size_t i=0; i<P.size(); i++){
      P[i]->update(t,dt/2.);
      Eigen::MatrixXcd M2 = M * de[i].get_Liouville(P[i]->M);
      P[i]->update(t+dt/2.,dt/2.);
      M = de[i].get_Liouville(P[i]->M) * M2;
    }

  }

  void set_fprop(std::shared_ptr<FreePropagator> fp, const DimensionExtender &de_){
    fprop=fp;
    f_de=de_;
    dim=f_de.dim_out();
    check_consistency();
  }
    
  void add(std::shared_ptr<ModePropagator> P_, const DimensionExtender &de_){
    P.push_back(P_);
    de.push_back(de_);
    check_consistency();
  }
 
  CombineFreeAndModePropagator(){
  }
  CombineFreeAndModePropagator(std::shared_ptr<FreePropagator> P_, const DimensionExtender &de_){
    set_fprop(P_,de_);
  }
  CombineFreeAndModePropagator(
          std::shared_ptr<FreePropagator> P1, const DimensionExtender &de1,
          std::shared_ptr<ModePropagator> P2, const DimensionExtender &de2){
    set_fprop(P1,de1);
    add(P2,de2);
  }

  CombineFreeAndModePropagator(std::shared_ptr<FreePropagator> P1, 
                               std::shared_ptr<ModePropagator> P2){
    int N=P1->get_dim();
    int M=P2->get_N_mode();
    if(P2->get_N_system()!=N){
      std::cerr<<"CombineFreeAndModePropagator: P2->get_N_system()!=N!"<<std::endl;   
      exit(1);
    }

    set_fprop(P1, DimensionExtender(1, N, M));
    add(P2, DimensionExtender(1, N*M, 1));
  }
  CombineFreeAndModePropagator(std::shared_ptr<FreePropagator> P1, 
                               std::shared_ptr<ModePropagator> P2,
                               std::shared_ptr<ModePropagator> P3){
    int N=P1->get_dim();
    int M1=P2->get_N_mode();
    int M2=P3->get_N_mode();
    if(P2->get_N_system()!=N){
      std::cerr<<"CombineFreeAndModePropagator: P2->get_N_system()!=N!"<<std::endl;   
      exit(1);
    }
    if(P3->get_N_system()!=N){
      std::cerr<<"CombineFreeAndModePropagator: P3->get_N_system()!=N!"<<std::endl;   
      exit(1);
    }

    set_fprop(P1, DimensionExtender(1, N, M1*M2));
    add(P2, DimensionExtender(1, N*M1, M2));
    add(P3, DimensionExtender(1, N, M1, M2, 1));
  }
 
  virtual ~CombineFreeAndModePropagator(){};
};
}//namespace
#endif
