#ifndef ACE_PT_UNIACE_DEFINED_H
#define ACE_PT_UNIACE_DEFINED_H

#include <memory>
#include "ProcessTensorForward.hpp"
#include "ProcessTensorRepeat.hpp"
#include "ModePropagatorGenerator.hpp"
#include "TimeGrid.hpp"
#include "TruncationLayout.hpp"
#include "CompressionTree.hpp"
#include "Parameters.hpp"
#include "MPG_Selector.hpp"
#include "ProcessTensorElement.hpp"

namespace ACE{

struct VectorSpace{
  Eigen::MatrixXcd V;

  operator Eigen::MatrixXcd() const{return V;}

  void add_ortho(Eigen::VectorXcd v_new, double thr=1e-14){
    for(int i=0; i<V.cols(); i++){
      v_new -= V.col(i).dot(v_new)*V.col(i);
    }
    double n=v_new.stableNorm();
//std::cout<<n<<std::endl;
    if(n<thr)return;
    v_new.stableNormalize();
    //second orthogonalization 
    for(int i=0; i<V.cols(); i++){
      v_new -= V.col(i).dot(v_new)*V.col(i);
    }
    v_new.stableNormalize();
    V.conservativeResize(v_new.rows(), V.cols()+1);
    V.col(V.cols()-1)=v_new;
  }

  inline void add_ortho_Mat(const Eigen::MatrixXcd & U,double thr=1e-14){
    for(int i=0; i<U.cols(); i++){
      add_ortho(U.col(i), thr);
    }	  
  }
  inline void add_ortho_Mat_normalize(const Eigen::MatrixXcd & U,double thr=1e-14){
    for(int i=0; i<U.cols(); i++){
      add_ortho(U.col(i)/U.col(i).stableNorm(), thr);
    }	  
  }
};

class PT_uniACE{
public:
  Eigen::VectorXcd bath_init;
  ProcessTensorElement e;
  std::shared_ptr<CompressionTree> TTree;
 
  void set_from_ModePropagator(ModePropagator &mpg, const TimeGrid & tgrid, double dict_zero);

  void join_symmetric_and_sweep_forward(PT_uniACE &other, const TruncatedSVD & trunc, int n_steps=2);
  
  void sweep_backward(const TruncatedSVD & trunc, int n_steps=2);

  void sweep_backward_and_eliminate(const TruncatedSVD & trunc, int n_steps=2, bool stop_if_min_dim=false);


  void add_modes(ModePropagatorGenerator &mpg, const TimeGrid & tgrid, 
                 const TruncationLayout & trunc, double dict_zero,
                 int verbosity=1);

  inline void add_from_generators(
                 std::vector<std::shared_ptr<ModePropagatorGenerator> > &mpgs,
                 const TimeGrid & tgrid, const TruncationLayout & trunc, 
                 double dict_zero){

    for(size_t i=0; i<mpgs.size(); i++){
      add_modes(*mpgs[i].get(), tgrid, trunc, dict_zero);
    }
  }


  std::shared_ptr<ProcessTensorRepeat> get_PTR(const std::string &write_PT="") const;
  operator std::shared_ptr<ProcessTensorRepeat>() const{
    return get_PTR();
  }


  PT_uniACE(){
  }
  PT_uniACE(Parameters &param){
    std::vector<std::shared_ptr<ModePropagatorGenerator> > mpgs=MPG_Selector(param);
    TimeGrid tgrid(param);
    TruncationLayout trunc(param);
    double dict_zero = param.get_as_double("dict_zero",-1.);
    add_from_generators(mpgs, tgrid, trunc, dict_zero);
  }
  PT_uniACE(std::vector<std::shared_ptr<ModePropagatorGenerator> > &mpgs,
                 const TimeGrid & tgrid, const TruncationLayout & trunc,
                 double dict_zero){
    add_from_generators(mpgs, tgrid, trunc, dict_zero);
  }
};

}//namespace
#endif
