#include "CompressedPropagator.hpp"
#include "CompressionTree.hpp"
#include "TimeGrid.hpp"
#include "ProcessTensorRepeat.hpp"
#include "LiouvilleTools.hpp"
#include "otimes.hpp"


namespace ACE{


void CompressedPropagator::compress(CompressedPropagator::PTE_rhoE0 &e_rhoE0,
                       std::shared_ptr<CompressionTree> TTree,
                       std::shared_ptr<CompressionTree> TTree_inv){
    if(TTree->T.cols()!=e_rhoE0.first.M.dim_d2){
      std::cerr<<"TTree->T.cols()!=e_rhoE0.first.M.dim_d2!"<<std::endl;
      throw DummyException();
    }
    if(TTree_inv->T.cols()!=e_rhoE0.first.M.dim_d1){
      std::cerr<<"TTree_inv->T.cols()!=e_rhoE0.first.M.dim_d1!"<<std::endl;
      throw DummyException();
    }
    if(TTree_inv->T.cols()!=e_rhoE0.first.closure.rows()){
      std::cerr<<"TTree_inv->T.cols()!=e_rhoE0.first.closure.rows()!"<<std::endl;
      throw DummyException();
    }
    if(TTree->T.cols()!=e_rhoE0.second.rows()){
      std::cerr<<"TTree->T.cols()!=e_rhoE0.second.rows()!"<<std::endl;
      throw DummyException();
    }
    e_rhoE0.first.M.inner_multiply_left(TTree_inv->T);
    e_rhoE0.first.M.inner_multiply_right(TTree->T.transpose());
    e_rhoE0.first.closure=TTree_inv->T*e_rhoE0.first.closure;
    e_rhoE0.second=TTree->T*e_rhoE0.second;
//For now:
    e_rhoE0.first.env_ops.clear();
}

CompressedPropagator::PTE_rhoE0 CompressedPropagator::calculate_single(
                       std::shared_ptr<CompressionTree> TTree,
                       std::shared_ptr<CompressionTree> TTree_inv,
//                       std::shared_ptr<ModePropagatorGenerator> mpg,
//                       int mode, 
                       std::shared_ptr<ModePropagator> mpp,
                       double ta, double dt, double dict_zero){
    PTE_rhoE0 e_rhoE0;
//    ModePropagatorPtr mpp=mpg->get_ModePropagator(mode);
    e_rhoE0.first.set_from_ModePropagator(*mpp.get(), ta, dt, dict_zero);
    e_rhoE0.second = H_Matrix_to_L_Vector(mpp->get_initial());

    compress(e_rhoE0, TTree, TTree_inv);

    return e_rhoE0;
}

CompressedPropagator::PTE_rhoE0 CompressedPropagator::calculate_ACE_sequential(
                       std::shared_ptr<CompressionTree> TTree, 
                       std::shared_ptr<CompressionTree> TTree_inv,
                       std::shared_ptr<ModePropagatorGenerator> mpg, 
                       int mode, double ta, double dt, double dict_zero){

  std::cout<<"CompressedPropagator::calculate_ACE_sequential: mode="<<mode<<std::endl;
  if(mode<0){
    std::cerr<<"CompressedPropagator::calculate_ACE_sequential: mode<0!"<<std::endl;
    throw DummyException();
  }else if(mode==0){
    return calculate_single(TTree, TTree_inv, mpg->get_ModePropagator(mode), ta, dt, dict_zero);
  }

  PTE_rhoE0 e_rhoE0=calculate_ACE_sequential(TTree->first, \
                     TTree_inv->first, mpg, mode-1, ta, dt, dict_zero);

  PTE_rhoE0 e_rhoE0_l=calculate_single(TTree->second, TTree_inv->second, \
                     mpg->get_ModePropagator(mode), ta, dt/2, dict_zero);

  PTE_rhoE0 e_rhoE0_r=calculate_single(TTree->second, TTree_inv->second, \
                     mpg->get_ModePropagator(mode), ta+dt/2, dt/2, dict_zero);

  std::cout<<"first:d1="<<e_rhoE0.first.M.dim_d1;
  std::cout<<" first:d2="<<e_rhoE0.first.M.dim_d2;
  std::cout<<" second:d1="<<e_rhoE0_l.first.M.dim_d1;
  std::cout<<" second:d_="<<e_rhoE0_l.first.M.dim_d2;
  std::cout<<" second:d2="<<e_rhoE0_r.first.M.dim_d2;
  std::cout<<std::endl;
  e_rhoE0.first.join_symmetric(e_rhoE0_l.first,e_rhoE0_r.first);
  e_rhoE0.second=Vector_otimes(e_rhoE0.second, e_rhoE0_l.second);
 
  compress(e_rhoE0, TTree, TTree_inv);

  std::cout<<"CompressedPropagator::calculate_sequential: done"<<std::endl;
  return e_rhoE0;
}


}//namespace
