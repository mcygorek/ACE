#include "AddPT.hpp"
#include "DummyException.hpp"

namespace ACE{
namespace AddPT{
bool split_alternating(ProcessTensorElement &E11, ProcessTensorElement &E12,
                       ProcessTensorElement &E21, ProcessTensorElement &E22){

  if(E11.M.dim_d2<2){
    int Nsys=sqrt(E11.M.dim_i);
    E21.set_trivial(E11.get_N());
    E22.set_trivial(E11.get_N());
    E21.M.set_zero();
    E22.M.set_zero();
    return false;
  }

  MPS_Matrix M1(E11.M.dim_i, E11.M.dim_d1, (E11.M.dim_d2+1)/2);
  MPS_Matrix M2(E21.M.dim_i, E21.M.dim_d1, E11.M.dim_d2-M1.dim_d2);
  
  //low-end elements: E11 and E21:
  { // M:
    for(int i=0; i<M1.dim_i; i++){  
      for(int d1=0; d1<M1.dim_d1; d1++){  
        for(int d2=0; d2<M1.dim_d2; d2++){
          M1(i, d1, d2)=E11.M(i, d1, 2*d2);
        }
      }  
    }
    for(int i=0; i<M2.dim_i; i++){  
      for(int d1=0; d1<M2.dim_d1; d1++){  
        for(int d2=0; d2<M2.dim_d2; d2++){
          M2(i, d1, d2)=E11.M(i, d1, 2*d2+1);
        }
      }  
    }
    E11.M.swap(M1);
    E21.M.swap(M2);
    E21.accessor=E11.accessor;
  } 
  { //closures

    Eigen::VectorXcd closure1(E11.M.dim_d2);
    Eigen::VectorXcd closure2(E21.M.dim_d2);
    for(int d2=0; d2<E11.M.dim_d2; d2++){
      closure1(d2)=E11.closure(2*d2); 
    }
    for(int d2=0; d2<E21.M.dim_d2; d2++){
      closure2(d2)=E11.closure(2*d2+1); 
    }
    E11.closure=closure1;
    E21.closure=closure2;

  }
  { //forwardNF:
    if(E11.forwardNF.rows()>0){
      Eigen::VectorXd forwardNF1(E11.M.dim_d2);
      Eigen::VectorXd forwardNF2(E21.M.dim_d2);
      for(int d2=0; d2<E11.M.dim_d2; d2++){
        forwardNF1(d2)=E11.forwardNF(2*d2); 
      }
      for(int d2=0; d2<E21.M.dim_d2; d2++){
        forwardNF2(d2)=E11.forwardNF(2*d2+1); 
      }
      E11.forwardNF=forwardNF1;
      E21.forwardNF=forwardNF2;
    }else{
      E21.forwardNF=Eigen::VectorXd(0);
    }
  }

  { //env_ops:
    for(int o=0; o<E11.env_ops.size(); o++){
      if(E11.env_ops[o].rows()!=E11.M.dim_d2+E21.M.dim_d2)continue;

      Eigen::VectorXcd env_ops1(E11.M.dim_d2);
      Eigen::VectorXcd env_ops2(E21.M.dim_d2);
      for(int d2=0; d2<E11.M.dim_d2; d2++){
        env_ops1(d2)=E11.env_ops[o](2*d2); 
      }
      for(int d2=0; d2<E21.M.dim_d2; d2++){
        env_ops2(d2)=E11.env_ops[o](2*d2+1); 
      }
      E11.env_ops.ops[o]=env_ops1;
      E21.env_ops.ops[o]=env_ops2;
    }
  } 

  //high-end elements: E12 and E22:
  { // M
    MPS_Matrix M1(E12.M.dim_i, (E12.M.dim_d1+1)/2, E12.M.dim_d2);
    MPS_Matrix M2(E12.M.dim_i, E12.M.dim_d1-M1.dim_d1, E22.M.dim_d2);
  
    for(int i=0; i<E12.M.dim_i; i++){  
      for(int d1=0; d1<E12.M.dim_d1; d1++){  
        for(int d2=0; d2<E12.M.dim_d2; d2++){
          M1(i, d1, d2)=E12.M(i, 2*d1, d2);
        }
      }  
    }
    for(int i=0; i<E22.M.dim_i; i++){  
      for(int d1=0; d1<E22.M.dim_d1; d1++){  
        for(int d2=0; d2<E22.M.dim_d2; d2++){
          M2(i, d1, d2)=E12.M(i, 2*d1+1, d2);
        }
      }  
    }
    E12.M.swap(M1);
    E22.M.swap(M2);
    E22.accessor=E12.accessor;
  } 
  {//backwardNF
    if(E12.backwardNF.rows()>0){
      Eigen::VectorXd backwardNF1(E12.M.dim_d1);
      Eigen::VectorXd backwardNF2(E22.M.dim_d1);
      for(int d1=0; d1<E12.M.dim_d1; d1++){
        backwardNF1(d1)=E12.backwardNF(2*d1); 
      }
      for(int d1=0; d1<E22.M.dim_d1; d1++){
        backwardNF2(d1)=E12.backwardNF(2*d1+1); 
      }
      E12.backwardNF=backwardNF1;
      E22.backwardNF=backwardNF2;
    }else{  
      E22.backwardNF=Eigen::VectorXd(0);
    }
  }
  return true;
}

void add(ProcessTensorElement &E11, 
         ProcessTensorElement &E21){

  ProcessTensorElement e;
  if(E11.M.dim_i!=E21.M.dim_i){
    std::cerr<<"AddPT::add: E11.M.dim_i!=E21.M.dim_i!"<<std::endl;
    throw DummyException();
  }
  if(E11.env_ops.size()!=E21.env_ops.size()){
    std::cerr<<"AddPT::add: E11.env_ops.size()!=E21.env_ops.size()!"<<std::endl;
    throw DummyException();
  }
 
  {
    MPS_Matrix M_tmp(E11.M.dim_i, E11.M.dim_d1+E21.M.dim_d1, E11.M.dim_d2+E21.M.dim_d2);
    M_tmp.set_zero();
    for(int i=0; i<M_tmp.dim_i; i++){
      for(int d1=0; d1<E11.M.dim_d1; d1++){
        for(int d2=0; d2<E11.M.dim_d2; d2++){
          M_tmp(i,d1,d2)=E11.M(i,d1,d2);
        }
      }
      for(int d1=0; d1<E21.M.dim_d1; d1++){
        for(int d2=0; d2<E21.M.dim_d2; d2++){
          M_tmp(i,d1+E11.M.dim_d1,d2+E21.M.dim_d2)=E21.M(i,d1,d2);
        }
      }
    }
    E11.M.swap(M_tmp);
  }

  { 
    Eigen::VectorXcd closure_tmp(E11.closure.rows()+E21.closure.rows());
    closure_tmp << E11.closure, E21.closure;
    E11.closure=closure_tmp;
  }

  for(size_t o=0; o<E11.env_ops.size(); o++){
    Eigen::VectorXcd obs_tmp(E11.env_ops[o].rows()+E21.env_ops[o].rows());
    obs_tmp << E11.env_ops[o], E21.env_ops[o];
    E11.env_ops.ops[o] = obs_tmp;
  }
 
  E11.clearNF(); //<- could be recovered by sorting
}

void add_head(ProcessTensorElement &E11, 
                ProcessTensorElement &E21){

  if(E11.M.dim_d1!=1 || E21.M.dim_d1!=1){
    std::cerr<<"AddPT::add_head: E11.M.dim_d1="<<E11.M.dim_d1<<"!=1 || E21.M.dim_d1"<<E21.M.dim_d1<<"!=1!"<<std::endl;
    throw DummyException();
  }
  
  add(E11, E21);

  MPS_Matrix M_tmp(E11.M.dim_i, 1, E11.M.dim_d2);
  M_tmp.set_zero();
  for(int i=0; i<E11.M.dim_i; i++){
    for(int d1=0; d1<E11.M.dim_d1; d1++){
      for(int d2=0; d2<E11.M.dim_d2; d2++){
        M_tmp(i, 0, d2)+=E11.M(i, d1, d2);
      }
    }
  }
  E11.M.swap(M_tmp);

  E11.clearNF();  //<- Could be recovered by sorting
}

void add_tail(ProcessTensorElement &E11,
              ProcessTensorElement &E21){

  if(E11.M.dim_d2!=1 || E21.M.dim_d2!=1){
    std::cerr<<"AddPT::add_tail: E11.M.dim_d2="<<E11.M.dim_d2<<"!=1 || E21.M.dim_d2"<<E21.M.dim_d2<<"!=1!"<<std::endl;
    throw DummyException();
  }
  
  add(E11, E21);

  MPS_Matrix M_tmp(E11.M.dim_i, E11.M.dim_d1, 1);
  M_tmp.set_zero();
  for(int i=0; i<E11.M.dim_i; i++){
    for(int d1=0; d1<E11.M.dim_d1; d1++){
      for(int d2=0; d2<E11.M.dim_d2; d2++){
        M_tmp(i, d1, 0)+=E11.M(i, d1, d2);
      }
    }
  }
  E11.M.swap(M_tmp);
  
  Eigen::VectorXcd closure_tmp(1);
  closure_tmp(0)=E11.closure(0)+E11.closure(1);
  E11.closure=closure_tmp;

  for(size_t o=0; o<E11.env_ops.size(); o++){
    Eigen::VectorXcd obs_tmp(1);
    obs_tmp(0)=E11.env_ops[o](0)+E11.env_ops[o](1);
    E11.env_ops.ops[o] = obs_tmp;
  }

  E11.clearNF();  //<- Could be recovered by sorting
}



}//namespace
}//namespace
