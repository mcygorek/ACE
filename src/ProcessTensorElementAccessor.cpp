#include "ProcessTensorElementAccessor.hpp"
#include <iostream>
#include <fstream>
#include "DummyException.hpp"
#include <exception>
#include <typeinfo>
#include <stdexcept>
#include "ReadPT_struct.hpp"

namespace ACE{

void ProcessTensorElementAccessor::check_N(int dim)const{
  if(dim!=get_N()){
    std::cerr<<"ProcessTensorElementAccessor: dim!=get_N() ("<<dim<<" vs. "<<get_N()<<")!"<<std::endl;
    throw DummyException();
  }
}
/*
void ProcessTensorElementAccessor::join_thisfirst_sameinner(ProcessTensorElement &e1, const ProcessTensorElementAccessor &e2){
  if(*this!=e1.accessor){
    std::cerr<<"ProcessTensorElementAccessor::join_thisfirst_sameinner: *this!=e1.accessor!"<<std::endl;
    throw DummyException();
  }
  if(e1.M.dim_d2 != e2.M.dim_d1){
    std::cerr<<"ProcessTensorElementAccessor::join_thisfirst_sameinner: e1.M.dim_d2 != e2.M.dim_d1!"<<std::endl;
    throw DummyException();
  }
  
  IF_OD_Dictionary dict_first=e1.accessor.dict;
  IF_OD_Dictionary dict_second=e2.accessor.dict;
  IF_OD_Dictionary dict_new=dict_first; dict_new.join(dict_second);
  int NL=dict_new.get_NL();
  std::vector<std::vector<int> > rev_new=dict_new.get_reverse_beta();

  MPS_Matrix M(dict_new.get_reduced_dim(), e1.M.dim_d1, e2.M.dim_d2);
  M.set_zero();
  for(int k=0; k<NL; k++){
    for(int j=0; j<NL; j++){
      int i_ind=dict_first.beta[j*NL+k]; if(i_ind<0)continue;
      for(int i=0; i<NL; i++){
        int i_ind2=dict_second.beta[i*NL+j]; if(i_ind2<0)continue;
        int i_ind3=dict_new.beta[i*NL+k]; if(i_ind3<0)continue;
        if(rev_new[i_ind3][0]!=i*NL+k)continue; //only modify first occurance
        for(int d1=0; d1<e1.M.dim_d1; d1++){
          for(int d2=0; d2<e1.M.dim_d2; d2++){
            for(int d3=0; d3<e2.M.dim_d2; d3++){
              M(i_ind3, d1, d3)+= e1.M(i_ind, d1, d2)*e2.M(i_ind2, d2, d3);
            }
          }
        }
      }
    }
  }
  e1.M.swap(M);
  e1.accessor.dict=dict_new;
  e1.closure=e.closure;
  e1.env_ops=e.env_ops;
  e1.forwardNF=e.forwardNF;
}
*/

void ProcessTensorElementAccessor::propagate(
          Eigen::MatrixXcd &state, const MPS_Matrix &M, int dim1_front, 
          const ReadPT_struct &expand)const{
 
  int dim1_back=state.cols()/dim1_front/M.dim_d1;
  if(state.cols() != dim1_front * M.dim_d1 * dim1_back){
    std::cerr<<"ProcessTensorElementAccessor::propagate: state.cols() != dim1_front * M.dim_d1 * dim1_back!"<<std::endl;
    throw DummyException();
  }

  IF_OD_Dictionary thisdict=dict;
  if(expand.have_to_expand()){
    thisdict.expand(expand,false);
  }

  int NL=thisdict.get_NL();
  if(state.rows()!=NL){
    std::cerr<<"ProcessTensorElementAccessor::propagate: state.rows()!=NL ("<<state.rows()<<" vs. "<<NL<<")!"<<std::endl; 
    throw DummyException();
  }
 
  Eigen::MatrixXcd state2 = Eigen::MatrixXcd::Zero(state.rows(), dim1_front * M.dim_d2 * dim1_back);

  for(int i=0; i<NL; i++){
    for(int j=0; j<NL; j++){
      int i_ind=thisdict.beta[i*NL+j]; if(i_ind<0)continue;
      for(int df=0; df<dim1_front; df++){
        for(int d1=0; d1<M.dim_d1; d1++){
          for(int d2=0; d2<M.dim_d2; d2++){
            for(int db=0; db<dim1_back; db++){
              state2(i, (df*M.dim_d2+d2)*dim1_back+db) += 
                M(i_ind, d1, d2) * state(j, (df*M.dim_d1+d1)*dim1_back+db);
            }
          }
        }
      }
    }
  }
  state.swap(state2);
}

void ProcessTensorElementAccessor::set_trivial(int N, MPS_Matrix &M){
  dict.set_trivial(N);
  M.resize(1, 1, 1); M.fill(1.);
}

void ProcessTensorElementAccessor::apply_HilbertSpaceRotation(
         const HilbertSpaceRotation &hs_rot, MPS_Matrix &M, double dict_zero){
  if(!hs_rot.used()){
    return;
  }
  int N=hs_rot.U.rows();
  if(dict.get_N()!=N){
    std::cerr<<"ProcessTensorElementAccessor::apply_HilbertSpaceRotation: dict.get_N()!=hs_rot.U.rows()!"<<std::endl;
    throw DummyException();
  }
  MPS_Matrix tmp(N*N*N*N, M.dim_d1, M.dim_d2);
  tmp.set_zero();

  Eigen::MatrixXcd Map=Eigen::MatrixXcd::Zero(N*N*N*N, M.dim_i);
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      for(int k=0; k<N; k++){
        for(int l=0; l<N; l++){
          for(int m1=0; m1<N; m1++){
            for(int o1=0; o1<N; o1++){
              for(int m2=0; m2<N; m2++){
                for(int o2=0; o2<N; o2++){
                  int i_ind=dict.beta[((m1*N+o1)*N+m2)*N+o2]; if(i_ind<0)continue;
                  Map(((i*N+j)*N+k)*N+l, i_ind) +=
                      hs_rot.U(i,m1) * hs_rot.U.adjoint()(m2,k) *
                      hs_rot.U(l,o2) * hs_rot.U.adjoint()(o1,j);
                }
              }
            }
          }
        }
      }
    }
  }
  for(int index=0; index<N*N*N*N; index++){
    for(int i=0; i<M.dim_i; i++){
      if(abs(Map(index, i))<1e-24)continue;
      for(int d1=0; d1<M.dim_d1; d1++){
        for(int d2=0; d2<M.dim_d2; d2++){
          tmp(index, d1, d2) += Map(index, i) * M(i, d1, d2);
        }
      }
    }
  }

  M.swap(tmp);
  dict.set_default(N);

  if(dict_zero>0.){
    dict.detect(M, dict_zero);
    dict.reduce_MPS_Matrix(M);
  }

//  calculate_dict(dict_zero);
//  reduce_to_dict();
}

std::pair<double, std::vector<int> > ProcessTensorElementAccessor::closure_indices()const{
  std::pair<double, std::vector<int> > ret;
  int N=get_N();
  ret.first=1./((double)N);
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      int i_ind=dict(((i*N+i)*N+j)*N+j);
      if(i_ind>=0)ret.second.push_back(i_ind);
    }
  }
  return ret;
}

ProcessTensorElementAccessor::VVPI ProcessTensorElementAccessor::
    join_thisfirst_indices(const ProcessTensorElementAccessor &other){

  int NL=dict.get_NL();
  IF_OD_Dictionary ndict=dict; ndict.join(other.dict);
  std::vector<std::vector<int> > newrev=ndict.get_reverse_beta();

  VVPI i_list(ndict.get_reduced_dim());

  for(int i=0; i<NL; i++){
    for(int j=0; j<NL; j++){
      int i_ind2=other.dict.beta[i*NL+j];   if(i_ind2<0)continue;
      for(int k=0; k<NL; k++){
        int i_ind_new=ndict.beta[i*NL+k];  if(i_ind_new<0)continue;
        if(newrev[i_ind_new][0]!=i*NL+k)continue;
        int i_ind=dict.beta[j*NL+k];     if(i_ind<0)continue;
        i_list[i_ind_new].push_back(std::pair<int,int>(i_ind, i_ind2));
      }
    } 
  }
  dict=ndict;
  return i_list;
}

ProcessTensorElementAccessor::VVPI ProcessTensorElementAccessor::
    join_thissecond_indices(const ProcessTensorElementAccessor &other){

  int NL=dict.get_NL();
  IF_OD_Dictionary ndict=dict; ndict.join(other.dict);
  std::vector<std::vector<int> > newrev=ndict.get_reverse_beta();

  VVPI i_list(ndict.get_reduced_dim());

  for(int i=0; i<NL; i++){
    for(int j=0; j<NL; j++){
      int i_ind=dict.beta[i*NL+j];         if(i_ind<0)continue;

      for(int k=0; k<NL; k++){
        int i_ind_new=ndict.beta[i*NL+k];  if(i_ind_new<0)continue;
        if(newrev[i_ind_new][0]!=i*NL+k)continue;

        int i_ind2=other.dict.beta[j*NL+k];  if(i_ind2<0)continue;
        i_list[i_ind_new].push_back(std::pair<int,int>(i_ind, i_ind2));
      }
    } 
  }
  dict=ndict;
  return i_list;
}


void ProcessTensorElementAccessor::join(const VVPI &i_list,  
                                   MPS_Matrix & M, const MPS_Matrix & M2){

  MPS_Matrix tmp(i_list.size(), M.dim_d1*M2.dim_d1, M.dim_d2*M2.dim_d2);
  tmp.fill(0.);
  for(int i=0; i<(int)i_list.size(); i++){
    for(int j=0; j<(int)i_list[i].size(); j++){
      for(int d1=0; d1<M.dim_d1; d1++){
        for(int d2=0; d2<M.dim_d2; d2++){
          for(int od1=0; od1<M2.dim_d1; od1++){
            for(int od2=0; od2<M2.dim_d2; od2++){
              tmp(i, d1*M2.dim_d1 + od1, d2*M2.dim_d2 + od2)+=
                M(i_list[i][j].first,d1,d2) * M2(i_list[i][j].second,od1,od2);
            }
          }
        }
      }
    }
  } 
  M.swap(tmp);   
}

/*
void ProcessTensorElementAccessor::join_meanfield(IF_OD_Dictionary other_dict, 
                                   MPS_Matrix & M, MPS_Matrix M2){

  if(dict.get_N()!=other_dict.get_N()){
    std::cerr<<"ProcessTensorElementAccessor::join_meanfield: dict.get_N()!=other_dict.get_N()!"<<std::endl;
    throw DummyException();
  }

  IF_OD_Dictionary old_dict=dict;
  dict.join(other_dict);
  dict.translate_MPS_Matrix(M, old_dict);
  dict.translate_MPS_Matrix(M2, other_dict);

  MPS_Matrix tmp(dict.get_reduced_dim(), M.dim_d1+M2.dim_d1, M.dim_d2+M2.dim_d2);
  tmp.fill(0.);
  for(int i=0; i<dict.get_reduced_dim(); i++){
    for(int d1=0; d1<M.dim_d1; d1++){
      for(int d2=0; d2<M.dim_d2; d2++){
        tmp(i, d1, d2)+=M(i,d1,d2); 
      }
    }
    for(int od1=0; od1<M2.dim_d1; od1++){
      for(int od2=0; od2<M2.dim_d2; od2++){
         tmp(i, M.dim_d1 + od1, M2.dim_d2 + od2)+=M2(i,od1,od2);
      }
    }
  } 
  M.swap(tmp);   
}
*/

void ProcessTensorElementAccessor::join_thisfirst(
                 const ProcessTensorElementAccessor &other,
                 MPS_Matrix & M, const MPS_Matrix & M2){

  VVPI i_list=join_thisfirst_indices(other);
  join(i_list, M, M2);  
}

void ProcessTensorElementAccessor::join_thissecond(
                 const ProcessTensorElementAccessor &other,
                 MPS_Matrix & M, const MPS_Matrix & M2){

  VVPI i_list=join_thissecond_indices(other);
  join(i_list, M, M2);  
}

void ProcessTensorElementAccessor::join_symmetric(
         const ProcessTensorElementAccessor & acc_L,
         const ProcessTensorElementAccessor & acc_R,
         MPS_Matrix & M, const MPS_Matrix & M_L, const MPS_Matrix & M_R){

  int dim_d1_orig=M.dim_d1;
  int dim_d2_orig=M.dim_d2;

  VVPI i_list=join_thissecond_indices(acc_L);
  join(i_list, M, M_L);  

  VVPI i_list2=join_thisfirst_indices(acc_R);
//  join(i_list2, M, M_R);  


  MPS_Matrix tmp2(i_list2.size(), dim_d1_orig*M_L.dim_d1, dim_d2_orig*M_R.dim_d2);
  tmp2.fill(0.);
  for(int i=0; i<(int)i_list2.size(); i++){
    for(int j=0; j<(int)i_list2[i].size(); j++){
      for(int d1=0; d1<dim_d1_orig; d1++){
        for(int d2=0; d2<dim_d2_orig; d2++){
          for(int od0=0; od0<M_L.dim_d1; od0++){
            for(int od1=0; od1<M_R.dim_d1; od1++){
              for(int od2=0; od2<M_R.dim_d2; od2++){
                tmp2(i, d1*M_L.dim_d1 + od0, d2*M_L.dim_d2 + od2) += 
                  M(i_list2[i][j].first, d1*M_L.dim_d1+od0, d2*M_L.dim_d2+od1) 
                * M_R(i_list2[i][j].second, od1, od2);
              }
            }
          }
        }
      }
    }
  }
  M.swap(tmp2);   
}


void ProcessTensorElementAccessor::join_select_indices(
               const VVPI &i_list, MPS_Matrix &M1, const MPS_Matrix &M2,
               const ProcessTensorElementAccessor & acc_other,
               const SelectIndices & k_left, const SelectIndices & k_right){

//std::cout<<"i_list.size()="<<i_list.size()<<" k_left.size()="<<k_left.size()<<" k_right.size()="<<k_right.size()<<std::endl;

  MPS_Matrix tmp;
  try{
    tmp.resize(i_list.size(), k_left.size(), k_right.size());
  }catch(const std::exception& e){
    std::cerr<<"Exception '"<<e.what()<<"' in ProcessTensorElementAccessor::join_select_indices!"<<std::endl;
    std::cerr<<"Requested memory for MPS_Matrix of dimensions "<<i_list.size()<<","<<k_left.size()<<","<<k_right.size()<<"."<<std::endl;
    throw; // DummyException();
//  }catch(...){
//    std::exception_ptr p = std::current_exception();
//    std::cerr<<"Exception '"<<(p ? p.__cxa_exception_type()->name() : "null")<<"' in ProcessTensorElementAccessor::join_select_indices!"<<std::endl;
//    std::cerr<<"Requested memory for MPS_Matrix of dimensions "<<i_list.size()<<","<<k_left.size()<<","<<k_right.size()<<"."<<std::endl;
//    throw; // DummyException();
  }
//std::cout<<"tmp.dim_i="<<tmp.dim_i<<" M1.dim_i="<<M1.dim_i<<" M2.dim_i="<<M2.dim_i<<std::endl;
  tmp.fill(0.);
  for(int i=0; i<(int)i_list.size(); i++){
//std::cout<<"i_list["<<i<<"].size()="<<(int)i_list[i].size()<<std::endl;
    for(int j=0; j<(int)i_list[i].size(); j++){
      for(int kl=0; kl<k_left.size(); kl++){
        for(int kr=0; kr<k_right.size(); kr++){
          tmp(i, kl, kr) +=
              M1(i_list[i][j].first, k_left[kl].first , k_right[kr].first) 
            * M2(i_list[i][j].second, k_left[kl].second, k_right[kr].second);
        }
      }
    }
  } 
  M1.swap(tmp);   
//std::cout<<"accessor: TEST2"<<std::endl;
}
void ProcessTensorElementAccessor::join_select_indices_alternate(
               int n, MPS_Matrix &M1, const MPS_Matrix &M2,
               const ProcessTensorElementAccessor & acc_other,
               const SelectIndices & k_left, const SelectIndices & k_right){
  VVPI i_list;
  if(n%2==0){
    i_list=join_thissecond_indices(acc_other);
  }else{
    i_list=join_thisfirst_indices(acc_other);
  }
  join_select_indices(i_list, M1, M2, acc_other, k_left, k_right);
}


void ProcessTensorElementAccessor::read_binary(std::istream &ifs){
  dict.read_binary(ifs);
}
void ProcessTensorElementAccessor::write_binary(std::ostream &ofs)const{
  dict.write_binary(ofs);
}


}//namspace
