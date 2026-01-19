#include "ProcessTensorElementAccessor.hpp"
#include <iostream>
#include <fstream>
#include "DummyException.hpp"
#include <exception>
#include <typeinfo>
#include <stdexcept>
#include "ReadPT_struct.hpp"
#include "Timings.hpp"

namespace ACE{

void ProcessTensorElementAccessor::check_N(int dim)const{
  if(dim!=get_N()){
    std::cerr<<"ProcessTensorElementAccessor: dim!=get_N() ("<<dim<<" vs. "<<get_N()<<")!"<<std::endl;
    throw DummyException();
  }
}

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

//  typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatriXcdRM;

  if(dim1_front==1 && dim1_back==1){
    for(int i=0; i<NL; i++){
      for(int j=0; j<NL; j++){
        int i_ind=thisdict.beta[i*NL+j]; if(i_ind<0)continue;
//        for(int d1=0; d1<M.dim_d1; d1++){
//          for(int d2=0; d2<M.dim_d2; d2++){
////              state2(i, d2) += M(i_ind, d1, d2) * state(j, d1);
//              state2(i, d2) +=  state(j, d1) * M(i_ind, d1, d2);
//          }
//        }
          state2.row(i).noalias() += state.row(j) * 
Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 0, Eigen::OuterStride<> >(M.mem+i_ind*M.dim_d2,M.dim_d1,M.dim_d2,Eigen::OuterStride<>(M.dim_i*M.dim_d2));
      }
    }
    state.swap(state2);
    return;
  }

/*
  if(dim1_front==1){
    for(int db=0; db<dim1_back; db++){
      Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>, 0 , Eigen::OuterStride<> > mp(&state(0,db), NL, M.dim_d1, Eigen::OuterStride<>(dim1_back*NL));
      Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>, 0 , Eigen::OuterStride<> > mp2(&state2(0,db), NL, M.dim_d2, Eigen::OuterStride<>(dim1_back*NL));
      for(int i=0; i<NL; i++){
        for(int j=0; j<NL; j++){
          int i_ind=thisdict.beta[i*NL+j]; if(i_ind<0)continue;
//          for(int d1=0; d1<M.dim_d1; d1++){
//            for(int d2=0; d2<M.dim_d2; d2++){
////                state2(i, d2*dim1_back+db) +=  state(j, d1*dim1_back+db) * M(i_ind, d1, d2);
//                mp2(i, d2) +=  mp(j, d1) * M(i_ind, d1, d2);
//            }
//          }
          mp2.row(i).noalias() += mp.row(j) * \
Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 0, Eigen::OuterStride<> >(M.mem+i_ind*M.dim_d2,M.dim_d1,M.dim_d2,Eigen::OuterStride<>(M.dim_i*M.dim_d2));
        }
      }
    }
    state.swap(state2);
    return;
  }
*/

  for(int df=0; df<dim1_front; df++){
    for(int db=0; db<dim1_back; db++){
      Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>, 0 , Eigen::OuterStride<> > mp(&state(0,df*M.dim_d1*dim1_back+db), NL, M.dim_d1, Eigen::OuterStride<>(dim1_back*NL));
      Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>, 0 , Eigen::OuterStride<> > mp2(&state2(0,df*M.dim_d2*dim1_back+db), NL, M.dim_d2, Eigen::OuterStride<>(dim1_back*NL));
      for(int i=0; i<NL; i++){
        for(int j=0; j<NL; j++){
          int i_ind=thisdict.beta[i*NL+j]; if(i_ind<0)continue;
//          for(int d1=0; d1<M.dim_d1; d1++){
//            for(int d2=0; d2<M.dim_d2; d2++){
//              state2(i, (df*M.dim_d2+d2)*dim1_back+db) +=  \
                M(i_ind, d1, d2) * state(j, (df*M.dim_d1+d1)*dim1_back+db);
//            }
//          }
          mp2.row(i).noalias() += mp.row(j) * \
Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 0, Eigen::OuterStride<> >(M.mem+i_ind*M.dim_d2,M.dim_d1,M.dim_d2,Eigen::OuterStride<>(M.dim_i*M.dim_d2));
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

  constexpr bool debug=false;
  if(debug){std::cout<<"join_thissecond_indices: Mark1"<<std::endl;}
  int NL=dict.get_NL();
  IF_OD_Dictionary ndict=dict; ndict.join(other.dict);
  std::vector<std::vector<int> > newrev=ndict.get_reverse_beta();

  if(debug){std::cout<<"join_thissecond_indices: Mark2"<<std::endl;}
  VVPI i_list(ndict.get_reduced_dim());

  if(debug){std::cout<<"join_thissecond_indices: Mark3"<<std::endl;}
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
  if(debug){std::cout<<"join_thissecond_indices: Mark4"<<std::endl;}
  dict=ndict;
  if(debug){std::cout<<"join_thissecond_indices: Mark5"<<std::endl;}
  return i_list;
} 

ProcessTensorElementAccessor::VVTI ProcessTensorElementAccessor::join_symmetric_indices(
                              const ProcessTensorElementAccessor &otherL,
                              const ProcessTensorElementAccessor &otherR){

  int NL=dict.get_NL();
  IF_OD_Dictionary idict=dict; idict.join(otherL.dict);
  IF_OD_Dictionary ndict=otherR.dict; ndict.join(idict);
  std::vector<std::vector<int> > newrev=ndict.get_reverse_beta();

  VVTI i_list(ndict.get_reduced_dim());

  for(int i=0; i<NL; i++){
    for(int l=0; l<NL; l++){
      int i_ind_new=ndict.beta[i*NL+l];  if(i_ind_new<0)continue;
      if(newrev[i_ind_new][0]!=i*NL+l)continue;
      for(int j=0; j<NL; j++){
        int i_ind_R=otherR.dict.beta[i*NL+j];    if(i_ind_R<0)continue;
        for(int k=0; k<NL; k++){
          int i_ind=dict.beta[j*NL+k];  if(i_ind<0)continue;
          int i_ind_L=otherL.dict.beta[k*NL+l];    if(i_ind_L<0)continue;
  i_list[i_ind_new].push_back(std::tuple<int,int,int>(i_ind_R, i_ind, i_ind_L));
        }
      }
    } 
  }
  dict=ndict;
  return i_list;
}

void ProcessTensorElementAccessor::join(const VVPI &i_list,  
                                   MPS_Matrix & M, const MPS_Matrix & M2){
//  time_point time1=now();
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
          //BLAS works but happens to be slower:
//          Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 0, Eigen::OuterStride<> >(      \
   tmp.mem+(d1*M2.dim_d1*tmp.dim_i+i)*tmp.dim_d2+d2*M2.dim_d2,M2.dim_d1,M2.dim_d2,Eigen::OuterStride<>(tmp.dim_i*tmp.dim_d2)).noalias()  +=  \
            M(i_list[i][j].first,d1,d2) *                    \
            Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 0, Eigen::OuterStride<> >(      \
   M2.mem+i_list[i][j].second*M2.dim_d2,M2.dim_d1,M2.dim_d2,Eigen::OuterStride<>(M2.dim_i*M2.dim_d2)); 
        }
      }
    }
  } 
//  time_point time2=now();
//  std::cout<<"runtime for join: "<<time_diff(time2-time1)<<"ms"<<std::endl;
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

  //contract B-A-B by first contracting C = B- -B:
  MPS_Matrix C(M_R.dim_i*M_L.dim_i, M_L.dim_d1, M_R.dim_d2);
  C.set_zero();
  for(int iR=0; iR<M_R.dim_i; iR++){
    for(int iL=0; iL<M_L.dim_i; iL++){
      for(int d1=0; d1<M_L.dim_d1; d1++){
        for(int d=0; d<M_L.dim_d2; d++){
          for(int d2=0; d2<M_R.dim_d2; d2++){
            C(iR*M_L.dim_i+iL, d1, d2)+=M_L(iL, d1, d)*M_R(iR, d, d2);
          }
        }
      }
    }
  }   

  VVTI i_list=join_symmetric_indices(acc_L, acc_R);

  MPS_Matrix tmp(i_list.size(), dim_d1_orig*M_L.dim_d1, dim_d2_orig*M_R.dim_d2);
  tmp.set_zero();
  for(int i=0; i<(int)i_list.size(); i++){
    for(int j=0; j<(int)i_list[i].size(); j++){
      for(int d1=0; d1<dim_d1_orig; d1++){
        for(int d2=0; d2<dim_d2_orig; d2++){
          for(int od1=0; od1<M_L.dim_d1; od1++){
            for(int od2=0; od2<M_R.dim_d2; od2++){
              tmp(i, d1*M_L.dim_d1 + od1, d2*M_R.dim_d2 + od2) +=
                M(std::get<1>(i_list[i][j]), d1, d2) * 
                C(std::get<0>(i_list[i][j])*M_L.dim_i+std::get<2>(i_list[i][j]), od1, od2);
            }
          }
        }
      }
    }
  }
      
  M.swap(tmp);   
/*
  VVPI i_list=join_thissecond_indices(acc_L);
  join(i_list, M, M_L);  

  VVPI i_list2=join_thisfirst_indices(acc_R);
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
*/
} 

void ProcessTensorElementAccessor::pass_on_before_join_symmetric(
         const ProcessTensorElementAccessor & acc_L,
         const ProcessTensorElementAccessor & acc_R,
         MPS_Matrix & M, const MPS_Matrix & M_L, const MPS_Matrix & M_R,
         PassOn &pass_on){
 
//  std::cout<<"TEST: P: "<<pass_on.P.rows()<<","<<pass_on.P.cols()<<std::endl; 

  if(pass_on.P.cols()!=M_L.dim_d1*M.dim_d1){
    std::cerr<<"PTEA::pass_on_before_join_symmetric : pass_on.P.cols()=!M_L.dim_d1*M.dim_d1!"<<std::endl;
    throw DummyException();
  }

  int dim_d1_orig=M.dim_d1;
  int dim_d2_orig=M.dim_d2;

  //contract B-A-B-P by first contracting C = B- -B and X = A-P-
  MPS_Matrix C(M_R.dim_i*M_L.dim_i, M_L.dim_d1, M_R.dim_d2);
  C.set_zero();
  for(int iR=0; iR<M_R.dim_i; iR++){
    for(int iL=0; iL<M_L.dim_i; iL++){
      for(int d1=0; d1<M_L.dim_d1; d1++){
        for(int d=0; d<M_L.dim_d2; d++){
          for(int d2=0; d2<M_R.dim_d2; d2++){
            C(iR*M_L.dim_i+iL, d1, d2)+=M_L(iL, d1, d)*M_R(iR, d, d2);
          }
        }
      }
    }
  }  

  MPS_Matrix X(M.dim_i*M_L.dim_d1, pass_on.P.rows(), M.dim_d2);
  X.set_zero();
  for(int i=0; i<M.dim_i; i++){
    for(int od1=0; od1<M_L.dim_d1; od1++){
      for(int d1=0; d1<M.dim_d1; d1++){
        for(int d2=0; d2<M.dim_d2; d2++){
          for(int k=0; k<pass_on.P.rows(); k++){ 
            X(i*M_L.dim_d1+od1, k, d2) += pass_on.P(k, d1*M_L.dim_d1+od1) \
                                        * M(i, d1, d2);
          }
        }
      }
    }
  }

  VVTI i_list=join_symmetric_indices(acc_L, acc_R);

  MPS_Matrix tmp(i_list.size(), pass_on.P.rows(), dim_d2_orig*M_R.dim_d2);
  tmp.set_zero();
  for(int i=0; i<(int)i_list.size(); i++){
    for(int j=0; j<(int)i_list[i].size(); j++){
      for(int k=0; k<pass_on.P.rows(); k++){
        for(int d2=0; d2<dim_d2_orig; d2++){
          for(int od1=0; od1<M_L.dim_d1; od1++){
            for(int od2=0; od2<M_R.dim_d2; od2++){
              tmp(i, k, d2*M_R.dim_d2 + od2) +=
                  X(std::get<1>(i_list[i][j])*M_L.dim_d1+od1, k, d2) * 
                  C(std::get<0>(i_list[i][j])*M_L.dim_i+std::get<2>(i_list[i][j]), od1, od2);
            }
          }
        }
      }
    }
  }
  M.swap(tmp);  

/*
  VVTI i_list=join_symmetric_indices(acc_L, acc_R);

  MPS_Matrix tmp(i_list.size(), dim_d1_orig*M_L.dim_d1, dim_d2_orig*M_R.dim_d2);
  tmp.set_zero();
  for(int i=0; i<(int)i_list.size(); i++){
    for(int j=0; j<(int)i_list[i].size(); j++){
      for(int d1=0; d1<dim_d1_orig; d1++){
        for(int d2=0; d2<dim_d2_orig; d2++){
          for(int od1=0; od1<M_L.dim_d1; od1++){
            for(int od2=0; od2<M_R.dim_d2; od2++){
              tmp(i, d1*M_L.dim_d1 + od1, d2*M_R.dim_d2 + od2) +=
                M(std::get<1>(i_list[i][j]), d1, d2) *
                C(std::get<0>(i_list[i][j])*M_L.dim_i+std::get<2>(i_list[i][j]), od1, od2);
            }
          }
        }
      }
    }
  }

  M.swap(tmp);
  M.inner_multiply_left(pass_on.P);
*/
}


void ProcessTensorElementAccessor::join_select_indices(
               const VVPI &i_list, MPS_Matrix &M1, const MPS_Matrix &M2,
               const ProcessTensorElementAccessor & acc_other,
               const SelectIndices & k_left, const SelectIndices & k_right){

  constexpr bool debug=false;
  if(debug){std::cout<<"join_select_indices: i_list.size()="<<i_list.size()<<" k_left.size()="<<k_left.size()<<" k_right.size()="<<k_right.size()<<std::endl;}

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
  if(debug){std::cout<<"join_select_indices: done."<<std::endl;}
}
void ProcessTensorElementAccessor::join_select_indices_alternate(
               int n, MPS_Matrix &M1, const MPS_Matrix &M2,
               const ProcessTensorElementAccessor & acc_other,
               const SelectIndices & k_left, const SelectIndices & k_right){

  constexpr bool debug=false;
  if(debug){std::cout<<"join_select_alternate: k_left.size()="<<k_left.size()<<" k_right.size()="<<k_right.size()<<std::endl;}

  VVPI i_list;
  if(n%2==0){
    i_list=join_thissecond_indices(acc_other);
  }else{
    i_list=join_thisfirst_indices(acc_other);
  }
  join_select_indices(i_list, M1, M2, acc_other, k_left, k_right);
  if(debug){std::cout<<"join_select_alternate: done."<<std::endl;}
}


void ProcessTensorElementAccessor::read_binary(std::istream &ifs){
  dict.read_binary(ifs);
}
void ProcessTensorElementAccessor::write_binary(std::ostream &ofs)const{
  dict.write_binary(ofs);
}


}//namespace
