#include "RandomizedCompression.hpp"
#include "otimes.hpp"
#include <random>

namespace ACE{

   
void Randomized_Combine(ProcessTensorBuffer & PTB, \
                        ProcessTensorBuffer & PTB2, \
                        int chi_new, const TruncatedSVD &trunc){
  constexpr bool debug=false; //true;
if(debug){std::cout<<"Randomize_Combine: begin"<<std::endl;}

  std::mt19937 mt;
  std::random_device rd;
  std::seed_seq seq{rd()};
  mt.seed(seq);
  std::normal_distribution<double> dist(0., 1.);

  if(PTB.n_tot!=PTB2.n_tot){
    std::cerr<<"Randomized_Combine: PTB.n_tot!=PTB2.n_tot!"<<std::endl;
    throw DummyException();
  }

  //Reuse PTB structre for C: d1: bond~PTB, d2: bond~PTB2, i: chi_new
  //-> Note that this introduces a mismatch between i and dict
  //-> Also, let's start at n=1 to make loop easier
  ProcessTensorBuffer C;
  C.set_new_temporary(PTB);
  C.resize(PTB.n_tot);
  Eigen::VectorXd kept=Eigen::VectorXd::Zero(PTB.n_tot); //weight of C kept in place
  {
   ProcessTensorElement &c=C.get(0, ForwardPreload);
   c.M.resize(chi_new, 1, 1);
   for(int p=0; p<chi_new; p++){ 
     c.M(p,0,0)=sqrt(1./chi_new);//std::complex<double>( dist(mt), dist(mt) ) ;
   }
   kept(0)=sqrt((double)chi_new);
  }
  for(int n=1; n<PTB.n_tot; n++){
if(debug){std::cout<<"Randomize_Combine: composing C: n="<<n<<std::endl;}
    bool this_second=((n-1)%2==1);
    ProcessTensorElement &e=PTB.get(n-1, ForwardPreload);
    ProcessTensorElement &e2=PTB2.get(n-1, ForwardPreload);
    
    ProcessTensorElementAccessor tmp_acc=e.accessor;
    ProcessTensorElementAccessor::VVPI i_list=
       this_second ? tmp_acc.join_thissecond_indices(e2.accessor) :
                     tmp_acc.join_thisfirst_indices(e2.accessor);
 
    if(n-1==0 && (PTB.get(0).M.dim_d1!=1 || PTB2.get(0).M.dim_d1!=1)){
      std::cerr<<"get(0).M.dim_d1!=1!"<<std::endl;
      throw DummyException();
    }   

    // contract c_last with e 
    ProcessTensorElement &c_last=C.get(n-1, ForwardPreload);
    MPS_Matrix c_last_e(chi_new*e.M.dim_i, e.M.dim_d2, e2.M.dim_d1);
    c_last_e.set_zero();
    for(int p=0; p<chi_new; p++){
      for(int i=0; i<e.M.dim_i; i++){
        for(int e2d1=0; e2d1<e2.M.dim_d1; e2d1++){
          for(int e1d1=0; e1d1<e.M.dim_d1; e1d1++){
            for(int e1d2=0; e1d2<e.M.dim_d2; e1d2++){
              c_last_e(p*e.M.dim_i+i, e1d2, e2d1) += \
                c_last.M(p, e1d1, e2d1) * e.M(i, e1d1, e1d2);
            }
          }
        }
      }
    }
    // contract with e2
    MPS_Matrix c_last_e_e2(chi_new*i_list.size(), e.M.dim_d2, e2.M.dim_d2);
    c_last_e_e2.set_zero();
    for(int p=0; p<chi_new; p++){
      for(int i=0; i<(int)i_list.size(); i++){ 
        for(int j=0; j<(int)i_list[i].size(); j++){
          for(int e1d2=0; e1d2<e.M.dim_d2; e1d2++){
            for(int e2d1=0; e2d1<e2.M.dim_d1; e2d1++){
              for(int e2d2=0; e2d2<e2.M.dim_d2; e2d2++){
                c_last_e_e2(p*i_list.size()+i, e1d2, e2d2) += \
                  c_last_e(p*e.M.dim_i+i_list[i][j].first, e1d2, e2d1) * \
                  e2.M(i_list[i][j].second, e2d1, e2d2);
              }
            }
          }
        }
      }
    }
    // Random matrix:
    Eigen::MatrixXcd Omega(i_list.size(), chi_new);
    for(int i=0; i<i_list.size(); i++){
      for(int p=0; p<chi_new; p++){
        Omega(i,p)=std::complex<double>( dist(mt), dist(mt) );
      }
    }

    ProcessTensorElement &c=C.get(n, ForwardPreload);
    c.M.resize(chi_new, e.M.dim_d2, e2.M.dim_d2);
    c.M.set_zero();
    for(int p=0; p<chi_new; p++){
      for(int i=0; i<(int)i_list.size(); i++){ 
        for(int e1d2=0; e1d2<e.M.dim_d2; e1d2++){
          for(int e2d2=0; e2d2<e2.M.dim_d2; e2d2++){
            c.M(p, e1d2, e2d2) += Omega(i,p) * \
                                   c_last_e_e2(p*i_list.size()+i, e1d2, e2d2);
          }
        }
      }
    }
    //update "kept"
    for(int i=0; i<c.M.dim_i; i++){for(int d1=0; d1<c.M.dim_d1; d1++){for(int d2=0; d2<c.M.dim_d2; d2++){
          kept(n)+=std::norm( c.M(i, d1, d2) ); //C++ function for SQUARED norm
    }}}
    kept(n)=sqrt(kept(n));
if(debug){std::cout<<"kept("<<n<<")="<<kept(n)<<std::endl;}
    double kept_inv=1./kept(n);
    for(int i=0; i<c.M.dim_i; i++){for(int d1=0; d1<c.M.dim_d1; d1++){for(int d2=0; d2<c.M.dim_d2; d2++){
      c.M(i, d1, d2)*=kept_inv;
    }}}
 
 
  }

  if(PTB.get(PTB.n_tot-1, BackwardPreload).M.dim_d2!=1 || PTB2.get(PTB.n_tot-1, BackwardPreload).M.dim_d2!=1){
    std::cerr<<"get(n, BackwardPreload).M.dim_d2!=1!"<<std::endl;
    throw DummyException();
  }


  MPS_Matrix last_trail; 
  last_trail.resize(1,1,1); last_trail.set_zero(); last_trail(0,0,0)=1.;

  for(int n=PTB.n_tot-1; n>0; n--){
if(debug){std::cout<<"Randomize_Combine: constructing e: n="<<n<<std::endl;}
    bool this_second=(n%2==1);
    ProcessTensorElement &e=PTB.get(n, BackwardPreload);
    ProcessTensorElement &e2=PTB2.get(n, BackwardPreload);
    
    //note: this also modifies dict:
    ProcessTensorElementAccessor::VVPI i_list=
       this_second ? e.accessor.join_thissecond_indices(e2.accessor) :
                     e.accessor.join_thisfirst_indices(e2.accessor);

    MPS_Matrix tmp1(e.M.dim_i*last_trail.dim_i, e.M.dim_d1, e2.M.dim_d2);
    tmp1.set_zero();
    for(int i=0; i<e.M.dim_i; i++){
      for(int p=0; p<last_trail.dim_i; p++){
        for(int e1d1=0; e1d1<e.M.dim_d1; e1d1++){
          for(int e1d2=0; e1d2<e.M.dim_d2; e1d2++){
            for(int e2d2=0; e2d2<e2.M.dim_d2; e2d2++){
              tmp1(i*last_trail.dim_i+p, e1d1, e2d2)+= \
                e.M(i, e1d1, e1d2) * last_trail(p, e1d2, e2d2);
            }
          }
        }
      }
    }
    Eigen::MatrixXcd B=Eigen::MatrixXcd::Zero(e.M.dim_d1*e2.M.dim_d1,\
                                              i_list.size()*last_trail.dim_i);
    for(int i=0; i<(int)i_list.size(); i++){ 
      for(int j=0; j<(int)i_list[i].size(); j++){
        for(int p=0; p<last_trail.dim_i; p++){
          for(int e1d1=0; e1d1<e.M.dim_d1; e1d1++){
            for(int e2d1=0; e2d1<e2.M.dim_d1; e2d1++){
              for(int e2d2=0; e2d2<e2.M.dim_d2; e2d2++){
                B(e1d1*e2.M.dim_d1+e2d1, i*last_trail.dim_i+p) +=
                    tmp1(i_list[i][j].first*last_trail.dim_i+p, e1d1, e2d2)*\
                    e2.M(i_list[i][j].second, e2d1, e2d2);
              }
            }
          }
        }
      }
    }

     //Eigen::MatrixXcd A = B;
     //multiply c*B to get A:
     ProcessTensorElement &c = C.get(n, BackwardPreload);

     Eigen::MatrixXcd A=Eigen::MatrixXcd::Zero(c.M.dim_i, \
                                               i_list.size()*last_trail.dim_i);
     for(int p=0; p<c.M.dim_i; p++){
       for(int i=0; i<i_list.size(); i++){
         for(int p2=0; p2<last_trail.dim_i; p2++){
           for(int e1d1=0; e1d1<e.M.dim_d1; e1d1++){
             for(int e2d1=0; e2d1<e2.M.dim_d1; e2d1++){
               A(p, i*last_trail.dim_i + p2) += c.M(p, e1d1, e2d1) * \
                         B(e1d1*e2.M.dim_d1+e2d1, i*last_trail.dim_i+p2);
             }
           }
         }
       }
     }

     PassOn pass_on=trunc.compress_backward(A, e.backwardNF);
if(debug){std::cout<<"singular values: "<<e.backwardNF.transpose()<<std::endl;}
    Eigen::MatrixXcd Ainv=A.adjoint();
    double keep=e.backwardNF(0);
    if(trunc.keep>0)keep=trunc.keep;
    if(fabs(keep-1.)>1e-6){
      pass_on.P/=keep;
      pass_on.Pinv*=keep;
      A*=keep;
      Ainv*=1./keep;
    }

    Eigen::MatrixXcd P = B*Ainv;

    last_trail.resize(P.cols(), e.M.dim_d1, e2.M.dim_d1);
    last_trail.set_zero();
    for(int k=0; k<P.cols(); k++){
      for(int d1=0; d1<e.M.dim_d1; d1++){
        for(int d2=0; d2<e2.M.dim_d1; d2++){
          last_trail(k, d1, d2) = P(d1*e2.M.dim_d1+d2, k); 
        }
      }
    }

    e.M.set_from_Matrix_d1_id2(A, i_list.size());
  }
  {//first element
if(debug){std::cout<<"Randomize_Combine: constructing e: n="<<0<<std::endl;}
    bool this_second=(0%2==1);
    ProcessTensorElement &e=PTB.get(0, BackwardPreload);
    ProcessTensorElement &e2=PTB2.get(0, BackwardPreload);
    
    //note: this also modifies dict:
    ProcessTensorElementAccessor::VVPI i_list=
       this_second ? e.accessor.join_thissecond_indices(e2.accessor) :
                     e.accessor.join_thisfirst_indices(e2.accessor);

    MPS_Matrix tmp1(e.M.dim_i*last_trail.dim_i, e.M.dim_d1, e2.M.dim_d2);
    tmp1.set_zero();
    for(int i=0; i<e.M.dim_i; i++){
      for(int p=0; p<last_trail.dim_i; p++){
        for(int e1d1=0; e1d1<e.M.dim_d1; e1d1++){
          for(int e1d2=0; e1d2<e.M.dim_d2; e1d2++){
            for(int e2d2=0; e2d2<e2.M.dim_d2; e2d2++){
              tmp1(i*last_trail.dim_i+p, e1d1, e2d2)+= \
                e.M(i, e1d1, e1d2) * last_trail(p, e1d2, e2d2);
            }
          }
        }
      }
    }
    e.M.resize(i_list.size(), 1, last_trail.dim_i);
    e.M.set_zero();
    for(int i=0; i<(int)i_list.size(); i++){
      for(int j=0; j<(int)i_list[i].size(); j++){
        for(int p=0; p<last_trail.dim_i; p++){
          for(int e2d2=0; e2d2<e2.M.dim_d2; e2d2++){
            e.M(i, 0, p)+= \
                    tmp1(i_list[i][j].first*last_trail.dim_i+p, 0, e2d2)*\
                    e2.M(i_list[i][j].second, 0, e2d2);
          }
        }
      }
    }
  } 

  PTB.calculate_closures();
 
if(debug){std::cout<<"Randomize_Combine: done"<<std::endl;}
}
}//namespace
