#ifndef FULL_INFLUENCE_FUNCTIONAL_COLUMN_DEFINED_H
#define FULL_INFLUENCE_FUNCTIONAL_COLUMN_DEFINED_H

#include "FullInfluenceFunctional_MPS.h"

class FullInfluenceFunctional_Column: public FullInfluenceFunctional_MPS{
public:

  void get_column(int n_max, double dt, DiagBB &diagBB, std::vector<Eigen::MatrixXcd> &B, bool expand){
    int NL=diagBB.get_dim()*diagBB.get_dim();

    std::vector<Eigen::MatrixXcd> b(n_max+1);
    for(int n=0; n<n_max+1; n++)b[n]=diagBB.calculate_expS(n,dt);

    if(n_max==0){
      B.clear();
      if(!expand){
        B.resize(NL,Eigen::MatrixXcd::Zero(1,1));
        for(int a=0; a<NL; a++){
          B[a](0,0)=b[0](a,a);
        } 
      }else{
        B.resize(NL,Eigen::MatrixXcd::Zero(NL,1));
        for(int a=0; a<NL; a++){
          B[a](a,0)=b[0](a,a);
        }
      } 
      return;
    }


    Tensor_Dimensions dims(n_max, NL);
    int BS=dims.get_total_size();
    int BS1=BS;
    if(expand)BS1*=NL;
    std::cout<<"n_max: "<<n_max<<" BS: "<<BS1<<" "<<BS<<std::endl;

    B.clear();
    B.resize(NL,Eigen::MatrixXcd::Zero(BS1,BS));
    for(int a=0; a<NL; a++){
      for(Tensor_Index ind(dims); !ind.done(dims); ind.increment_from_back(dims)){
        std::complex<double> fac=b[0](a,a);
        for(int n=0; n<n_max; n++)fac*=b[n+1](a, ind[n]);
        Tensor_Index ind2(dims);
        ind2[0]=a;
        for(int n=1; n<n_max; n++){
          ind2[n]=ind[n-1];
        }
        if(expand){
          B[a](dims.get_block_index(ind2)*NL+ind.back(), dims.get_block_index(ind))=fac;
        }else{
          B[a](dims.get_block_index(ind2), dims.get_block_index(ind))=fac;
        }
      }
    }
  }


  virtual void calculate(int n_max, double dt, double t_tot, DiagBB &diagBB, RankCompressor &compressor){
    int NL=diagBB.get_dim()*diagBB.get_dim();
    int n_tot=t_tot/dt;
    if(n_tot<n_max)n_max=n_tot;
  
    MPS::resize_fill_one(n_tot+1, NL);
     
    for(int n=0; n<n_max; n++){
      std::vector<Eigen::MatrixXcd> B;
      get_column(n, dt, diagBB, B, true);
      a[n_tot-n].resize(NL, B[0].rows(), B[0].cols());
      a[n_tot-n].fill(0);
      for(int i=0; i<NL; i++){
        for(int d1=0; d1<B[i].rows(); d1++){
          for(int d2=0; d2<B[i].cols(); d2++){
//            a[n_tot-n](i,i*B[i].rows()+d1,d2)=B[i](d1,d2);
            a[n_tot-n](i,d1,d2)=B[i](d1,d2);
          }
        }
      }
    }
    {//n_max
      int n=n_max;
      std::vector<Eigen::MatrixXcd> B;
      get_column(n, dt, diagBB, B, false);
      a[n_tot-n].resize(NL, B[0].rows(), B[0].cols());
      a[n_tot-n].fill(0);
      for(int i=0; i<NL; i++){
        for(int d1=0; d1<B[i].rows(); d1++){
          for(int d2=0; d2<B[i].cols(); d2++){
            a[n_tot-n](i,d1,d2)=B[i](d1,d2);
          }
        }
      }
    }  
    for(int n=n_max+1; n<n_tot+1; n++){
      a[n_tot-n].copy(a[n_tot-n_max]);
    }
    MPS_Matrix c(NL, 1, a[0].dim_d2);
    c.fill(0.);
    for(int i=0; i<NL; i++){
      for(int d1=0; d1<a[0].dim_d1; d1++){
        for(int d2=0; d2<a[0].dim_d2; d2++){
          c(i,0,d2)+=a[0](i,d1,d2);
        }
      }
    }
    std::cout<<"c: "; c.print_dims(); std::cout<<std::endl;
    a[0].swap(c);

    print_dims();

    calculate_closures();
  }
  FullInfluenceFunctional_Column(int n_max, double dt, double t_tot, DiagBB &diagBB, RankCompressor &compressor){
    calculate(n_max, dt, t_tot, diagBB, compressor);
  }
  FullInfluenceFunctional_Column(){
  }
};


#endif
