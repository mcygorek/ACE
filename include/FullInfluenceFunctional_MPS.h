#ifndef FULL_INFLUENCE_FUNCTIONAL_MPS_DEFINED_H
#define FULL_INFLUENCE_FUNCTIONAL_MPS_DEFINED_H

#include "IF_Line_MPS.h"


class FullInfluenceFunctional_MPS: public MPS{
public:
  //closures:
  std::vector<Eigen::MatrixXcd> c;
 

  void calculate_closures(){
    int ZI=0; //zero index, must be a sojourn
    if(a.size()<1){
      std::cerr<<"FullInfluenceFunctional_MPS::calculate_closures: a.size()<1!"<<std::endl;
      exit(1);
    }

    c.resize(a.size());
    c[0].resize(a[0].dim_i, a[0].dim_d2);
    for(int i=0; i<a[0].dim_i; i++){
      for(int d2=0; d2<a[0].dim_d2; d2++){
        c[0](i,d2)=a[0](i,0,d2);
      }
    }
    for(size_t n=1; n<a.size(); n++){
      c[n]=Eigen::MatrixXcd::Zero(a[n].dim_i, a[n].dim_d2);
      for(int i=0; i<a[n].dim_i; i++){
        for(int d1=0; d1<a[n].dim_d1; d1++){
          for(int d2=0; d2<a[n].dim_d2; d2++){
            c[n](i,d2)+=c[n-1](ZI,d1)*a[n](i,d1,d2);
          }
        }
      }
    }
  }

  virtual void calculate(int n_max, double dt, double t_tot, DiagBB &diagBB, RankCompressor &compressor){
    int NL=diagBB.get_dim()*diagBB.get_dim();
    int n_tot=t_tot/dt;

    IF_Line_MPS line(n_max, dt, diagBB, n_tot);
    MPS::copy(line);
   
//    print_dims(); 
    for(int i=0; i<get_rank(); i++){
      std::cout<<a[i].dim_d2<<" ";
    }std::cout<<std::endl;

    for(int n=0; n<n_tot; n++){
      int n_cut=n_max;
      if(n_tot-1-n < n_max)n_cut=n_tot-1-n;

      IF_Line_MPS line2(n_cut, dt, diagBB, n_tot-1-n);

      MPS::multiply_front_SVD(line2, compressor);

//      print_dims(); 
      for(int i=0; i<get_rank(); i++){
        std::cout<<a[i].dim_d2<<" ";
      }std::cout<<std::endl;
    }
/*
    for(int n=n_max-1; n>=0; n--){
      IF_Line_MPS line2(n, dt, diagBB, n_tot-(n_max-n));
//      MPS::multiply_front(line2);
      MPS::multiply_front_SVD(line2, eps);

      print_dims(); 
    }
*/
    calculate_closures();
  }
  virtual void read_binary(const std::string &filename){
    MPS::read_binary(filename);
    calculate_closures();
  }
  void set_none(int n_max, int NL){
    a.resize(n_max+1);
    for(int i=0; i<n_max+1; i++){  
      a[i].resize_fill_one(NL, 1, 1);
    }
    calculate_closures();
  }

  FullInfluenceFunctional_MPS(int n_max, double dt, double t_tot, DiagBB &diagBB, RankCompressor &compressor){
    calculate(n_max, dt, t_tot, diagBB, compressor);
  }
  FullInfluenceFunctional_MPS(const std::string &filename){
    read_binary(filename);
  }
  FullInfluenceFunctional_MPS(){
  }
};


#endif
