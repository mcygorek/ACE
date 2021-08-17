#ifndef INFLUENCE_FUNCTIONAL_D_DEFINED_H
#define INFLUENCE_FUNCTIONAL_D_DEFINED_H

#include "MPS.h"
#include "IF_Line_MPS_low_to_high.h"

class InfluenceFunctional_D: public MPS{
public:
  IF_TimeGrid tgrid;

  //closures:
  std::vector<Eigen::VectorXcd> c;
 
  //groups
  Coupling_Groups_Liouville lgroups; 


  void calculate_closures(){
    int ZI=0; //zero index, must be a sojourn
    if(a.size()<1){
      std::cerr<<"InfluenceFunctional_OD::calculate_closures: a.size()<1!"<<std::endl;
      exit(1);
    }

    int NL=a[0].dim_i;
    int N=sqrt(NL);

    if(NL!=a[0].dim_i){
      std::cerr<<"InfluenceFunctional_OD::calculate_closures: NL!=a[0].dim_i!"<<std::endl;
      exit(1);
    }
    if(N*N!=NL){
      std::cerr<<"InfluenceFunctional_OD::calculate_closures: N*N!=NL!"<<std::endl;
      exit(1);
    }

    c.resize(a.size());
    c[c.size()-1].resize(1);
    c[c.size()-1](0)=1;

    for(int n=(int)c.size()-2; n>=0; n--){
      c[n]=Eigen::VectorXcd::Zero(a[n].dim_d2);
      for(int d1=0; d1<a[n+1].dim_d1; d1++){
        for(int d2=0; d2<a[n+1].dim_d2; d2++){
          c[n](d1)+=c[n+1](d2)*a[n+1]( ZI ,d1, d2);
        }
      }
    }
  }


  void calculate_diagBB(DiagBB &diagBB, RankCompressor &compressor){
    lgroups=Coupling_Groups_Liouville(diagBB.groups);
    int NL=diagBB.get_dim()*diagBB.get_dim();
    double dt=tgrid.dt;
    int n_tot=tgrid.n_tot;
    int n_max=tgrid.n_mem; //tgrid.n_calc;

    std::cout<<"InfluenceFunctional_D: Calculating single line: "<<n_max<<" "<<n_tot<<std::endl;
    IF_Line_MPS_low_to_high line(n_max, dt, diagBB, n_tot);
    std::cout<<"Done."<<std::endl;
    MPS::copy(line);
   

    for(int n=1; n<n_tot; n++){
      std::cout<<"calculate_diagBB: "<<n<<"/"<<n_tot<<std::endl;
      std::cout<<"max_dim: "<<get_max_dim()<<std::endl;
//      int n_cut=n_max;
//      if(n_tot-n < n_max)n_cut=n_tot-n;

      line.a.pop_back();
      MPS_Matrix ma(line.a.back().dim_i, line.a.back().dim_d1, 1);
      for(int i=0; i<ma.dim_i; i++){
        for(int d1=0; d1<ma.dim_d1; d1++){
          ma(i,d1,0)=0.;
          for(int d2=0; d2<line.a.back().dim_d2; d2++){
            ma(i,d1,0)+=line.a.back()(i,d1,d2);
          }
        }
      }
      line.a.back().swap(ma);
//std::cout<<"a.size(): "<<a.size()<<" line.a.size(): "<<line.a.size()<<std::endl;
//      IF_Line_MPS_low_to_high line2(n_cut, dt, diagBB, n_tot-n);

      MPS::high_end_multiply_and_compress(line, compressor, sqrt(NL));

//      print_dims(); 
    }

    calculate_closures();
  }

  virtual void read_binary(const std::string &filename){
    MPS::read_binary(filename);
    tgrid.set_default(MPS::get_rank()-1);
    calculate_closures();
  }

  void set_none(int n_max, int N){
    int NL=N*N;
    a.resize(n_max);
    for(int i=0; i<n_max; i++){  
      a[i].resize(NL, 1, 1);
      for(int j=0; j<NL; j++){
        a[i](j,0,0)=1;
      }
    }
    tgrid.set_default(n_max);
    calculate_closures();
  }


  InfluenceFunctional_D(const IF_TimeGrid &tgr, DiagBB &diagBB, RankCompressor &compressor){
    tgrid=tgr;
    calculate_diagBB(diagBB, compressor);
  }


  InfluenceFunctional_D(const std::string &filename){
    read_binary(filename);
  }

  InfluenceFunctional_D(){
  }
};


#endif
