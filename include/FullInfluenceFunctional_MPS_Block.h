#ifndef FULL_INFLUENCE_FUNCTIONAL_MPS_BLOCK_DEFINED_H
#define FULL_INFLUENCE_FUNCTIONAL_MPS_BLOCK_DEFINED_H

#include "FullInfluenceFunctional_MPS.h"


class FullInfluenceFunctional_MPS_Block: public FullInfluenceFunctional_MPS{
public:

  virtual void calculate(int n_max, double dt, double t_tot, DiagBB &diagBB, RankCompressor &compressor){
    int NL=diagBB.get_dim()*diagBB.get_dim();
    int n_tot=t_tot/dt;
    if(n_tot<n_max)n_tot=n_max;

    //find out smallest exponent so that n_tot<=2^ex-1
    int ex=0;
    int neff=1;
    while(n_tot> neff-1){
      ex++;
      neff*=2;
    }
    std::cout<<"input n_tot+1: "<<n_tot+1<<" next power of two: "<<neff<<" = 2^"<<ex<<std::endl;
    n_tot=neff-1;



    IF_Line_MPS line(n_max, dt, diagBB, n_tot);
    MPS::copy(line);
    calculate_closures();
   
//    print_dims();
    std::cout<<"rank: "<<get_rank()<<std::endl; 
    for(int i=0; i<get_rank(); i++){
      std::cout<<a[i].dim_d2<<" ";
    }std::cout<<std::endl;

    
    int shift=1;
    for(int step=0; step<ex; step++){
      std::cout<<"step: "<<step<<" shift: "<<shift<<std::endl;
      MPS line2;
      line2.a.resize(get_rank()-shift);
      line2.a[0].resize(a[shift].dim_i, 1, a[shift].dim_d2);
      for(int i=0; i<a[shift].dim_i; i++){
        for(int d2=0; d2<a[shift].dim_d2; d2++){
          line2.a[0](i,0,d2)=c[shift](i,d2);
        }
      }
      for(int i=1; i<line2.get_rank(); i++){
        line2.a[i].copy(a[shift+i]); 
      }

      
      if(false && line2.a[0].dim_d2>line2.a[0].dim_i){ //otherwise it doesnt help
        std::cout<<"inner dimensions of line2 before sweep:"<<std::endl;
        for(int i=0; i<line2.get_rank(); i++){
          std::cout<<line2.a[i].dim_d2<<" ";
        }std::cout<<std::endl;
        for(size_t i=1; i<line2.a.size(); i++){
          line2.sweep_block_left(i-1, compressor);
        }
        for(int i=line2.a.size()-1; i>1; i--){
          sweep_block_right(i, compressor);
        }
      }

      std::cout<<"inner dimensions of line2:"<<std::endl;
      for(int i=0; i<line2.get_rank(); i++){
        std::cout<<line2.a[i].dim_d2<<" ";
      }std::cout<<std::endl;


      MPS::multiply_front_SVD(line2, compressor);


      shift*=2;
      calculate_closures();
 
//      print_dims(); 
      std::cout<<"inner dimensions:"<<std::endl;
      for(int i=0; i<get_rank(); i++){
        std::cout<<a[i].dim_d2<<" ";
      }std::cout<<std::endl;
    }


    calculate_closures();
  }

  FullInfluenceFunctional_MPS_Block(int n_max, double dt, double t_tot, DiagBB &diagBB, RankCompressor &compressor){
    calculate(n_max, dt, t_tot, diagBB, compressor);
  }
  FullInfluenceFunctional_MPS_Block(const std::string &filename){
    read_binary(filename);
  }
  FullInfluenceFunctional_MPS_Block(){
  }
};


#endif
