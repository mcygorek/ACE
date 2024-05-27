#ifndef CAVITYFY_DEFINED_H
#define CAVITYFY_DEFINED_H

#include "otimes.hpp"
#include "Operators_Boson.hpp"

namespace ACE{
/** Purpose: add a cavity to the system. Expand the IFs and the initia_rho
             correspondingly
*/

void cavityfy(FreePropagator &fprop, 
              std::vector<Smart_Ptr<InfluenceFunctional_OD> > &IFV, 
              Eigen::MatrixXcd &initial_rho, 
              Output_Ops &output_Op,
              Parameters &param){

  int n_max=param.get_as_size_t("add_system_cavity",0);
  if(n_max<1)return;
  
  int n_dim=n_max+1;
  int N=fprop.get_dim();
  int NL=N*N;
  int ML=n_dim*n_dim;
  Eigen::MatrixXcd id=Eigen::MatrixXcd::Identity(n_dim, n_dim);
  //FreePropagator:
  {
    fprop.const_H=otimes(fprop.const_H, id);

    for(size_t i=0; i<fprop.timedep_H.size(); i++){
      fprop.timedep_H[i].second=otimes(fprop.timedep_H[i].second,id);
    }

    if(fprop.nonH.rows()>0){
      Eigen::MatrixXcd nonH(NL*ML, NL*ML);
      for(int i=0; i<NL; i++){
        for(int j=0; j<NL; j++){
          for(int k=0; k<ML; k++){
            nonH(i*ML+k, j*ML+k)=fprop.nonH(i,j);
          }
        }
      }
      fprop.nonH=nonH;
    }
    fprop.precalculated=false;
    
    if(param.is_specified("system_cavity_coupling")){
      double g=param.get_as_double("system_cavity_coupling");
      if(N!=2){
        std::cerr<<"system_cavity_coupling only implemented for 2LSs!"<<std::endl; 
        exit(1);
      }
      fprop.const_H+=Constants::*g*
 (otimes(Operators2x2::sigma_plus(), Operators_Boson::a(n_dim))
 +otimes(Operators2x2::sigma_minus(), Operators_Boson::adagger(n_dim)));
    }
  }

  //Influence Functionals:
  for(size_t ifi=0; ifi<IFV.size(); ifi++){
    if(IFV[ifi]->dict.beta.size()!=NL*NL){
      std::cerr<<"Mismatch in dimensions between influence functional and system!"<<std::endl;
      exit(1);
    }
    IF_OD_Dictionary dict;
    dict.N=N*n_dim;
    dict.beta.clear();
    dict.beta.resize(NL*ML*NL*ML, -1);
    for(int i1=0; i1<N; i1++){
      for(int k1=0; k1<n_dim; k1++){
        for(int i2=0; i2<N; i2++){
          for(int k2=0; k2<n_dim; k2++){
            int alpha=((i1*n_dim+k1)*N+i2)*n_dim+k2;
            for(int j1=0; j1<N; j1++){
              for(int j2=0; j2<N; j2++){
                int tildealpha=((j1*n_dim+k1)*N+j2)*n_dim+k2;
                dict.beta[alpha*NL*ML+tildealpha]=
                  IFV[ifi]->dict.beta[((i1*N+i2)*N+j1)*N+j2];
              }
            }
          }
        }
      }
    }

    /*for(int i=0; i<NL; i++){
      for(int j=0; j<NL; j++){
        for(int k=0; k<ML; k++){
          dict.beta[(i*ML+k)*NL*ML+(j*ML+k)]=IFV[ifi]->dict.beta[i*NL+j];
        }
      }
    }*/
    dict.calculate_reduced_dim();
std::cout<<"cavity test dict before: ";IFV[ifi]->dict.print_beta(); std::cout<<std::endl;
    IFV[ifi]->dict=dict;
std::cout<<"cavity test dict after: ";IFV[ifi]->dict.print_beta(); std::cout<<std::endl;
    //IFV[ifi]->calculate_closures();
  }

  //Initial state:
  initial_rho=otimes(initial_rho, Operators_Boson::vacuum(n_dim));

  //Output:
  for(size_t o=0; o<output_Op.size(); o++){
    Eigen::MatrixXcd op=Eigen::MatrixXcd::Zero(N*n_dim, N*n_dim);
    for(int i=0; i<output_Op[o].rows(); i++){
      for(int j=0; j<output_Op[o].cols(); j++){
        for(int k=0; k<n_dim; k++){
          op(i*n_dim+k, j*n_dim+k)=output_Op[o](i,j);
        }
      } 
    }
    output_Op[o]=op;
  }


/*
  TODO: -system-cavity coupling
        -radiative decay
        -output: photon number
*/
}



}//namespace
#endif
