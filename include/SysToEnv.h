#ifndef SYS_TO_ENV_DEFINED_H
#define SYS_TO_ENV_DEFINED_H
#include "InfluenceFunctional_OD.hpp"

namespace ACE{

class SysToEnv{
public:

Smart_Ptr<InfluenceFunctional_OD> IF;

MPS_Matrix Z_to_MPS_Matrix(int NS2, int NS1, const Eigen::MatrixXcd &M){
    int NL1=NS1*NS1;
    int NL2=NS2*NS2;

    MPS_Matrix a(NL2*NL2,NL1,NL1);
    for(int nu1=0; nu1<NS2; nu1++){
      for(int mu1=0; mu1<NS2; mu1++){
        for(int nu2=0; nu2<NS2; nu2++){
          for(int mu2=0; mu2<NS2; mu2++){
            for(int xi1=0; xi1<NS1; xi1++){
              for(int xi_1=0; xi_1<NS1; xi_1++){
                for(int xi2=0; xi2<NS1; xi2++){
                  for(int xi_2=0; xi_2<NS1; xi_2++){
a((nu2*NS2+mu2)*NL2+(nu1*NS2+mu1), xi1*NS1+xi_1, xi2*NS1+xi_2)=
//M((nu1*NS1+xi1)*NS1*NS2+(mu1*NS1+xi_1),(nu2*NS1+xi2)*NS1*NS2+(mu2*NS1+xi_2));
M((nu2*NS1+xi2)*NS1*NS2+(mu2*NS1+xi_2),(nu1*NS1+xi1)*NS1*NS2+(mu1*NS1+xi_1));
                  }
                }
              }
            }
          }
        }
      }
    }
    return a;
}


void calculate(Parameters &param){

  param.complain_if_not_specified("SysToEnv_Op");
  Eigen::MatrixXcd SysToEnv_Op=param.get_as_operator("SysToEnv_Op");  


  Eigen::MatrixXcd initial_rho=InitialState(param);
  FreePropagator fprop(param);
  FreePropagator Zprop; Zprop.add_Hamiltonian(SysToEnv_Op);

 
  //get and check dimensions
  int NS1=fprop.get_dim();
  if(NS1!=initial_rho.cols() || NS1!=initial_rho.rows()){
    std::cerr<<"NS1!=initial_rho.cols() || NS1!=initial_rho.rows()!"<<std::endl;
    exit(1);
  }
  if(SysToEnv_Op.rows()<1){std::cerr<<"SysToEnv_op: dimension < 1!"<<std::endl;exit(1);}
  int NS2=SysToEnv_Op.rows()/NS1;
  if(NS2*NS1!=SysToEnv_Op.rows()){
    std::cerr<<"NS2*NS1!=SysToEnv_Op.rows()!"<<std::endl; 
    exit(1);
  }
 // int NL1=NS1*NS1;
  int NL2=NS2*NS2;
  


  std::vector<Smart_Ptr<InfluenceFunctional_OD> > IF_E=IF_from_Parameters(param);
  RankCompressor_Ptr compressor=RankCompressor_Selector(param);
param.get_as_bool("IF_print_timesteps",false);

  TimeGrid tgrid(param);
  IF->tgrid=tgrid;
  if(IF->tgrid.n_tot!=IF_E[0]->tgrid.n_tot){
    std::cerr<<"IF->tgrid.n_tot!=IF_E[0]->tgrid.n_tot!"<<std::endl;
    exit(1);
  }

  bool print_timesteps=param.get_as_bool("IF_print_timesteps",false);
  bool do_sweep=param.get_as_bool("do_sweep",true);
  int keep_weight=0;



  // Build PT as:  Tr_S ( Z/2 M/2 Q M/2 Z/2 )^n   (*rho_0)
  //First, initialize without environment (use only Z)  (for fixing dictionary)
  IF->a.resize(IF->tgrid.n_tot);
  for(int n=0; n<(int)IF->a.size(); n++){ 
//    fprop.update(IF->tgrid.get_t(n), IF->tgrid.dt);
    Zprop.update(IF->tgrid.get_t(n), IF->tgrid.dt);
    IF->a[n]=Z_to_MPS_Matrix(NS2, NS1, Zprop.M);
  }
  {
    MPS_Matrix M(IF->a[0].dim_i, 1, IF->a[0].dim_d2);
    M.fill(0.);
    for(int beta=0; beta<IF->a[0].dim_i; beta++){
      for(int d2=0; d2<IF->a[0].dim_d2; d2++){
        for(int xi1=0; xi1<NS1; xi1++){
          for(int xi2=0; xi2<NS1; xi2++){
            M(beta, 0, d2)+=initial_rho(xi1,xi2)*IF->a[0](beta,xi1*NS1+xi2,d2);
          }
        }
      }
    }
    IF->a[0].swap(M);
  }
  {
    MPS_Matrix M(IF->a.back().dim_i, IF->a.back().dim_d1, 1);
    M.fill(0.);
    for(int beta=0; beta<IF->a.back().dim_i; beta++){
      for(int d1=0; d1<IF->a.back().dim_d1; d1++){
        for(int xi1=0; xi1<NS1; xi1++){
          M(beta, d1, 0)+=IF->a[0](beta,d1,xi1*NS1+xi1);
        }
      }
    }
    IF->a.back().swap(M);
  }

  //fixing dictionary 
  IF->calculate_dict(param.get_as_double("dict_zero",1e-20));
  IF->reduce_to_dict();
  std::vector<std::vector<int> > rev=IF->dict.get_reverse_beta();


  const InfluenceFunctional_OD & OIF=IF_E[0].ref();
  std::cout<<"SysToEnv: inner dict: ";OIF.dict.print_beta();std::cout<<std::endl;
  std::cout<<"SysToEnv: outer dict: ";IF->dict.print_beta();std::cout<<std::endl;

  //now:  Tr_S ( M/2 Z/2 Q Z/2 M/2 )^n   (*rho_0):
//  Zprop.precalculated=false;
  for(int n=0; n<(int)IF->a.size(); n++){
    fprop.update(IF->tgrid.get_t(n), IF->tgrid.dt/2);   
    Zprop.update(IF->tgrid.get_t(n), IF->tgrid.dt/2);  
    MPS_Matrix Zmat=Z_to_MPS_Matrix(NS2, NS1, Zprop.M);
    IF->dict.reduce_MPS_Matrix(Zmat);

    MPS_Matrix ZM(Zmat.dim_i, fprop.M.cols(), Zmat.dim_d2); ZM.fill(0.);
    for(int i=0; i<Zmat.dim_i; i++){
      for(int d1=0; d1<Zmat.dim_d1; d1++){
        for(int d2=0; d2<Zmat.dim_d2; d2++){
          for(int d0=0; d0<fprop.M.cols(); d0++){
            ZM(i, d0, d2)+=Zmat(i, d1, d2)*fprop.M(d1,d0);
          }
        }
      }
    }

// tmp = Q * Z/2 * M/2
    MPS_Matrix tmp(ZM.dim_i, ZM.dim_d1*OIF.a[n].dim_d1, ZM.dim_d2*OIF.a[n].dim_d2);
    tmp.fill(0.);
    for(int i=0; i<ZM.dim_i; i++){
      for(int d1=0; d1<ZM.dim_d1; d1++){
        for(int d2=0; d2<ZM.dim_d2; d2++){
          for(int d3=0; d3<ZM.dim_d2; d3++){
            int beta=OIF.dict.beta[d3*ZM.dim_d2+d2]; if(beta<0)continue;
            for(int od1=0; od1<OIF.a[n].dim_d1; od1++){
              for(int od2=0; od2<OIF.a[n].dim_d2; od2++){
                tmp(i, d1*OIF.a[n].dim_d1+od1, d3*OIF.a[n].dim_d2+od2)+=
                                     OIF.a[n](beta, od1, od2)*ZM(i, d1,d2);
              }
            }
          }
        }
      }
    }
   
//    IF->a[n].swap(tmp);

// M/2 * Z/2 * tmp: 
    fprop.update(IF->tgrid.get_t(n)+IF->tgrid.dt/2, IF->tgrid.dt/2);   
    Zprop.update(IF->tgrid.get_t(n)+IF->tgrid.dt/2, IF->tgrid.dt/2);  
    Zmat=Z_to_MPS_Matrix(NS2, NS1, Zprop.M);
    IF->dict.reduce_MPS_Matrix(Zmat);

    MPS_Matrix MZ(Zmat.dim_i, Zmat.dim_d1, fprop.M.rows()); MZ.fill(0.);
    for(int i=0; i<Zmat.dim_i; i++){
      for(int d1=0; d1<Zmat.dim_d1; d1++){
        for(int d2=0; d2<Zmat.dim_d2; d2++){
          for(int d3=0; d3<fprop.M.rows(); d3++){
            MZ(i, d1, d3)+=fprop.M(d3,d2)*Zmat(i, d1, d2);
          }
        }
      }
    }


    IF->a[n].resize(ZM.dim_i, MZ.dim_d1*OIF.a[n].dim_d1, MZ.dim_d2*OIF.a[n].dim_d2);
    IF->a[n].fill(0.);
    for(int l1=0; l1<NL2; l1++){
      for(int l2=0; l2<NL2; l2++){
        int beta1=IF->dict.beta[l2*NL2+l1]; if(beta1<0)continue;
        for(int l3=0; l3<NL2; l3++){
          int beta2=IF->dict.beta[l3*NL2+l2]; if(beta2<0)continue;
          int beta_new=IF->dict.beta[l3*NL2+l1]; if(beta_new<0)continue;
//don't add to same beta_new multiple times:
          if(rev[beta_new][0]!=l3*NL2+l1)continue;

          for(int d1=0; d1<Zmat.dim_d1; d1++){
            for(int d2=0; d2<Zmat.dim_d2; d2++){
              for(int d3=0; d3<Zmat.dim_d1; d3++){
                for(int od1=0; od1<OIF.a[n].dim_d1; od1++){
                  for(int od2=0; od2<OIF.a[n].dim_d2; od2++){
IF->a[n](beta_new, d1*OIF.a[n].dim_d1+od1, d3*OIF.a[n].dim_d2+od2)+=
MZ(beta2, d2, d3)*tmp(beta1, d1*OIF.a[n].dim_d1+od1, d2*OIF.a[n].dim_d2+od2);
                  }
                }
              }
            }
          }
        }
      }
    }
    if(n==0){
      MPS_Matrix M(IF->a[0].dim_i, 1, IF->a[0].dim_d2);
      M.fill(0.);
      for(int beta=0; beta<IF->a[0].dim_i; beta++){
        for(int d2=0; d2<IF->a[0].dim_d2; d2++){
          for(int xi1=0; xi1<NS1; xi1++){
            for(int xi2=0; xi2<NS1; xi2++){
              M(beta, 0, d2)+=initial_rho(xi1,xi2)*IF->a[0](beta,xi1*NS1+xi2,d2);
            }
          }  
        }
      }
      IF->a[0].swap(M);
    }
    if(do_sweep && n>0){
      if(print_timesteps){
        if(n==1){std::cout<<"Sweep forward: "<<n<<std::flush;}
        else{
          std::stringstream ss_last; ss_last<<n-1;
          for(int i=ss_last.str().length(); i>0; i--)std::cout<<'\b';
          std::cout<<n<<std::flush;
        }
        if(n==(int)IF->a.size()-1)std::cout<<std::endl;
      }

      Eigen::MatrixXcd R;
      IF->sweep_block_low_to_high(n-1, compressor.ref(), keep_weight, &R);
//      if(tgrid.use_rep)cta.set_R_if_correct_n(R,n-1);
//      for(size_t i=0; i<env_ops[n-1].size(); i++){
//        env_ops[n-1][i]=R*env_ops[n-1][i];
//      }
    }
  }
  {
    MPS_Matrix M(IF->a.back().dim_i, IF->a.back().dim_d1, 1);
    M.fill(0.);
    for(int beta=0; beta<IF->a.back().dim_i; beta++){
      for(int d1=0; d1<IF->a.back().dim_d1; d1++){
        for(int xi1=0; xi1<NS1; xi1++){
          M(beta, d1, 0)+=IF->a.back()(beta,d1,xi1*NS1+xi1);
        }
      }
    }
    IF->a.back().swap(M);
  }
  if(do_sweep){
   for(int n=(int)IF->a.size()-1; n>1; n--){
      if(print_timesteps){
        if(n==(int)IF->a.size()-1){std::cout<<"Sweep backward: "<<n<<std::flush;}
        else{
          std::stringstream ss_last; ss_last<<n+1;
          for(int i=ss_last.str().length(); i>0; i--)std::cout<<"\b \b";
          std::cout<<n<<std::flush;
        }
        if(n==2)std::cout<<std::endl;
      }

      Eigen::MatrixXcd L;
      IF->sweep_block_high_to_low(n, compressor.ref(), keep_weight, &L);
//      if(tgrid.use_rep)cta.set_L_if_correct_n(L,n);
//      Eigen::VectorXcd vL(L.cols());
//      for(int c=0; c<L.cols(); c++)vL(c)=L.col(c).norm();
//      Eigen::MatrixXcd Linv=L.adjoint();
//      for(int c=0; c<Linv.rows(); c++)Linv.row(c)/=(vL(c)*vL(c));

//      for(size_t i=0; i<env_ops[n-1].size(); i++){
//        env_ops[n-1][i]=Linv*env_ops[n-1][i];
//      }
    }
  }
}

operator Smart_Ptr<InfluenceFunctional_OD> (){
  return IF;
}
SysToEnv(Parameters &params){
  IF=new InfluenceFunctional_OD();
  calculate(params);
}

};

}//namespace
#endif
