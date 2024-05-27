#include "Rep_GIF.hpp"
#include "MPS_Matrix.hpp"
#include "IF_OD_Abstract.hpp"
#include "otimes.hpp"
#include "LiouvilleTools.hpp"
#include "ModePropagator.hpp"
#include "Compress_Trafo_At.hpp"

namespace ACE{

  void Rep_GIF::set_default(int N){
    init=Eigen::VectorXcd(1); init(0)=1.;
    M.resize(N*N*N*N, 1,1); M.fill(1.);
    dict.set_default(N);
    env_ops.resize(1);
    env_ops[0].resize(1); env_ops[0](0)=1;
  }

  void Rep_GIF::expand_ops(const ModePropagator &mprop){
    int dim_before=init.size();
    init=Vector_otimes(init, H_Matrix_to_L_Vector(mprop.get_bath_init()));

#ifdef DEBUG_REP
std::cout<<"mprop.env_ops.size(): "<<mprop.env_ops.size()<<std::endl;
std::cout<<"env_ops.size(): "<<env_ops.size()<<std::endl;
#endif
    env_ops[0]=Vector_otimes(env_ops[0], H_Matrix_to_L_Vector(Eigen::MatrixXcd::Identity(mprop.get_N_mode(), mprop.get_N_mode())));
    for(size_t i=0; i<mprop.env_ops.size(); i++){
      if(i+1<env_ops.size()){
        env_ops[i+1]=Vector_otimes(env_ops[i+1], H_Matrix_to_L_Vector(mprop.env_ops[i]));
        if(init.size()!=env_ops[i+1].size()){
          std::cerr<<"Error: Rep_GIF::expand_ops: init.size()!=env_ops["<<i+1<<"].size()!"<<std::endl;
          exit(1);
        }
      }else{
        if(dim_before==1){ 
          env_ops.push_back(H_Matrix_to_L_Vector(mprop.env_ops[i]));
          if(init.size()!=env_ops[i+1].size()){
            std::cerr<<"Error: Rep_GIF::expand_ops: init.size()!=env_ops["<<i+1<<"].size()!"<<std::endl;
            exit(1);
          }
        }else{
          std::cerr<<"Error: Rep_GIF::expand_ops: too many env_ops in mprop!"<<std::endl;
          exit(1);
        }
      }
    }
#ifdef DEBUG_REP
std::cout<<"env_ops.size(): "<<env_ops.size()<<std::endl;
#endif
//    ident=Vector_otimes(ident, H_Matrix_to_L_Vector(Eigen::MatrixXcd::Identity(mprop.get_N_mode(), mprop.get_N_mode())));
  }

  void Rep_GIF::apply_compress_trafo(Compress_Trafo_At &cta, bool low_high_low){
#ifdef DEBUG_COMPRESS_TRAFO
std::cout<<"Apply trafo to M!"<<std::endl;
std::cout<<"M: "<<M.dim_d1<<" "<<M.dim_d2;
std::cout<<" R: "<<cta.R.rows()<<" "<<cta.R.cols();
std::cout<<" L: "<<cta.L.rows()<<" "<<cta.L.cols()<<std::endl;
#endif
    cta.apply_trafo(M, low_high_low);
#ifdef DEBUG_COMPRESS_TRAFO
    std::cout<<"Apply to trafo successful"<<std::endl;
#endif
    cta.apply_to_init(init, low_high_low);
    for(size_t e=0; e<env_ops.size(); e++){
#ifdef DEBUG_COMPRESS_TRAFO
      std::cout<<"e: "<<e<<std::endl;
#endif
      cta.apply_to_op(env_ops[e], low_high_low);
#ifdef DEBUG_COMPRESS_TRAFO
      std::cout<<"done: "<<e<<std::endl;
#endif
    }
  }

  void Rep_GIF::regularize(){
//    int N=dict.get_N();
    int NL=dict.get_NL();
   
    //Trace after first step:
    Eigen::MatrixXcd A=Eigen::MatrixXcd::Zero(NL, NL);
    for(int d1=0; d1<M.dim_d1; d1++){
      for(int d2=0; d2<M.dim_d2; d2++){
        for(int i=0; i<NL; i++){
          for(int j=0; j<NL; j++){
            A(i,j)+=init(d1)*M(dict.beta[i*NL+j],d1,d2)*env_ops[0](d2);
          }
        }
      }
    }
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd( A, Eigen::ComputeFullU | Eigen::ComputeFullV);
    std::cout<<"SVDs: ";
    for(int s=0; s<svd.singularValues().size(); s++){
       std::cout<<svd.singularValues()(s)<<" ";
    }
    std::cout<<std::endl;
    std::cout<<"U.adjoint: "<<std::endl<<svd.matrixU().adjoint()<<std::endl;
    std::cout<<"V.adjoint: "<<std::endl<<svd.matrixV().adjoint()<<std::endl;
  }


}//namespace
