#include "ACE.hpp"
#include "InfluenceFunctional.hpp"
#include "ProcessTensorBuffer.hpp"
#include "ProcessTensorForwardList.hpp"
#include "OutputPrinter.hpp"
#include "InitialState.hpp"

using namespace ACE;

Tensor_Dense get_dense(ProcessTensorBuffer &PTB){
  int n_tot=PTB.n_tot;
  if(n_tot<1){
    std::cerr<<"get_dense:: n_tot<1"<<std::endl;
    exit(1);
  }
  int NL=PTB.get(0).M.dim_i;
  if(n_tot==1){
    Tensor_Dense ten(n_tot, NL);
    ten.fill(0.);
      for(int i=0; i<NL; i++){ten[i]=PTB.get(0).M(i,0,0);}
    return ten;
  }

  std::vector<Tensor_Dense> ten(PTB.get(0).M.dim_d2,Tensor_Dense(1, NL));
  for(int d2=0; d2<ten.size(); d2++){
    for(int i=0; i<NL; i++){
      ten[d2][i]=PTB.get(0).M(i,0,d2);
    }
  }

  for(int n=1; n<n_tot-1; n++){
    std::vector<Tensor_Dense> ten2(PTB.get(n).M.dim_d2,Tensor_Dense(n+1, NL));
    for(int d2=0; d2<PTB.get(n).M.dim_d2; d2++){
      ten2[d2].fill(0.);
    }
    int BS=ten[0].get_total_size();
    for(int d1=0; d1<PTB.get(n).M.dim_d1; d1++){
      for(int d2=0; d2<PTB.get(n).M.dim_d2; d2++){
        for(int i=0; i<NL; i++){
          for(int b=0; b<BS; b++){
            ten2[d2][i*BS+b]+=ten[d1][b]*PTB.get(n).M(i,d1,d2);
          }
        }
      }
    }
    ten.swap(ten2);
  }

  Tensor_Dense result(n_tot, NL); result.fill(0.); 
  for(int d2=0; d2<ten.size(); d2++){
    int BS=ten[0].get_total_size();
    for(int i=0; i<NL; i++){
      for(int b=0; b<BS; b++){
        result[i*BS+b]+=ten[d2][b]*PTB.get(n_tot-1).M(i,d2,0);
//        result[b*NL+i]+=ten[d2][b]*PTB.get(n_tot-1).M(i,d2,0);
      }
    }
  }
  return result;
}

int main(int args, char** argv){

  Parameters param(args, argv, true);

  TimeGrid tgrid(param);
  int n_max=param.get_as_double("n_max", tgrid.n_mem);

  DiagBB diagBB(param, param.get_as_string("Gaussian_prefix","Boson"));

  InfluenceFunctional IF(n_max, tgrid.dt, diagBB);
  int N=IF.get_Ngrps();
  int NL=N*N;
  int BS=IF.ten.back().get_total_size();
  std::cout<<"n_max="<<n_max<<" NL="<<NL<<" BS="<<BS<<std::endl;
  ProcessTensorBuffer PTB;
  PTB.set_from_DiagBB_single_line(diagBB, tgrid.dt, n_max);

  ProcessTensorForwardList PTlist(param);
  
  Tensor_Dense IFbig(IF.ten[n_max-1]);
  for(int n=n_max-2; n>=0; n--){
    int this_tot=1; for(int i=0; i<=n; i++)this_tot*=NL;
    int other_tot=1; for(int i=n+1; i<n_max; i++)other_tot*=NL;
std::cout<<"n="<<n<<" this_tot="<<this_tot<<" other_tot="<<other_tot<<std::endl;
    for(int o=0; o<other_tot; o++){
      for(int b=0; b<this_tot; b++){
        IFbig[o*this_tot+b]*=IF.ten[n][b];
      }
    }
  }

//  Tensor_Dense PTten=get_dense(PTB);
  ProcessTensorBuffer *ptr=dynamic_cast<ProcessTensorBuffer*>(PTlist.list[0].get());
  Tensor_Dense PTten=get_dense(*ptr);

  double err=0.;
  for(size_t b=0; b<BS; b++){
//   std::cout<<IF.ten.back()[b].real()<<" "<<IF.ten.back()[b].imag();
//    std::cout<<"   "<<PTten[b].real()<<" "<<PTten[b].imag();
//    std::cout<<std::endl;
    std::cout<<IFbig[b].real()<<" "<<IFbig[b].imag();
    std::cout<<"   "<<PTten[b].real()<<" "<<PTten[b].imag();
    std::cout<<"   "<<std::abs(IFbig[b]-PTten[b]);
    std::cout<<std::endl;
    err+=std::abs(IFbig[b]-PTten[b])*std::abs(IFbig[b]-PTten[b]);
  }
  std::cout<<"Accumulated error: "<<sqrt(err)<<std::endl;

  FreePropagator prop(param);
  OutputPrinter printer(param);
  Eigen::MatrixXcd rho;
  {
    Tensor_Dense U(IFbig);
    prop.update(tgrid.get_t(0), tgrid.dt/2.);
    Eigen::VectorXcd rho1=prop.M*H_Matrix_to_L_Vector(InitialState(param));
    int BS=U.get_total_size()/NL;
    for(int b=0; b<BS; b++){
      for(int i=0; i<NL; i++){
        U[b*NL+i]*=rho1(i);
      }
    }
    for(int n=0; n<n_max-1; n++){
      std::cout<<"TEST: n="<<n<<std::endl;
      prop.update(tgrid.get_t(n)+tgrid.dt/2, tgrid.dt);
      Tensor_Dense U2(U.get_rank()-1, NL); U2.fill(0.);
      int BS=U2.get_total_size()/NL;
      for(int b=0; b<BS; b++){
        for(int i=0; i<NL; i++){
          for(int j=0; j<NL; j++){
            U2[b*NL+i]+=prop.M(i,j)*U[(b*NL+i)*NL+j];
          }
        }
      }
      U=U2;
    }
    prop.update(tgrid.get_t(n_max)-tgrid.dt/2, tgrid.dt/2);
    Eigen::VectorXcd rhoF=Eigen::VectorXcd::Zero(NL);
    for(int i=0; i<NL; i++){ 
      for(int j=0; j<NL; j++){ 
        rhoF(i)+=prop.M(i,j)*U[j];
      }
    }
    rho=L_Vector_to_H_Matrix(rhoF);
  }
  std::cout<<"rho:"<<std::endl<<rho<<std::endl;
  //printer.print(0, tgrid.get_t(0), init);
  

  return 0;
}
