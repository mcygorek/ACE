#include "ACE.hpp"
#include "Simulation_PT.hpp"
#include "ProcessTensorForwardList.hpp"
#include "ProcessTensorBuffer.hpp"
#include "Timings.hpp"

using namespace ACE;

//Setup: E1-S1-S2 -> E2-S2 
//via Qtilde ~ Q M_{S1} e^{L_{S1S2}\Delta t} 

int main(int args, char **argv){
  Parameters param(args, argv, true);
 
  Eigen::MatrixXcd initial=InitialState(param);
  int D1=initial.rows(); int DL1=D1*D1;
  if(D1<2){
    std::cerr<<"initial.rows()<2!"<<std::endl;
    exit(1);
  }

  TimeGrid tgrid(param);

  //Process coupling term:
  Parameters param_coupling;   // define coupling propagator 
  param_coupling.add_from_prefix("coupling", param);
  FreePropagator coupling_prop(param_coupling);

  //For dictionary reduction: assume time-independent coupling H:
//  std::cout<<"Coupling: const_H:"<<std::endl<<coupling_prop.const_H<<std::endl;
  coupling_prop.update(tgrid.ta,tgrid.dt);
  const Eigen::MatrixXcd & ES1S2 = coupling_prop.M;
  int DL2=ES1S2.rows()/DL1; int D2=sqrt(DL2);
  std::cout<<"S1 dim: "<<D1<<", S2 dim: "<<D2<<std::endl;
  //How many outer bonds?
  IF_OD_Dictionary dict(D2);
  std::vector<Eigen::MatrixXcd> coupling_mats;
  double dict_zero=param.get_as_double("dict_zero",0);
  for(int i1=0; i1<D2; i1++)for(int i2=0; i2<D2; i2++){ int i=i2*D2+i1;
    for(int j1=0; j1<D2; j1++)for(int j2=0; j2<D2; j2++){ int j=j2*D2+j1;
      Eigen::MatrixXcd thismat=Eigen::MatrixXcd::Zero(DL1,DL1);
      for(int k1=0; k1<D1; k1++)for(int k2=0; k2<D1; k2++){ int k=k2*D1+k1;
        for(int l1=0; l1<D1; l1++)for(int l2=0; l2<D1; l2++){ int l=l2*D1+l1;
          thismat(k,l)=ES1S2(((k1*D2+i1)*D1+k2)*D2+i2,((l1*D2+j1)*D1+l2)*D2+j2);
        }
      }
      if(dict_zero<=0){
        coupling_mats.push_back(thismat);
      }else{
        if(thismat.norm()<dict_zero){
          dict.beta[i*DL2+j]=-1;
          continue;
        }
        int oselect=-1;
        for(int o=0; o<(int)coupling_mats.size(); o++){
          if( (coupling_mats[o]-thismat).norm()<dict_zero ){
            oselect=o; break;
          }
        }
        if(oselect<0){
          dict.beta[i*DL2+j]=(int)coupling_mats.size();
          coupling_mats.push_back(thismat);
        }else{
          dict.beta[i*DL2+j]=oselect;
        }
      }
    }
  }
  dict.calculate_reduced_dim();

  //check:
  Eigen::MatrixXcd check=Eigen::MatrixXcd::Zero(DL1*DL2,DL1*DL2);
  for(int i1=0; i1<D2; i1++){for(int i2=0; i2<D2; i2++){ int i=i2*D2+i1;
    for(int j1=0; j1<D2; j1++){for(int j2=0; j2<D2; j2++){ int j=j2*D2+j1;
      int I=dict.beta[i*DL2+j];
      if(I<0)continue;
      for(int k1=0; k1<D1; k1++)for(int k2=0; k2<D1; k2++){ int k=k2*D1+k1;
        for(int l1=0; l1<D1; l1++)for(int l2=0; l2<D1; l2++){ int l=l2*D1+l1;
          check(((k1*D2+i1)*D1+k2)*D2+i2,((l1*D2+j1)*D1+l2)*D2+j2)+=coupling_mats[I](k,l);
        }
      }
    }}
  }}
  std::cout<<"(check-ES1S2).norm()="<<(check-ES1S2).norm()<<std::endl;
  std::cout<<"dict.reduced_dim="<<dict.reduced_dim<<" (total: "<<DL2*DL2<<")"<<std::endl;

 try{
  //Now: put S1 dynamics into PT-MPO form. Reuse functions of Simulation_PT
  FreePropagator fprop(param,D1);
  Parameters param_PT=param; param_PT.erase("write_PT");
  ProcessTensorForwardList PT(param_PT, D1);
  Simulation_PT sim(param);
  Output_Ops ops(param,false);

  std::string write_PT=param.get_as_string_check("write_PT");
  int buffer_blocksize=param.get_as_int("buffer_blocksize",-1);
//  ProcessTensorBufferSpec PTBspec(write_PT, buffer_blocksize);
  ProcessTensorBuffer PTB(ProcessTensorBufferSpec(write_PT, buffer_blocksize));
  PTB.resize(tgrid.n_tot);
  double forward_threshold=param.get_as_double("forward_threshold", -1.);
  double backward_threshold=param.get_as_double("backward_threshold", -1.);
  TruncatedSVD forward_trunc(forward_threshold);
  TruncatedSVD backward_trunc(backward_threshold);

  int DI2=DL1;
  PassOn pass_on(1);
  PT.reset();
  for(int n=0; n<tgrid.n_tot; n++){
    std::cout<<"n="<<n<<"/"<<tgrid.n_tot<<std::endl;
    ProcessTensorElement &e=PTB.get(n);
    e.accessor.dict=dict;
    int DI1=DI2;
    int DE1=DI1/DL1; 
    //for every beta: propagate S1 using fuctions of sim; 
    //then multiply with the corresponding coupling operator
    time_point time1=now();
    for(int o=0; o<(int)dict.reduced_dim; o++){
      for(int d1=0; d1<DI1; d1++){
        Eigen::MatrixXcd state=Eigen::MatrixXcd::Zero(DL1, DE1);
        state(d1/DE1, d1%DE1)=1.;
        if(sim.propagate_alternate && n%2==1){
          sim.propagate_state(state, n, tgrid, fprop, PT);
          state=coupling_mats[o]*state;
        }else{
          state=coupling_mats[o]*state;
          sim.propagate_state(state, n, tgrid, fprop, PT);
        }
        int DE2=state.cols();
        DI2=DL1*DE2;  
        if(o==0&&d1==0){
          e.M=MPS_Matrix(dict.reduced_dim, DI1, DI2);
        }
        for(int d2=0; d2<DI2; d2++){
          e.M(o, d1, d2)=state(d2/DE2, d2%DE2);
        }
      }
    }
    time_point time2=now();
    std::cout<<"runtime for construction of M: "<<time_diff(time2-time1)<<"ms"<<std::endl;
    
    e.closure=Eigen::VectorXcd::Zero(DI2);
    int DE2=DI2/DL1;
std::cout<<"DI2="<<DI2<<" DL1="<<DL1<<" DE2="<<DE2<<std::endl;
    for(int d2=0; d2<DE2; d2++){
      Eigen::MatrixXcd closure_state=Eigen::MatrixXcd::Zero(1, DE2);
      closure_state(0, d2)=1.;
      Eigen::VectorXcd rho_reduced=PT.get_rho_reduced(closure_state);
      for(int i=0; i<D1; i++){
        e.closure((i*D1+i)*DE2+d2)=rho_reduced(0);
      }
    }
    for(size_t oo=0; oo<ops.size(); oo++){
      Eigen::VectorXcd ops_vec=Eigen::VectorXcd::Zero(DL1*DE2);
      for(int i=0; i<D1; i++){
        for(int j=0; j<D1; j++){
          for(int d2=0; d2<DE2; d2++){
            ops_vec((i*D1+j)*DE2+d2)+=ops[oo](j,i)*e.closure(d2);
          }
        }
      }
      e.env_ops.ops.push_back(ops_vec);
    }

    if(n==0){
      Eigen::MatrixXcd ini(1,DL1);
      ini.row(0)=H_Matrix_to_L_Vector(initial);
      e.M.inner_multiply_left(ini);
    } 
    if(n==tgrid.n_tot-1){
      e.close_off();
    }
//    e.M.print_HR(std::string("M")+int_to_string(n)+std::string(".dat"),param.get_as_double("print_thr"));

    if(forward_threshold>=0.){
      time1=now();
      std::cout<<"dim_d2: "<<e.M.dim_d2<<std::flush;
      e.sweep_forward(forward_trunc, pass_on, (n==tgrid.n_tot-1));
      time_point time2=now();
      std::cout<<" -> "<<e.M.dim_d2;
      std::cout<<" (runtime: "<<time_diff(time2-time1)<<"ms)"<<std::endl;
    }

    PT.load_next();
  }
  
  if(backward_threshold>=0.){
    PTB.sweep_backward(backward_trunc,1);
  }
 
 }catch (DummyException &e){
  return 1;
 }
#ifdef EIGEN_USE_MKL_ALL
  mkl_free_buffers();
#endif
  return 0;
}
