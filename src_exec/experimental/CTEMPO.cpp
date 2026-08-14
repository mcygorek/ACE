#include "DummyException.hpp"
#include "Parameters.hpp"
#include "ProcessTensorBuffer.hpp"
#include "FreePropagator.hpp"
#include "InitialState.hpp"
#include "OutputPrinter.hpp"
#include "Timings.hpp"
#include "TruncationLayout.hpp"
#include "Coupling_Groups.hpp"
#include "otimes.hpp"
#include "ProcessTensorForwardList.hpp"



using namespace ACE;

class Simulation_CTEMPO{
public:
  bool print_timesteps;
  std::string print_dims_file;
  int print_dims_step;

  Eigen::MatrixXcd run(Propagator &prop, DiagBB &diagBB,
           const Eigen::MatrixXcd & initial_rho, const TimeGrid &tgrid,
           OutputPrinter &printer, const TruncationLayout &trunc_layout){

    int N=initial_rho.rows();
    int NL=N*N;
    if(N<2){
      std::cerr<<"Simulation_CTEMPO::run: N<2!"<<std::endl;
      exit(1);
    }
    if(prop.get_dim() != N){
      std::cerr<<"Simulation_CTEMPO: prop.get_dim() != N!"<<std::endl;
      throw DummyException();
    }

    int n_tot=tgrid.n_tot;
    if(n_tot<1){
      std::cerr<<"Simulation_CTEMPO::run: tgrid.n_tot<1!"<<std::endl;
      exit(1);
    }
    int n_mem=tgrid.n_mem;
    if(n_mem<2){
      std::cerr<<"Simulation_CTEMPO::run: tgrid.n_mem<1!"<<std::endl;
      exit(1);
    }

    int NdiagBB=diagBB.get_dim();
    int NLdiagBB=NdiagBB*NdiagBB;
    std::cout<<"diagBB.sys_dim()="<<diagBB.sys_dim()<<" ";
    std::cout<<"diagBB.get_dim()="<<diagBB.get_dim()<<std::endl;
    if(diagBB.sys_dim() != N){
      std::cerr<<"Simulation_CTEMPO::run: diagBB.sys_dim() != N!"<<std::endl;
      exit(1);
    }


    //row of Influence Functional
    std::vector<Eigen::MatrixXcd> b(n_mem);
    for(int n=0; n<n_mem; n++){
      b[n]=diagBB.calculate_expS(n, tgrid.dt);
    }
    Coupling_Groups_Liouville lgroups(diagBB.groups);
    int Ngrps2=lgroups.get_Ngrps();
    const std::vector<int> & grp2=lgroups.grp;


    Eigen::VectorXcd rho_reduced=H_Matrix_to_L_Vector(initial_rho);
    printer.print(0, tgrid.get_t(0), rho_reduced);

    ProcessTensorBuffer PTB; 
    //PTB.set_trivial(n_mem, N);
    {
      // work with dimension diagBB.get_dim(). 
      // Only at the end, expand using dict.expand_DiagBB(diagBB);
      // Expansion may only be needed just before merging with sys. prop.      
      ProcessTensorElement e;
      e.set_trivial(NdiagBB);
      e.accessor.dict.set_default_diag(NdiagBB); 
      e.M.set_from_Matrix_d1i_d2(Eigen::VectorXcd::Ones(NLdiagBB), NLdiagBB);
      e.env_ops.clear();
      PTB.resize(n_mem, e); 
    }
    { // set initial state
      ProcessTensorElement & e = PTB.get(0);
      e.accessor.dict.set_default_diag(N); // full set of outer bonds
      e.M.set_from_Matrix_d1i_d2(rho_reduced, NL); 
    }

    for(int n=0; n<tgrid.n_tot; n++){
      if(print_timesteps){
        std::cout<<"step: "<<n<<"/"<<tgrid.n_tot<<std::endl;
      } 
      int cur_length=n_mem;
      int new_length=n_mem;   //length of new MPO without PT
      if(n+n_mem > tgrid.n_tot){ 
        cur_length = tgrid.n_tot-n+1; 
        new_length = tgrid.n_tot-n; 
      }
      //std::cout<<"n="<<n<<" n_tot="<<tgrid.n_tot<<" cur_length="<<cur_length<<" new_length="<<new_length<<std::endl;

      //TruncatedSVD trunc = trunc_layout.get_base();
      TruncatedSVD trunc_fwd=trunc_layout.get_forward(n, tgrid.n_tot);
      trunc_fwd.keep = diagBB.get_dim();
      trunc_fwd.print_info(); std::cout<<std::endl;

      ProcessTensorElement & e0 = PTB.get(0, ForwardPreload);
      const ProcessTensorElement & e1 = PTB.peek(1);
      // propagate first element:
      ProcessTensorElement Q0;
      Q0.set_trivial(N);
      Q0.accessor.dict.set_default_diag(N);
      Q0.M.resize(NL,1, e0.M.dim_d2);
      Q0.closure = e0.closure;
      Q0.env_ops.clear();// = e0.env_ops;
      Q0.M.set_zero();
      
      prop.update(tgrid.get_t(n), tgrid.dt/2);
      for(int d0=0; d0<e0.M.dim_d2; d0++){
        for(int a1=0; a1<NL; a1++){ 
          for(int a0=0; a0<NL; a0++){ 
            Q0.M(a1, 0, d0) += prop.M(a1, a0) * e0.M(a0, 0, d0);
          }
        }
      }
      
      ProcessTensorElement Qtmp;
      Qtmp.set_trivial(N);
      Qtmp.accessor.dict.set_default_diag(N);
      Qtmp.M.resize(NL,1, e1.M.dim_d2);
      Qtmp.closure = e1.closure;
      Qtmp.env_ops.clear(); // = e1.env_ops;
      Qtmp.M.set_zero();
      for(int a1=0; a1<NL; a1++){ 
        int g1=grp2[a1];
        for(int d1=0; d1<e1.M.dim_d2; d1++){
          for(int d0=0; d0<e1.M.dim_d1; d0++){
            Qtmp.M(a1, 0, d1) += e1.M(g1, d0, d1) * Q0.M(a1, 0, d0);
          }
        }
      }      
      
      prop.update(tgrid.get_t(n)+tgrid.dt/2, tgrid.dt/2);

      Q0.M.resize(NL, 1, Qtmp.M.dim_d2*Ngrps2);
      Q0.M.set_zero();
      for(int d1=0; d1<e1.M.dim_d2; d1++){
        for(int a1=0; a1<NL; a1++){ 
          int g1=grp2[a1];
          for(int a_=0; a_<NL; a_++){ 
            Q0.M(a_, 0, d1*Ngrps2+g1) += prop.M(a_, a1) * b[0](g1,g1) * Qtmp.M(a1, 0, d1);
          } 
        }
      }
      Q0.closure = Vector_otimes(Qtmp.closure, Eigen::VectorXcd::Ones(Ngrps2));
      e0.swap(Q0);
     
      for(int l=1; l<cur_length-1; l++){
        ProcessTensorElement & e = PTB.get(l, ForwardPreload);
        const ProcessTensorElement & e1 = PTB.peek(l+1);
        e.M.resize(e1.M.dim_i, e1.M.dim_d1*Ngrps2, e1.M.dim_d2*Ngrps2);
        e.M.set_zero();
        for(int i=0; i<Ngrps2; i++){
          for(int beta=0; beta<Ngrps2; beta++){
            for(int d1=0; d1<e1.M.dim_d1; d1++){
              for(int d2=0; d2<e1.M.dim_d2; d2++){
                e.M(i, d1*Ngrps2+beta, d2*Ngrps2+beta) = b[l](i,beta)*e1.M(i, d1, d2);
              }
            }
          } 
        }
        e.closure = Vector_otimes(e1.closure, Eigen::VectorXcd::Ones(Ngrps2));
      }

      if(new_length==cur_length){ //have to expand
        ProcessTensorElement & e = PTB.get(new_length-1);
        e.M.resize(e.M.dim_i, Ngrps2, 1);
        e.M.set_zero();
        for(int i=0; i<Ngrps2; i++){
          for(int beta=0; beta<Ngrps2; beta++){
            e.M(i, beta, 0) = b[new_length-1](i, beta);
          }
        }
      }else{
        PTB.get(new_length-1).close_off();
      }

      //std::cout<<"Dims before compression (i, d1, d2): "; for(int l=0; l<new_length; l++){PTB.get(l, ForwardPreload).M.print_dims();std::cout<<"; ";if(l>=9){std::cout<<"..."; break;}} std::cout<<std::endl;
      
      PTB.sweep_forward(trunc_fwd, 1, 0, new_length);

      TruncatedSVD trunc_bwd=trunc_layout.get_backward(n, tgrid.n_tot);
      trunc_bwd.keep = diagBB.get_dim();
      trunc_bwd.print_info(); std::cout<<std::endl;
      PTB.sweep_backward(trunc_bwd, 1, 0, new_length);

      //std::cout<<"Dims after compression (i, d1, d2): "; for(int l=0; l<new_length; l++){PTB.get(l, ForwardPreload).M.print_dims();std::cout<<"; "; if(l>=9){std::cout<<"..."; break;}} std::cout<<std::endl;
      if(print_dims_file!="" && n==print_dims_step){
        std::ofstream ofs(print_dims_file);
        for(int l=0; l<new_length; l++){
          ofs<<PTB.get(l, ForwardPreload).M.dim_d1<<" " \
             <<PTB.get(l, ForwardPreload).M.dim_d2<<std::endl;
        }
      }

      //extract rho_reduced:
      rho_reduced = e0.M.get_Matrix_d1i_d2()*e0.closure;

      //print output:
      double t_next=tgrid.get_t(n+1);
      printer.print(n+1, t_next, rho_reduced);

    }

    return L_Vector_to_H_Matrix(rho_reduced);
  }


  void setup(Parameters &param){
    print_timesteps=param.get_as_bool("print_timesteps",true);
    print_dims_file=param.get_as_string("print_dims_file","");
    print_dims_step=param.get_as_int("print_dims_step",-1);
  }

  Simulation_CTEMPO(Parameters &param){
    setup(param);
  }
  Simulation_CTEMPO(){
    Parameters param;
    setup(param);
  }

};

int main(int args, char** argv){
 try{
  Parameters param(args, argv, true);

  TimeGrid tgrid(param);
  InitialState initial(param);
  FreePropagator fprop(param);
  if(!fprop.dim_fixed){
    fprop.set_dim(initial.rho.rows());
  }else{
    if(fprop.get_dim()!=initial.rho.rows()){
      std::cerr<<"Mismatch in dimensions between initial state and propagator!"<<std::endl;
      exit(1);
    }
  }

  std::string prefix=param.get_as_string("prefix","Boson");
  DiagBB diagBB; diagBB.setup(param,prefix);  
  OutputPrinter printer(param);
  TruncationLayout trunc_layout(param);
  Simulation_CTEMPO sim(param);
  sim.run(fprop, diagBB, initial.rho, tgrid, printer, trunc_layout);

/*
  GenericSimulation sim(param);

  std::cout<<"Setting up Process Tensor..."<<std::endl;
  Parameters paramPT=TransferTensor::get_paramPT(param);
 
  ProcessTensorForwardList PT;
  std::cout<<"Propagating..."<<std::endl;
  time_point time1=now();
  sim.run(tgrid, PT);
  time_point time2=now();
  std::cout<<"runtime for propagation: "<<time_diff(time2-time1)<<"ms"<<std::endl;
*/

 }catch (DummyException &e){
  return 1;
 }

#ifdef EIGEN_USE_MKL_ALL
  mkl_free_buffers();
#endif

  return 0;
}

