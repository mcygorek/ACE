#ifndef SIMULATION_OD_DEFINED_H
#define SIMULATION_OD_DEFINED_H

#include "IF_from_Parameters.h"
#include "Propagator.h"
#include "Simulation_Results.h"
#include "Simulation.h"
#include "InitialState.h"
#include "FT_Output.h"
#include "Which_Env_Ops.h"
#include "Cavityfy.h"

#include "ModePropagatorGenerator.h"
#include "Pulse_Printer.h"


class Simulation_OD{
public:

  Eigen::MatrixXcd rho;
  bool print_timestep;


  Simulation_Results results;
  Output_Ops output_Op;
  Which_Env_Ops_List which_env_ops;
  FT_Output FTO;
  

  void add_output_Op(const Eigen::MatrixXcd &op){
    output_Op.add(op);
  }
  void setup_output(Parameters &param){
    output_Op.setup(param);
    FTO.setup(param);
    which_env_ops.setup(param);
  }

  void run_nobath(Propagator &prop,
    double ta, double dt, double te, const Eigen::MatrixXcd &initial_rho){

std::cout<<"RUN_NOBATH!"<<std::endl;

    int N=initial_rho.rows();
    int NL=N*N;
    int n_max=(te-ta)/dt;
    if(n_max<1){
      rho=initial_rho;
      return;
    }

    Eigen::MatrixXcd state(NL, 1);
    for(int i=0; i<N; i++){
      for(int j=0; j<N; j++){
        state(i*N+j,0)=initial_rho(i,j);
      }
    }

    results.clear();
    results.resize(n_max+1);
    results.set(0, ta, output_Op, initial_rho);

    for(int n=0; n<n_max; n++){
      double t=ta+n*dt;
      if(print_timestep)std::cout<<"Step: "<<n<<"/"<<n_max<<" ("<<t<<")"<<std::endl;
      prop.update(t,dt);

      //do the actual propagation: first only system-free
      {
        Eigen::MatrixXcd state2=Eigen::MatrixXcd::Zero(NL,state.cols());
        for(int d1=0; d1<state.cols(); d1++){
          for(int j=0; j<NL; j++){
            for(int k=0; k<NL; k++){
              state2(j,d1)+=prop.M(j,k)*state(k,d1);
            }
          }
        }
        state=state2;
      }
      rho=Eigen::MatrixXcd::Zero(N,N);
      for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
          rho(i,j) =  state(i*N+j, 0);
        }
      }
      results.set(n+1, ta+(n+1)*dt, output_Op, rho);
    }
  }

  void run(Propagator &prop, 
    std::vector<Smart_Ptr<InfluenceFunctional_OD> > & IFV,
    double ta, double dt, double te, const Eigen::MatrixXcd &initial_rho){
 
    if(initial_rho.rows()!=initial_rho.cols()){
      std::cerr<<"Simulation_OD::calculate: initial_rho.rows()!=initial_rho.cols()!"<<std::endl;
      exit(1);
    } 
    if(initial_rho.rows()!=prop.get_dim()){
      std::cerr<<"Simulation_OD::calculate: initial_rho.rows()!=prop.get_dim()! ["<<initial_rho.rows()<<" vs. "<<prop.get_dim()<<"]"<<std::endl;
      exit(1);
    }

    if(IFV.size()<1){
      std::cout<<"Simulation_OD::run: No InfluenceFunctional specified!"<<std::endl;
      run_nobath(prop, ta, dt, te, initial_rho);
      return;
    }   


    int N=initial_rho.rows();
    int NL=N*N;
    int Ngrps2=NL;

    int n_max=(te-ta)/dt;
    if(n_max<1){
      rho=initial_rho;
      return; 
    }


    for(size_t ifi=0; ifi<IFV.size(); ifi++){
      std::cout<<"ifi: "<<ifi<<" "<<IFV[ifi]->get_rank()<<std::endl;
      IFV[ifi]->check_within_limits(n_max);
      IFV[ifi]->check_consistency();

      if(Ngrps2*Ngrps2!=IFV[ifi]->dict.beta.size()){
        std::cerr<<"Simulation_OD::calculate: IF.dict.beta.size()!=NL*NL!"<<std::endl;
        std::cerr<<"IFV["<<ifi<<"]->dict.beta.size(): "<<IFV[ifi]->dict.beta.size()<<" NL: "<<NL<<std::endl;
        std::cerr<<"Please check Influence Functional or initial rho!"<<std::endl;
        exit(1);
      }
    }


    //initialize
    Eigen::MatrixXcd state(NL, 1);
    for(int i=0; i<N; i++){
      for(int j=0; j<N; j++){
        state(i*N+j,0)=initial_rho(i,j);
      }
    }
//std::cout<<"Initial state: "<<std::endl<<initial_rho<<std::endl;

    results.clear();
    results.resize(n_max+1);
    results.set(0, ta, output_Op, initial_rho);

    for(int n=0; n<n_max; n++){
      double t=ta+n*dt;
      if(print_timestep)std::cout<<"Step: "<<n<<"/"<<n_max<<" ("<<t<<")"<<std::endl;
      prop.update(t,dt);

      //do the actual propagation: first only system-free
      {
        Eigen::MatrixXcd state2=Eigen::MatrixXcd::Zero(NL,state.cols());
        for(int d1=0; d1<state.cols(); d1++){
          for(int j=0; j<NL; j++){
            for(int k=0; k<NL; k++){
              state2(j,d1)+=prop.M(j,k)*state(k,d1);
            }
          }
        }
        state=state2;
      }


      //next: account for influence functionals:  
      //strategy: take global d1, decompose into lower, middle and upper block
      //where the middle block belongs to IVF[ifi].
      //Then, multiply the middle block an pack the lower and upper blocks back.       
      int d1_tot=1;
      for(size_t ifi=0; ifi<IFV.size(); ifi++)d1_tot*=IFV[ifi]->get_a(n).dim_d1;

      int d1_lower=d1_tot;
      for(size_t ifi=0; ifi<IFV.size(); ifi++){
        int d1_this=IFV[ifi]->get_a(n).dim_d1;
        int d2_this=IFV[ifi]->get_a(n).dim_d2;
        d1_lower/=d1_this;

        int d2_tot=(d1_tot/d1_this)*d2_this;
        Eigen::MatrixXcd state2=Eigen::MatrixXcd::Zero(NL,d2_tot);
        
        for(int i=0; i<NL; i++){
          for(int j=0; j<NL; j++){
            int i_ind=IFV[ifi]->dict.beta[i*Ngrps2+j];
            if(i_ind<0)continue;

            for(int glob_d1=0; glob_d1<d1_tot; glob_d1++){
              for(int d2=0; d2<IFV[ifi]->get_a(n).dim_d2; d2++){
                int d1=(glob_d1/d1_lower)%d1_this;
                int up=(glob_d1/d1_lower)/d1_this;
                int glob_d2=(up*d2_this+d2)*d1_lower+glob_d1%d1_lower;

                state2(i,glob_d2)+=IFV[ifi]->get_a(n)(i_ind,d1,d2)*state(j,glob_d1);
              }
            }
          }
        }
        state=state2;
        d1_tot/=d1_this; d1_tot*=d2_this;
      }


      //extract reduced system density matrix
      int d2_tot=1;
      for(size_t ifi=0; ifi<IFV.size(); ifi++)d2_tot*=IFV[ifi]->get_a(n).dim_d2;

      Eigen::MatrixXcd c_state=state;
      for(size_t ifi=0; ifi<IFV.size(); ifi++){
        int d2_this=IFV[ifi]->get_a(n).dim_d2;
        int d2_factor=c_state.cols()/d2_this;

        Eigen::MatrixXcd c_state2=Eigen::MatrixXcd::Zero(NL, d2_factor);
        for(int j=0; j<NL; j++){
          for(int d2=0; d2<d2_this; d2++){
            for(int o=0; o<d2_factor; o++){
              c_state2(j, o)+=IFV[ifi]->get_c(n)(d2)*c_state(j, d2*d2_factor+o);
            }
          }
        }

        c_state=c_state2;       
      }
// * / Eigen::MatrixXcd c_state=state;

      rho=Eigen::MatrixXcd::Zero(N,N);
      for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
          rho(i,j) =  c_state(i*N+j, 0);
        }
      }
      results.set(n+1, ta+(n+1)*dt, output_Op, rho);


      //Environment operators
      for(size_t w=0; w<which_env_ops.size(); w++){
        if(which_env_ops[w].i>=IFV.size()){
          std::cerr<<"Error printing environment operators: which_env_ops[w].i>=IFV.size()!"<<std::endl;
          exit(1);
        }
        if(which_env_ops[w].o>=IFV[which_env_ops[w].i]->get_env_ops(n).size()){
          std::cerr<<"Error printing environment operators: which_env_ops[w].o>=IFV[which_env_ops[w].i].env_ops[n].size()!"<<std::endl;
          exit(1);
        }
        if(which_env_ops[w].A.rows()<2){
          which_env_ops[w].A=Eigen::MatrixXcd::Identity(N,N);
        }else if(which_env_ops[w].A.rows()!=N || which_env_ops[w].A.cols()!=N ){
          std::cerr<<"Error printing environment operators: which_env_ops[w].A.rows()!=N || which_env_ops[w].A.cols()!=N !"<<std::endl;
          exit(1);
        }
         
        if(n==0){
          results.add_back_nan(n);
        }else if(n<n_max-1){
          const Eigen::VectorXcd & env_op=
                IFV[which_env_ops[w].i]->get_env_ops(n)[which_env_ops[w].o];

          if(state.cols()!=env_op.size()){
            std::cerr<<"run: n="<<n<<": state.cols()!=env_op.size(): "<<state.cols()<<" vs. "<<env_op.size()<<std::endl;
            exit(1);
          }
          std::complex<double> res=0.;
          for(int d=0; d<env_op.size(); d++){
            for(int nu=0; nu<N; nu++){
              for(int mu=0; mu<N; mu++){
//                res+=env_op(d)*which_env_ops[w].A(nu,mu)*state(nu*N+mu, d);
                res+=env_op(d)*which_env_ops[w].A(mu,nu)*state(nu*N+mu, d);
              }
            }
          }
          results.add_back(n+1, res);
        }else if(n==n_max-1){
          results.add_back_nan(n+1);
        }
      }

    }// loop over n
  }

  void run(Propagator &prop, Smart_Ptr<InfluenceFunctional_OD> IF,
    double ta, double dt, double te, const Eigen::MatrixXcd &initial_rho){
  
    std::vector<Smart_Ptr<InfluenceFunctional_OD> > IFV(1, IF);
    
    run(prop, IFV, ta, dt, te, initial_rho);
  }

  void print_results(const std::string &fname)const{
    results.print(fname);
    FTO.print(results);
  }

  void initialize(){
    print_timestep=false;
  }
  Simulation_OD(Propagator &prop, Smart_Ptr<InfluenceFunctional_OD> IF, 
    double ta, double dt, double te, const Eigen::MatrixXcd &initial_rho){
     
    initialize();
    run(prop, IF, ta, dt, te, initial_rho);
  }
  Simulation_OD(){
    initialize();
  }
};

#endif
