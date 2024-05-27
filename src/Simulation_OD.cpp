#include "Simulation_OD.hpp"
#include "IF_from_Parameters.hpp"
#include "Propagator.hpp"
#include "Simulation_Results.hpp"
#include "InitialState.hpp"
#include "FT_Output.hpp"
#include "Which_Env_Ops.hpp"

#include "ModePropagatorGenerator.hpp"

namespace ACE{

  void Simulation_OD::add_output_Op(const Eigen::MatrixXcd &op){
    output_Op.add(op);
  }
  void Simulation_OD::setup_output(Parameters &param){
    output_Op.setup(param);
    FTO.setup(param);
    which_env_ops.setup(param);
  }

  void Simulation_OD::run_nobath_dt0_nmax(Propagator &prop,
                           double ta, double dt, double dt0, int n_max, 
                           const Eigen::MatrixXcd &initial_rho){

std::cout<<"RUN_NOBATH!"<<std::endl;

    int N=initial_rho.rows();
    int NL=N*N;

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
      double t=ta;
      if(n>0) t+=dt0+(n-1)*dt;
      if(print_timestep)std::cout<<"Step: "<<n<<"/"<<n_max<<" ("<<t<<")"<<std::endl;
      if(n==0){
        prop.update(t,dt0);
      }else{
        prop.update(t,dt);
      }

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
      results.set(n+1, ta+dt0+n*dt, output_Op, rho);
    }
  }

  void Simulation_OD::run_nobath(Propagator &prop,
    double ta, double dt, double te, const Eigen::MatrixXcd &initial_rho){
    
    run_nobath_dt0_nmax(prop, ta, dt, dt, (te-ta)/dt+1e-12, initial_rho);

  }

  void Simulation_OD::run_dt0_nmax(Propagator &prop, 
                    std::vector<std::shared_ptr<IF_OD_Abstract> > & IFV,
                    double ta, double dt, double dt0, int n_max, 
                    const Eigen::MatrixXcd &initial_rho){
 
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
      run_nobath_dt0_nmax(prop, ta, dt, dt0, n_max, initial_rho);
      return;
    }   


    int N=initial_rho.rows();
    int NL=N*N;
    int Ngrps2=NL;


    for(size_t ifi=0; ifi<IFV.size(); ifi++){
//      std::cout<<"ifi: "<<ifi<<" "<<IFV[ifi]->get_rank()<<std::endl;
      IFV[ifi]->check_within_limits(n_max);
//      IFV[ifi]->check_consistency();

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
    results.add_back_nan(0,which_env_ops.size());

    for(int n=0; n<n_max; n++){
      double t = n<=0 ? ta : ta+dt0+(n-1)*dt;
      if(print_timestep)std::cout<<"Step: "<<n<<"/"<<n_max<<" ("<<t<<")"<<std::endl;

      //do the actual propagation: first only system-free
      double dt1=dt; if(n==0)dt1=dt0;
      if(use_symmetric_Trotter){
        prop.update(t,dt1/2.);
      }else{
        prop.update(t,dt1);
      }
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


      //apply again system propagator for symmetric Trotter decomposition
      if(use_symmetric_Trotter){
        prop.update(t+dt1/2.,dt1/2.);

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
      results.set(n+1, ta+dt0+n*dt, output_Op, rho);


      //Environment operators
      for(size_t w=0; w<which_env_ops.size(); w++){
        if(which_env_ops[w].i>=(int)IFV.size()){
          std::cerr<<"Error printing environment operators: which_env_ops[w].i>=IFV.size()!"<<std::endl;
          exit(1);
        }
        if(which_env_ops[w].o>=(int)IFV[which_env_ops[w].i]->get_env_ops(n).size()){
          std::cerr<<"Error printing environment operators: which_env_ops[w].o>=IFV[which_env_ops[w].i].env_ops[n].size() ";
          std::cerr<<"("<<which_env_ops[w].o<<" vs. ";
          std::cerr<<(int)IFV[which_env_ops[w].i]->get_env_ops(n).size();
          std::cerr<<")!"<<std::endl;
          exit(1);
        }
        if(which_env_ops[w].A.rows()<2){
          which_env_ops[w].A=Eigen::MatrixXcd::Identity(N,N);
        }else if(which_env_ops[w].A.rows()!=N || which_env_ops[w].A.cols()!=N ){
          std::cerr<<"Error printing environment operators: which_env_ops[w].A.rows()!=N || which_env_ops[w].A.cols()!=N !"<<std::endl;
          exit(1);
        }
         
        if(n<n_max-1){
          Eigen::VectorXcd env_op(1); env_op<<1;
          for(int ifi=0; ifi<(int)IFV.size(); ++ifi){
            if(ifi==which_env_ops[w].i){
              env_op=Vector_otimes(env_op, IFV[which_env_ops[w].i]->get_env_ops(n)[which_env_ops[w].o]);
            }else{
              env_op=Vector_otimes(env_op, IFV[ifi]->get_c(n));
            }
          }

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

#include "Simulation_OD_sym_multiPT.hpp"
 
  void Simulation_OD::run_sym_multiPT(Propagator &prop, 
    std::vector<std::shared_ptr<IF_OD_Abstract> > & IFV,
    const TimeGrid &tgrid, const Eigen::MatrixXcd &initial_rho){
    
    run_dt0_nmax_sym_multiPT(prop, IFV, tgrid.ta, tgrid.dt, tgrid.dt0, tgrid.n_tot, initial_rho);
  }

  void Simulation_OD::run_sym_multiPT(Propagator &prop,
    std::vector<std::shared_ptr<InfluenceFunctional_OD> > & IFV,
    const TimeGrid &tgrid, const Eigen::MatrixXcd &initial_rho){

    std::vector<std::shared_ptr<IF_OD_Abstract> > IFabs(IFV.size());
    for(size_t i=0; i<IFV.size(); i++){
      IFabs[i]=std::make_shared<InfluenceFunctional_OD>(*IFV[i].get());
    } 
    run_sym_multiPT(prop, IFabs, tgrid, initial_rho);
  }



  void Simulation_OD::run(Propagator &prop, 
    std::vector<std::shared_ptr<IF_OD_Abstract> > & IFV,
    double ta, double dt, double te, const Eigen::MatrixXcd &initial_rho){
    
    run_dt0_nmax(prop, IFV, ta, dt, dt, (te-ta)/dt+1e-12, initial_rho);
  }

  void Simulation_OD::run(Propagator &prop, 
    std::vector<std::shared_ptr<InfluenceFunctional_OD> > & IFV,
    double ta, double dt, double te, const Eigen::MatrixXcd &initial_rho){
    
    std::vector<std::shared_ptr<IF_OD_Abstract> > IFabs(IFV.size());
    for(size_t i=0; i<IFV.size(); i++){
      IFabs[i]=std::make_shared<InfluenceFunctional_OD>(*IFV[i].get());
    } 
    run(prop, IFabs, ta, dt, te, initial_rho);
  }


  void Simulation_OD::run(Propagator &prop, 
    std::vector<std::shared_ptr<IF_OD_Abstract> > & IFV,
    const TimeGrid &tgrid, const Eigen::MatrixXcd &initial_rho){
    
    run_dt0_nmax(prop, IFV, tgrid.ta, tgrid.dt, tgrid.dt0, tgrid.n_tot, initial_rho);
  }
 
  void Simulation_OD::run(Propagator &prop,
    std::vector<std::shared_ptr<InfluenceFunctional_OD> > & IFV,
    const TimeGrid &tgrid, const Eigen::MatrixXcd &initial_rho){

    std::vector<std::shared_ptr<IF_OD_Abstract> > IFabs(IFV.size());
    for(size_t i=0; i<IFV.size(); i++){
      IFabs[i]=std::make_shared<InfluenceFunctional_OD>(*IFV[i].get());
    } 
    run(prop, IFabs, tgrid, initial_rho);
  }


  void Simulation_OD::run(Propagator &prop, std::shared_ptr<IF_OD_Abstract> IF,
    double ta, double dt, double te, const Eigen::MatrixXcd &initial_rho){
  
    std::vector<std::shared_ptr<IF_OD_Abstract> > IFV(1, IF);
    run(prop, IFV, ta, dt, te, initial_rho);
  }

  void Simulation_OD::print_results(const std::string &fname)const{
    results.print(fname);
    FTO.print(results);
  }

  void Simulation_OD::initialize(){
    print_timestep=false;
    use_symmetric_Trotter=true;
  }

}//namespace
