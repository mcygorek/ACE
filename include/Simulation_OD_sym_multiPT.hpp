#ifndef SIMULATION_OD_SYM_MULTIPT_DEFINED_H
#define SIMULATION_OD_SYM_MULTIPT_DEFINED_H

//included by:   src/Simulation_OD.cpp

  void Simulation_OD::run_dt0_nmax_sym_multiPT(Propagator &prop, 
                    std::vector<std::shared_ptr<IF_OD_Abstract> > & IFV,
                    double ta, double dt, double dt0, int n_max, 
                    const Eigen::MatrixXcd &initial_rho){

    if(n_max%2==1)n_max--;
 
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
    results.resize(n_max/2+1);
    results.set(0, ta, output_Op, initial_rho);
    results.add_back_nan(0,which_env_ops.size());

    for(int n=0; n<n_max/2; n++){
      double t = n<=0 ? ta : ta+dt0+(2*n-1)*dt;
      if(print_timestep)std::cout<<"Step: "<<2*n<<"/"<<n_max<<" ("<<t<<")"<<std::endl;

      //do the actual propagation: first only system-free
      prop.update(t, n==0?dt0:dt);
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
//For symmetric Trotter: multiply IFs twice, first is ascending then in descending order. 
//ascending:

      int d1_tot=1;
      for(size_t ifi=0; ifi<IFV.size(); ifi++)d1_tot*=IFV[ifi]->get_a(2*n).dim_d1;

      int d1_lower=d1_tot;
      for(size_t ifi=0; ifi<IFV.size(); ifi++){
        int d1_this=IFV[ifi]->get_a(2*n).dim_d1;
        int d2_this=IFV[ifi]->get_a(2*n).dim_d2;
        d1_lower/=d1_this;

        int d2_tot=(d1_tot/d1_this)*d2_this;
        Eigen::MatrixXcd state2=Eigen::MatrixXcd::Zero(NL,d2_tot);
        
        for(int i=0; i<NL; i++){
          for(int j=0; j<NL; j++){
            int i_ind=IFV[ifi]->dict.beta[i*Ngrps2+j];
            if(i_ind<0)continue;

            for(int glob_d1=0; glob_d1<d1_tot; glob_d1++){
              for(int d2=0; d2<IFV[ifi]->get_a(2*n).dim_d2; d2++){
                int d1=(glob_d1/d1_lower)%d1_this;
                int up=(glob_d1/d1_lower)/d1_this;
                int glob_d2=(up*d2_this+d2)*d1_lower+glob_d1%d1_lower;

                state2(i,glob_d2)+=IFV[ifi]->get_a(2*n)(i_ind,d1,d2)*state(j,glob_d1);
              }
            }
          }
        }
        state=state2;
        d1_tot/=d1_this; d1_tot*=d2_this;
      }

//descending
      d1_tot=1;
      for(size_t ifi=0; ifi<IFV.size(); ifi++)d1_tot*=IFV[ifi]->get_a(2*n+1).dim_d1;

      d1_lower=1;
      for(int ifi=(int)IFV.size()-1; ifi>=0; ifi--){
        if(ifi<IFV.size()-1){
          d1_lower*=IFV[ifi+1]->get_a(2*n+1).dim_d2;
        }
        int d1_this=IFV[ifi]->get_a(2*n+1).dim_d1;
        int d2_this=IFV[ifi]->get_a(2*n+1).dim_d2;

        int d2_tot=(d1_tot/d1_this)*d2_this;
        Eigen::MatrixXcd state2=Eigen::MatrixXcd::Zero(NL,d2_tot);
        
        for(int i=0; i<NL; i++){
          for(int j=0; j<NL; j++){
            int i_ind=IFV[ifi]->dict.beta[i*Ngrps2+j];
            if(i_ind<0)continue;

            for(int glob_d1=0; glob_d1<d1_tot; glob_d1++){
              for(int d2=0; d2<IFV[ifi]->get_a(2*n+1).dim_d2; d2++){
                int d1=(glob_d1/d1_lower)%d1_this;
                int up=(glob_d1/d1_lower)/d1_this;
                int glob_d2=(up*d2_this+d2)*d1_lower+glob_d1%d1_lower;

                state2(i,glob_d2)+=IFV[ifi]->get_a(2*n+1)(i_ind,d1,d2)*state(j,glob_d1);
              }
            }
          }
        }
        state=state2;
        d1_tot=d2_tot;
      }



      //apply again system propagator for symmetric Trotter decomposition
      {
        prop.update(t+(n==0?dt0:dt), dt);

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
      for(size_t ifi=0; ifi<IFV.size(); ifi++)d2_tot*=IFV[ifi]->get_a(2*n+1).dim_d2;

      Eigen::MatrixXcd c_state=state;
      for(size_t ifi=0; ifi<IFV.size(); ifi++){
        int d2_this=IFV[ifi]->get_a(2*n+1).dim_d2;
        int d2_factor=c_state.cols()/d2_this;

        Eigen::MatrixXcd c_state2=Eigen::MatrixXcd::Zero(NL, d2_factor);
        for(int j=0; j<NL; j++){
          for(int d2=0; d2<d2_this; d2++){
            for(int o=0; o<d2_factor; o++){
              c_state2(j, o)+=IFV[ifi]->get_c(2*n+1)(d2)*c_state(j, d2*d2_factor+o);
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
      results.set(n+1, ta+dt0+(2*n+1)*dt, output_Op, rho);


      //Environment operators
      for(size_t w=0; w<which_env_ops.size(); w++){
        if(which_env_ops[w].i>=(int)IFV.size()){
          std::cerr<<"Error printing environment operators: which_env_ops[w].i>=IFV.size()!"<<std::endl;
          exit(1);
        }
        if(which_env_ops[w].o>=(int)IFV[which_env_ops[w].i]->get_env_ops(2*n+1).size()){
          std::cerr<<"Error printing environment operators: which_env_ops[w].o>=IFV[which_env_ops[w].i].env_ops[n].size()!"<<std::endl;
          exit(1);
        }
        if(which_env_ops[w].A.rows()<2){
          which_env_ops[w].A=Eigen::MatrixXcd::Identity(N,N);
        }else if(which_env_ops[w].A.rows()!=N || which_env_ops[w].A.cols()!=N ){
          std::cerr<<"Error printing environment operators: which_env_ops[w].A.rows()!=N || which_env_ops[w].A.cols()!=N !"<<std::endl;
          exit(1);
        }
         
        if(2*n+1<n_max-1){
          Eigen::VectorXcd env_op(1); env_op<<1;
          for(int ifi=0; ifi<(int)IFV.size(); ++ifi){
            if(ifi==which_env_ops[w].i){
              env_op=Vector_otimes(env_op, IFV[which_env_ops[w].i]->get_env_ops(2*n+1)[which_env_ops[w].o]);
            }else{
              env_op=Vector_otimes(env_op, IFV[ifi]->get_c(2*n+1));
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
        }else if(2*n+1==n_max-1){
          results.add_back_nan(n+1);
        }
      }

    }// loop over n
  }

#endif
