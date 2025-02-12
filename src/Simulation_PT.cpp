#include "Simulation_PT.hpp"
#include "LiouvilleTools.hpp"
#include "DummyException.hpp"
#include "TransferTensor.hpp"
#include "LindbladMasterEquation.hpp"
#include "ReaderBasics.hpp"

namespace ACE{

void Simulation_PT::propagate_system(
         Eigen::MatrixXcd & state, Propagator &prop, double t, double dt){

  prop.update(t, dt);
  if(prop.M.rows()!=prop.M.cols() || prop.M.rows()!=state.rows()){
    std::cerr<<"Simulation_PT::propagate_system: prop.M.rows()!=prop.M.cols() || prop.M.rows()!=state.rows()!"<<std::endl;
    exit(1);
  }

//  Eigen::MatrixXcd state2=Eigen::MatrixXcd::Zero(state.rows(),state.cols());
//  for(int d1=0; d1<state.cols(); d1++){
//    for(int j=0; j<state.rows(); j++){
//      for(int k=0; k<state.rows(); k++){
//        state2(j,d1)+=prop.M(j,k)*state(k,d1);
//      }
//    }
//  }
//  state.swap(state2);
  state=prop.M*state;
}

void Simulation_PT::propagate_state(Eigen::MatrixXcd &state, int n, const TimeGrid &tgrid, Propagator &prop, ProcessTensorForwardList &PT)const{

  double t=tgrid.get_t(n);
  double dt=tgrid.get_dt(n);
  if(use_symmetric_Trotter && (!propagate_alternate)){  
    propagate_system(state, prop, t, dt/2.);
    PT.propagate(state, false);
    propagate_system(state, prop, t+dt/2., dt/2.);

  }else if(n%2==0 || !propagate_alternate){
    propagate_system(state, prop, t, dt);
    PT.propagate(state, false);
  }else{
    PT.propagate(state, true);
    propagate_system(state, prop, t, dt);
  }
}

Eigen::MatrixXcd Simulation_PT::run_std(
           Propagator &prop, ProcessTensorForwardList &PT,
           const Eigen::MatrixXcd & initial_rho, const TimeGrid &tgrid,
           OutputPrinter &printer){

  int N=initial_rho.rows();
  int NL=N*N;
  if(N<2){
    std::cerr<<"Simulation_PT::run: N<2!"<<std::endl;
    exit(1);
  }
  if(tgrid.n_tot<1){
    std::cerr<<"Simulation_PT::run: tgrid.n_tot<1!"<<std::endl;
    exit(1);
  }
  for(int i=0; i<PT.size(); i++){
    if( (!PT.list[i]) || PT.list[i]->get_n_tot()<tgrid.n_tot){
      std::cerr<<"Simulation_PT::run: PT["<<i<<"]->get_n_tot()="<<PT.list[i]->get_n_tot()<<"<tgrid.n_tot="<<tgrid.n_tot<<"!"<<std::endl;
      exit(1);
    }
  }

  Eigen::MatrixXcd state(NL, 1);
  state.col(0)=H_Matrix_to_L_Vector(initial_rho);

  Eigen::VectorXcd rho_reduced = state.col(0);
  std::vector<Eigen::MatrixXcd> pair_reduced(PT.size(), state);

  printer.print(0, tgrid.get_t(0), state.col(0), 
    std::vector<std::complex<double> >(printer.which_env_ops.size(), 0.)); 

  printer.print_eigenstate_occupations(tgrid.get_t(0), prop.get_Htot(tgrid.get_t(0)), L_Vector_to_H_Matrix(rho_reduced));
 
  int maxdim=1;
  PT.reset();
  for(int n=0; n<tgrid.n_tot; n++){
    if(print_timesteps){ 
      std::cout<<"step: "<<n<<"/"<<tgrid.n_tot<<std::endl;
    }
    propagate_state(state, n, tgrid, prop, PT);
    if(maxdim<state.cols())maxdim=state.cols();

    rho_reduced=PT.get_rho_reduced(state);

    double t_next=tgrid.get_t(n+1);
    printer.print(n+1, t_next, rho_reduced, 
                     PT.get_env_reduced(state, printer.which_env_ops));

    printer.print_eigenstate_occupations(t_next, prop.get_Htot(t_next), L_Vector_to_H_Matrix(rho_reduced));

    PT.load_next();
  }
  PT.reset();
  if(print_final_maxdim){std::cout<<"Final Maxdim: "<<maxdim<<std::endl;}

  return L_Vector_to_H_Matrix(rho_reduced);
}


/*
Eigen::MatrixXcd Simulation_PT::run_TT(
           Propagator &prop, ProcessTensorForwardList &PT,
           const Eigen::MatrixXcd & initial_rho, const TimeGrid &tgrid,
           OutputPrinter &printer){
      std::cout<<"Calculating Transfer Tensor over "<<TT_n_mem<<" time steps..."<<std::endl;
      
      TransferTensor TT;
      std::vector<Eigen::MatrixXcd> E=TT.calculate_E(prop, PT, *this, tgridPT);
      TT.set_T_from_E(E);
      if(param.get_as_bool("TT_use_LT")){
        TT.set_LT_from_E(E);
        TT.use_LT=true;
      }


//std::cout<<"TT[0]:"<<std::endl<<TT.T[0]<<std::endl;
//std::cout<<"TT norms:"<<std::endl; TT.print_norms();
      std::string TT_print_norms=param.get_as_string("TT_print_norms");
      TT.print_norms(TT_print_norms, tgrid.dt);
      TT.propagate(sim.initial_rho, tgrid, sim.printer);
}
*/

Eigen::MatrixXcd Simulation_PT::run(
           Propagator &prop, ProcessTensorForwardList &PT,
           const Eigen::MatrixXcd & initial_rho, const TimeGrid &tgrid,
           OutputPrinter &printer){

  if(!use_TT || TT_n_from>=tgrid.n_tot){
    return run_std(prop, PT, initial_rho, tgrid, printer);

  }else{
    TransferTensor TT;
    TimeGrid tgridTT=tgrid; 
    tgridTT.ta=tgrid.get_t(TT_n_from); 
    tgridTT.n_tot=TT_n_mem;
    std::cout<<"TT: calculating propagators for "<<TT_n_mem<<" steps..."<<std::endl;
    std::vector<Eigen::MatrixXcd> E=TT.calculate_E(prop, PT, *this, tgridTT, 1);
    
    int NL=initial_rho.rows()*initial_rho.cols();

    if(ME_print_rates!="" || ME_print_L!=""){
      std::ofstream ofs_rates;
      if(ME_print_rates!=""){
        ofs_rates.open(ME_print_rates.c_str());
      }
      std::vector<std::ofstream> ofs_L(NL-1);
      if(ME_print_L!=""){
        for(int k=0; k<NL-1; k++){
          std::string str=ME_print_L+"_"+double_to_string(k);
          ofs_L[k].open(str.c_str());
        }
      }

      Eigen::MatrixXcd log_dE_last;
      for(size_t l=0; l<E.size(); l++){
        Eigen::MatrixXcd dE=E[l];
        Eigen::MatrixXcd log_dE=dE.log();
        Eigen::MatrixXcd Liou=log_dE/tgrid.get_dt(l);
        if(l>0){
          Eigen::JacobiSVD<Eigen::MatrixXcd> svd(E[l-1], Eigen::ComputeFullU | Eigen::ComputeFullV);
//std::cout<<"singular values: "<<svd.singularValues().transpose()<<std::endl;
          Eigen::VectorXcd inv_sigma(NL);
          for(int i=0; i<NL; i++)inv_sigma(i)=1./svd.singularValues()(i);
          dE=E[l]*svd.matrixV()*(inv_sigma.asDiagonal())*svd.matrixU().adjoint();
          Eigen::MatrixXcd log_dE=dE.log();
std::cout<<"(log_dE.exp()*E[l-1]-E[l]).norm(): "<<(log_dE.exp()*E[l-1]-E[l]).norm()<<std::endl;
          Liou=(log_dE+log_dE_last)/(2.*tgrid.get_dt(l));
        }
        log_dE_last=log_dE;

        LindbladMasterEquation LME;
        LME.set_from_Liouvillian(Liou, 0, 0);
std::cout<<"t="<<tgrid.get_t(l)<<": |Liou-reconstructed|="<<(Liou-LME.construct_Liouvillian()).norm()<<std::endl;

        if(ME_print_rates!=""){
          ofs_rates<<tgrid.get_t(l);
          for(size_t k=0; k<LME.L.size(); k++){
            ofs_rates<<" "<<LME.L[k].first;
          }
          ofs_rates<<std::endl;
        }

        if(ME_print_L!=""){
          for(int k=0; k<LME.L.size(); k++){
            ofs_L[k]<<tgrid.get_t(l);
            for(int r=0; r<LME.L[k].second.rows(); r++){
              for(int c=0; c<LME.L[k].second.cols(); c++){
                std::complex<double> value=LME.L[k].second(r,c);
                ofs_L[k]<<" "<<value.real()<<" "<<value.imag();
              }
            }
            ofs_L[k]<<std::endl;
          }
/*
          std::cout<<"Traces: t="<<tgrid.get_t(l)<<std::endl;
          for(int k=0; k<LME.L.size(); k++){
            std::cout<<" "<<LME.L[k].second.trace();
          }
          std::cout<<std::endl;
*/
        }      
      }
    }

    TT.set_T_from_E(E);
    if(TT_print_norms!=""){TT.print_norms(TT_print_norms, tgrid.dt);}
    if(use_LT){
      TT.use_LT=true;
      TT.set_LT_from_E(E);
      std::cout<<"Using LT"<<std::endl;
    }

    if(TT_n_from<=0){
      std::cout<<"Scheduled: Use TTs over "<<tgrid.n_tot<<" steps"<<std::endl; 
      TT.propagate(H_Matrix_to_L_Vector(initial_rho), tgrid, printer);

    }else{
      std::cout<<"Scheduled: PT to n="<<TT_n_from<<" then TT to "<<tgrid.n_tot<<std::endl; 
   
      TimeGrid tgrid1=tgrid;
      tgrid1.n_tot=TT_n_from;
      
      if(printer.do_extract){
        std::cerr<<"Simulation_PT::run: printer.do_extract with use_TT NOT IMPLEMENTED YET!"<<std::endl;
        throw DummyException();
      }
      printer.do_extract=true;
      printer.start_extract=TT_n_from-TT_n_mem;
      
      run_std(prop, PT, initial_rho, tgrid1, printer);
      printer.do_extract=false;
      std::vector<Eigen::VectorXcd> rho_t; rho_t.swap(printer.rho_t);
      std::reverse(rho_t.begin(), rho_t.end());
      rho_t.resize(TT_n_mem,Eigen::VectorXcd::Zero(rho_t[0].rows()));

      TimeGrid tgrid_bulk=tgrid;
      tgrid_bulk.ta=tgrid.get_t(TT_n_from); tgrid_bulk.n_tot-=TT_n_from;
      if(TT.use_LT){
        TT.propagate_LT(rho_t, tgrid_bulk, printer);
      }else{
        TT.propagate_bulk(rho_t, tgrid_bulk, printer);
      }
    } 
    return Eigen::MatrixXcd();
  }
}

void Simulation_PT::setup(Parameters & param){
  print_timesteps=param.get_as_bool("print_timesteps",false);
  propagate_alternate=param.get_as_bool("propagate_alternate",false);
  print_final_maxdim=param.get_as_bool("print_final_maxdim",false);
  use_symmetric_Trotter=param.get_as_bool("use_symmetric_Trotter",true);


  //control use of Transfer Tensors
  TimeGrid tgrid(param);
  use_LT=param.get_as_bool("use_LT",false);
  use_TT=TransferTensor::get_use_TT(param);
if(use_TT){std::cout<<"use_TT=true"<<std::endl;}//else{std::cout<<"use_TT=false"<<std::endl;}
  TT_n_from=TransferTensor::get_TT_n_from(param);
  TT_n_mem=TransferTensor::get_TT_n_mem(param, use_TT);
  TT_print_norms=param.get_as_string("TT_print_norms");
  ME_print_rates=param.get_as_string("ME_print_rates");
  ME_print_L=param.get_as_string("ME_print_L");
}


}//namespace
