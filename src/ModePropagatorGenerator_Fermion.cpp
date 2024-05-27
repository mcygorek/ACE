#include "ModePropagatorGenerator_Fermion.hpp"
#include "ModePropagatorGenerator.hpp"
#include "Parameters.hpp"
#include "ReadTemperature.hpp"
#include "Operators.hpp"
#include "otimes.hpp"

namespace ACE{

  double ModePropagatorGenerator_Fermion::k_label(int k)const{ 
    return get_E(k)/hbar_in_meV_ps; 
  }

  std::vector<Eigen::MatrixXcd> ModePropagatorGenerator_Fermion::get_env_ops(int k) const{
    Operators op(2); 
    std::vector<Eigen::MatrixXcd> mats;
    if(!no_env_ops){
      mats.push_back(op.ketbra(1,1));
      mats.push_back(op.ketbra(1,0));
      mats.push_back(2.*E_g.get_g(k)* op.ketbra(1,0) );
    }
    return mats;
  }

  void ModePropagatorGenerator_Fermion::setup(Parameters &param){
    E_g.setup(param, name());
    set_N_modes(E_g.N);
    setup_skip(param);

    EFermi=param.get_as_double("EFermi", -1e6);
    EFermi=param.get_as_double(add_name("EFermi"), EFermi);

    temperature=readTemperature(param, name());

    //to introduce line broading in the hope of reducing inner dimensions
    //double value corresponds to thermalization rate
    thermalize=param.get_as_double(add_name("thermalize"),0);

    Operators op(2);
    sysop=op.ketbra(0,1);  
    if(param.is_specified(add_name("SysOp"))){
      sysop=param.get_as_operator(add_name("SysOp"));
    }

    no_env_ops=param.get_as_bool(add_name("no_env_ops"),false);

    double low_pass_cutoff=param.get_as_double(add_name("low_pass_cutoff"),0);
    if(low_pass_cutoff>1e-16){
      low_pass.use=true;
      low_pass.cutoff=low_pass_cutoff;
      low_pass.factor=param.get_as_double(add_name("low_pass_cutoff"),1.,0,1);
    }

    continuum_subdiv_N=param.get_as_double(add_name("continuum_subdiv"),0);
    

    gparam.add_from_prefix(add_name("Propagator"),param);

    std::string print_E_g=param.get_as_string(add_name("print_E_g"));
    if(print_E_g!=""){
      std::ofstream ofs(print_E_g.c_str());
      for(int i=0; i<E_g.N; i++){
        ofs<<E_g.get_E(i)<<" "<<E_g.get_g(i)<<std::endl;
      }
    } 
//std::cout<<"FERMION: setup ended"<<std::endl;
  }

  ModePropagatorPtr ModePropagatorGenerator_Fermion::getModePropagator(int k)const{
    if(k<0||k>=get_N_modes()){
      std::cerr<<"ModePropagatorGenerator_Fermion: k<0||k>=get_N_modes()!"<<std::endl; 
      exit(1);
    }

    int sysdim=sysop.rows();
    Operators op(sysdim);
    Operators op2(2);
    Eigen::MatrixXcd HB_diag=otimes(op.id(), op2.ketbra(1,1));
    Eigen::MatrixXcd HB_base=otimes(sysop, op2.ketbra(1,0)) + \
                             otimes(sysop.adjoint(), op2.ketbra(0,1));


    Eigen::MatrixXcd HB = get_E(k)*HB_diag
           +hbar_in_meV_ps*get_g(k)*HB_base;

    
    Parameters kparam=gparam;
    ModePropagatorPtr ptr=std::make_shared<ModePropagator>(sysdim, get_bath_init(k));
    ptr->FreePropagator::setup(kparam);
    ptr->env_ops=get_env_ops(k);
    ptr->add_Hamiltonian(HB);
  

    //damping: 
    // n = gamma_up/gamma_down , gamma=gamma_up+gamma_down 
    // -> gamma_up=gamma-gamma_down
    // n = gamma/gamma_down -1 -> gamma_down = gamma/(n+1)
    if(thermalize>0){
      double n_eq=get_n_eq(k);
      if((1.-n_eq)<1e-8){  //only rate ->up
        ptr->add_Lindblad(thermalize, 
          otimes(Eigen::MatrixXcd::Identity(sysdim,sysdim),op2.ketbra(1,0)) );
      }else{
        //double gamma_down=thermalize/(n_eq+1.);
        //double gamma_up=thermalize-gamma_down;
        double gamma_up=thermalize*(n_eq);
        double gamma_down=thermalize*(1.-n_eq);

        ptr->add_Lindblad(gamma_up,
          otimes(Eigen::MatrixXcd::Identity(sysdim,sysdim),op2.ketbra(1,0)) );
        ptr->add_Lindblad(gamma_down,
          otimes(Eigen::MatrixXcd::Identity(sysdim,sysdim),op2.ketbra(0,1)) );
      }
    }

    ptr->low_pass=low_pass;

    if(continuum_subdiv_N>0){  
      ptr->continuum_subdiv.N=continuum_subdiv_N;
      ptr->continuum_subdiv.dH=get_dE(k)*HB_diag;
    }

/* Test: reduced basis:
    ptr->rBasis->use_reduce=true;
    ptr->rBasis->U=Eigen::MatrixXcd::Zero(4,4);
    ptr->rBasis->U(0,0)=1;
    ptr->rBasis->U(1,2)=1;
    ptr->rBasis->U(2,1)=0.99;
    ptr->rBasis->U(3,3)=1;
*/
    return ptr;

//    return ModePropagatorPtr(new ModePropagator(sysdim,get_bath_init(k),HB,get_env_ops()));
  }

  double ModePropagatorGenerator_Fermion::get_n_eq(int k)const{
    if(temperature<1e-6){
      if(get_E(k)>EFermi)return 0;
      else return 1;
    }
    double x=(get_E(k)-EFermi)/(kB_in_meV_by_K*temperature);
    return fermi(x);
  }
  Eigen::MatrixXcd ModePropagatorGenerator_Fermion::get_bath_init(int k)const{
    double n_eq=get_n_eq(k);
    return (1.-n_eq)*Operators(2).ketbra(0,0)+n_eq*Operators(2).ketbra(1,1);
  }

}//namespace
