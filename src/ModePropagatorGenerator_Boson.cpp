#include "ModePropagatorGenerator_Boson.hpp"
#include "ModePropagatorGenerator.hpp"
#include "Parameters.hpp"
#include "Operators.hpp"
#include "Operators_Boson.hpp"
#include "Equilibrium.hpp"
#include "ReadTable.hpp"
#include "Reader.hpp"
#include "ReadTemperature.hpp"
#include "ReducedLiouvilleBasis_Boson.hpp"
#include "ReducedLiouvilleBasis_Boson_FB.hpp"
#include "RealFunction_Interpolate.hpp"
#include "LiouvilleTools.hpp"


namespace ACE{

  double ModePropagatorGenerator_Boson::k_label(int k)const{ 
    return get_E(k)/hbar_in_meV_ps;
  }

  //print initial boson number per mode
  void ModePropagatorGenerator_Boson::print_initial_n(const std::string &fname){
    if(fname=="")return;
    std::ofstream ofs(fname.c_str());
    for(int k=0; k<get_N_modes(); k++){
      Eigen::MatrixXcd init=get_bath_init(k);
      double n=0; for(int r=0; r<init.rows(); r++)n+=r*init(r,r).real();
      ofs<<get_E(k)<<" "<<n<<" "<<k<<std::endl;
    }
  }

  double ModePropagatorGenerator_Boson::env_ops_filter(int k)const{
    if(use_env_filter){
      return env_filter->f(E_g.get_omega(k));
    }else{
      return 1.;
    }
  }

  EnvironmentOperators ModePropagatorGenerator_Boson::get_env_ops(int k) const{
    std::vector<Eigen::MatrixXcd> mats;
    double filter=env_ops_filter(k);
    mats.push_back(Eigen::MatrixXcd::Identity(M,M));
    mats.push_back(filter* Operators_Boson::n(M) );
    mats.push_back(filter* get_HE_diag(k) );
    mats.push_back(filter* 2.*get_g(k) * Operators_Boson::adagger(M) );
    mats.push_back(filter* 2.*get_E(k)*get_g(k) * Operators_Boson::adagger(M) );
    return EnvironmentOperators(mats);
  }


  Eigen::MatrixXcd ModePropagatorGenerator_Boson::get_HE_diag(int k)const{
    int sysdim=sysop.rows();
    Eigen::MatrixXcd HE_diag=Operators_Boson::n(M);

    if(use_anharmonic){
      HE_diag-=anharmonic_chi*otimes(
               Eigen::MatrixXcd::Identity(sysdim, sysdim),
               Operators_Boson::n(M)*Operators_Boson::n(M));
    }
    HE_diag*=get_E(k);

    return HE_diag;
  }

  Eigen::MatrixXcd ModePropagatorGenerator_Boson::get_HE(int k)const{
    int sysdim=sysop.rows();
    Eigen::MatrixXcd HE_diag=otimes(
         Eigen::MatrixXcd::Identity(sysdim, sysdim), get_HE_diag(k));

    Eigen::MatrixXcd coupling_op=get_pos_phase_op(k);

    Eigen::MatrixXcd HE_coupling= hbar_in_meV_ps*get_g(k)*
                 otimes(coupling_op, Operators_Boson::adagger(M));


    Eigen::MatrixXcd HE = HE_diag + HE_coupling + HE_coupling.adjoint();
    
    if(use_polaron_shift){
      double pshift=0;
      if(fabs(get_E(k))>1e-12){
        double hg=hbar_in_meV_ps*get_g(k);
        pshift+=hg*hg/get_E(k);
      }
      HE += pshift*otimes(sysop.adjoint()*sysop, Operators_Boson::id(M));
    }
    return HE;
  }

  void ModePropagatorGenerator_Boson::setup(Parameters &param){
    E_g.setup(param, name());
    set_N_modes(E_g.N);
    setup_skip(param);

    M=param.get_as_size_t(add_name("M"), 2);

    Operators op(2);
    sysop=op.ketbra(1,1);   //independent boson model
//    sysop=op.ketbra(0,1);  //Jaynes-Cummings
    if(param.is_specified(add_name("SysOp"))){
      sysop=param.get_as_operator(add_name("SysOp"));
    }
 
    use_pos_phase=param.is_specified(add_name("pos_phase"));
    if(use_pos_phase){
      pos_phase_pos=param.get_as_double(add_name("pos_phase"));
      pos_phase_offset=param.get_as_double(add_name("pos_phase"),0.,0,1);
      pos_phase_c=param.get_as_double(add_name("pos_phase"),c_in_nm_by_ps,0,2);
    }
    
    use_anharmonic=param.is_specified(add_name("anharmonic_chi"));
    anharmonic_chi=param.get_as_double(add_name("anharmonic_chi"),0.);

    use_polaron_shift=param.get_as_bool(add_name("subtract_polaron_shift"),true);
    interaction_picture=param.get_as_bool(add_name("interaction_picture"),false);
    interaction_picture_dt=param.get_as_double(add_name("interaction_picture_dt"),param.get_as_double("dt"));

    std::string print_E_g=param.get_as_string(add_name("print_E_g"));
    if(print_E_g!=""){
      std::ofstream ofs(print_E_g.c_str());
      for(int i=0; i<E_g.N; i++){
        ofs<<E_g.get_E(i)<<" "<<E_g.get_g(i)<<std::endl;
      }
    }
 
    gparam.add_from_prefix(add_name("Propagator"),param);


    use_initial_coherent=false; use_initial_thermal=false;
    if(param.is_specified(add_name("initial_coherent"))){
      use_initial_coherent=true;
      initial_coherent.real(param.get_as_double(add_name("initial_coherent")));
      initial_coherent.imag(param.get_as_double(add_name("initial_coherent"),0.,0,1));
    }else if(param.is_specified(add_name("temperature")) || 
             param.is_specified("temperature") ||
             param.is_specified("temperature_unitless") ||
             param.is_specified(add_name("temperature_unitless")) ||
             param.is_specified(add_name("E_shift_init"))){
     
      use_initial_thermal=true;
      temperature=readTemperature(param,name());

      E_shift_init=param.get_as_double(add_name("E_shift_init"), 
                    param.get_as_double(add_name("omega_shift_init"),0)/
                    hbar_in_meV_ps);
    }

    print_initial_n(param.get_as_string(add_name("print_initial_n")));
    reduced=param.get_as_size_t(add_name("reduced"),0);
    reduced_fb=param.get_as_size_t(add_name("reduced_fb"),0);

    //to introduce line broading in the hope of reducing inner dimensions
    //double value corresponds to thermalization rate
    thermalize=param.get_as_double(add_name("thermalize"),0);

    std::string env_filter_file=param.get_as_string(add_name("env_filter"));
    if(env_filter_file==""){
      use_env_filter=false;
    }else{
      check_file_exists(env_filter_file);
      use_env_filter=true;
      env_filter=std::make_shared<RealFunction_Interpolate>(env_filter_file);
    }

  }

  

  /* multiply position-dependent phase factor exp(i k*r) to coupling
     Note: here we need to distinguish between phyiscal wave vectors k 
     and mode indices, which we now call "j". Actually, the indices correspond
     to an energy discretization (1d):
     hbar*omega= offset +  get_E 
     k=omega/c
 
  */
  Eigen::MatrixXcd ModePropagatorGenerator_Boson::get_pos_phase_op(int j)const{
    Eigen::MatrixXcd pos_phase_op=sysop;
    if(!use_pos_phase){
      return pos_phase_op;
    }
    double k=(pos_phase_offset+get_E(j))/(hbar_in_meV_ps*pos_phase_c);
    std::complex<double> ep=exp(std::complex<double>(0., k*pos_phase_pos));
    pos_phase_op=ep*sysop;

    return pos_phase_op;
  };

  ModePropagatorPtr ModePropagatorGenerator_Boson::get_ModePropagator(int k)const{
    if(k<0||k>=get_N_modes()){
      std::cerr<<"ModePropagatorGenerator_Boson: k<0||k>=get_N_modes()!"<<std::endl; 
      exit(1);
    }

    int sysdim=sysop.rows();

    Parameters kparam=gparam;
    ModePropagatorPtr ptr=std::make_shared<ModePropagator>(sysdim, get_bath_init(k));
    ptr->FreePropagator::setup(kparam);
    ptr->env_ops=get_env_ops(k);
    if(interaction_picture){
      ptr->add_Hamiltonian(Eigen::MatrixXcd::Zero(M*get_N(),M*get_N()));
      double pshift=0;
      if(fabs(get_E(k))>1e-12){
        double hg=hbar_in_meV_ps*get_g(k);
        double w=get_E(k)/hbar_in_meV_ps;
        if(!use_polaron_shift){
          pshift-=hg*hg/get_E(k);
        }
        pshift+=hg*hg/get_E(k)*sin(w*interaction_picture_dt)/(w*interaction_picture_dt);
        ptr->add_Hamiltonian(pshift*otimes(sysop.adjoint()*sysop, Operators_Boson::id(M)));
      }
       
      ComplexFunctionPtr tmpfct=std::make_shared<InteractionPictureFunction>(get_E(k)/hbar_in_meV_ps, interaction_picture_dt);
      Eigen::MatrixXcd tmpop=hbar_in_meV_ps*get_g(k)*otimes(get_pos_phase_op(k), Operators_Boson::adagger(M));

      ptr->add_Pulse(tmpfct, tmpop);

    }else{
      Eigen::MatrixXcd HE=get_HE(k);
      ptr->add_Hamiltonian(HE);
    }
 
    if(reduced_fb>0){
      ptr->rBasis=std::make_shared<ReducedLiouvilleBasis_Boson_FB>(get_bath_init(k),reduced_fb);
   
      if(k==0){
        std::cout<<"Using reduced_fb="<<reduced_fb<<" (dim=";
        std::cout<<ptr->rBasis->U.cols()<<"): "<<std::endl;
print_diff_from_ortho(ptr->rBasis->U);
      }
    }else if(reduced>0){
      ptr->rBasis=std::make_shared<ReducedLiouvilleBasis_Boson>(get_bath_init(k),reduced);
   
      if(k==0){
        std::cout<<"Using reduced="<<reduced<<" (dim=";
        std::cout<<ptr->rBasis->U.cols()<<"): "<<std::endl;
print_diff_from_ortho(ptr->rBasis->U);
      }
    }

    //damping: 
    // e^{-\beta E} = gamma_up/gamma_down , gamma=gamma_up+gamma_down 
    // -> gamma_up=gamma-gamma_down
    // e^{-\beta E} = gamma/gamma_down -1 -> gamma_down = gamma/(e^{-\beta E}+1)
    if(thermalize>0){
      double boltz=1.;
      if(temperature<1e-5){
        boltz=0.;
      }else{
        double x=get_E(k)/(kB_in_meV_by_K*temperature);
        if(x>1e-8){
          boltz=exp(-x);
        }
      }

      if(boltz<1e-8){  //only rate -> down
        ptr->add_Lindblad(thermalize, 
          otimes(Eigen::MatrixXcd::Identity(sysdim,sysdim),Operators_Boson::a(M)) );
      }else{
        double gamma_down=thermalize/(boltz+1.);
        double gamma_up=thermalize-gamma_down;

        ptr->add_Lindblad(gamma_up,
          otimes(Eigen::MatrixXcd::Identity(sysdim,sysdim),Operators_Boson::adagger(M)) );
        ptr->add_Lindblad(gamma_down,
          otimes(Eigen::MatrixXcd::Identity(sysdim,sysdim),Operators_Boson::a(M)) );
      }
    }
 
    return ptr;
  }

  Eigen::MatrixXcd ModePropagatorGenerator_Boson::get_bath_init(int k)const{
    if(use_initial_coherent){
      Eigen::VectorXcd v(M);
      std::complex<double> fac=1.;
      double anorm=exp(-abs(initial_coherent)*abs(initial_coherent)/2.);
      double N2=0.;
      for(int n=0; n<v.size(); n++){
        v(n)=anorm*fac;
        N2+=abs(anorm*fac)*abs(anorm*fac);
        fac*=initial_coherent/sqrt(n+1);
      }
      std::cout<<"Initialize with coherent state "<<initial_coherent<<": Loss of norm due to truncation: "<<1.-N2<<std::endl;
      
      //for(int n=0; n<v.size(); n++)v(n)/=sqrt(N2);

      Eigen::MatrixXcd initial(M,M);
      for(int i=0; i<initial.rows(); i++){
        for(int j=0; j<initial.cols(); j++){
          initial(i,j)=std::conj(v(i))*v(j);
        }
      }
      return initial;
    }else if(use_initial_thermal){
      Eigen::MatrixXcd HB_diag=Operators_Boson::n(M);
      if(use_anharmonic){
        HB_diag-=anharmonic_chi*Operators_Boson::n(M)*Operators_Boson::n(M);
      }
      HB_diag*=(get_E(k)+E_shift_init);
      return Boson_Equilibrium(HB_diag, temperature);
    }else{
      return Operators_Boson::vacuum(M);
    }
  }

}//namespace
