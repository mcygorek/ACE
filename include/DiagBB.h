#ifndef DIAG_BB_DEFINED_H
#define DIAG_BB_DEFINED_H

#include "Function.h"
#include "Constants.h"
#include "Coupling_Groups.h"
#include "Parameters.h"
#include "ReadTemperature.h"
#include "EnergyRange.h"
#include "SpectralDensity_Selector.h"
#include "HilbertSpaceRotation.h"



class DiagBB{
public:
  RealFunctionPtr J;
  double temperature;
  double E_shift_init;
  double omega_min;
  double omega_max;
  size_t Ndiscr;
  bool noSubPS;
  bool separate_freq;
  HilbertSpaceRotation hs_rot;
  std::vector<double> couplings;
  Coupling_Groups groups;

  virtual int get_dim()const{ return couplings.size(); }
  virtual int sys_dim()const{ return groups.sys_dim(); }

  
  double get_beta() {
    if(temperature<1e-10)return 1e12;
    return 1./(Constants::kB_in_meV_by_K*temperature);
  }
  static double get_coth(double beta, double E_shift, double w){
    if(beta>=0.999e12)return 1.;
    else return Constants::coth( 0.5*beta*(Constants::hbar_in_meV_ps*w+E_shift) );
  }

  #include "DiagBB_K_integrands.h"


  //obtain K_XX for diagonally coupled bath
  std::complex<double> calculate_K(int n, double dt){
    double beta=get_beta();
    if(n==0){
      if(noSubPS){
        K0_integrand_class_noSubPS integrand(J, beta, E_shift_init, dt);
        return integrand.integrate(omega_min, omega_max, Ndiscr);
      }else{
        K0_integrand_class integrand(J, beta, E_shift_init, dt);
        return integrand.integrate(omega_min, omega_max, Ndiscr);
      }
    }else{
/*
        std::complex<double> res=0.; 
        K0_posdt_integrand_class integrand_posdt(J, beta, E_shift_init);
        res+=integrand_posdt.integrate_times_expI(omega_min, omega_max, Ndiscr, dt);
        K0_negdt_integrand_class integrand_negdt(J, beta, E_shift_init);
        res+=integrand_negdt.integrate_times_expI(omega_min, omega_max, Ndiscr, -dt);
        if(noSubPS){
          K0_const_integrand_class_noSubPS integrand(J, beta, E_shift_init, dt);
          return res+integrand.integrate(omega_min, omega_max, Ndiscr);
        }else{
          K0_const_integrand_class integrand(J, beta, E_shift_init, dt);
          return res+integrand.integrate(omega_min, omega_max, Ndiscr);
        }
*/
      if(separate_freq){
        std::complex<double> res=0.; 
        Kn_posdt_integrand_class integrand_posdt(J, beta, E_shift_init, dt);
        res+=integrand_posdt.integrate_times_expI(omega_min, omega_max, Ndiscr, n*dt);
        Kn_negdt_integrand_class integrand_negdt(J, beta, E_shift_init, dt);
        res+=integrand_negdt.integrate_times_expI(omega_min, omega_max, Ndiscr, -n*dt);
        return res;

      }else{
        Kn_integrand_class integrand(J, beta, E_shift_init, dt, n*dt);
        return integrand.integrate(omega_min, omega_max, Ndiscr);
      }
    }
  }

  Eigen::MatrixXcd calculate_expS(int n, double dt){

    int dim=get_dim();
    int NL=dim*dim;
    Eigen::MatrixXcd expS=Eigen::MatrixXcd::Zero(NL, NL);

    std::complex<double> K=calculate_K(n,dt);
    for(int nu_k=0; nu_k<dim; nu_k++){
      for(int mu_k=0; mu_k<dim; mu_k++){
        for(int nu_l=0; nu_l<dim; nu_l++){
          for(int mu_l=0; mu_l<dim; mu_l++){
            std::complex<double> S=
              (-couplings[nu_l]*couplings[nu_k]
               +couplings[nu_l]*couplings[mu_k])*K
             +(-couplings[mu_k]*couplings[mu_l]
               +couplings[nu_k]*couplings[mu_l])*std::conj(K);

            expS(nu_k*dim+mu_k,nu_l*dim+mu_l)=exp(S);
          }
        }
      }
    }
    if(n==0){
      for(int i=0; i<NL; i++){
        for(int j=0; j<NL; j++){  
          if(i!=j)expS(i,j)=0.;
        }
      }
    }
    return expS;
  }

  void print_K(const std::string &fname, int n_max, double dt){
    std::ofstream ofs(fname.c_str());
    for(int n=0; n<n_max; n++){
      double t=n*dt;
      std::complex<double> K=calculate_K(n, dt);
      ofs<<t<<" "<<K.real()<<" "<<K.imag()<<std::endl;
    }
  }

  int estimate_memory_length(int n_max, double dt, double threshold, bool verbose){
    // Find memory time: relative to maximal (max) value around K(t=0).
    // K can be oscillatory. So, after finding first |K| < threshold * max, go at least 2 twice as long to see if values larger than that are found
  
    if(threshold<=0.){
      return n_max;
    }

    double max=0.;
    int below_thr=-1;
    int first_below_thr=-1;
    int n_break=-1;

    if(verbose)std::cout<<"Estimating memory time"<<std::endl;
    for(int n=0; n<n_max; n++){
      double t=n*dt;
      std::complex<double> K=calculate_K(n, dt);

      if(n<=10){ 
        if(abs(K)>max)max=abs(K);
      }else{
        if(abs(K)>threshold*max){
          below_thr=-1.;
        }else{
          if(below_thr<0){
            below_thr=n;
            if(first_below_thr<0)first_below_thr=n;
          }
        }
        //test break criterion:
        if(below_thr>0 && n>below_thr+2*first_below_thr>0){
          n_break=n;
          break;
        }
      }
    }

    if(n_break>=0){
      if(verbose){
        std::cout<<"First time below threshold: "<<first_below_thr*dt<<std::endl;
        std::cout<<"Estimated memory time: "<<n_break*dt<<std::endl;
      }
      return n_break;
    }else{
      if(verbose){
        std::cout<<"No break before n_max*dt: "<<n_max*dt<<std::endl;
      }
      return n_max;
    }
  }


  void setup_groups_and_couplings(const Eigen::MatrixXcd &Op){
    groups.setup_from_matrix(Op);

    couplings.resize(groups.Ngrps);
    for(int i=0; i<Op.rows(); i++){
      couplings[groups[i]]=Op(i,i).real();
    }
    if(couplings.size()<2){
      std::cerr<<"DiagBB:set_defaults: couplings.size()<2!"<<std::endl;
      exit(1);
    }
  }
  void setup(const Coupling_Groups & grp, 
                const Eigen::MatrixXcd &couplings_, 
                RealFunctionPtr J_, double temperature_=0., 
                bool noSubPS_=false,
                double omega_min_ = 0., double omega_max_ = 50., 
                double E_shift_init_ =0. ) {
    J=J_;
    temperature=temperature_;
    groups=grp;
    Ndiscr=1e6;
    noSubPS=noSubPS_;
    separate_freq=false;
    omega_min=omega_min_;
    omega_max=omega_max_;
    E_shift_init=E_shift_init_;

    std::cout<<"DiagBB NGRPS: "<<grp.Ngrps<<std::endl;
    std::cout<<"Couplings Matrix: "<<std::endl<<couplings_<<std::endl;
//    J->print("TEST.dat",omega_min_, omega_max_, 10000);
    
    couplings.resize(grp.Ngrps);
    for(int i=0; i<couplings_.rows(); i++){
      couplings[grp[i]]=couplings_(i,i).real();
    }
    if(couplings.size()<2){
      std::cerr<<"DiagBB:set_defaults: couplings.size()<2!"<<std::endl;
      exit(1);
    }    

  }
  static void complain_if_not_Hermitian(const Eigen::MatrixXcd & sysop){
    if(sysop.rows()!=sysop.cols()){
      std::cerr<<"DiagBB: SysOp is not Hermitian!"<<std::endl;
      std::cerr<<sysop<<std::endl;
      exit(1);
    }
    for(int i=0; i<sysop.rows(); i++){
      for(int j=0; j<sysop.cols(); j++){
        if(std::abs(sysop(i,j)-std::conj(sysop(j,i)))>1e-16){
          std::cerr<<"DiagBB: SysOp is not Hermitian!"<<std::endl;
          std::cerr<<sysop<<std::endl;
          exit(1);
        }
      }
    }
  }
  static bool is_offdiagonal(const Eigen::MatrixXcd & sysop){
    for(int i=0; i<sysop.rows(); i++){
      for(int j=0; j<sysop.cols(); j++){
        if(std::abs(sysop(i,j))>1e-16){
          return true;
        } 
      }
    }
    return false;
  }
  void setup(Parameters &param, const std::string &prefix){

    Operators2x2  op;
    Eigen::MatrixXcd sysop=op.ketbra(1,1);
    bool SysOp_norotate=false;
    if(param.is_specified(add_prefix(prefix,"SysOp"))){
      sysop=param.get_as_operator(add_prefix(prefix,"SysOp"));
    }
    if(param.is_specified(add_prefix(prefix,"SysOp_norotate"))){
      sysop=param.get_as_operator(add_prefix(prefix,"SysOp_norotate"));
      SysOp_norotate=true;
    }
    complain_if_not_Hermitian(sysop);

    if( (!SysOp_norotate) && is_offdiagonal(sysop) ){ 
      hs_rot.setup_by_diagonalizing(sysop);
      std::cout<<"Rotate system part of coupling operator "<<std::endl;
      std::cout<<sysop<<std::endl<<"to "<<std::endl;
      sysop=hs_rot.apply(sysop);
      std::cout<<sysop<<std::endl;
    }
    setup_groups_and_couplings(sysop);

    J=SpectralDensity_Selector(param, add_prefix(prefix,"J"));
    param.add_to(add_prefix(prefix,"omega_max"),"50.");
    EnergyRange Erange(param,prefix);
    omega_min=Erange.omega_min();
    omega_max=Erange.omega_max();
    Ndiscr=param.get_as_size_t(add_prefix(prefix,"override_Ndiscr"),1e6);

    noSubPS= !param.get_as_bool(add_prefix(prefix,"subtract_polaron_shift"), false);
    noSubPS = !param.get_as_bool("subtract_polaron_shift", !noSubPS);
  
    separate_freq = param.get_as_bool("DiagBB_integrate_separate_osc", false);

    temperature=readTemperature(param, prefix);
    E_shift_init=param.get_as_double(add_prefix(prefix,"E_shift_init"),0.);

  }

  DiagBB(const Coupling_Groups & grp, 
                const Eigen::MatrixXcd &couplings_, 
                RealFunctionPtr J_, double temperature_=0., 
                bool noSubPS_=false,
                double omega_min_ = 0., double omega_max_ = 50., 
                double E_shift_init_ =0. ) : J(J_){

    setup(grp, couplings_, J_, temperature_, noSubPS_, omega_min_, omega_max_, E_shift_init_);
  }
  DiagBB(Parameters &param, const std::string &prefix){
    setup(param, prefix);
  }
  DiagBB() : J (RealFunctionPtr_Zero) {}
};

#endif
