#include "DiagBB.hpp"

#include "TimeGrid.hpp"
#include "discreteFT.hpp"
#include "Function.hpp"
#include "Constants.hpp"
#include "Coupling_Groups.hpp"
#include "Parameters.hpp"
#include "ReadTemperature.hpp"
#include "EnergyRange.hpp"
#include "SpectralDensity_Selector.hpp"
#include "HilbertSpaceRotation.hpp"
#include "Operators.hpp"
#include "Reader.hpp"
#include "ReadTable.hpp"
#include "DummyException.hpp"
#include <cmath>

namespace ACE{

  int DiagBB::sys_dim()const{ return groups.sys_dim(); }

  double DiagBB::get_beta()const {
    if(temperature<1e-10)return 1e12;
    return 1./(kB_in_meV_by_K*temperature);
  }
  double DiagBB::get_coth(double beta, double E_shift, double w){
    if(beta>=0.999e12)return 1.;
    else return coth( 0.5*beta*(hbar_in_meV_ps*w+E_shift) );
  }

  std::complex<double> DiagBB::calculate_K_explicit(int n, double dt){
    double beta=get_beta();
    if(n==0){
      if(noSubPS){
        K0_integrand_class_noSubPS integrand(J.get(), beta, E_shift_init, dt);
        return integrand.integrate(omega_min, omega_max, Ndiscr);
      }else{
        K0_integrand_class integrand(J.get(), beta, E_shift_init, dt);
        return integrand.integrate(omega_min, omega_max, Ndiscr);
      }
    }else{
      if(separate_freq){
        std::complex<double> res=0.; 
        Kn_posdt_integrand_class integrand_posdt(J.get(), beta, E_shift_init, dt);
        res+=integrand_posdt.integrate_times_expI(omega_min, omega_max, Ndiscr, n*dt);
        Kn_negdt_integrand_class integrand_negdt(J.get(), beta, E_shift_init, dt);
        res+=integrand_negdt.integrate_times_expI(omega_min, omega_max, Ndiscr, -n*dt);
        return res;

      }else{
        Kn_integrand_class integrand(J.get(), beta, E_shift_init, dt, n*dt);
        return integrand.integrate(omega_min, omega_max, Ndiscr);
      }
    }
  }

  std::complex<double> DiagBB::calculate_K(int n, double dt){
    double beta=get_beta();
    double damp=1.;
    if(damping>0.){
      damp*=exp(-n*dt/damping);
    }
    if(damping_Gaussian>0.){
      damp*=exp(-pow(n*dt/damping_Gaussian,2));
    }
    
    if(n<K_precalc.size()){
      return K_precalc[n]*damp;
    }else{
      return calculate_K_explicit(n, dt)*damp;
    }
  }

  Eigen::MatrixXcd DiagBB::calculate_expS(int n, double dt){

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

  void DiagBB::print_K(const std::string &fname, int n_max, double dt){
    std::ofstream ofs(fname.c_str());
    for(int n=0; n<n_max; n++){
      double t=n*dt;
      std::complex<double> K=calculate_K(n, dt); ///(dt*dt);
      ofs<<t<<" "<<K.real()<<" "<<K.imag()<<std::endl;
    }
  }

  void DiagBB::read_K_int(const std::string &fname, int n_max, double dt){
std::cout<<"DiagBB::read_K_int called with '"<<fname<<"', "<<n_max<<", "<<dt<<std::endl;
    ReadTable table(fname, 0, 1, 2);
    if(table.size()<2){
      std::cerr<<"DiagBB::read_K_int: file '"<<fname<<"' contains fewer than 2 lines!"<<std::endl;
      throw DummyException();
    }
    double dt2=table[1][0]-table[0][0];
    if(dt2<1e-16){
      std::cerr<<"DiagBB::read_K_int: dt<1e-16!"<<std::endl;
      throw DummyException();
    }
    if( fabs(dt-dt2) > 0.01*dt2 ){
      std::cerr<<"DiagBB::read_K_int: file '"<<fname<<"': dt="<<dt<<" vs. dt2="<<dt2<<"!"<<std::endl;
      throw DummyException();
    }
    if(table.size()<n_max){  
      std::cerr<<"DiagBB::read_K_int: file '"<<fname<<"': rows="<<table.size()<<" < n_max="<<n_max<<"!"<<std::endl;
      throw DummyException();
    }
    K_precalc.clear();
    K_precalc.resize(n_max);
    for(int j=0; j<n_max; j++){
      K_precalc[j].real(table[j][1]);
      K_precalc[j].imag(table[j][2]);
    }
  }


  int DiagBB::estimate_memory_length(int n_max, double dt, double threshold, bool verbose){
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
      //double t=n*dt;
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

  void DiagBB::setup_groups_and_couplings(const Eigen::MatrixXcd &Op){
    groups.setup_from_matrix(Op);

    couplings.resize(groups.Ngrps);
    for(int i=0; i<Op.rows(); i++){
      couplings[groups[i]]=Op(i,i).real();
    }
    if(couplings.size()<2){
      std::cerr<<"DiagBB:set_defaults: couplings.size()<2!"<<std::endl;
      throw DummyException();
    }

    std::cout<<"Coupling groups identified: "<<groups.Ngrps<<":";
    for(size_t i=0; i<groups.grp.size(); ++i){
      std::cout<<" "<<groups.grp[i];
    }
    std::cout<<std::endl;
  }

  void DiagBB::setup(const Coupling_Groups & grp, 
                const Eigen::MatrixXcd &couplings_, 
                RealFunctionPtr J_, double temperature_, 
                bool noSubPS_,
                double omega_min_, double omega_max_, 
                double E_shift_init_  ) {
    J=J_;
    temperature=temperature_;
    groups=grp;
    Ndiscr=1e6;
    noSubPS=noSubPS_;
    separate_freq=false;
    omega_min=omega_min_;
    omega_max=omega_max_;
    damping=0;
    damping_Gaussian=0;
    E_shift_init=E_shift_init_;

//    {std::vector<std::complex<double> > tmp;
//    K_precalc.swap(tmp);}
    K_precalc.clear();

    std::cout<<"DiagBB NGRPS: "<<grp.Ngrps<<std::endl;
    std::cout<<"Couplings Matrix: "<<std::endl<<couplings_<<std::endl;
//    J->print("TEST.dat",omega_min_, omega_max_, 10000);
    
    couplings.resize(grp.Ngrps);
    for(int i=0; i<couplings_.rows(); i++){
      couplings[grp[i]]=couplings_(i,i).real();
    }
    if(couplings.size()<2){
      std::cerr<<"DiagBB:set_defaults: couplings.size()<2!"<<std::endl;
      throw DummyException();
    }    
  }

  void DiagBB::complain_if_not_Hermitian(const Eigen::MatrixXcd & sysop){
    if(sysop.rows()!=sysop.cols()){
      std::cerr<<"DiagBB: SysOp is not Hermitian!"<<std::endl;
      std::cerr<<sysop<<std::endl;
      throw DummyException();
    }
    for(int i=0; i<sysop.rows(); i++){
      for(int j=0; j<sysop.cols(); j++){
        if(std::abs(sysop(i,j)-std::conj(sysop(j,i)))>1e-16){
          std::cerr<<"DiagBB: SysOp is not Hermitian!"<<std::endl;
          std::cerr<<sysop<<std::endl;
          throw DummyException();
        }
      }
    }
  }

  bool DiagBB::is_offdiagonal(const Eigen::MatrixXcd & sysop){
    for(int i=0; i<sysop.rows(); i++){
      for(int j=0; j<sysop.cols(); j++){
        if(std::abs(sysop(i,j))>1e-16){
          return true;
        } 
      }
    }
    return false;
  }

  void DiagBB::setup(Parameters &param, const std::string &prefix){

    Operators op(2);
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

    K_precalc.clear();
    std::string read_K_int_file=param.get_as_string(add_prefix(prefix,"read_K_int"),"");
    if(read_K_int_file!=""){
      TimeGrid tgrid(param);
      read_K_int(read_K_int_file, tgrid.n_tot, tgrid.dt);
      J=RealFunctionPtr_Zero;
    }else{   
      J=SpectralDensity_Selector(param, add_prefix(prefix,"J"));
    }
    param.add_to(add_prefix(prefix,"omega_max"),"50.");
    EnergyRange Erange(param,prefix);
    omega_min=Erange.omega_min();
    omega_max=Erange.omega_max();
    Ndiscr=param.get_as_size_t(add_prefix(prefix,"override_Ndiscr"),1e6);

    noSubPS= !param.get_as_bool(add_prefix(prefix,"subtract_polaron_shift"), true);
    noSubPS = !param.get_as_bool("subtract_polaron_shift", !noSubPS);
  
    separate_freq = param.get_as_bool("DiagBB_integrate_separate_osc", false);

    temperature=readTemperature(param, prefix);
    E_shift_init=param.get_as_double(add_prefix(prefix,"E_shift_init"),0.);

    //additional damping (manual multiplication: K(tau)*e^{-damping tau}:
    damping=param.get_as_double(add_prefix(prefix,"damping"),0.);
    damping_Gaussian=param.get_as_double(add_prefix(prefix,"damping_Gaussian"),0.);

    //precalculate K using FFT
   if(K_precalc.size()<1){
    if(param.get_as_bool("Gaussian_precalc_FFT", true)){
      TimeGrid tgrid(param);
      int n_mem=tgrid.n_mem; if(n_mem<=0)n_mem=tgrid.n_tot;
//std::cout<<"TEST: tgrid.n_mem="<<tgrid.n_mem<<" tgrid.n_tot="<<tgrid.n_tot<<" DiagBB: n_mem="<<n_mem<<std::endl; 
      if(!param.is_specified("override_Ndiscr")){
        Ndiscr=tgrid.n_mem*4;
        if(Ndiscr<1e6)Ndiscr=1e6;
      }
      precalc_FFT(n_mem, tgrid.dt);
    }
   }

  }


// Conventional integrands:
    std::complex<double> DiagBB::K0_integrand_class::f(double w)const{
      if(w*w<1e-20)return 0.;
      return J->f(w)/(w*w)*
         std::complex<double>(DiagBB::get_coth(beta,E_shift,w)*(1.-cos(w*dt)) , sin(w*dt));
    }

    std::complex<double> DiagBB::K0_integrand_class_noSubPS::f(double w)const{
      if(w*w<1e-20)return 0.;
      return J->f(w)/(w*w) *
std::complex<double>(DiagBB::get_coth(beta,E_shift,w)*(1.-cos(w*dt)) , sin(w*dt)-w*dt );
    }

    std::complex<double> DiagBB::Kn_integrand_class::f(double w)const{
      if(w*w<1e-20)return 0.;
      return 2.*J->f(w)/(w*w)*(1.-cos(w*dt)) *
std::complex<double>(DiagBB::get_coth(beta,E_shift,w)*cos(w*tau) , -sin(w*tau) );
    }

// new ones: integrate cos and sin parts separately:

    std::complex<double> DiagBB::K0_const_integrand_class::f(double w)const{
      if(w*w<1e-20)return 0.;
      return J->f(w)/(w*w)*DiagBB::get_coth(beta,E_shift,w);
    }

    std::complex<double> DiagBB::K0_const_integrand_class_noSubPS::f(double w)const{
      if(w*w<1e-20)return 0.;
      return J->f(w)/(w*w)*std::complex<double>(DiagBB::get_coth(beta,E_shift,w), -w*dt );
    }

    std::complex<double> DiagBB::K0_posdt_integrand_class::f(double w)const{
      if(w*w<1e-20)return 0.;
      return J->f(w)/(w*w)* 1./2.*( -DiagBB::get_coth(beta,E_shift,w) + 1.);
    }
  
    std::complex<double> DiagBB::K0_negdt_integrand_class::f(double w)const{
      if(w*w<1e-20)return 0.;
      return J->f(w)/(w*w)* 1./2.*( -DiagBB::get_coth(beta,E_shift,w) - 1.);
    }

    std::complex<double> DiagBB::Kn_posdt_integrand_class::f(double w)const{
      if(w*w<1e-20)return 0.;
      return 2.*J->f(w)/(w*w)*(1.-cos(w*dt))*1./2.*(DiagBB::get_coth(beta,E_shift,w)-1. );
    }

    std::complex<double> DiagBB::Kn_negdt_integrand_class::f(double w)const{
      if(w*w<1e-20)return 0.;
      return 2.*J->f(w)/(w*w)*(1.-cos(w*dt))*1./2.*(DiagBB::get_coth(beta,E_shift,w)+1. );
    }

//Fast fourier routine:
  
  void DiagBB::precalc_FFT(int n_max, double dt){

    if(n_max<1){
      std::cerr<<"DiagBB::precalc_FFT: n_max="<<n_max<<"<1!"<<std::endl;
      throw DummyException();
    }

    double omega_range=omega_max-omega_min;

    double my_omega_range=(2.*M_PI)/dt;
    int N_new=pow(2, ceil(log((double)Ndiscr)/log(2.)) );

    std::cout<<"DiagBB::precalc_FFT: n_max="<<n_max;
    std::cout<<" dt="<<dt;
    std::cout<<" Ndiscr="<<Ndiscr<<" -> "<<N_new;
    std::cout<<" orig_omega_range="<<omega_range;
    std::cout<<" my_omega_range="<<my_omega_range;
    std::cout<<std::endl;
  
    double domega=my_omega_range/N_new;

    std::complex<double> *in_array=new std::complex<double>[N_new];

    for(size_t i=0; i<N_new; i++){
      double w=omega_min+i*domega;
      if(w<=omega_max){
        in_array[i]=(w*w<1e-20) ? 0. : J->f(w)/(w*w)*DiagBB::get_coth(get_beta(),E_shift_init,w);
      }else{
        in_array[i]=0.;
      }
    }
    std::complex<double> *out1_array=new std::complex<double>[N_new];
    discreteFT(in_array, out1_array, N_new, 1, +1);
    for(size_t i=0; i<N_new; i++){
      out1_array[i]*=exp(std::complex<double>(0,(+1)*i*omega_min*dt));
    }


    for(size_t i=0; i<N_new; i++){
      double w=omega_min+i*domega;
      if(w<=omega_max){
        in_array[i]=(w*w<1e-20) ? 0. : J->f(w)/(w*w);
      }else{
        in_array[i]=0.;
      }
    }
    std::complex<double> *out2_array=new std::complex<double>[N_new];
    discreteFT(in_array, out2_array, N_new, 1, +1);
 
    for(size_t i=0; i<N_new; i++){
#ifdef DIAGBB_FFT_CORRECTION
      double theta=2.*M_PI*i/N_new;
      out2_array[i]*=FFT_trapezoidal_correction_W(theta); 
      out2_array[i]+=FFT_trapezoidal_correction_a0(theta)*in_array[0];
      out2_array[i]+=exp(std::complex<double>(0,((+1)*i*N_new*domega*dt))) 
           *std::conj(FFT_trapezoidal_correction_a0(theta))*in_array[N_new-1];
#endif
      out2_array[i]*=exp(std::complex<double>(0,(+1)*i*omega_min*dt));
    }
    delete[] in_array;

    K_precalc=std::vector<std::complex<double> >(n_max,0.);

    for(int n=1; n<n_max && n<N_new-1; n++){
      K_precalc[n].real(domega*(2.*out1_array[n]-out1_array[n-1]-out1_array[n+1]).real());
      K_precalc[n].imag(-domega*(2.*out2_array[n]-out2_array[n-1]-out2_array[n+1]).imag());
    }
    K_precalc[0]=calculate_K_explicit(0, dt);

    delete[] out1_array;
    delete[] out2_array;
  }

}//namespace
