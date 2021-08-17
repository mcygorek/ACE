#ifndef MODE_PROPAGATOR_GENERATOR_BOSON_DEFINED_H
#define MODE_PROPAGATOR_GENERATOR_BOSON_DEFINED_H

#include "ModePropagatorGenerator.h"
#include "Parameters.h"
#include "Operators.h"
#include "Operators_Boson.h"
#include "ReadTable.h"
#include "ReadTemperature.h"
#include "ReducedLiouvilleBasis_Boson.h"
#include "ReducedLiouvilleBasis_Boson_FB.h"

class ModePropagatorGenerator_Boson: public ModePropagatorGenerator{
public:

  MPG_Discretization_E_g E_g;
  int M;
  int reduced, reduced_fb;
//  int M_base;
  Eigen::MatrixXcd sysop;
  Parameters gparam;
  bool use_anharmonic; double anharmonic_chi;

  bool use_initial_thermal; double E_shift_init, temperature;
  bool use_initial_coherent; std::complex<double> initial_coherent;
  bool use_polaron_shift;

  bool use_pos_phase; double pos_phase_offset, pos_phase_pos, pos_phase_c;

  virtual std::string name()const{return std::string("Boson");}

  double get_E(int k)const{ return E_g.get_E(k); }
  double get_g(int k)const{ return E_g.get_g(k); }

  virtual int get_N()const{ 
    if(sysop.rows()<=2) return 2; 
    else return sysop.rows();
  }
  virtual double k_label(int k)const{ 
    return get_E(k)/Constants::hbar_in_meV_ps; 
  }

  //print initial boson number per mode
  void print_initial_n(const std::string &fname){
    if(fname=="")return;
    std::ofstream ofs(fname.c_str());
    for(int k=0; k<get_N_modes(); k++){
      Eigen::MatrixXcd init=get_bath_init(k);
      double n=0; for(int r=0; r<init.rows(); r++)n+=r*init(r,r).real();
      ofs<<get_E(k)<<" "<<n<<" "<<k<<std::endl;
    }
  }
  virtual std::vector<Eigen::MatrixXcd> get_env_ops() const{
    std::vector<Eigen::MatrixXcd> mats;
    mats.push_back(Operators_Boson::n(M));
//    mats.push_back(Operators_Boson::n(M)*Operators_Boson::n(M));
    return mats;
  }
  virtual void setup(Parameters &param){
    E_g.setup(param, name());
    set_N_modes(E_g.N);
    setup_skip(param);

    M=param.get_as_size_t(add_name("M"), 2);

    Operators2x2 op;
    sysop=op.ketbra(1,1);   //independent boson model
//    sysop=op.ketbra(0,1);  //Jaynes-Cummings
    if(param.is_specified(add_name("SysOp"))){
      sysop=param.get_as_operator(add_name("SysOp"));
    }
 
    use_pos_phase=param.is_specified(add_name("pos_phase"));
    if(use_pos_phase){
      pos_phase_pos=param.get_as_double(add_name("pos_phase"));
      pos_phase_offset=param.get_as_double(add_name("pos_phase"),0.,0,1);
      pos_phase_c=param.get_as_double(add_name("pos_phase"),Constants::c_in_nm_by_ps,0,2);
    }
    
    use_anharmonic=param.is_specified(add_name("anharmonic_chi"));
    anharmonic_chi=param.get_as_double(add_name("anharmonic_chi"),0.);

    use_polaron_shift=param.get_as_bool(add_name("subtract_polaron_shift"),false);

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
                    Constants::hbar_in_meV_ps);
    }

    print_initial_n(param.get_as_string(add_name("print_initial_n")));
    reduced=param.get_as_size_t(add_name("reduced"),0);
    reduced_fb=param.get_as_size_t(add_name("reduced_fb"),0);
  }

  

  /* multiply position-dependent phase factor exp(i k*r) to coupling
     Note: here we need to distinguish between phyiscal wave vectors k 
     and mode indices, which we now call "j". Actually, the indices correspond
     to an energy discretization (1d):
     hbar*omega= offset +  get_E 
     k=omega/c
 
  */
  Eigen::MatrixXcd get_pos_phase_op(int j)const{
    Eigen::MatrixXcd pos_phase_op=sysop;
    if(!use_pos_phase){
      return pos_phase_op;
    }
    double k=(pos_phase_offset+get_E(j))/(Constants::hbar_in_meV_ps*pos_phase_c);
    std::complex<double> ep=exp(std::complex<double>(0., k*pos_phase_pos));
    pos_phase_op=ep*sysop;

    return pos_phase_op;
  };

  virtual ModePropagatorPtr getModePropagator(int k)const{
    if(k<0||k>=get_N_modes()){
      std::cerr<<"ModePropagatorGenerator_Boson: k<0||k>=get_N_modes()!"<<std::endl; 
      exit(1);
    }

    int sysdim=sysop.rows();
    Eigen::MatrixXcd HB_diag=otimes(
         Eigen::MatrixXcd::Identity(sysdim, sysdim), Operators_Boson::n(M));

    if(use_anharmonic){
      HB_diag-=anharmonic_chi*otimes(
               Eigen::MatrixXcd::Identity(sysdim, sysdim), 
               Operators_Boson::n(M)*Operators_Boson::n(M));
    }

    Eigen::MatrixXcd coupling_op=get_pos_phase_op(k);

    Eigen::MatrixXcd HB_base=
         otimes(coupling_op, Operators_Boson::adagger(M)) +
         otimes(coupling_op.adjoint(), Operators_Boson::a(M));


    Eigen::MatrixXcd HB = get_E(k)*HB_diag 
                         + Constants::hbar_in_meV_ps*get_g(k)*HB_base;
    
    if(use_polaron_shift){
      double pshift=0;
      if(fabs(get_E(k))>1e-12){
        double hg=Constants::hbar_in_meV_ps*get_g(k);
        pshift+=hg*hg/get_E(k);
      }
      HB -= pshift*otimes(sysop.adjoint()*sysop, Operators_Boson::id(M));
    }


    Parameters kparam=gparam;
    ModePropagatorPtr ptr=new ModePropagator(sysdim, get_bath_init(k));
    ptr->FreePropagator::setup(kparam);
    ptr->env_ops=get_env_ops();
    ptr->add_Hamiltonian(HB);
 
    if(reduced_fb>0){
      ptr->rBasis=new ReducedLiouvilleBasis_Boson_FB(get_bath_init(k),reduced_fb);
   
      if(k==0){
        std::cout<<"Using reduced_fb="<<reduced_fb<<" (dim=";
        std::cout<<ptr->rBasis->U.cols()<<"): "<<std::endl;
print_diff_from_ortho(ptr->rBasis->U);
      }
    }else if(reduced>0){
      ptr->rBasis=new ReducedLiouvilleBasis_Boson(get_bath_init(k),reduced);
   
      if(k==0){
        std::cout<<"Using reduced="<<reduced<<" (dim=";
        std::cout<<ptr->rBasis->U.cols()<<"): "<<std::endl;
print_diff_from_ortho(ptr->rBasis->U);
      }
    }
 
    return ptr;

//    return ModePropagatorPtr(new ModePropagator(sysdim,get_bath_init(k),HB,get_env_ops()));
  }

  virtual Eigen::MatrixXcd get_bath_init(int k)const{

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
      return Operators_Boson::equilibrium(HB_diag, temperature);
    }else{
      return Operators_Boson::vacuum(M);
    }
  }

  ModePropagatorGenerator_Boson(Parameters &param){
    setup(param);
  }
};

#endif
