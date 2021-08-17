#ifndef MODE_PROPAGATOR_GENERATOR_POTENTIAL1D_DEFINED_H
#define MODE_PROPAGATOR_GENERATOR_POTENTIAL1D_DEFINED_H

#include "ModePropagatorGenerator.h"
#include "Operators.h"
#include "Operators_Boson.h"
#include "Potential1D.h"
#include "Equilibrium.h"
#include "ReadTemperature.h"


class ModePropagatorGenerator_Potential1D: public ModePropagatorGenerator{
public:

  MPG_Discretization_E_g E_g;

  Potential1D p1d;
  int M;
  int reduced, reduced_fb;
//  int M_base;
  Eigen::MatrixXcd sysop;
  Parameters gparam;

  double temperature;
//  bool use_polaron_shift; 
  bool use_renorm_x;


  virtual std::string name()const{return std::string("Potential1D");}

  double get_E(int k)const{ return E_g.get_E(k); }
  double get_g(int k)const{ return E_g.get_g(k); }

  virtual int get_N()const{ 
    if(sysop.rows()<=2) return 2; 
    else return sysop.rows();
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

    p1d.setup(param, name());

    M=p1d.M;

    Operators2x2 op;
    sysop=op.ketbra(1,1);   //independent boson model
//    sysop=op.ketbra(0,1);  //Jaynes-Cummings
    if(param.is_specified(add_name("SysOp"))){
      sysop=param.get_as_operator(add_name("SysOp"));
    }
 

 //   use_polaron_shift=param.get_as_bool(add_name("subtract_polaron_shift"),false);
    use_renorm_x=param.get_as_bool(add_name("renorm_x"),false);

    std::string print_E_g=param.get_as_string(add_name("print_E_g"));
    if(print_E_g!=""){
      std::ofstream ofs(print_E_g.c_str());
      for(int i=0; i<E_g.N; i++){
        ofs<<E_g.get_E(i)<<" "<<E_g.get_g(i)<<std::endl;
      }
    }
 
    gparam.add_from_prefix(add_name("Propagator"),param);


    temperature=readTemperature(param,name());


    print_initial_n(param.get_as_string(add_name("print_initial_n")));

    reduced=param.get_as_size_t(add_name("reduced"),0);
    reduced_fb=param.get_as_size_t(add_name("reduced_fb"),0);
  }

  

  virtual ModePropagatorPtr getModePropagator(int k)const{
    if(k<0||k>=get_N_modes()){
      std::cerr<<"ModePropagatorGenerator_Potential1D: k<0||k>=get_N_modes()!"<<std::endl; 
      exit(1);
    }

    int sysdim=sysop.rows();
    Eigen::MatrixXcd HB_diag=otimes(
         Eigen::MatrixXcd::Identity(sysdim, sysdim), p1d.E.asDiagonal());


    Eigen::MatrixXcd HB_base=
         otimes(sysop+sysop.adjoint(), 0.5*p1d.X);


    Eigen::MatrixXcd HB = get_E(k)*HB_diag 
                         + Constants::hbar_in_meV_ps*get_g(k)*HB_base;
    

    if(use_renorm_x){
      Eigen::MatrixXcd b_init=get_bath_init(k);
      std::complex<double> Xinit=(p1d.X*b_init).trace();

      Eigen::MatrixXcd polaron_op = Constants::hbar_in_meV_ps * get_g(k)* sysop * 0.5 * Xinit;
      HB -= otimes(polaron_op+polaron_op.adjoint(), Eigen::MatrixXcd::Identity(M,M));
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
//      Eigen::MatrixXcd HB_diag=get_E(k)*p1d.E.asDiagonal();
    return Equilibrium(get_E(k)*p1d.E, temperature);
  }

  ModePropagatorGenerator_Potential1D(Parameters &param){
    setup(param);
  }
};

#endif
