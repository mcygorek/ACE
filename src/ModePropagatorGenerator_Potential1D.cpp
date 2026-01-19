#include "ModePropagatorGenerator_Potential1D.hpp"
#include "ModePropagatorGenerator.hpp"
#include "Operators.hpp"
#include "Operators_Boson.hpp"
#include "Potential1D.hpp"
#include "Equilibrium.hpp"
#include "ReadTemperature.hpp"
#include "otimes.hpp"
#include "LiouvilleTools.hpp"
#include "ReducedLiouvilleBasis_Boson.hpp"
#include "ReducedLiouvilleBasis_Boson_FB.hpp"

namespace ACE{

  //print initial boson number per mode
  void ModePropagatorGenerator_Potential1D::print_initial_n(const std::string &fname)const{
    if(fname=="")return;
    std::ofstream ofs(fname.c_str());
    for(int k=0; k<get_N_modes(); k++){
      Eigen::MatrixXcd init=get_bath_init(k);
      double n=0; for(int r=0; r<init.rows(); r++)n+=r*init(r,r).real();
      ofs<<get_E(k)<<" "<<n<<" "<<k<<std::endl;
    }
  }

  EnvironmentOperators ModePropagatorGenerator_Potential1D::get_env_ops(int k)const{
    std::vector<Eigen::MatrixXcd> mats;
    mats.push_back(Eigen::MatrixXcd::Identity(M,M));
    mats.push_back(Operators_Boson::n(M));
//    mats.push_back(Operators_Boson::n(M)*Operators_Boson::n(M));
    return EnvironmentOperators(mats);
  }

  void ModePropagatorGenerator_Potential1D::setup(Parameters &param){
    E_g.setup(param, name());
    set_N_modes(E_g.N);
    setup_skip(param);

    p1d.setup(param, name());

    M=p1d.M;

    Operators op(2);
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

  ModePropagatorPtr ModePropagatorGenerator_Potential1D::get_ModePropagator(int k)const{
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
                         + hbar_in_meV_ps*get_g(k)*HB_base;
    

    if(use_renorm_x){
      Eigen::MatrixXcd b_init=get_bath_init(k);
      std::complex<double> Xinit=(p1d.X*b_init).trace();

      Eigen::MatrixXcd polaron_op = hbar_in_meV_ps * get_g(k)* sysop * 0.5 * Xinit;
      HB -= otimes(polaron_op+polaron_op.adjoint(), Eigen::MatrixXcd::Identity(M,M));
    }


    Parameters kparam=gparam;
    ModePropagatorPtr ptr=std::make_shared<ModePropagator>(sysdim, get_bath_init(k));
    ptr->FreePropagator::setup(kparam);
    ptr->env_ops=get_env_ops(k);
    ptr->add_Hamiltonian(HB);
 
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
 
    return ptr;

//    return ModePropagatorPtr(new ModePropagator(sysdim,get_bath_init(k),HB,get_env_ops()));
  }

  Eigen::MatrixXcd ModePropagatorGenerator_Potential1D::get_bath_init(int k)const{
    return Nlevel_Equilibrium(get_E(k)*p1d.E, temperature);
  }


}//namespace
