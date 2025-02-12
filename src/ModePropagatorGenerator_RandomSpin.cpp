#include "ModePropagatorGenerator_RandomSpin.hpp"
#include "ModePropagatorGenerator.hpp"
#include "Parameters.hpp"
#include "Operators.hpp"
#include "otimes.hpp"
#include "RandomDirection.hpp"

namespace ACE{

  bool ModePropagatorGenerator_RandomSpin::compare_abs_smaller(
                          const double &p1, const double &p2){
    return fabs(p1)<fabs(p2);
  } 
  bool ModePropagatorGenerator_RandomSpin::compare_abs_larger(
                          const double &p1, const double &p2){
    return fabs(p1)>fabs(p2);
  } 
  bool ModePropagatorGenerator_RandomSpin::compare_dir_x(
                 const Eigen::Vector3d &v1, const Eigen::Vector3d &v2){
    return v1(0)<v2(0);
  }
  bool ModePropagatorGenerator_RandomSpin::compare_dir_z(
                 const Eigen::Vector3d &v1, const Eigen::Vector3d &v2){
    return v1(2)<v2(2);
  }

  EnvironmentOperators ModePropagatorGenerator_RandomSpin::get_env_ops(int k) const{
    std::vector<Eigen::MatrixXcd> mats(4);
    mats[0]=Eigen::MatrixXcd::Identity(2,2);
    mats[1]=0.5*sigma_x();
    mats[2]=0.5*sigma_y();
    mats[3]=0.5*sigma_z();
    return EnvironmentOperators(mats);
  }

  void ModePropagatorGenerator_RandomSpin::setup(int Nmod, double J_max, double J_min, size_t seed){

    set_N_modes(Nmod);

    std::srand(seed);

    J.clear();
    J.resize(Nmod, J_max);
    
    for(int i=0; i<Nmod; i++){
      J[i]=J_min+random_double()*(J_max-J_min);
    }

    initialDirs.clear();
    initialDirs.resize(Nmod);
    for(int i=0; i<Nmod; i++){
      initialDirs[i]=RandomDirection(); 
    }
  }

  void ModePropagatorGenerator_RandomSpin::setup(Parameters &param){
    setup_default(param);
    int N_modes=param.get_as_size_t("RandomSpin_N_modes", 0);

    if(!param.is_specified("RandomSpin_J_max")){
      std::cerr<<"Please specify 'RandomSpin_J_max'!"<<std::endl;  
      exit(1);
    }
    double J_max=param.get_as_double("RandomSpin_J_max",0);
    double J_min=param.get_as_double("RandomSpin_J_min",0);
    size_t seed=param.get_as_size_t("RandomSpin_seed",1);


    setup(N_modes, J_max, J_min, seed);

 
    std::vector<Eigen::Vector3d> B_eff_init(N_modes,Eigen::Vector3d::Zero());
    { //B_eff
      std::vector<double> dvec(3,0.);
      if(param.is_specified("RandomSpin_B_eff")){
        dvec=param.get_row_doubles("RandomSpin_B_eff", 0, 3);
      }
      Eigen::Vector3d evec; for(int i=0; i<3; i++)evec(i)=dvec[i];

      B_eff=std::vector<Eigen::Vector3d>(N_modes, evec);
      B_eff_init=std::vector<Eigen::Vector3d>(N_modes, evec);
    }
    { //B_eff_init
      std::vector<double> dvec(3,0.);
      if(param.is_specified("RandomSpin_B_init")){
        dvec=param.get_row_doubles("RandomSpin_B_init", 0, 3);
      }
      Eigen::Vector3d evec; for(int i=0; i<3; i++)evec(i)=dvec[i];

      B_eff_init=std::vector<Eigen::Vector3d>(N_modes, evec);
    }




    double T=param.get_as_double("RandomSpin_T",0.);
    double beta=1./(kB_in_meV_by_K*T);

    if(param.is_specified("RandomSpin_set_initial_dir")){
      std::vector<double> dvec=param.get_row_doubles("RandomSpin_set_initial_dir", 0, 3);
      Eigen::Vector3d vvec; for(int i=0; i<3; i++)vvec(i)=dvec[i];
      vvec.normalize();

      initialDirs.clear();
      initialDirs.resize(get_N_modes());
      for(size_t i=0; i<initialDirs.size(); i++){
        initialDirs[i]=vvec;
      }
    }else if(T>1e-10){
      initialDirs.clear();
      initialDirs.resize(get_N_modes());

      for(int i=0; i<get_N_modes(); i++){
        initialDirs[i]=RandomDirection(); 

        //Energy wrt. maximal possible Energy:
        double E=0.5*(B_eff_init[i].dot(initialDirs[i]) - B_eff_init[i].norm());  
        if(random_double()>exp(beta*E)){
          i--;
          continue;
        }
      }

    }
 
//sort modes:
    if(param.get_as_bool("RandomSpin_sort_dir_x")){
      sort(initialDirs.begin(), initialDirs.end(), compare_dir_x);
    }else if(param.get_as_bool("RandomSpin_sort_dir_z")){
      sort(initialDirs.begin(), initialDirs.end(), compare_dir_z);
    }else if(param.get_as_bool("RandomSpin_sort_J_smaller")){
      sort(J.begin(), J.end(), compare_abs_smaller);
    }else if(param.get_as_bool("RandomSpin_sort_J_larger")){
      sort(J.begin(), J.end(), compare_abs_larger);
    }


// Print initial spins:
    std::string print_initial=param.get_as_string("RandomSpin_print_initial");
    std::ofstream print_ofs;
    if(print_initial!="")print_ofs.open(print_initial.c_str());

    std::cout<<"RandomSpin Modes: (J, dirx, diry, dirz):"<<std::endl;
    for(int i=0; i<get_N_modes(); i++){
      std::cout<<J[i]<<": "<<initialDirs[i](0)<<" ";
      std::cout<<initialDirs[i](1)<<" "<<initialDirs[i](2)<<std::endl;
      
      if(print_ofs.good()){
        print_ofs<<J[i]<<": "<<initialDirs[i](0)<<" ";
        print_ofs<<initialDirs[i](1)<<" "<<initialDirs[i](2)<<std::endl;
      }
    }
 
    Eigen::Vector3d totalenvspin=Eigen::Vector3d::Zero();
    for(size_t i=0; i<initialDirs.size(); i++)totalenvspin+=0.5*initialDirs[i];
    std::cout<<"Total environment spin: ";
    std::cout<<totalenvspin(0)<<" "<<totalenvspin(1)<<" "<<totalenvspin(2)<<std::endl;
  }

  ModePropagatorPtr ModePropagatorGenerator_RandomSpin::getModePropagator(int k)const{
    if(k<0||k>=get_N_modes()){
      std::cerr<<"ModePropagatorGenerator_RandomSpin: k<0||k>=get_N_modes()!"<<std::endl; 
      exit(1);
    }

   
    Operators op(2);
    Eigen::MatrixXcd HB_base = 0.25 * J[k] * (
        otimes(sigma_x(), sigma_x())
      + otimes(sigma_y(), sigma_y())
      + otimes(sigma_z(), sigma_z()));

    Eigen::MatrixXcd HB_B = 0.5 *  (
        B_eff[k](0) * otimes(op.id(), sigma_x())
      + B_eff[k](1) * otimes(op.id(), sigma_y())
      + B_eff[k](2) * otimes(op.id(), sigma_z()));

//    std::cout<<"HB_base["<<k<<"]:"<<std::endl<<HB_base<<std::endl;

    return ModePropagatorPtr(new ModePropagator(2,get_bath_init(k),HB_base+HB_B,get_env_ops(k)));
  }

  Eigen::MatrixXcd ModePropagatorGenerator_RandomSpin::get_bath_init(int k)const{
    Operators op(2);
    Eigen::MatrixXcd mat= 0.5*( op.id() + initialDirs[k](0) * sigma_x() +
          initialDirs[k](1) * sigma_y() + initialDirs[k](2) * sigma_z() );
 
//    mat=0.5*(op.id()+sigma_x());
//    mat=0.5*(op.id()+sigma_y());
 
    std::cout<<"k: "<<k;
    std::cout<<" <sigma_x>: "<<(sigma_x()*mat).trace();
    std::cout<<" <sigma_y>: "<<(sigma_y()*mat).trace();
    std::cout<<" <sigma_z>: "<<(sigma_z()*mat).trace();
    std::cout<<std::endl;
    return mat;
  }

}//namespace
