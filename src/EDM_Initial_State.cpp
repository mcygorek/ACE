#include "EDM_Initial_State.hpp"
#include "Operators.hpp"
#include "LiouvilleTools.hpp"
#include "DummyException.hpp"
#include "ReadExpression.hpp"
#include "ReaderBasics.hpp"

namespace ACE {

void EDM_Initial_State::set_default_closure(int r){
  const EDM_Index &Ldim=rhos.Ldim;
  Ldim.check_in_range(r);

  if(rhos.closures.size()<Ldim.size()){
    rhos.closures.resize(Ldim.size());
  }
  Eigen::VectorXcd &q=rhos.closures[r];
  q=Eigen::VectorXcd::Zero(Ldim[r]);

  int dimH=sqrt(Ldim[r]);
  if(dimH*dimH==Ldim[r]){ //site can be system
    for(int i=0; i<dimH; i++){
      q(i*dimH+i)=1.;
    }
  }else{
    if(Ldim[r]>0)q(0)=1.;
  }
}

void EDM_Initial_State::setup(Parameters &param){
  rhos.clear();
  double thr=param.get_as_double("initial_threshold",1e-16);

  int N_sites=param.get_as_size_t("N_sites",1);
  rhos.Ldim.list=std::vector<int>(N_sites,1);

//std::cout<<"N_sites="<<N_sites<<std::endl;

  std::vector<Eigen::VectorXcd> rho_initial(N_sites);
  for(int r=0; r<N_sites; r++){
    std::string key="S"+int_to_string(r)+"_initial";
//std::cout<<"key: '"<<key<<"'"<<std::endl;
    rho_initial[r]=H_Matrix_to_L_Vector(param.get_as_operator(key,Eigen::MatrixXcd::Identity(1,1)));
  }
  for(int r=0; r<N_sites; r++){
    std::cout<<rho_initial[r].transpose()<<std::endl;
  }

//  rho_initial.push_back(H_Matrix_to_L_Vector(Operators(2).ketbra(1,1)));
//  rho_initial.push_back(H_Matrix_to_L_Vector(Operators(2).ketbra(0,0)));
//  Eigen::MatrixXcd M(2,2); M<<0.5, 0.5 , 0.5 ,0.5;
//  rho_initial.push_back(H_Matrix_to_L_Vector(M));
  rhos.set_from_product(rho_initial, thr);

  for(int r=0; r<rhos.Ldim.size(); r++){
    set_default_closure(r);
  } 


 
}

}
