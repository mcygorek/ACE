
#include "Potential1D.hpp"
#include "RealFunction_Interpolate.hpp"
#include "Parameters.hpp"

namespace ACE{

  void Potential1D::calculate(const std::string &printEV){
    int N=vx.size(); 
    if(N<M){std::cerr<<"Potential1D::calculate: N<M!"<<std::endl;exit(1);}

    double fac=-ddx2_scale/(dx*dx);

    Eigen::VectorXd h_d(N);
    for(int i=0; i<N; i++){
      h_d(i)=vx[i] - 2.*fac;
    }
    Eigen::VectorXd h_od(N-1);
    for(int i=0; i<N-1; i++){
      h_od(i)=fac;
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd > solver;
    solver.computeFromTridiagonal(h_d, h_od);
    
    E=E_scale*solver.eigenvalues().head(M);

    X=Eigen::MatrixXd::Zero(M,M);
    X2=Eigen::MatrixXd::Zero(M,M);
    for(int n=0; n<M; n++){
      for(int m=0; m<M; m++){
        for(int i=0; i<N; i++){
          double x=x_scale*(xa+i*dx);
          X(n,m)+=solver.eigenvectors()(i,n)*x*solver.eigenvectors()(i,m);
          X2(n,m)+=solver.eigenvectors()(i,n)*x*x*solver.eigenvectors()(i,m);
        }
      }
    }
   
    X*=x_scale*x_coupling;
    X2*=x_scale*x_scale*x_coupling*x_coupling;


    if(printEV!=""){
      std::ofstream ofsEV(printEV.c_str());
      for(int i=0; i<N; i++){
        ofsEV<<xa+i*dx<<" "<<vx[i];
        for(int n=0; n<M; n++){
          ofsEV<<" "<<solver.eigenvectors()(i,n);
        }
        ofsEV<<std::endl;
      }
    }
  }

  std::string Potential1D::add_name(const std::string &name, const std::string &str){
    if(name=="")return str;
    return std::string(name + "_" + str);
  }

  void Potential1D::setup(Parameters &param, const std::string &name){
    
    xa=param.get_as_double(add_name(name,"xa"),0.);  
    double xe=param.get_as_double(add_name(name,"xe"),10.);  
    int Ndiscr=param.get_as_int(add_name(name,"Ndiscr"),1000);  
    dx=(xe-xa)/((double)Ndiscr-1.);

    M=param.get_as_int(add_name(name,"M"),0);  
    if(param.is_empty())M=2;

    E_scale=param.get_as_double(add_name(name,"E_scale"),1.);  
    x_scale=param.get_as_double(add_name(name,"x_scale"),1.);  
    ddx2_scale=1.;
    x_coupling=sqrt(2.);
  
    std::string type=param.get_as_string(add_name(name,"type"));
    if(type==""||type=="harmonic"||type=="HO"){
      double E=param.get_as_double(add_name(name,"E"),1.);
      vx.clear(); vx.resize(Ndiscr);
      for(int i=0; i<Ndiscr; i++){
        double x=xa+i*dx;
        vx[i]=0.5*E*E*x*x;
      }
      ddx2_scale=0.5;
      x_coupling=sqrt(2.*E);
    }else if(type=="morse"||type=="morse_scaleHO"||type=="morse_matchLowestHO"){
      double Lambda=param.get_as_double_check(add_name(name,"Lambda"));
      if(Lambda<=1e-20){
        std::cerr<<"Potential1D: Lambda must be strictly positive!"<<std::endl;
        exit(1);
      }
      if(type=="morse_scaleHO")E_scale/=(2.*Lambda);
      if(type=="morse_matchLowestHO")E_scale/=(2.*Lambda*(1.-1./Lambda));

      x_coupling=sqrt(2.*Lambda);
      vx.clear(); vx.resize(Ndiscr);
      for(int i=0; i<Ndiscr; i++){
        double x=xa+i*dx;
        vx[i]=Lambda*Lambda*(exp(-2.*x)-2.*exp(-x));
      }
      if(M<1){
        int M2=Lambda+0.5;
        while(!((double)M2<(Lambda+0.5))){
          M2--;
        }
        M=M2;
      }
    }else{
      std::cerr<<"Potential1D: type '"<<type<<"' not known!"<<std::endl;
      exit(1);
    }



    if(M<1){
      std::cerr<<"Potential1D 'M' must be larger than 0!"<<std::endl;
      exit(1);
    }
    if(M>Ndiscr){
      std::cerr<<"Potential1D: M>Ndiscr!"<<std::endl;
      exit(1);
    }

    calculate(param.get_as_string(add_name(name,"print_EV")));


    if(param.get_as_bool(add_name(name,"clear_x_diag"),false)){
      for(int i=0; i<X.rows(); i++){
        X(i,i)=0.;
      }
    }


    std::string print_E=param.get_as_string(add_name(name,"print_E"));
    if(print_E!=""){
      std::ofstream ofs(print_E.c_str());
      for(int i=0; i<E.size(); i++){ 
        ofs<<E(i)<<std::endl;
      }
    }
    
    std::string print_X=param.get_as_string(add_name(name,"print_X"));
    if(print_X!=""){
      std::ofstream ofs(print_X.c_str());
      ofs<<X<<std::endl;
    }

    std::string print_X2=param.get_as_string(add_name(name,"print_X2"));
    if(print_X2!=""){
      std::ofstream ofs(print_X2.c_str());
      ofs<<X2<<std::endl;
    }
  }
  Potential1D::Potential1D(){
    Parameters param;
    setup(param);
  }

}//namespace
