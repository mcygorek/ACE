#include "PCH.hpp"
#include "Equilibrium.hpp"
#include <Eigen/Core>
#include <iostream>
#include "Constants.hpp"
#include "Operators_Boson.hpp"

namespace ACE{

Eigen::MatrixXcd Nlevel_Equilibrium(const Eigen::VectorXcd &E, double T){ 
  int dim=E.size();
  if(T<1e-12){  //zero temperature:
    //determine degeneracy:
    int n=1;
    for(int i=1; i<dim; i++){
      if(abs(E(i)-E(0))<1e-12)n++;
    }
    Eigen::MatrixXcd mat=Eigen::MatrixXcd::Zero(dim,dim);
    for(int i=0; i<n; i++){
      mat(i,i)=1./( (double)n);
    }
    return mat;
  }else if(T>=1e6){ // infinite temperature
    Eigen::MatrixXcd mat=1./((double)dim)*Eigen::MatrixXcd::Identity(dim,dim);
    return mat;
  }else{ // finite temperature
    double beta=1./(kB_in_meV_by_K*T);
    Eigen::MatrixXcd mat=Eigen::MatrixXcd::Zero(dim,dim);
    double norm=0.;
    for(int i=0; i<dim; i++){
      double contrib=exp(-beta*(E(i).real()-E(0).real()));
      mat(i,i)=contrib;
      norm+=contrib;
    }
    for(int i=0; i<dim; i++){
      mat(i,i)/=norm;
    }
    return mat;
  }
}


Eigen::MatrixXcd Boson_Equilibrium(int n_max, double dE, double T){ //x=(E-mu)
    if(T<1e-12){  //effectively: step function
      if(dE<0){
        return 1./(double)n_max*Eigen::MatrixXcd::Identity(n_max, n_max);
      }else{
        return Operators_Boson::vacuum(n_max);
      }
    }   

    double x=dE/(kB_in_meV_by_K*T);
    if(x<1e-8){
      return 1./(double)n_max*Eigen::MatrixXcd::Identity(n_max, n_max);
    }
    Eigen::VectorXd v(n_max);
    double norm=0;
    for(int i=0; i<n_max; i++){
      double e=exp(-x*i);
      norm+=e;
      v(i)=e;
    }
    for(int i=0; i<n_max; i++){
      v(i)/=norm;
    }
    return v.asDiagonal();
}

// Equilibrium w.r.t. to Hamiltonian matrix
Eigen::MatrixXcd Boson_Equilibrium(Eigen::MatrixXcd H, double T, double E_shift){ 
    int n_max=H.rows();
    if(H.rows()!=H.cols()){
      std::cerr<<"Operators_Boson::equilibrium: H.rows()!=H.cols()!"<<std::endl;
      exit(1);
    }
    if(H.rows()<1){
      std::cerr<<"Operators_Boson::equilibrium: H.rows()<1!"<<std::endl;
      exit(1);
    }

    if(T<1e-20){
      Eigen::MatrixXcd ret=Eigen::MatrixXcd::Zero(H.rows(),H.rows());
      ret(0,0)=1.;
      return ret;
    }

    for(int i=0; i<H.rows(); i++){
      H(i,i)+=E_shift;
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(H);
    Eigen::MatrixXcd V=solver.eigenvectors();

    int nrltz=0;
    for(int i=0; i<n_max; i++){
      if(solver.eigenvalues()(i)<0.)nrltz++;
    }
    if(nrltz>0){
      std::cerr<<"Operators_Boson::equilibrium: "<<nrltz<<" of "<<n_max<<" eigenvalues smaller than 0!"<<std::endl;
      std::cerr<<solver.eigenvalues()<<std::endl;
      exit(1);
    }

    if(T<1e-12){  //effectively: step function
      return Operators_Boson::vacuum(n_max);
    }   

    double beta=1./(kB_in_meV_by_K*T);

    Eigen::VectorXd v(n_max);
    double norm=0;
    for(int i=0; i<n_max; i++){
      double e=exp(-beta*solver.eigenvalues()(i));
      norm+=e;
      v(i)=e;
    }
    for(int i=0; i<n_max; i++){
      v(i)/=norm;
    }
 
    Eigen::MatrixXcd ret=Eigen::MatrixXcd::Zero(n_max,n_max);
    for(int i=0; i<n_max; i++){
       ret += v(i) * V.col(i)*( V.col(i).adjoint());
    }
    return ret;
}

double Boson_Equilibrium_invert(int nmax, double nav){
  auto n = [nmax] (double x) -> double {
      double tmp=0, norm=0; 
      for(int i=0; i<nmax; i++){
        double p=exp(-x*i); tmp+=i*p; norm+=p;
      }
      return tmp/norm;
    };

  double a=100, b=1e-4;
  double na=n(a); double nb=n(b);

  if(nav<=na){
    std::cerr<<"Boson_Equilibrium_invert: nav="<<nav<<" <= "<<na<<"!"<<std::endl;
    exit(1);
  }
  if(nav>nb){
    std::cerr<<"Boson_Equilibrium_invert: nav="<<nav<<" > "<<nb<<"!"<<std::endl;
    exit(1);
  }
  
  for(int loop=0; loop<100; loop++){
    double mid=(a+b)/2.;
    double nmid=n(mid);
//    std::cout<<"loop="<<loop<<" mid="<<mid<<" nmid="<<nmid<<std::endl;
    if(fabs(nmid-nav)<1e-10)break;
    if(nav>nmid){ 
      a=mid;
    }else{
      b=mid;
    }
  }

  return (a+b)/2.;
}

}//namespace
