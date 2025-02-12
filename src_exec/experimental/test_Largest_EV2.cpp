#include "Largest_EV.hpp"
#include "Parameters.hpp"
#include "Timings.hpp"
#include <iostream>
#include <fstream>
#include <Eigen/Dense>

using namespace ACE;

int main(int args, char **argv){

//OBSERVATION: Yes, the power iteration is accelerated
//             BUT: it typically converges to a different EV!!!

  Parameters param(args, argv, true);

  int dim=param.get_as_int("dim", 100);
  int iter=param.get_as_int("iter", 1000);
  double epsilon=param.get_as_double("epsilon", 1e-15);
  bool reortho=param.get_as_bool("reortho",false);

//  Eigen::MatrixXd B = Eigen::MatrixXd::Random(dim,dim);
//  Eigen::MatrixXcd A = B;
  Eigen::MatrixXcd A = Eigen::MatrixXcd::Random(dim,dim);

  //Exact:
  time_point time_exact1=now();
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solve(A);
  time_point time_exact2=now();
//  int imax=0; for(int i=1; i<dim; i++){if(std::abs(solve.eigenvalues()(i))>std::abs(solve.eigenvalues()(imax))){imax=i;}}
//  int imax2=0; for(int i=1; i<dim; i++){if((std::abs(solve.eigenvalues()(i))>std::abs(solve.eigenvalues()(imax2)) && (i!=imax))){imax2=i;}}
  std::vector<std::pair<std::complex<double>,int> > srt(dim);
  for(int i=0; i<dim; i++){
    srt[i]=std::make_pair(solve.eigenvalues()(i), i);
  }
  std::sort(srt.begin(),srt.end(), [](const std::pair<std::complex<double>,int> &a, const std::pair<std::complex<double>,int> &b){ return std::real(a.first)>std::real(b.first); });
//  std::sort(srt.begin(),srt.end(), [](const std::pair<std::complex<double>,int> &a, const std::pair<std::complex<double>,int> &b){ return std::norm(a.first)>std::norm(b.first); });

  std::cout<<"Exact eigenvalues: ";
  for(int i=0; i<dim; i++){
    std::cout<<srt[i].first<<" ";
  } std::cout<<std::endl;

  std::cout<<std::endl<<"Gap: "<<std::abs(srt[0].first)-std::abs(srt[1].first)<<std::endl;

  Eigen::VectorXcd exact0=solve.eigenvectors().col(srt[0].second);

  //Power iteration:
 {
  std::ofstream ofs("test_Largest_EV2.power");
  time_point time_power1=now();
  Eigen::VectorXcd v0=Eigen::VectorXcd::Ones(dim);
  Eigen::VectorXcd v1(v0.size());

  for(int i=0; i<iter; i++){
    double h0=v0.norm();
    v0/=h0; 
    v1=A*v0;
    std::complex<double> nu0=v0.dot(v1);
    double diff0=(v1-nu0*v0).norm();
    v0.swap(v1); 

    time_point time_power2=now();
    ofs<<i<<" "<<diff0<<" "<<1.-std::norm(exact0.dot(v0)/v0.norm())<<" "<<time_diff(time_power2-time_power1)<<" "<<nu0<<std::endl;
  }
 }
 {
  std::ofstream ofs("test_Largest_EV2.opt");
  time_point time_power1=now();

//  std::complex<double> beta=srt[0].first*srt[0].first/4.;
  std::complex<double> beta=srt[1].first*srt[1].first/4.;
  Eigen::VectorXcd v0=Eigen::VectorXcd::Ones(dim);
  Eigen::VectorXcd vm1=Eigen::VectorXcd::Zero(v0.size());
  Eigen::VectorXcd v1(v0.size());

  for(int i=0; i<iter; i++){
    v1=A*v0;
    std::complex<double> nu0=v0.dot(v1);
    double diff0=(v1-nu0*v0).norm();
//    if(i>=iter/1){ v1-=beta*vm1*(i/(iter/1)); }
    v1-=beta*vm1;
    double h1=v1.norm();
    v1/=h1; v0/=h1;
    vm1.swap(v0); v0.swap(v1); 

    time_point time_power2=now();
    ofs<<i<<" "<<diff0<<" "<<1.-std::norm(exact0.dot(v0)/v0.norm())<<" "<<time_diff(time_power2-time_power1)<<" "<<nu0<<std::endl;
  }
 }
 {
  std::ofstream ofs("test_Largest_EV2.dyn");
  time_point time_power1=now();

//  std::complex<double> beta=srt[1].first*srt[1].first/4.;
  Eigen::VectorXcd v0=Eigen::VectorXcd::Ones(dim);
  Eigen::VectorXcd v1(v0.size());
  Eigen::VectorXcd v2(v0.size());
  Eigen::VectorXcd v3(v0.size());

  double h0=v0.norm(); v0/=h0;
  v1=A*v0;
  std::complex<double> nu0=v0.dot(v1);
  double diff0=(v1-nu0*v0).norm();

  double h1=v1.norm(); v1/=h1;
  v2=A*v1;
  std::complex<double> nu1=v1.dot(v2);
  double diff1=(v2-nu1*v1).norm();

  std::complex<double> r=1.;
  for(int i=2; i<iter; i++){
    std::complex<double> beta=nu1*nu1*r*r/4.;
//    std::complex<double> beta=srt[1].first*srt[1].first/4.;
    v2-=(beta/h1)*v0;
    double h2=v2.norm(); v2/=h2;
    v3=A*v2; 
    std::complex<double> nu2=v2.dot(v3);
    double diff2=(v3-nu2*v2).norm();
    double rho=1.; if(diff2<diff1)rho=diff2/diff1;
    r=2.*rho/(1.+rho*rho);

    v0.swap(v1); v1.swap(v2); v2.swap(v3); h1=h2; diff1=diff2; nu1=nu2;

    time_point time_power2=now();
    ofs<<i<<" "<<diff2<<" "<<1.-std::norm(exact0.dot(v2)/v2.norm())<<" "<<time_diff(time_power2-time_power1)<<" "<<nu1*nu1*r*r/4.<<" "<<srt[1].first<<" "<<nu1<<std::endl;

//    if(diff2<epsilon){
//      std::cout<<"dyn mom pow iter converged at iteration "<<i<<std::endl;
//      break;
//    }

  }
 }

/* 
  //Arnoldi:
  time_point time_Arnoldi1=now();
  v=Eigen::VectorXcd::Ones(dim);
  v.normalize();
  std::complex<double> Arnoldi=Largest_EV_Arnoldi(v, iter, epsilon, [&A](const Eigen::VectorXcd &x){return A*x;}, 2, reortho);
  Eigen::VectorXcd ArnoldiEV=v;
  time_point time_Arnoldi2=now();
*/

/*
 {
  time_point time_mpower1=now();
  std::ofstream ofs("test_Largest_EV2.opt");
  Eigen::VectorXcd v0=Eigen::VectorXcd::Ones(dim);

  Eigen::VectorXcd v1(v0.size());
  Eigen::VectorXcd v2(v0.size());
  Eigen::VectorXcd v3(v0.size());

  double h0=v0.norm();
  v0/=h0; 
      v1=A*v0;
      std::complex<double> nu0=v1.dot(v0);
      double diff0=(v1-nu0*v0).norm();
  
      double h1=v1.norm();
      v1/=h1; 
      v2=A*v1;
      std::complex<double> nu1=v2.dot(v1);
      double diff1=(v2-nu1*v1).norm();
 
      
      double r1=1.; if(diff1<diff0)r1=diff1/diff0;
      for(int it=12; it<iter2; it++){
        std::complex<double> beta=nu1*nu1*r1*r1/4.;
        v2-=(beta/h1)*v0;
        double h2=v2.norm();
        v2/=h2;
        v3=A*v2;
        std::complex<double> nu2=v3.dot(v2);
        double diff2=(v3-nu2*v2).norm();
        double rho1=1.; if(diff2<diff1)rho1=diff2/diff1;
        double r2=2.*rho1/(1.+rho1*rho1);

        v0.swap(v1); v1.swap(v2); v2.swap(v3);
        r1=r2; h0=h1; h1=h2; nu1=nu2; diff1=diff2;
        if(diff2<epsilon){
          break;
        }
      }
      Eigen::VectorXcd mpowerEV=v2;
      mpowerEV.normalize();
      std::complex<double> mpower=mpowerEV.dot(A*mpowerEV);
      time_point time_mpower2=now();
     
      ofs<<" "<<std::abs(1.-std::abs(mpowerEV.dot(exactEV)));
      ofs<<" "<<time_diff(time_mpower2-time_mpower1);
      ofs<<std::endl;
    }
  }
*/
  return 0;
}
