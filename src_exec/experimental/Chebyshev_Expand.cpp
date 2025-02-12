#include "Parameters.hpp"
#include "ReadTable.hpp"
#include "Constants.hpp"
#include <Eigen/Eigen>
#include <iostream>
#include <fstream>

using namespace ACE;

/**
Expand a function: f(x)=\sum_n a_n T_n(x)
in terms of Chebyshev Polynomials.
Both, the input and the output is given on regular grids, where the output
grid is finer, effectively constituting an interpolation. 

This is, however, very different from 'Chebyshev interpolation', which uses
special, unevenly distributed points. Here, we assume we can't choose our
sample points.

Given the fact that the n-th order polynomials interpolating n+1 data points
is unique, our expansion has to use more than n+1 expansion coefficients to
be beyond conventional polynomial interpolation. We use NT coefficients instead,
arriving at an undetermined system, from which we then demand the norm 
of the expansion coefficient vector to be minimal. 

In effect, compared to conventional polynomial interpolation, we add 
higher-order terms, such that Runge's phenomenon is suppressed.


We write the equation: y_j  = \sum_n a_n T_n(x_j)
as:                    y = A a
and construct the underdetermined solution for a with minimal norm:
                   a_min = A^T (A A^T)^-1 y,
Which is given by QR decomposition [https://see.stanford.edu/materials/lsoeldsee263/08-min-norm.pdf]:
     A^T = Q*R  =>  a_min = Q*R^{-T} y


 !!! NOTE !!!: 
 Turns out not to work. Produces bed-of-nails approximation for NT>>N
 BUT: Can be partially mitigated by setting 'scale', which adds heigher weights       to higher T_n. Try: 1, 1.5, and 2

*/

int main(int args, char ** argv){
  Parameters param(args, argv);

  std::string outfile=param.get_as_string_check("outfile");
  std::string infile=param.get_as_string_check("infile");
  double scale=param.get_as_double("scale", 1.);
  int col=param.get_as_size_t("col",2);


//Read data points from file:
  ReadTable readtab(infile, 0, col-1);
  std::vector<std::vector<double> > &tab=readtab.table;

  int N=tab.size();
  if(N<3){
    std::cerr<<"Please provide at least 3 data point!"<<std::endl;
    exit(1);
  }

  int NT=param.get_as_size_t("NT",N);
  std::cout<<"N="<<N<<" NT="<<NT<<std::endl;
  if(NT<N){
    std::cerr<<"NT has to be at least as large as N!"<<std::endl;
    exit(1);
  }
 
  double lim_a=tab[0][0];
  double lim_b=tab.back()[0];
  for(int i=1; i<tab.size(); i++){
    if(tab[i][0]<tab[i-1][0]){
      std::cerr<<"data not ordered!"<<std::endl;
      exit(1);
    }
  }

  for(size_t i=0; i<tab.size(); i++){
    std::cout<<tab[i][0]<<" "<<tab[i][1]<<std::endl;
  }
  
//map input interval to [-1:1]
  for(size_t i=0; i<tab.size(); i++){
    tab[i][0]-=lim_a; tab[i][0]*=2./(lim_b-lim_a); tab[i][0]-=1.;
//    std::cout<<tab[i][0]<<" "<<tab[i][1]<<std::endl;
  }

//Construct matrix A^T:
  Eigen::MatrixXd AT=Eigen::MatrixXd::Zero(NT,N);
  for(int j=0; j<N; j++){  
    double T1=1;
    double T2=tab[j][0]/scale;
    AT(0, j)=T1;
    AT(1, j)=T2;
    for(int n=2; n<NT; n++){
      double T=(2.*tab[j][0]*T2-T1/scale)/scale;
      AT(n, j)=T;
      T1=T2; T2=T;
    }
  }
  
  std::cout<<"A^T:"<<std::endl<<AT<<std::endl;

//Calculate least-norm coefficient vector
  Eigen::VectorXd y(N);
  for(size_t i=0; i<N; i++){
    y(i)=tab[i][1];
  }
  Eigen::MatrixXd AAT = AT.transpose()*AT;
  
  Eigen::VectorXd coeff = AT * (AAT.householderQr().solve(y));
std::cout<<"coeff: "<<coeff.transpose()<<std::endl;

//Print sampled function:
  double dt=param.get_as_double("dt", (lim_b-lim_a)/(N-1));
  int N_sample=(lim_b-lim_a)/dt+0.5;

  std::ofstream ofs(outfile.c_str());
  for(int sample=0; sample<=N_sample; sample++){
    double t=lim_a+sample*dt;
    double x=(2.*sample)/N_sample-1.;

    double val=0.;
    double T1=1;
    double T2=x/scale;
    val+=coeff(0)*T1;
    val+=coeff(1)*T2;
    for(int n=2; n<NT; n++){
      double T=(2.*x*T2-T1/scale)/scale;
      val+=coeff(n)*T;
      T1=T2; T2=T;
    }
    std::cout<<t<<" "<<val<<std::endl;
    ofs<<t<<" "<<val<<std::endl;
  }
}
