#include "OutputExtractor.hpp"
#include "Parameters.hpp"
#include "LiouvilleTools.hpp"
#include "DummyException.hpp"

namespace ACE{

void OutputExtractor::setup(Parameters & param, int setdim){
  { std::vector<Eigen::VectorXcd> tmp; rho_t.swap(tmp); }
}

void OutputExtractor::print(int n, double t,const Eigen::VectorXcd &rho_reduced,
                       const std::vector<std::complex<double> > & env_reduced){ 

  rho_t.push_back(rho_reduced);

}

void OutputExtractor::print_eigenstate_occupations(double t, const Eigen::MatrixXcd &H, const Eigen::MatrixXcd &rho){

}

void OutputExtractor::finish(){
}


}//namespace
