#ifndef ACE_OUTPUT_EXTRACTOR_DEFINED_H
#define ACE_OUTPUT_EXTRACTOR_DEFINED_H

#include "OutputPrinter.hpp"

namespace ACE{

class OutputExtractor: public OutputPrinter{
public:  
  
//  std::vector<Eigen::VectorXcd> rho_t;


  virtual void setup(Parameters & param, int setdim=-1);

  virtual void print(int n, double t, const Eigen::VectorXcd & rho_reduced,
                 const std::vector<std::complex<double> > & env_reduced);

  virtual void print_eigenstate_occupations(double t, const Eigen::MatrixXcd &H, const Eigen::MatrixXcd &prop);

  virtual void finish();

  OutputExtractor(Parameters & param, int setdim=-1){
    setup(param, setdim);
  }
  OutputExtractor(){
  }
};
}//namespace
#endif
