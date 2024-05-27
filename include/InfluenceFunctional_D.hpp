#ifndef INFLUENCE_FUNCTIONAL_D_DEFINED_H
#define INFLUENCE_FUNCTIONAL_D_DEFINED_H

#include <vector>
#include <complex>
#include "Eigen_fwd.hpp"

#include "TimeGrid.hpp"
#include "DiagBB.hpp"
#include "MPS.hpp"

namespace ACE{
template<typename T> class RankCompressor_ScalarType;
typedef RankCompressor_ScalarType<std::complex<double> > RankCompressor;
class Coupling_Groups_Liouville;

class InfluenceFunctional_D: public MPS{
public:
  TimeGrid tgrid;

  //closures:
  std::vector<Eigen::VectorXcd> c;
 
  //groups
  Coupling_Groups_Liouville lgroups; 

  void calculate_closures();

  void calculate_diagBB(DiagBB &diagBB, RankCompressor &compressor, int n_extra=0);

  virtual void read_binary(const std::string &filename);

  void set_none(int n_max, int N);

  inline InfluenceFunctional_D(const TimeGrid &tgr, DiagBB &diagBB, RankCompressor &compressor, int n_extra=0){
    tgrid=tgr;
    calculate_diagBB(diagBB, compressor, n_extra);
  }

  inline InfluenceFunctional_D(const std::string &filename){
    read_binary(filename);
  }

  inline InfluenceFunctional_D(){
  }
};

}//namespace
#endif
