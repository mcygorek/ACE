#ifndef ANALYZE_PT_DEFINED_H_
#define ANALYZE_PT_DEFINED_H_

#include <iosfwd>
#include <complex>
//#include "MPS_Matrix.hpp"

namespace ACE{

class PPM_Canvas;
template<typename T> class MPS_Matrix_ScalarType;
typedef MPS_Matrix_ScalarType<std::complex<double> > MPS_Matrix;
class InfluenceFunctional_OD;


// Set of functions to extract (debug) information from a process tensor
class AnalyzePT{
public:

  static void canvas_add_MPS_Matrix(PPM_Canvas &canv, const MPS_Matrix &m, int a, int posx, int posy);
  
  static void print_single_pgm(const std::string &file, const InfluenceFunctional_OD &IF, int n, int a);

  static void print_pgm(const std::string &file, const InfluenceFunctional_OD &IF, int margins=10);

  static void print_SVD(const std::string &file, const InfluenceFunctional_OD &IF, const std::string &print_PT="");

  static void print_summary(std::ostream &ofs, const InfluenceFunctional_OD &IF);

  static void print_summary(const std::string &file, const InfluenceFunctional_OD &IF);
};

}
#endif
