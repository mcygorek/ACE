#ifndef ACE_EDM_OUTPUT_PRINTER
#define ACE_EDM_OUTPUT_PRINTER

#include "Parameters.hpp"
#include <fstream>
#include <memory>
#include <complex>

namespace ACE{

class EDM_OutputPrinter{
public:
  std::unique_ptr<std::ofstream> ofs;

  inline operator bool()const{
    return (bool)ofs;
  }

  void print_time(double t);
  void print_value(std::complex<double> val);
  void print_endl();

  void set_file(const std::string &str);
  void setup(Parameters &param, const std::string &outfile_default="ACE.out", const std::string &outfile_key="outfile");

  EDM_OutputPrinter(const std::string &str=""){
    set_file(str);
  }
  EDM_OutputPrinter(Parameters &param, const std::string &outfile_default="ACE.out", const std::string &outfile_key="outfile"){
    setup(param, outfile_default, outfile_key);
  }
};

}//namespace

#endif
