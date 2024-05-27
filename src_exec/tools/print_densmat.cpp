#include "BinaryReader.hpp"
#include "Parameters.hpp"
#include <iomanip>

using namespace ACE;

int main(int args, char **argv){
  Parameters param(args, argv, true);
  std::string read_densmat=param.get_as_string_check("read_densmat");
  int prec=param.get_as_size_t("precision",4);

  Eigen::MatrixXcd M=binary_read_EigenMatrixXcd(read_densmat);

  std::cout<<std::setprecision(4);
  std::cout<<M<<std::endl;

  return 0;
}
