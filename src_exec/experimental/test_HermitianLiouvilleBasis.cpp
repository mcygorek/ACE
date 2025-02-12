#include "ACE.hpp"
#include "HermitianLiouvilleBasis.hpp"

using namespace ACE;

int main(int args, char **argv){
  Parameters param(args, argv, true);
  int N=param.get_as_int("N",2);

  Eigen::MatrixXcd M = HermitianLiouvilleBasis(N).get_Matrix();

  std::cout<<round_Matrix(M,1e-15)<<std::endl;
 
  std::cout<<std::endl;

  Eigen::MatrixXcd O(N*N,N*N);
  for(int i=0; i<N*N; i++){
    for(int j=0; j<N*N; j++){
      O(i,j)=M.col(i).dot(M.col(j));
    }
  }
  std::cout<<round_Matrix(O,1e-15)<<std::endl;

  return 0;
}
