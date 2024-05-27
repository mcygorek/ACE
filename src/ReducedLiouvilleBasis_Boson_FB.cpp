#include "ReducedLiouvilleBasis_Boson_FB.hpp"
#include <iostream>
#include <vector>
#include "otimes.hpp"
#include "LiouvilleTools.hpp"
#include "Operators_Boson.hpp"

/**  Reduce to the following basis:
     - start with initial state vector (initial mode density matrix)
     - construct new vectors by applying b^dagger and b at most 'nr_climb' times

HERE: Forward-backward version: 
     Add b/b^dagger once to forward and once backward direction but then also 
     include application of b/b^dagger at both, forward and backward direction.
     This matches more closely the effect of application of H_JC on the env.


*/
namespace ACE{

  void ReducedLiouvilleBasis_Boson_FB::setup_from_initial(const Eigen::MatrixXcd &initial, int nr_climb, double threshold){
    int M=initial.rows();
    if(M<2){
      std::cerr<<"ReducedLiouvilleBasis_Boson::setup_from_initial: initial.rows()<2!"<<std::endl;
      exit(1);
    }

    std::vector<Eigen::VectorXcd> v;
    v.push_back(H_Matrix_to_L_Vector(initial));
    v[0].normalize();

    std::vector<Eigen::MatrixXcd> OPS;

    OPS.push_back(otimes(
      Operators_Boson::adagger(M), Eigen::MatrixXcd::Identity(M,M))); 

    OPS.push_back(otimes(
      Eigen::MatrixXcd::Identity(M,M), Operators_Boson::adagger(M))); 

    OPS.push_back(otimes(
      Operators_Boson::a(M), Eigen::MatrixXcd::Identity(M,M))); 

    OPS.push_back(otimes(
      Eigen::MatrixXcd::Identity(M,M), Operators_Boson::a(M))); 

    OPS.push_back(otimes(
      Operators_Boson::adagger(M), Operators_Boson::adagger(M))); 

    OPS.push_back(otimes(
      Operators_Boson::adagger(M), Operators_Boson::a(M))); 

    OPS.push_back(otimes(
      Operators_Boson::a(M), Operators_Boson::adagger(M))); 

    OPS.push_back(otimes(
      Operators_Boson::a(M), Operators_Boson::a(M))); 

    int lastbase=0; 
    for(int n=0; n<nr_climb; n++){
      int added=0;
      int last_vsize=v.size();
      for(int i=lastbase; i<last_vsize; i++){
        for(size_t o=0; o<OPS.size(); o++){
//std::cout<<"TeSt: "<<n<<" "<<i<<" "<<o<<std::endl;
          Eigen::VectorXcd v2=OPS[o]*v[i];

          for(int ortho=0; ortho<2; ortho++){ //orthogonalize multiple times
            for(int j=0; j<(int)v.size(); j++){
              std::complex<double> proj=v[j].dot(v2);
              v2-=proj*v[j];
            }
          }
          double norm=v2.norm();
//std::cout<<"norm: "<<norm<<std::endl;
          if(norm>threshold){
//std::cout<<"accepted as element "<<v.size()<<std::endl;
            v.push_back(v2);
            v.back().normalize();
            added++;
          }
        }
      }
      if(added==0)break;
      lastbase=last_vsize;
    }
    
    U=Eigen::MatrixXcd(v[0].size(), v.size());
    for(size_t i=0; i<v.size(); i++){
      U.col(i)=v[i];
    }
    use_reduce=true;
  }

}//namespace
