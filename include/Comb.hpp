#ifndef ACE_COMB_DEFINED_H
#define ACE_COMB_DEFINED_H

#include "CombElement.hpp"
#include "BufferedContainer.hpp"
#include "Parameters.hpp"
#include "MPO.hpp"
#include "CCTEMPO_IF_Nodes.hpp"

namespace ACE{

class Comb: public BufferedContainer<CombElement>{
public:


//Operations on MPO
std::vector<Eigen::VectorXcd> get_cap_d2();
std::vector<Eigen::VectorXcd> get_cap_d4();
std::complex<double> evaluate(const std::vector<int> &list_d3);


void join_d3_d1_compress_d2_d4(const MPO & other, const TruncatedSVD &trunc, int verbosity=1);
void sweep_d4_d2(const TruncatedSVD &trunc, int verbosity=1);


void set_teeth(const CCTEMPO_IF_Nodes &nodes);
//void set_teeth(Parameters & param);
void set_initial_state(const MPO & mpo);

void setup(const MPO & initial, const CCTEMPO_IF_Nodes &nodes, int buffer_blocksize_, std::string buffer_filename_="");
//void setup(Parameters & param);

void print_dims(std::ostream &os=std::cout);

Comb(const MPO & initial, const CCTEMPO_IF_Nodes &nodes, int buffer_blocksize_, std::string buffer_filename_=""){
  magicString="COMB0";
  setup(initial, nodes, buffer_blocksize_, buffer_filename_);
}
Comb(){
  magicString="COMB0";
}
~Comb(){}
};
}//namespace
#endif
