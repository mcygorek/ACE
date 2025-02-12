#pragma once
#ifndef OUTPUT_OPS_DEFINED_H
#define OUTPUT_OPS_DEFINED_H

#include <vector>
#include <Eigen/Core>

namespace ACE{
class HilbertSpaceRotation;
class Parameters;

class Output_Ops{
public:
  std::vector<Eigen::MatrixXcd> ops;
  std::vector<int> proj;

  bool use_IP;
  Eigen::MatrixXcd H_IP;
  
  inline size_t size()const{return ops.size();}
  inline int get_dim()const{if(ops.size()>0)return ops[0].rows(); else return 0;}
  inline Eigen::MatrixXcd & operator[](size_t i){return ops[i];}
  inline const Eigen::MatrixXcd & operator[](size_t i)const {return ops[i];}
  inline std::vector<Eigen::MatrixXcd>::iterator begin(){return ops.begin();}
  inline std::vector<Eigen::MatrixXcd>::iterator insert(
                     std::vector<Eigen::MatrixXcd>::iterator position, 
                                          const Eigen::MatrixXcd& val){
    return ops.insert(position, val);
  }

  void add(const Eigen::MatrixXcd &op);
  
  void add_proj(int i);

  void rotate(const HilbertSpaceRotation &hs_rot);
 
  Eigen::MatrixXcd trafoIP(const Eigen::MatrixXcd &A, double t) const;

  void setup(Parameters &param, bool set_default_output=true);

  inline Output_Ops(const std::vector<Eigen::MatrixXcd> & list){
    ops=list;
  }
  inline Output_Ops(Parameters &param, bool set_default_output=true){
    setup(param, set_default_output);
  }
  inline Output_Ops(){
    use_IP=false;
  }
};

}//namespace
#endif
