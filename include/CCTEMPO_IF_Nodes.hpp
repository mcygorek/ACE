#pragma once
#ifndef ACE_CCTEMPO_IF_DEINFED_H
#define ACE_CCTEMPO_IF_DEINFED_H

#include "DiagBB.hpp"
#include "TimeGrid.hpp"
#include "PCH.hpp"
#include "ProcessTensorForwardList.hpp"

namespace ACE{

class CCTEMPO_IF_Nodes{
public:
  struct Node{
    std::vector<Eigen::MatrixXcd> b;
    std::shared_ptr<ProcessTensorForwardList> PT;
    DiagBB diagBB;
  };

  std::vector< Node > list;

  inline Node & operator[](int i){ return list[i]; }
  inline const Node & operator[](int i)const{ return list[i]; }
  inline size_t size()const{return list.size();}

  void setup(Parameters &param);
  CCTEMPO_IF_Nodes(Parameters &param){
    setup(param);
  }
  CCTEMPO_IF_Nodes(){}
  ~CCTEMPO_IF_Nodes(){}
};
}//namespace
#endif
