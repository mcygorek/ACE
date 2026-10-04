#ifndef ACE_COMB_ELEMENT_DEFINED_H
#define ACE_COMB_ELEMENT_DEFINED_H

#include "BufferedElement.hpp"
#include "FourLeg.hpp"
#include "ProcessTensorBuffer.hpp"
#include "ProcessTensorForwardList.hpp"

namespace ACE{
class CombElement: public BufferedElement{
public:
  FourLeg f;
  Eigen::VectorXcd closure_d1;
  std::unique_ptr<ProcessTensorBuffer> tooth;
//  std::vector<Eigen::MatrixXcd> b;
//  std::shared_ptr<ProcessTensorForwardList> PT;


  void sweep_d2_d4(const TruncatedSVD &trunc, Eigen::MatrixXcd &P, bool is_last, int verbosity=1);
  void sweep_d4_d2(const TruncatedSVD &trunc, Eigen::MatrixXcd &P, bool is_last, int verbosity=1);


  virtual void read_binary(std::istream &is);
  virtual void write_binary(std::ostream &os)const;

  inline void copy(const CombElement & other){
    f=other.f;
    closure_d1=other.closure_d1;
    if(other.tooth){
      tooth=std::make_unique<ProcessTensorBuffer>();
      tooth->copy_content(*other.tooth);
    }else{
      tooth.reset(nullptr);
    }
//    b=other.b;
//    PT=other.PT;
  }
  CombElement & operator=(const CombElement & other){
    copy( other ); 
    return (*this);
  }
  CombElement(const CombElement & other){ copy( other ); }
  CombElement(){};
  virtual ~CombElement();
};

}//namespace
#endif
