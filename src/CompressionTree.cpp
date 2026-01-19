#include "CompressionTree.hpp"
#include <memory>

namespace ACE{
void CompressionTree::combine_select(std::shared_ptr<CompressionTree> & other, 
                               const SelectIndices &select){
  if(!other){
    std::cerr<<"CompressionTree::combine_select: shared_ptr 'other' not set!"<<std::endl;
    throw DummyException();
  }
//  int dim1=get_dim();
//  int dim2=other->get_dim();
  std::shared_ptr<CompressionTree> newfirst=std::make_shared<CompressionTree>();
  newfirst->first=first;
  newfirst->second=second;
  newfirst->T=T;
  newfirst->S=S;
  
  S=select;
  initialize_T(S.size());
  first=newfirst;
  second=other;
 
  print_info(); 
}
void CompressionTree::combine(std::shared_ptr<CompressionTree> & other){
  if(!other){
    std::cerr<<"CompressionTree::combine: shared_ptr 'other' not set!"<<std::endl;
    throw DummyException();
  }
  SelectIndices sel; sel.set_full(get_dim(), other->get_dim());
  combine_select(other, sel);
}
void CompressionTree::read(std::istream &ifs, const std::string &context){
  T=binary_read_EigenMatrixXcd(ifs, context);
  S=SelectIndices(ifs, context);
  int has_first = binary_read_int(ifs, context);
  if(has_first){
    first=std::make_shared<CompressionTree>(ifs, context);
  }
  int has_second = binary_read_int(ifs, context);
  if(has_second){
    second=std::make_shared<CompressionTree>(ifs, context);
  }
}

void CompressionTree::write(std::ostream &ofs)const{
  binary_write_EigenMatrixXcd(ofs, T);
  S.write(ofs);
  if(first){
    binary_write_int(ofs,1);
    first->write(ofs);
  }else{
    binary_write_int(ofs,0);
  }
  if(second){
    binary_write_int(ofs,1);
    second->write(ofs);
  }else{
    binary_write_int(ofs,0);
  }
}

void sweep_forward_pre(int n, int TTree_at, \
                   std::shared_ptr<CompressionTree> & TTree, \
                   std::shared_ptr<CompressionTree> & TTree_inv, \
                   const PassOn &pass_on){
  if(n==TTree_at+1 && TTree_inv){
    std::cout<<"sweep_forward: TTree_inv"<<std::endl;
    std::cout<<"P.rows()="<<pass_on.P.rows()<<" ";
    std::cout<<"P.cols()="<<pass_on.P.cols()<<" ";
    std::cout<<"T.rows()="<<TTree_inv->T.rows()<<std::endl;
    //update Tinv:
    TTree_inv->T = (pass_on.P)*TTree_inv->T; 
  }
}
void sweep_forward_post(int n, int TTree_at, \
                   std::shared_ptr<CompressionTree> & TTree, \
                   std::shared_ptr<CompressionTree> & TTree_inv, \
                   const PassOn &pass_on){
  if(n==TTree_at && TTree){
    std::cout<<"sweep_forward: TTree"<<std::endl;
    std::cout<<"P.rows()="<<pass_on.P.rows()<<" ";
    std::cout<<"P.cols()="<<pass_on.P.cols()<<" ";
    std::cout<<"T.rows()="<<TTree->T.rows()<<std::endl;
    //update Tinv:
    TTree->T = (pass_on.Pinv.transpose())*TTree->T; 
  }
}
void sweep_backward_pre(int n, int TTree_at, \
                   std::shared_ptr<CompressionTree> & TTree, \
                   std::shared_ptr<CompressionTree> & TTree_inv, \
                   const PassOn &pass_on){
  if(n==TTree_at && TTree){
    std::cout<<"sweep_backward: TTree"<<std::endl;
    std::cout<<"P.rows()="<<pass_on.P.rows()<<" ";
    std::cout<<"P.cols()="<<pass_on.P.cols()<<" ";
    std::cout<<"T.rows()="<<TTree->T.rows()<<std::endl;
    //update Tinv:
    TTree->T = (pass_on.P.transpose())*TTree->T; 
  }
}
void sweep_backward_post(int n, int TTree_at, \
                   std::shared_ptr<CompressionTree> & TTree, \
                   std::shared_ptr<CompressionTree> & TTree_inv, \
                   const PassOn &pass_on){
  if(n==TTree_at+1 && TTree_inv){
    std::cout<<"sweep_backward: TTree_inv"<<std::endl;
    std::cout<<"P.rows()="<<pass_on.P.rows()<<" ";
    std::cout<<"P.cols()="<<pass_on.P.cols()<<" ";
    std::cout<<"T.rows()="<<TTree_inv->T.rows()<<std::endl;
    //update Tinv:
    TTree_inv->T = (pass_on.Pinv)*TTree_inv->T; 
  }
}

}//namespace
