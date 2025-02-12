#include "TrafoTree.hpp"
#include <memory>

namespace ACE{

void TrafoTree::combine_select(std::unique_ptr<TrafoTree> other, 
                               const SelectIndices &select){
  if(!other){
    std::cerr<<"TrafoTree::combine_select: unique_ptr 'other' not set!"<<std::endl;
    throw DummyException();
  }
  int dim1=get_dim();
  int dim2=other->get_dim();
  std::unique_ptr<TrafoTree> newfirst=std::unique_ptr<TrafoTree>(new TrafoTree());
  newfirst->first=std::move(first);
  newfirst->second=std::move(second);
  newfirst->T=T;
  newfirst->Tinv=Tinv;
  newfirst->S=S;
  
  S=select;
  initialize_T(S.size());
  first=std::move(newfirst);
  second=std::move(other);
  
}

void TrafoTree::combine(std::unique_ptr<TrafoTree> other){
  if(!other){
    std::cerr<<"TrafoTree::combine: unique_ptr 'other' not set!"<<std::endl;
    throw DummyException();
  }
  SelectIndices sel; sel.set_full(get_dim(), other->get_dim());
  combine_select(std::move(other), sel);
}
void TrafoTree::read(const std::string &filename){
  std::cerr<<"TrafoTree::read: NOT IMPLEMENTED YET!"<<std::endl;
  throw DummyException();
}

void TrafoTree::write(const std::string &filename)const{
  std::cerr<<"TrafoTree::write: NOT IMPLEMENTED YET!"<<std::endl;
  throw DummyException();
}

}//namespace
