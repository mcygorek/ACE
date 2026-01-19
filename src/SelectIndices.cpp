#include "SelectIndices.hpp"
#include <Eigen/Dense>

namespace ACE{

Eigen::MatrixXcd SelectIndices::otimes(const Eigen::MatrixXcd & M1, 
                          const Eigen::MatrixXcd & M2, 
                          const SelectIndices & other) const{

  Eigen::MatrixXcd res(size(), other.size());
  for(int i=0; i<(int)list.size(); i++){
    for(int j=0; j<(int)other.list.size(); j++){
      res(i,j) = M1(list[i].first, other.list[j].first) 
               * M2(list[i].second, other.list[j].second);
    }
  }
  return res;
}

Eigen::VectorXcd SelectIndices::Vector_otimes(const Eigen::VectorXcd & M1, 
                          const Eigen::VectorXcd & M2) const{

  Eigen::VectorXcd res(size());
  for(int i=0; i<(int)list.size(); i++){
    res(i) = M1(list[i].first) * M2(list[i].second);
  }
  return res;
}

void SelectIndices::read(std::istream &ifs, const std::string &context){
  int sz=binary_read_int(ifs, context);
  {std::vector<std::pair<int, int> > tmp(sz);
  list.swap(tmp);}
  for(int i=0; i<sz; i++){
    list[i].first=binary_read_int(ifs, context);
    list[i].second=binary_read_int(ifs, context);
  }
}
void SelectIndices::write(std::ostream &ofs)const{
  int sz=list.size();
  binary_write_int(ofs, sz);
  for(int i=0; i<sz; i++){
    binary_write_int(ofs, list[i].first);
    binary_write_int(ofs, list[i].second);
  }
}
void SelectIndices::print_info(std::ostream &ofs)const{
  ofs<<"size()="<<size()<<" min_dim_first()="<<min_dim_first()<<" min_dim_second()="<<min_dim_second();
}

}//namespace 
