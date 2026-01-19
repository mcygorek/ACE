#ifndef ACE_SELECT_INDICES_DEFINED_H
#define ACE_SELECT_INDICES_DEFINED_H
#include "BinaryReader.hpp"
#include <vector>

namespace ACE{

//used to project "otimes" of two matrices on combination of indices
struct SelectIndices{
  std::vector<std::pair<int, int> > list;

  inline size_t size()const{return list.size();}
  inline void resize(size_t sz){list.resize(sz);}
  inline void clear(){list.clear();}
  inline void push_back(const std::pair<int, int> &p){ list.push_back(p); }
  inline void push_back(int i1, int i2){ 
    list.push_back(std::pair<int,int>(i1,i2)); 
  }
  inline std::pair<int, int> & operator[](int i){ 
    return list[i];
  }
  inline const std::pair<int, int> & operator[](int i)const{ 
    return list[i];
  }
  inline operator const std::vector<std::pair<int, int> >&() const{
    return list;
  }
  inline void set_full(int dim1, int dim2){
    clear();
    for(int k1=0; k1<dim1; k1++){
      for(int k2=0; k2<dim2; k2++){
        push_back(k1,k2);
      }
    }
  }
  inline int min_dim_first()const{
    int dim=0; 
    for(int i=0; i<(int)list.size(); i++){
      if(list[i].first>=dim)dim=list[i].first+1;
    }
    return dim;
  }
  inline int min_dim_second()const{
    int dim=0; 
    for(int i=0; i<(int)list.size(); i++){
      if(list[i].second>=dim)dim=list[i].second+1;
    }
    return dim;
  }

  Eigen::MatrixXcd otimes(const Eigen::MatrixXcd & M1, 
                          const Eigen::MatrixXcd & M2, 
                          const SelectIndices & other) const;

  Eigen::VectorXcd Vector_otimes(const Eigen::VectorXcd & M1, 
                          const Eigen::VectorXcd & M2) const;


  void read(std::istream &ifs, const std::string &context="");
  void write(std::ostream &ofs)const;

  inline void read(const std::string &filename){
    std::unique_ptr<std::ifstream> ifs = open_file_check(filename);
    read(*(ifs.get()),filename);
  }
  inline void write(const std::string &filename)const{
    std::ofstream ofs(filename);
    write(ofs);
  }
  void print_info(std::ostream &ofs=std::cout)const;

  SelectIndices(){}
  SelectIndices(const std::string &filename){read(filename);}
  SelectIndices(std::istream &ifs, const std::string &context=""){
    read(ifs,context);
  }
  ~SelectIndices(){}
};


}//namespace
#endif
