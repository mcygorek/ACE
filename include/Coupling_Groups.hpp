#ifndef COUPLING_GROUPS_DEFINED_H
#define COUPLING_GROUPS_DEFINED_H

//#include <Eigen/Core>
#include "Eigen_fwd.hpp"
#include <vector>
#include "Printable.hpp"

namespace ACE{

class Coupling_Groups: public Printable{
  public:

  int Ngrps;
  std::vector<int> grp;
  
  inline int operator[](int i)const{ return grp[i]; }
  inline int get_Ngrps()const{ return Ngrps; }
  inline int sys_dim()const{ return grp.size(); }
  inline bool is_set_up()const{ return (grp.size()>0); }

  void update_Ngrps();
  
  void set_default(int N);
  
  virtual std::ostream &print(std::ostream &os=std::cout)const;
  
  void autodetect_groups(const Eigen::MatrixXcd &couplings);
  
  void setup_from_matrix(const Eigen::MatrixXcd &couplings_);
  
  inline Coupling_Groups(const Eigen::MatrixXcd &couplings_){
    setup_from_matrix(couplings_);
  }
  inline Coupling_Groups(const std::vector<int> &gr) : grp(gr){
    update_Ngrps();
  }
  inline Coupling_Groups(int N=2){
    set_default(N);
  }
};

class Coupling_Groups_Liouville: public Printable{
public:
  int Ngrps;
  std::vector<int> grp;

  inline int operator[](int i)const{ return grp[i]; }
  inline int get_Ngrps()const{ return Ngrps; }
  inline int sys_dim()const{ return grp.size(); }

  void calculate(const Coupling_Groups &groups);
  
  void set_default(int N);
  
  virtual std::ostream &print(std::ostream &os=std::cout)const;
  
  inline Coupling_Groups_Liouville(const Coupling_Groups &groups){
    calculate(groups);
  }
  inline Coupling_Groups_Liouville(int N=2){
    set_default(N);
  }
};

}
#endif
