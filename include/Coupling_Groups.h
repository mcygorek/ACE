#ifndef COUPLING_GROUPS_DEFINED_H
#define COUPLING_GROUPS_DEFINED_H

#include <Eigen/Core>
#include <vector>
#include "Printable.h"

class Coupling_Groups: public Printable{
  public:

  int Ngrps;
  std::vector<int> grp;
  
  int operator[](int i)const{ return grp[i]; }
  int get_Ngrps()const{ return Ngrps; }
  int sys_dim()const{ return grp.size(); }

  void update_Ngrps(){
    Ngrps=0;
    for(size_t i=0; i<grp.size(); i++)if(Ngrps<=grp[i])Ngrps=grp[i]+1;
  }
  void set_default(int N){
    grp.resize(N);
    for(int i=0; i<N; i++)grp[i]=i;
    Ngrps=N;
  }
  virtual std::ostream &print(std::ostream &os=std::cout)const{
    os<<Ngrps<<":";
    for(size_t i=0; i<grp.size(); i++){
      os<<" "<<grp[i];
    }
    return os;
  }
  void autodetect_groups(const Eigen::MatrixXcd &couplings){
    grp.clear(); grp.resize(couplings.rows(),-1);
    Ngrps=0;
   
    for(int i=0; i<couplings.rows(); i++){
      int found=-1;
      for(int j=0; j<i; j++){
        if(abs(couplings(j,j)-couplings(i,i))<1e-8){
          found=j;
          break;
        }
      }
      if(found<0){
        grp[i]=Ngrps++;
      }else{
        grp[i]=grp[found];
      }
    }

#ifdef DEBUGF
    print();
#endif
  }
  void setup_from_matrix(const Eigen::MatrixXcd &couplings_){
    autodetect_groups(couplings_);
    update_Ngrps();
  }
  Coupling_Groups(const Eigen::MatrixXcd &couplings_){
    setup_from_matrix(couplings_);
  }
  Coupling_Groups(const std::vector<int> &gr) : grp(gr){
    update_Ngrps();
  }
  Coupling_Groups(int N=2){
    set_default(N);
  }
};



class Coupling_Groups_Liouville: public Printable{
public:
  int Ngrps;
  std::vector<int> grp;

  int operator[](int i)const{ return grp[i]; }
  int get_Ngrps()const{ return Ngrps; }
  int sys_dim()const{ return grp.size(); }


  void calculate(const Coupling_Groups &groups){
    size_t N=groups.grp.size();
    size_t NL=N*N;
    grp.resize(NL);
    for(size_t i=0; i<N; i++){
      for(size_t j=0; j<N; j++){
        grp[i*N+j]=groups[i]*groups.Ngrps+groups[j];
      }
    }
    Ngrps=groups.Ngrps*groups.Ngrps;
  }
  void set_default(int N){
    calculate(Coupling_Groups(N));
  }
  virtual std::ostream &print(std::ostream &os=std::cout)const{
    os<<Ngrps<<":";
    for(size_t i=0; i<grp.size(); i++){
      os<<" "<<grp[i];
    }
    return os;
  }
  Coupling_Groups_Liouville(const Coupling_Groups &groups){
    calculate(groups);
  }
  Coupling_Groups_Liouville(int N=2){
    set_default(N);
  }
};

#endif
