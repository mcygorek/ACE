#include "Coupling_Groups.hpp"
#include <Eigen/Core>
#include <vector>
#include "Printable.hpp"

namespace ACE{

  void Coupling_Groups::update_Ngrps(){
    Ngrps=0;
    for(size_t i=0; i<grp.size(); i++)if(Ngrps<=grp[i])Ngrps=grp[i]+1;
  }
  void Coupling_Groups::set_default(int N){
    grp.resize(N);
    for(int i=0; i<N; i++)grp[i]=i;
    Ngrps=N;
  }
  std::ostream & Coupling_Groups::print(std::ostream &os)const{
    os<<Ngrps<<":";
    for(size_t i=0; i<grp.size(); i++){
      os<<" "<<grp[i];
    }
    return os;
  }
  void Coupling_Groups::autodetect_groups(const Eigen::MatrixXcd &couplings){
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

  void Coupling_Groups::setup_from_matrix(const Eigen::MatrixXcd &couplings_){
    autodetect_groups(couplings_);
    update_Ngrps();
  }


//Coupling_Groups_Liouville:
  void Coupling_Groups_Liouville::calculate(const Coupling_Groups &groups){
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

  void Coupling_Groups_Liouville::set_default(int N){
    calculate(Coupling_Groups(N));
  }

  std::ostream & Coupling_Groups_Liouville::print(std::ostream &os)const{
    os<<Ngrps<<":";
    for(size_t i=0; i<grp.size(); i++){
      os<<" "<<grp[i];
    }
    return os;
  }

}
