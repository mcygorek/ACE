#include "IF_OD_Dictionary.hpp"
#include "Coupling_Groups.hpp"
#include "MPS.hpp"
#include "DummyException.hpp"

namespace ACE{

  void IF_OD_Dictionary::calculate_reduced_dim(){  
    reduced_dim=0;
    for(size_t i=0; i<beta.size(); i++){
      if(beta[i]>=0 && beta[i]>=reduced_dim)reduced_dim=beta[i]+1;
    }
  }
  bool IF_OD_Dictionary::is_diagonal()const{
    int NL=get_NL();
    for(int aa=0; aa<NL; aa++){
      for(int bb=0; bb<NL; bb++){
        if(aa!=bb && beta[bb*NL+aa]>=0)return false;
      }
    }
    return true;
  }
  
  bool IF_OD_Dictionary::operator==(const IF_OD_Dictionary &other)const{ 
    if(beta.size()!=other.beta.size())return false;
    for(size_t i=0; i<beta.size(); i++)if(beta[i]!=other.beta[i])return false;
    return true;
  }

  void IF_OD_Dictionary::print_beta(std::ostream &os)const{
    if(beta.size()<1)return;
    os<<beta[0];
    for(size_t i=1; i<beta.size(); i++){
      os<<" "<<beta[i];
    }
  }

  void IF_OD_Dictionary::set_default(int n){
    N=n;
    beta.resize(get_NL2());
    for(size_t aa=0; aa<beta.size(); aa++)beta[aa]=aa;
    reduced_dim=beta.size();
  }

  void IF_OD_Dictionary::set_default_diag(int n){
    N=n;
    beta.clear();
    beta.resize(get_NL2(), -1);
    for(int aa=0; aa<get_NL(); aa++){
      beta[aa*get_NL()+aa]=aa;
    }
    reduced_dim=get_NL();
  }

  void IF_OD_Dictionary::set_trivial(int n){
    N=n;
    beta.clear();
    beta.resize(get_NL2(), -1);
    for(int aa=0; aa<get_NL(); aa++){
      beta[aa*get_NL()+aa]=0;
    }
    reduced_dim=1;
  }

  std::vector<std::vector<int> > IF_OD_Dictionary::get_reverse_beta()const{
    //calculate_reduced_dim();
    std::vector<std::vector<int> > rev(reduced_dim);
    for(size_t i=0; i<beta.size(); i++){
      if(beta[i]>=0)rev[beta[i]].push_back(i);
    }
    return rev;
  }

  /** 
    For joining two IF_ODs with possibly different coupling structure:
    Find smallest common dictionary. 
    Returns true if it has changed
  */
  bool IF_OD_Dictionary::join(const IF_OD_Dictionary &other){
    if(is_diagonal() && *this==other)return false;

    if(beta.size()!=other.beta.size()){
      std::cerr<<"IF_OD_Dictionary::join: beta.size()!=other.beta.size() ("<<beta.size()<<" vs. "<<other.beta.size()<<")!"<<std::endl;
      exit(1);
    }

    std::vector<std::vector<int> > rev=get_reverse_beta();
    std::vector<std::vector<int> > rev2=other.get_reverse_beta();

    /*  Strategy: 
        start from new full dictionary, then reduce
    */

    std::vector<int> newbeta(beta.size());
    for(size_t i=0; i<newbeta.size(); i++)newbeta[i]=i;

    //remove simultaneous zeros
    for(size_t i=0; i<newbeta.size(); i++){
      if(beta[i]<0 && other.beta[i]<0){
        newbeta[i]=-1;
        for(size_t j=i+1; j<newbeta.size(); j++)newbeta[j]--;
      }
    }
//std::cout<<">";for(size_t i=0; i<newbeta.size(); i++)std::cout<<" "<<newbeta[i];std::cout<<std::endl;
    //remove simultaneous dubilcates
    for(int i=0; i<(int)newbeta.size(); i++){
      if(beta[i]<0){ 
        if(other.beta[i]<0)continue;
        for(size_t k=0; k<rev2[other.beta[i]].size(); k++){
          int i3=rev2[other.beta[i]][k];
          if(i==i3)continue;
          newbeta[i3]=newbeta[i];
        } 
      }else for(size_t j=0; j<rev[beta[i]].size(); j++){
        int i2=rev[beta[i]][j]; //point to same state as i
        if(i==i2)continue;

        if(other.beta[i]<0){
          if(other.beta[i2]<0){
            newbeta[i2]=newbeta[i];
          }
        }else{
          for(size_t k=0; k<rev2[other.beta[i]].size(); k++){
            int i3=rev2[other.beta[i]][k];
            if(i==i3)continue;
            if(i3==i2){
              newbeta[i2]=newbeta[i];
              break;
            }
          }
        }
      }
    }

//std::cout<<">";for(size_t i=0; i<newbeta.size(); i++)std::cout<<" "<<newbeta[i];std::cout<<std::endl; //check missing entries
    int maxentry=-1;
    for(size_t i=0; i<newbeta.size(); i++){
      if(newbeta[i]>maxentry)maxentry=newbeta[i];
    } 
    for(int m=0; m<maxentry; m++){
      bool found=false;
      for(size_t i=0; i<newbeta.size(); i++){
        if(newbeta[i]==m){ found=true; break; }
      }
      if(!found){
        for(size_t i=0; i<newbeta.size(); i++){
          if(newbeta[i]>m)newbeta[i]--;
        }
        maxentry--; m--;
      }
    }


    bool was_expanded=false;
    for(size_t i=0; i<beta.size(); i++){
      if(newbeta[i]!=beta[i])was_expanded=true;
    }
   
    beta.swap(newbeta);
    calculate_reduced_dim();
    return was_expanded;
  }
  
  template <typename T> 
  void IF_OD_Dictionary::detect(const MPS_Matrix_ScalarType<T> &a, double zero){
    int N=sqrt(sqrt(a.dim_i));
    if(N*N*N*N!=a.dim_i){
      std::cerr<<"IF_OD_Dictionary::detect: N*N*N*N!=a.dim_i (N="<<N<<" a.dim_i="<<a.dim_i<<")!"<<std::endl;
      throw DummyException();
    }
    set_default(N);   
    if(zero<0.)return;

    //identify zeros:
    for(size_t aa=0; aa<beta.size(); aa++){
      bool nonzero=false;
      for(int d1=0; d1<a.dim_d1; d1++){
        for(int d2=0; d2<a.dim_d2; d2++){
          if(abs(a(aa, d1, d2))>zero){
            nonzero=true;
            break;
          }
        }
        if(nonzero)break;
      }
      if(!nonzero){
        beta[aa]=-1;
        for(int aa3=aa+1; aa3<(int)beta.size(); aa3++)beta[aa3]--;
      }
    }

    for(size_t aa=1; aa<beta.size(); aa++){
      if(beta[aa]<0)continue;
      for(size_t aa2=0; aa2<aa; aa2++){
        if(beta[aa2]<0)continue;
        bool different=false;
        for(int d1=0; d1<a.dim_d1; d1++){
          for(int d2=0; d2<a.dim_d2; d2++){
            if(abs(a(aa, d1, d2)-a(aa2, d1, d2))>zero){
              different=true;
              break;
            }
          }
          if(different)break;
        }
        if(!different){
          beta[aa]=beta[aa2];
          for(int aa3=aa+1; aa3<(int)beta.size(); aa3++){
            if(beta[aa3]>=0)beta[aa3]--;
          }
          break;
        }
      }
    }
    calculate_reduced_dim();
  }
template void IF_OD_Dictionary::detect(const MPS_Matrix_ScalarType<double> &a, double zero);
template void IF_OD_Dictionary::detect(const MPS_Matrix_ScalarType<std::complex<double>> &a, double zero);

  template <typename T>
  void IF_OD_Dictionary::detect(const MPS_ScalarType<T> &m, double zero){
    if(m.a.size()<1){
      set_default(2);
      return;
    }
    detect(m.a[0], zero);
    for(size_t i=1; i<m.a.size(); i++){
      IF_OD_Dictionary dict2;
      dict2.detect(m.a[i],zero);
      join(dict2);
    }
  }
template void IF_OD_Dictionary::detect(const MPS_ScalarType<double> &m, double zero);
template void IF_OD_Dictionary::detect(const MPS_ScalarType<std::complex<double>> &m, double zero);


  double IF_OD_Dictionary::get_keep_weight()const{
    int NL=get_NL();
   if(true){
    return sqrt(NL);
   }else{
    std::vector<std::vector<int> > rev=get_reverse_beta();
    //how many diagonals are kept?
    int diags=0;
    for(size_t i=0; i<rev.size(); i++){
      for(size_t j=0; j<rev[i].size(); j++){
        if(rev[i][j]/NL==rev[i][j]%NL){
          diags++;
          break;
        }
      }
    }
    return sqrt(diags);
   }
  }

  void IF_OD_Dictionary::read_binary(std::istream &ifs){
    int size;
    ifs.read((char*)&size, sizeof(int));
//std::cout<<"TEST: IF_OD_Dictionary::read_binary: size: "<<size<<std::endl;
    beta.resize(size);
    for(size_t i=0; i<beta.size(); i++)ifs.read((char*)&beta[i], sizeof(int));
    N=sqrt(sqrt(size));
    if(size!=N*N*N*N){
      std::cerr<<"dict::read_binary: size!=N*N*N*N!"<<std::endl;
      exit(1);
    }
    calculate_reduced_dim();
  }

  void IF_OD_Dictionary::write_binary(std::ostream &ofs)const{
    int size=beta.size();
//std::cout<<"TEST: IF_OD_Dictionary::write_binary: size: "<<size<<std::endl;
    ofs.write((char*)&size, sizeof(int));
    for(size_t i=0; i<beta.size(); i++)ofs.write((char*)&beta[i], sizeof(int));
  }

  //Assume $\mathcal{H_E}\to \mathcal{H_X}\otimes \mathcal{H_E}$
  void IF_OD_Dictionary::expand_space_front(int Nfront){
    IF_OD_Dictionary IF_bck(*this);
    N*=Nfront;
    beta.clear();
    beta.resize(get_NL2(), -1);
    for(int nu=0; nu<Nfront; nu++){
      for(int mu=0; mu<Nfront; mu++){
        for(int nu1=0; nu1<IF_bck.get_N(); nu1++){
          for(int mu1=0; mu1<IF_bck.get_N(); mu1++){
            for(int nu2=0; nu2<IF_bck.get_N(); nu2++){
              for(int mu2=0; mu2<IF_bck.get_N(); mu2++){
                int alpha1=(nu*IF_bck.get_N()+nu1)*get_N()
                          +(mu*IF_bck.get_N()+mu1);
                int alpha2=(nu*IF_bck.get_N()+nu2)*get_N()
                          +(mu*IF_bck.get_N()+mu2);
                beta[alpha1*get_NL()+alpha2]=
IF_bck.beta[(nu1*IF_bck.get_N()+mu1)*IF_bck.get_NL()+(nu2*IF_bck.get_N()+mu2)];
              }
            }
          }
        }
      }
    }
    calculate_reduced_dim();
  }
  //Assume $\mathcal{H_E}\to \mathcal{H_E}\otimes\mathcal{H_X}$
  void IF_OD_Dictionary::expand_space_back(int Nback){
    IF_OD_Dictionary IF_bck(*this);
    N*=Nback;
    beta.clear();
    beta.resize(get_NL2(), -1);
    for(int nu=0; nu<Nback; nu++){
      for(int mu=0; mu<Nback; mu++){
        for(int nu1=0; nu1<IF_bck.get_N(); nu1++){
          for(int mu1=0; mu1<IF_bck.get_N(); mu1++){
            for(int nu2=0; nu2<IF_bck.get_N(); nu2++){
              for(int mu2=0; mu2<IF_bck.get_N(); mu2++){
                int alpha1=(nu+nu1*Nback)*get_N()
                          +(mu+mu1*Nback);
                int alpha2=(nu+nu2*Nback)*get_N()
                          +(mu+mu2*Nback);
                beta[alpha1*get_NL()+alpha2]=
IF_bck.beta[(nu1*IF_bck.get_N()+mu1)*IF_bck.get_NL()+(nu2*IF_bck.get_N()+mu2)];
              }
            }
          }
        }
      }
    }
    calculate_reduced_dim();
  }
  
  void IF_OD_Dictionary::expand(const ReadPT_struct &read_PT, bool verbose){
    if(read_PT.expand_front>1){
      if(verbose)std::cout<<"Expanding system dimension by "<<read_PT.expand_front<<" (front)"<<std::endl;
      expand_space_front(read_PT.expand_front);
    }
    if(read_PT.expand_back>1){
      if(verbose)std::cout<<"Expanding system dimension by "<<read_PT.expand_back<<" (back)"<<std::endl;
      expand_space_back(read_PT.expand_back);
    }
  }

  void IF_OD_Dictionary::expand_DiagBB(const DiagBB &diagBB){
    Coupling_Groups_Liouville lgroups(diagBB.groups);
    int NL=lgroups.sys_dim();
    beta.clear();
    N=sqrt(NL);
    beta.resize(NL*NL,-1);
    for(int i=0; i<NL; i++){
      beta[i*NL+i]=lgroups[i];
    }
    calculate_reduced_dim();

    //std::cerr<<"IF_OD_Dictionary::expand_DiagBB: Not implemented yet!"<<std::endl;
    //exit(1);
//  throw DummyException();
  }

  template <typename T>
  void IF_OD_Dictionary::reduce_MPS_Matrix(MPS_Matrix_ScalarType<T> &M)const{
    MPS_Matrix_ScalarType<T> m(get_reduced_dim(), M.dim_d1, M.dim_d2);
    std::vector<std::vector<int> > rev=get_reverse_beta();
    for(int i=0; i<m.dim_i; i++){
      if(rev[i].size()<1){ 
        std::cerr<<"IF_OD_Dictionary::reduce_MPS_Matrix: rev[i].size()<1!"<<std::endl;
        exit(1);
      }
      for(int d1=0; d1<m.dim_d1; d1++){
        for(int d2=0; d2<m.dim_d2; d2++){
          m(i, d1, d2)=M(rev[i][0], d1, d2);
        }
      }
    } 
    M.swap(m);
  }

template void IF_OD_Dictionary::reduce_MPS_Matrix(MPS_Matrix_ScalarType<double> &M) const;
template void IF_OD_Dictionary::reduce_MPS_Matrix(MPS_Matrix_ScalarType<std::complex<double>> &M) const;
 

  template <typename T>
  void IF_OD_Dictionary::expand_MPS_Matrix(MPS_Matrix_ScalarType<T> &M)const{
    if(M.dim_i!=get_reduced_dim()){
      std::cerr<<"IF_OD_Dictionary::expand_MPS_Matrix: M.dim_i!=get_reduced_dim()!"<<std::endl;
      exit(1);
    }
    MPS_Matrix_ScalarType<T> m(beta.size(), M.dim_d1, M.dim_d2);
    m.fill(0.);
    for(int i=0; i<m.dim_i; i++){
      if(beta[i]<0)continue;
      for(int d1=0; d1<m.dim_d1; d1++){
        for(int d2=0; d2<m.dim_d2; d2++){
          m(i, d1, d2)=M(beta[i], d1, d2);
        }
      }
    }
    M.swap(m);
  }
template void IF_OD_Dictionary::expand_MPS_Matrix(MPS_Matrix_ScalarType<double> &M) const;
template void IF_OD_Dictionary::expand_MPS_Matrix(MPS_Matrix_ScalarType<std::complex<double>> &M) const;
 

  template <typename T>
  void IF_OD_Dictionary::translate_MPS_Matrix(MPS_Matrix_ScalarType<T> &M,
                                         const IF_OD_Dictionary &other)const{

    if(get_N()!=other.get_N()){
      std::cerr<<"IF_OD_Dictionary::translate_MPS_Matrix: get_N()!=other.get_N()!"<<std::endl;
      exit(1);
    }
    if(M.dim_i!=other.get_reduced_dim()){
      std::cerr<<"IF_OD_Dictionary::translate_MPS_Matrix: M.dim_i!=get_reduced_dim()!"<<std::endl;
      exit(1);
    }

    std::vector<std::vector<int> > rev=get_reverse_beta();


    MPS_Matrix_ScalarType<T> m(get_reduced_dim(), M.dim_d1, M.dim_d2);
    m.fill(0.);
    for(int i=0; i<m.dim_i; i++){
      if(other.beta[rev[i][0]]<0)continue;
      for(int d1=0; d1<m.dim_d1; d1++){
        for(int d2=0; d2<m.dim_d2; d2++){
          m(i, d1, d2)=M(other.beta[rev[i][0]], d1, d2);
        }
      }
    }
    M.swap(m);
  }
template void IF_OD_Dictionary::translate_MPS_Matrix(MPS_Matrix_ScalarType<double> &M, const IF_OD_Dictionary &other) const;
template void IF_OD_Dictionary::translate_MPS_Matrix(MPS_Matrix_ScalarType<std::complex<double>> &M, const IF_OD_Dictionary &other) const;
 
  template <typename T>
  void IF_OD_Dictionary::expand_all_diagonal(MPS_Matrix_ScalarType<T> &M){
    if(M.dim_i!=get_reduced_dim()){
      std::cerr<<"IF_OD_Dictionary::expand_add_diagonal: M.dim_i!=get_reduced_dim()!"<<std::endl;
      exit(1);
    }
    
    if(!is_diagonal()){
      std::cerr<<"IF_OD_Dictionary::expand_all_diagonal: is_diagonal() is false!"<<std::endl;
      exit(1);
    }
    IF_OD_Dictionary other; other.set_default_diag(N);
    if(*this==other)return;

    int NL=get_NL();

    MPS_Matrix_ScalarType<T> m(NL, M.dim_d1, M.dim_d2);
    m.fill(0.);
    for(int i=0; i<NL; i++){
      if(beta[i*NL+i]<0)continue;
      for(int d1=0; d1<m.dim_d1; d1++){
        for(int d2=0; d2<m.dim_d2; d2++){
          m(i, d1, d2)=M(beta[i*NL+i], d1, d2);
        }
      }
    }
    M.swap(m);
    *this=other;
  }

template void IF_OD_Dictionary::expand_all_diagonal(MPS_Matrix_ScalarType<double> &M);
template void IF_OD_Dictionary::expand_all_diagonal(MPS_Matrix_ScalarType<std::complex<double>> &M);
 
void IF_OD_Dictionary::copy(const IF_OD_Dictionary &other){
  N=other.N;
  reduced_dim=other.reduced_dim;
  beta=other.beta;
//  beta.resize(other.beta.size());
//  for(size_t i=0; i<beta.size(); i++)beta[i]=other.beta[i];
}

}//namespace
