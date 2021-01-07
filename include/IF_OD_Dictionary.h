#ifndef IF_OD_DICTIONARY_DEFINED_H
#define IF_OD_DICTIONARY_DEFINED_H


/**  maps (alpha, tilde{alpha}) \to beta; store identical only once */
class IF_OD_Dictionary{
public:
  int N;  
  int reduced_dim;
  //Mapping. Note: if -1: contribution is zero
  std::vector<int> beta;

  int get_N()const{ return N; }
  int get_NL()const{ return get_N()*get_N(); }
  int get_NL2()const{ return get_NL()*get_NL(); }
  int get_reduced_dim()const{ return reduced_dim; }
  void calculate_reduced_dim(){  
    reduced_dim=0;
    for(size_t i=0; i<beta.size(); i++){
      if(beta[i]>=0 && beta[i]>=reduced_dim)reduced_dim=beta[i]+1;
    }
  }
  
  int operator()(int i)const{ return beta[i];}
  bool operator==(const IF_OD_Dictionary &other)const{ 
    if(beta.size()!=other.beta.size())return false;
    for(size_t i=0; i<beta.size(); i++)if(beta[i]!=other.beta[i])return false;
    return true;
  }

  void print_beta(std::ostream &os=std::cout)const{
    if(beta.size()<1)return;
    os<<beta[0];
    for(size_t i=1; i<beta.size(); i++){
      os<<" "<<beta[i];
    }
  }

  void set_default(int n){
    N=n;
    beta.resize(get_NL2());
    for(size_t aa=0; aa<beta.size(); aa++)beta[aa]=aa;
    reduced_dim=beta.size();
  }
  void set_default_diag(int n){
    N=n;
    beta.clear();
    beta.resize(get_NL2(), -1);
    for(size_t aa=0; aa<get_NL(); aa++){
      beta[aa*get_NL()+aa]=aa;
    }
    reduced_dim=get_NL();
  }




  std::vector<std::vector<int> > get_reverse_beta()const{
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
  bool join(const IF_OD_Dictionary &other){
//std::cout<<"dict join: first: "; print_beta(); std::cout<<" second: ";other.print_beta(); std::cout<<std::endl;
    if(*this==other)return false;

    if(beta.size()!=other.beta.size()){
      std::cerr<<"IF_OD_Dictionary::compare_expand: beta.size()!=other.beta.size()!"<<std::endl;
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
    for(size_t i=0; i<newbeta.size(); i++){
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
  
  void detect(const MPS_Matrix &a, double zero=1e-12){
//std::cout<<"DETECT: zero: "<<zero<<std::endl;
    int N=sqrt(sqrt(a.dim_i));
    if(N*N*N*N!=a.dim_i){
      std::cerr<<"IF_OD_Dictionary::detect: N*N*N*N!=a.dim_i!"<<std::endl;
      exit(1);
    }
    set_default(N);   
    if(zero<0.)return;
//std::cout<<"TEST1: "; print_beta(); std::cout<<std::endl;

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
        for(int aa3=aa+1; aa3<beta.size(); aa3++)beta[aa3]--;
      }
    }

//std::cout<<"TEST2: "; print_beta(); std::cout<<std::endl;
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
          for(int aa3=aa+1; aa3<beta.size(); aa3++){
            if(beta[aa3]>=0)beta[aa3]--;
          }
          break;
        }
      }
    }
//std::cout<<"TEST3: "; print_beta(); std::cout<<std::endl;
    calculate_reduced_dim();
  }
  void detect(const MPS &m, double zero=1e-12){
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
    
/*
    if(m.a.size()<1){
      set_default(2);
    }else if(m.a.size()<2){
      detect(m.a[0], zero);
    }else{
      detect(m.a[1], zero);
    }
*/
  }
  double get_keep_weight()const{
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

  void read_binary(std::istream &ifs){
    int size;
    ifs.read((char*)&size, sizeof(int));
    beta.resize(size);
    for(size_t i=0; i<beta.size(); i++)ifs.read((char*)&beta[i], sizeof(int));
    N=sqrt(sqrt(size));
    if(size!=N*N*N*N){
      std::cerr<<"dict::read_binary: size!=N*N*N*N!"<<std::endl;
      exit(1);
    }
    calculate_reduced_dim();
  }
  void write_binary(std::ostream &ofs)const{
    int size=beta.size();
    ofs.write((char*)&size, sizeof(int));
    for(size_t i=0; i<beta.size(); i++)ofs.write((char*)&beta[i], sizeof(int));
  }
  //Assume $\mathcal{H_E}\to \mathcal{H_X}\otimes \mathcal{H_E}$
  void expand_space_front(int Nfront){
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
  void expand_space_back(int Nback){
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
 
  void copy(const IF_OD_Dictionary &other){
    N=other.N;
    reduced_dim=other.reduced_dim;
    beta.resize(other.beta.size());
    for(size_t i=0; i<beta.size(); i++)beta[i]=other.beta[i];
  }
  IF_OD_Dictionary & operator=(const IF_OD_Dictionary &other){
    copy(other);
    return *this;
  }
  IF_OD_Dictionary (const IF_OD_Dictionary &other){
    copy(other);
  }
  IF_OD_Dictionary(const MPS_Matrix &a, double zero=1e-12){
    detect(a, zero);
  }
  IF_OD_Dictionary(const MPS &m, double zero=1e-12){
    detect(m, zero);
  }
  IF_OD_Dictionary(int def=2){
    set_default(def);
  }

};

#endif
