#include "ProcessTensor_real.hpp"
#include "InfluenceFunctional_OD.hpp"
#include "HermitianLiouvilleBasis.hpp"
#include "Closure_Ops.hpp"
#include "Env_State_Filter.hpp"
#include "Trafo_Chain.hpp"
#include "LiouvilleTools.hpp"
#include "MPS_Matrix.hpp"

namespace ACE{

  MPS_Matrix ProcessTensor_real::get_a_phys(int n)const{  
    if(n<0||n>=a.size()){
      std::cerr<<"ProcessTensor_real::get_a_phys: n<0||n>=a.size()!"<<std::endl;
      exit(1);
    }
    int N=dict.get_N();
    int NL=N*N;
    if(HLU.rows()!=NL){
      std::cerr<<"ProcessTensor_real::get_a_phys: HLU.rows()!=NL ("<<HLU.rows()<<"/"<<NL<<")!"<<std::endl;
      exit(1);
    }
    MPS_Matrix M(NL*NL, a[n].dim_d1, a[n].dim_d2);
    M.set_zero();
    for(int i=0; i<NL; i++){
      for(int j=0; j<NL; j++){
        for(int l=0; l<NL; l++){
          if(dict.beta[i*NL+l]<0)continue;
          for(int d1=0; d1<a[n].dim_d1; d1++){
            for(int d2=0; d2<a[n].dim_d2; d2++){
              M(i*NL+j, d1, d2) += a[n](dict.beta[i*NL+l], d1, d2)*HLU.adjoint()(l,j);
            } 
          }
        }
      } 
    }
    MPS_Matrix M2(NL*NL, a[n].dim_d1, a[n].dim_d2);
    M2.set_zero();
    for(int i=0; i<NL; i++){
      for(int j=0; j<NL; j++){
        for(int l=0; l<NL; l++){
          for(int d1=0; d1<a[n].dim_d1; d1++){
            for(int d2=0; d2<a[n].dim_d2; d2++){
              M2(i*NL+j, d1, d2) += HLU(i,l) * M(l*NL+j, d1, d2);
            } 
          }
        }
      } 
    }
    return M2;
  }

  void ProcessTensor_real::trafo_d1(const Eigen::MatrixXd &T, int n){
    if(n<0 || n>=a.size()){
      std::cerr<<"ProcessTensor_real::trafo_d1: n<0 || n>=a.size()!"<<std::endl;
      exit(1);
    }
    if(a[n].dim_d1!=T.cols()){
      std::cerr<<"ProcessTensor_real::trafo_d1: a[n].dim_d1!=T.cols()!"<<std::endl;
      exit(1);
    }
    MPS_Matrix_real M(a[n].dim_i, T.rows(), a[n].dim_d2);
    M.set_zero();
    for(int i=0; i<a[n].dim_i; i++){
      for(int d1=0; d1<T.rows(); d1++){
        for(int d2=0; d2<a[n].dim_d2; d2++){
          for(int d=0; d<a[n].dim_d1; d++){
            M(i, d1, d2) += T(d1, d) * a[n](i, d, d2);
          }
        }
      }
    }
    a[n].swap(M);   
  }

  void ProcessTensor_real::trafo_d2(const Eigen::MatrixXd &T, const Eigen::MatrixXd &Tback, int n){
    if(n<0 || n>=a.size()){
      std::cerr<<"ProcessTensor_real::trafo_d2: n<0 || n>=a.size()!"<<std::endl;
      exit(1);
    }
    if(a[n].dim_d2!=T.rows()){
      std::cerr<<"ProcessTensor_real::trafo_d2: a[n].dim_d2!=T.rows()!"<<std::endl;
      exit(1);
    }
    MPS_Matrix_real M(a[n].dim_i, a[n].dim_d1, T.cols());
    M.set_zero();
    for(int i=0; i<a[n].dim_i; i++){
      for(int d1=0; d1<a[n].dim_d1; d1++){
        for(int d2=0; d2<T.cols(); d2++){
          for(int d=0; d<a[n].dim_d2; d++){
            M(i, d1, d2) += a[n](i, d1, d) * T(d, d2);
          }
        }
      }
    }
    a[n].swap(M);   

    if(n>=c.size()){
      std::cerr<<"ProcessTensor_real::trafo_d2: n>=c.size()!"<<std::endl;
      exit(1);
    }
    if(c[n].rows()!=Tback.cols()){
      std::cerr<<"ProcessTensor_real::trafo_d2: c[n].rows()!=Tback.cols()!"<<std::endl;
      exit(1);
    }
    Eigen::VectorXcd c_tmp=Eigen::VectorXcd::Zero(Tback.rows());
    for(int d1=0; d1<Tback.rows(); d1++){ 
      for(int d=0; d<Tback.cols(); d++){ 
        c_tmp(d1)+=Tback(d1,d)*c[n](d);
      }
    }
    c[n]=c_tmp;
   
    if(n<env_ops.size()){
      std::vector<Eigen::VectorXcd> env_ops_tmp(env_ops[n].size(),
                                        Eigen::VectorXcd::Zero(Tback.rows()));
      for(size_t o=0; o<env_ops[n].size(); o++){
        for(int d1=0; d1<Tback.rows(); d1++){ 
          for(int d=0; d<Tback.cols(); d++){ 
             env_ops_tmp[o](d1)+=Tback(d1,d)*env_ops[n][o](d);
          }
        }
      }
      env_ops[n]=env_ops_tmp;
    }
  }
 
  std::vector<std::vector<Eigen::MatrixXd> > ProcessTensor_real::propA_to_real(
        const std::vector<std::vector<Eigen::MatrixXcd> > & propA, 
        const Eigen::MatrixXcd &HLBasisN, const Eigen::MatrixXcd &HLBasisM,
        double *res){
 
    int NL=HLBasisN.rows();  
    int ML=HLBasisM.rows();  

    std::vector<std::vector<Eigen::MatrixXcd> > tmp(NL, 
      std::vector<Eigen::MatrixXcd>(NL, Eigen::MatrixXcd::Zero(ML,ML)));
    for(int a1=0; a1<NL; a1++){
      for(int a2=0; a2<NL; a2++){
        tmp[a2][a1]=
          HLBasisM.adjoint()*(propA[a2][a1])*HLBasisM;
      }
    }

    std::vector<std::vector<Eigen::MatrixXcd> > tmp2(NL,
      std::vector<Eigen::MatrixXcd>(NL, Eigen::MatrixXcd::Zero(ML,ML)));
    for(int a1=0; a1<NL; a1++){
      for(int a2=0; a2<NL; a2++){
        for(int a=0; a<NL; a++){
          tmp2[a2][a1]+=HLBasisN.adjoint()(a2,a)*tmp[a][a1];
        }
      }
    }

    double im=0.;
    std::vector<std::vector<Eigen::MatrixXd> > B(NL,
      std::vector<Eigen::MatrixXd>(NL, Eigen::MatrixXd::Zero(ML,ML)));
    for(int a1=0; a1<NL; a1++){
      for(int a2=0; a2<NL; a2++){
        Eigen::MatrixXcd c=Eigen::MatrixXcd::Zero(ML,ML);
        for(int a=0; a<NL; a++){
          c+=tmp2[a2][a]*HLBasisN(a,a1);
        }
        B[a2][a1]+=c.real();
        im+=(c.imag().squaredNorm());
      }
    }
    if(res!=NULL)*res=sqrt(im);
    return B;
  }


  void ProcessTensor_real::calculate_closures(){
    c.clear();
    if(a.size()<1){
      return;
    }

    int N=dict.get_N();

    c.resize(a.size());
    c.back()=Eigen::VectorXcd::Zero(a.back().dim_d2);
    c.back()(0)=1;

    for(int n=(int)c.size()-2; n>=0; n--){
      MPS_Matrix aphys=get_a_phys(n+1);

      c[n]=Eigen::VectorXcd::Zero(a[n].dim_d2);
      for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
          int i_ind=dict(((i*N+i)*N+j)*N+j);
          if(i_ind<0)continue;
          for(int d1=0; d1<a[n+1].dim_d1; d1++){
            for(int d2=0; d2<a[n+1].dim_d2; d2++){
              c[n](d1)+=c[n+1](d2)*aphys(i_ind, d1, d2)/((double)N);
            }
          }
        }
      }
    }
  }

  void ProcessTensor_real::calculate_dict(double zero, bool verbose){
    dict.detect(*this, zero);
    if(verbose){
      std::cout<<"InfluenceFunctional_OD: Dictionary(";
      std::cout<<dict.get_reduced_dim()<<"): ";
      dict.print_beta(); 
      std::cout<<std::endl;
    }
  }

  //keep only non-redundent terms with respect to outer indices
  void ProcessTensor_real::reduce_to_dict(){
    if(a.size()<1)return;
    if(a[0].dim_i!=(int)dict.beta.size()){
      std::cerr<<"InfluenceFunctional_OD::reduce_to_dict: a[0].dim_i!=dict.beta.size()!"<<std::endl;
      exit(1);
    }
    for(size_t n=0; n<a.size(); n++){
      dict.reduce_MPS_Matrix(a[n]);
    }
//    calculate_closures();
  }

  //build full MPS matrices from reduced (dictionary)
  void ProcessTensor_real::expand_from_dict(){
    if(a.size()<1)return;
    if(a[0].dim_i!=dict.get_reduced_dim()){
      std::cerr<<"InfluenceFunctional_OD::expand_from_dict: a[0].dim_i!=dict.get_reduced_dim()!"<<std::endl;
      exit(1);
    }
    for(size_t n=0; n<a.size(); n++){
      MPS_Matrix_real m(dict.beta.size(), a[n].dim_d1, a[n].dim_d2);
      m.fill(0.);
      for(int i=0; i<m.dim_i; i++){
        if(dict.beta[i]<0)continue;
        for(int d1=0; d1<m.dim_d1; d1++){
          for(int d2=0; d2<m.dim_d2; d2++){
            m(i, d1, d2)=a[n](dict.beta[i], d1, d2);
          }
        }
      } 
      m.swap(a[n]);
    }
    calculate_closures();
  }

  void ProcessTensor_real::make_longer(int n_template, int how_many){
    if(n_template<0||n_template>=a.size()){
      std::cerr<<"ProcessTensor_real::make_longer: n_template<0||n_template>=a.size()!"<<std::endl;
      exit(1);
    }
    if(a[n_template].dim_d1!=a[n_template].dim_d2){
      std::cerr<<"ProcessTensor_real::make_longer: a[n_template].dim_d1!=a[n_template].dim_d2!"<<std::endl;
      exit(1);
    }
   
    size_t oldsize=a.size();
    MPS_Matrix_real a_n=a[n_template];
    a.insert(a.begin()+n_template, (size_t)how_many, a_n);
//    a.insert((std::vector<MPS_Matrix_real>::iterator)(a.begin()+n_template), (size_t)how_many, a_n);
    
    if(c.size()!=oldsize){
      std::cerr<<"ProcessTensor_real::make_longer: c.size()!=a.size() ("<<c.size()<<" vs. "<<a.size()<<")!"<<std::endl;
      exit(1);
    }
    Eigen::VectorXcd c_n=c[n_template];
    c.insert(c.begin()+n_template, how_many, c_n);
  
    if(env_ops.size()==oldsize){
      std::vector<Eigen::VectorXcd>  env_ops_n=env_ops[n_template];
      env_ops.insert(env_ops.begin()+n_template, how_many, env_ops_n);
    }
  }

  void ProcessTensor_real::calculate_dt0_ndt0(ModePropagator &mprop, int n_max, double ta, double dt, double dt0, int ndt0){

//std::cout<<"MYTEST"<<std::endl;
    if(closure_ops)closure_ops.print_info();

    if(mprop.rBasis->use()){
      std::cerr<<"Using reduced environment Liouville basis not implemented for real-valued PT!"<<std::endl;
      exit(1);
    }

    a.clear();
    if(n_max<1)return;
    a.resize(n_max);
    int N=mprop.get_N_system();
    int M=mprop.get_N_mode();
    int NL=N*N;
    int ML=M*M;

    HLU = HermitianLiouvilleBasis(N).get_Matrix();
    Eigen::MatrixXcd HLU2= HermitianLiouvilleBasis(M).get_Matrix();

    const Eigen::VectorXcd c_bath_init=
         (HLU2.adjoint()*H_Matrix_to_L_Vector(mprop.get_bath_init()));
    const Eigen::VectorXd bath_init=c_bath_init.real();
//std::cout<<"c_bath_init: "<<c_bath_init.transpose()<<std::endl;


    //Environment operators: 
    std::vector<Eigen::VectorXcd> env_ops_template(1+mprop.env_ops.size(), 
                                                   Eigen::VectorXcd::Zero(ML));
    for(int m=0; m<M; m++){
      env_ops_template[0](m*M+m)=1;
    }
    env_ops_template[0]=(env_ops_template[0].transpose()*HLU2).transpose();

    for(size_t i=0; i<mprop.env_ops.size(); i++){
      if(mprop.env_ops[i].cols()!=M || mprop.env_ops[i].rows()!=M){
        std::cerr<<"mprop.env_ops[i].cols()!=M || mprop.env_ops[i].rows()!=M"<<std::endl;
        exit(1);
      }
      for(int m1=0; m1<M; m1++){
        for(int m2=0; m2<M; m2++){
          env_ops_template[i+1](m1*M+m2)=mprop.env_ops[i](m2,m1);
        }
      } 
      env_ops_template[i+1]=(env_ops_template[i+1].transpose()*HLU2).transpose();
    }

    env_ops.resize(n_max-1); 
    for(int n=0; n<n_max-1; n++){ 
      env_ops[n].resize(env_ops_template.size());
      for(size_t i=0; i<env_ops_template.size(); i++){
        env_ops[n][i]=env_ops_template[i];
      }
    }


    std::vector<Eigen::VectorXcd> endv(env_ops[0].size(), Eigen::VectorXcd::Zero(1));
    env_ops.push_back(endv);


    //actual MPR:
    double max_res=0.;
    for(int n=0; n<n_max; n++){
      a[n].resize(NL*NL, ML, ML); 
      if(n<ndt0){
        mprop.update(ta+n*dt0, dt0); 
      }else{
        mprop.update(ta+ndt0*dt0+(n-ndt0)*dt, dt);
      }
      double res=0;
      std::vector<std::vector<Eigen::MatrixXd> > B=propA_to_real(mprop.A, HLU, HLU2, &res);
      if(res>max_res)max_res=res;
      for(int i=0; i<NL; i++){
        for(int j=0; j<NL; j++){
          for(int d1=0; d1<ML; d1++){
            for(int d2=0; d2<ML; d2++){
              a[n](i*NL+j, d1, d2)=B[i][j](d2,d1);
            }
          }
        }
      }
    }
    std::cout<<"Imaginary remainder: "<<max_res<<std::endl;
    

    MPS_Matrix_real ainit(NL*NL, 1, ML);
    ainit.fill(0.);
    for(int i=0; i<NL; i++){
      for(int j=0; j<NL; j++){
        for(int d1=0; d1<ML; d1++){
          for(int d2=0; d2<ML; d2++){
            ainit(i*NL+j, 0, d2) += a[0](i*NL+j, d1, d2)*bath_init(d1);
          }
        }
      }
    }
    a[0].swap(ainit);


    if(!closure_ops){
      MPS_Matrix_real &a_ref=a.back();
      MPS_Matrix_real afinal(NL*NL, a_ref.dim_d1, 1);
      afinal.fill(0.);
      for(int i=0; i<NL; i++){
        for(int j=0; j<NL; j++){
          for(int d1=0; d1<a_ref.dim_d1; d1++){
            for(int xi0=0; xi0<M; xi0++){

#ifndef HERMITIAN_BASIS_USE_TRACE
              afinal(i*NL+j, d1, 0) += a_ref(i*NL+j, d1, xi0); 
#else
              afinal(i*NL+j, d1, 0) += a_ref(i*NL+j, d1, 0)*sqrt(N)/M; 
#endif
            }
          }
        }
      }
      a_ref.swap(afinal);
    }else{
      MPS_Matrix_real &a_ref=a.back();
      std::vector<Eigen::VectorXcd> closures;

      //only boolean "use_env_closures" set 
      if(closure_ops.use_env_ops && closure_ops.env_ops_nr.size()<1){
        closures=env_ops_template;
      }else{
        closures.push_back(env_ops_template[0]);
        if(!closure_ops.use_env_ops){
          for(size_t i=0; i<closure_ops.ops.size(); i++){
           closures.push_back((closure_ops.ops[i].transpose()*HLU2).transpose());
          }
        } 
        for(size_t i=0; i<closure_ops.env_ops_nr.size(); i++){
          if( (closure_ops.env_ops_nr[i].first>=env_ops_template.size())
              || (closure_ops.env_ops_nr[i].first<0)){
            std::cerr<<"ProcessTensor_real: closure_ops.env_ops_nr[i]>=env_ops_template.size()!"<<std::endl;
            exit(1);
          }
          closures.push_back(closure_ops.env_ops_nr[i].second
                        * env_ops_template[closure_ops.env_ops_nr[i].first]);
        }
        if(closure_ops.H_scale>0){
          Eigen::MatrixXcd Htot=mprop.get_Htot(ta+n_max*dt);
          for(int i=0; i<N; i++){
            for(int j=0; j<N; j++){
              Eigen::VectorXcd v(ML);
              for(int m1=0; m1<M; m1++){
                for(int m2=0; m2<M; m2++){
                  v(m1*M+m2)=Htot(i*M+m1,j*M+m2);
                }
              }
              closures.push_back(closure_ops.H_scale*(v.transpose()*HLU2).transpose());
            }
          }
        }
      }
std::cout<<"closures.size()="<<closures.size()<<std::endl;
         
      MPS_Matrix_real afinal(NL*NL, a_ref.dim_d1, closures.size());
      afinal.fill(0.);
      for(int i=0; i<NL; i++){
        for(int j=0; j<NL; j++){
          for(int d1=0; d1<a_ref.dim_d1; d1++){
            for(int o=0; o<closures.size(); o++){
              for(int d2=0; d2<ML; d2++){
                afinal(i*NL+j, d1, o) += closures[o](d2).real() 
                                              * a_ref(i*NL+j, d1, d2); 
              }
            }
          }
        }
      }
      a_ref.swap(afinal);
    }
    
    dict.set_default(N);
    calculate_closures();
  }

  void ProcessTensor_real::calculate(ModePropagator &mprop, int n_max, 
                    double ta, double dt, const Closure_Ops &closure_ops_){
    closure_ops=closure_ops_;
    calculate_dt0_ndt0(mprop, n_max, ta, dt, dt, 0);
  }

  void ProcessTensor_real::calculate(ModePropagator &mprop, 
                  const TimeGrid &tgrid, const Closure_Ops &closure_ops_){
    closure_ops=closure_ops_;
    calculate_dt0_ndt0(mprop, tgrid.n_calc, tgrid.ta, tgrid.dt, tgrid.dt0, tgrid.ndt0);
  }

  //second order combination of two MPS matrices
  void ProcessTensor_real::single_join_halfdt(
      MPS_Matrix_real &M, const MPS_Matrix_real &M1, const MPS_Matrix_real &M2, 
      const IF_OD_Dictionary &thisdict, const IF_OD_Dictionary &odict, 
      const IF_OD_Dictionary &ndict){

      int NL=thisdict.get_NL();
      std::vector<std::vector<int> > rev=thisdict.get_reverse_beta();
      std::vector<std::vector<int> > rev2=odict.get_reverse_beta();
      std::vector<std::vector<int> > newrev=ndict.get_reverse_beta();


      MPS_Matrix_real btmp(ndict.get_reduced_dim(), M.dim_d1*M1.dim_d1, M.dim_d2*M1.dim_d2);
      btmp.fill(0.);

      MPS_Matrix_real b(ndict.get_reduced_dim(), M.dim_d1*M1.dim_d1, M.dim_d2*M2.dim_d2);
      b.fill(0.);


      //first: btmp <- e^{this} e^{other}
      for(int i=0; i<NL; i++){
        for(int j=0; j<NL; j++){
          int i_ind=thisdict.beta[i*NL+j];         if(i_ind<0)continue;

          for(int k=0; k<NL; k++){
            int i_ind_new=ndict.beta[i*NL+k];  if(i_ind_new<0)continue;
            if(newrev[i_ind_new][0]!=i*NL+k)continue;

            int i_ind2=odict.beta[j*NL+k];  if(i_ind2<0)continue;
  
            for(int d1=0; d1<M.dim_d1; d1++){
              for(int d2=0; d2<M.dim_d2; d2++){
                for(int od1=0; od1<M1.dim_d1; od1++){
                  for(int od2=0; od2<M1.dim_d2; od2++){
  btmp(i_ind_new, d1*M1.dim_d1 + od1, d2*M1.dim_d2 + od2)+=
                    M(i_ind,d1,d2) * M1(i_ind2,od1,od2);
                  }
                }
              }
            }
          }
        }
      }

      //second: b <- e^{other} btmp 
      for(int i=0; i<NL; i++){
        for(int j=0; j<NL; j++){
          int i_ind=odict.beta[i*NL+j];   if(i_ind<0)continue;

          for(int k=0; k<NL; k++){
            int i_ind_new=ndict.beta[i*NL+k];  if(i_ind_new<0)continue;
            if(newrev[i_ind_new][0]!=i*NL+k)continue;

            int i_ind2=ndict.beta[j*NL+k];     if(i_ind2<0)continue;
 
  
            for(int d1=0; d1<M.dim_d1; d1++){
              for(int d2=0; d2<M.dim_d2; d2++){
                for(int od0=0; od0<M1.dim_d1; od0++){
                  for(int od1=0; od1<M2.dim_d1; od1++){
                    for(int od2=0; od2<M2.dim_d2; od2++){
  b(i_ind_new, d1*M1.dim_d1 + od0, d2*M2.dim_d2 + od2)
    += M2(i_ind,od1,od2) * btmp(i_ind2, d1*M1.dim_d1 + od0, d2*M1.dim_d2 + od1);
                    }
                  }
                }
              }
            }
          }
        }
      }
      M.swap(b);
  }

  void ProcessTensor_real::join_and_sweep_halfdt(
                             ProcessTensor_real &other, 
                             RankCompressor_real &compressor,
                             bool do_sweep, double keep_weight,
                             bool do_join){
  #ifdef NOSWEEP 
    do_sweep=false;
  #endif

    int aback_dim_d2=a.back().dim_d2;
    if(closure_ops){
      std::cout<<"closure_ops=true aback_dim_d2="<<aback_dim_d2<<std::endl;     
      if(!other.closure_ops){
        std::cerr<<"ProcessTensor_real::join_and_sweep_halfdt: other.use_env_closure==false!"<<std::endl;
        exit(1);
      }
    }else{
      //std::cout<<"use_env_closure=false"<<std::endl;     
      if(other.closure_ops){
        std::cerr<<"ProcessTensor_real::join_and_sweep_halfdt: other.use_env_closure==true!"<<std::endl;
        exit(1);
      }
    }
     
    if(a.size()<1)return;
//    int NL=sqrt(dict.beta.size());
    IF_OD_Dictionary newdict(dict);


    if(do_join){
      if(2*a.size()!=other.a.size()){
        std::cerr<<"ProcessTensor_real::join: 2*a.size()!=other.size()!"<<std::endl;
        exit(1);
      }
      if(dict.beta.size() != other.dict.beta.size()){
        std::cerr<<"ProcessTensor_real::join: dict.beta.size() != other.dict.beta.size()!"<<std::endl;
        exit(1);
      }
      newdict.join(other.dict);
    }


    for(size_t n=0; n<a.size(); n++){
     
      if(do_join){

        if(env_state_filter.use_count()!=0){
          env_state_filter->preprocess(a, env_ops, n,
                                       other.a, other.env_ops, 2*n+1);
        }

        single_join_halfdt(a[n], other.a[2*n], other.a[2*n+1], dict, other.dict, newdict);

       //environment operators:
       if(n<env_ops.size() && 2*n+1 < other.env_ops.size()){

        if(env_ops[n].size()<=1){
          Eigen::VectorXcd thisid=env_ops[n][0];
          env_ops[n].clear();
          env_ops[n].resize(other.env_ops[2*n+1].size());
          for(size_t i=0; i<other.env_ops[2*n+1].size(); i++){
            env_ops[n][i]=Vector_otimes(thisid, other.env_ops[2*n+1][i]);
          }
        }else if(env_ops[n].size()==other.env_ops[2*n+1].size()){
          Eigen::VectorXcd thisid=env_ops[n][0];
          Eigen::VectorXcd otherid=other.env_ops[2*n+1][0];
          env_ops[n][0]=Vector_otimes(thisid, otherid);
          for(size_t i=1; i<env_ops[n].size(); i++){
            env_ops[n][i]=Vector_otimes(env_ops[n][i], otherid)
                         +Vector_otimes(thisid, other.env_ops[2*n+1][i]);
          }
        }else if(env_ops[n].size()>other.env_ops[2*n+1].size() && 
                 other.env_ops[2*n+1].size()>0){
          Eigen::VectorXcd thisid=env_ops[n][0];
          Eigen::VectorXcd otherid=other.env_ops[2*n+1][0];
          env_ops[n][0]=Vector_otimes(thisid, otherid);
          for(size_t i=1; i<env_ops[n].size(); i++){
            env_ops[n][i]=Vector_otimes(env_ops[n][i], otherid);
            if(i<other.env_ops[2*n+1].size()){
              env_ops[n][i]+=Vector_otimes(thisid, other.env_ops[2*n+1][i]);
            }
          }
        }else{
          std::cerr<<"Problem joining environment operators: env_ops[n].size()!=other.env_ops[n].size()!"<<std::endl;
          exit(1);
        }
       }

        if(env_state_filter.use_count()!=0){
          env_state_filter->filter(a, env_ops,n);
        }

      }
      
      if(do_sweep && n>0){
        if(print_timesteps){
          if(n==1){std::cout<<"Sweep forward: "<<n<<std::flush;}
          else{
            std::stringstream ss_last; ss_last<<n-1;
            for(int i=ss_last.str().length(); i>0; i--)std::cout<<'\b';
            std::cout<<n<<std::flush;
          }
          if(n==a.size()-1)std::cout<<std::endl;
        }

        compressor.sweep_block_low_to_high(n-1, *this, keep_weight, this);
      }
    }

    //backward sweep
    if(do_sweep){
      for(int n=(int)a.size()-1; n>1; n--){
        if(print_timesteps){
          if(n==(int)a.size()-1){std::cout<<"Sweep backward: "<<n<<std::flush;}
          else{
            std::stringstream ss_last; ss_last<<n+1;
            for(int i=ss_last.str().length(); i>0; i--)std::cout<<"\b \b";
            std::cout<<n<<std::flush;
          }
          if(n==2)std::cout<<std::endl;
        }

        compressor.sweep_block_high_to_low(n, *this, keep_weight, this);
      }
    }

    if(do_join){
      dict=newdict;
    }

    calculate_closures();

    if(closure_ops && a.size()>0){
      int nmax=a.size();
      if(aback_dim_d2<2){
      }else{
        if(aback_dim_d2!=other.a[2*nmax-1].dim_d2){
          std::cerr<<"ProcessTensor_real: sweep: use_env_closure: a[nmax-1].dim_d2!=other.a[2*nmax-1].dim_d2 ("<<aback_dim_d2<<" vs. "<<other.a[2*nmax-1].dim_d2<<")!"<<std::endl;
          exit(1);
        }
        MPS_Matrix_real M(a[nmax-1].dim_i, a[nmax-1].dim_d1, aback_dim_d2);
        M.set_zero();
        for(int i=0; i<a[nmax-1].dim_i; i++){
          for(int d1=0; d1<a[nmax-1].dim_d1; d1++){
            M(i, d1, 0)=a[nmax-1](i, d1, 0); 
            for(int d2=1; d2<aback_dim_d2; d2++){
              M(i, d1, d2)+=a[nmax-1](i, d1, 0*aback_dim_d2+d2);
              M(i, d1, d2)+=a[nmax-1](i, d1, d2*aback_dim_d2+0);
            }
          }
        }
        a[nmax-1].swap(M);
      }
    }
  }

  void ProcessTensor_real::sweep_low_high_low(RankCompressor_real &compressor, double keep_weight){
    for(size_t n=1; n<a.size(); n++){
      compressor.sweep_block_low_to_high(n-1, *this, keep_weight, this);
    }
    for(int n=(int)a.size()-1; n>1; n--){
      compressor.sweep_block_high_to_low(n, *this, keep_weight, this);
    }
    calculate_closures();
  }

  //Define tracking of env_ops with sweeps:
  void ProcessTensor_real::process_low_to_high(int n, const Eigen::MatrixXd &R){
//    if(tgrid.use_rep)cta.set_R_if_correct_n(R,n);
    if(n<env_ops.size())for(size_t i=0; i<env_ops[n].size(); i++){
      env_ops[n][i]=R*env_ops[n][i];
    }

    for(size_t i=0; i<trafo_chains.size(); i++){
      if(trafo_chains[i].second==n){
        trafo_chains[i].first.add_low_to_high(R);
      }
    }
  }

  void ProcessTensor_real::process_high_to_low(int n, const Eigen::MatrixXd &L){
//    if(tgrid.use_rep)cta.set_L_if_correct_n(L,n);

    //What we need is the inverse of L.
    //Note: L also contains (parts of) the singular values, which 
    //can be obtained from the norms of the columns.
    //L is orthogonal, but not orthonormal. 
    //To get the inverse, we must transpose+conjugate L 
    //and devide by the norms of the columns twice.
    Eigen::VectorXd vL(L.cols());
    for(int c=0; c<L.cols(); c++)vL(c)=L.col(c).norm();
    Eigen::MatrixXd Linv=L.adjoint();
    for(int c=0; c<Linv.rows(); c++)Linv.row(c)/=(vL(c)*vL(c));

    if(n-1<env_ops.size())for(size_t i=0; i<env_ops[n-1].size(); i++){
      env_ops[n-1][i]=Linv*env_ops[n-1][i];
    }

    for(size_t i=0; i<trafo_chains.size(); i++){
      if(trafo_chains[i].second==n-1){
        trafo_chains[i].first.add_high_to_low(L);
//        trafo_chains[i].first.T.back() = trafo_chains[i].first.T.back() * Linv.adjoint();
      }
    }
  }


  void ProcessTensor_real::sweep(RankCompressor_real &compressor, bool do_sweep, double keep_weight){
    ProcessTensor_real other;
    join_and_sweep_halfdt(other, compressor, do_sweep, keep_weight, false);  
  }
  
  InfluenceFunctional_OD ProcessTensor_real::get_IF_OD()const{
    InfluenceFunctional_OD IF;

    int N=dict.get_N(); 
//    int NL=N*N;
//    HLU = HermitianLiouvilleBasis(N).get_Matrix();
    
    IF.a.clear();
    IF.a.resize(a.size());
    for(int n=0; n<(int)a.size(); n++){
      IF.a[n]=get_a_phys(n);
    }
    IF.env_ops=env_ops;
    
    IF.dict.set_default(N);

    //IF.calculate_closures();
    IF.c.clear();
    IF.c.resize(a.size());
    for(int n=0; n<(int)a.size(); n++){
      IF.c[n]=c[n];
    }

std::cout<<"IF.env_ops.size(): "<<IF.env_ops.size()<<std::endl;
    return IF;
  }

//read - write 
  void ProcessTensor_real::read_env_ops_binary(std::istream &ifs){
    int itmp;
    std::complex<double> c;

    env_ops.clear();
    ifs.read((char*)&itmp, sizeof(int));
    env_ops.resize(itmp);  
   
    for(size_t i=0; i<env_ops.size(); i++){
      ifs.read((char*)&itmp, sizeof(int));
      env_ops[i].clear();
      env_ops[i].resize(itmp); 
      for(size_t j=0; j<env_ops[i].size(); j++){
        ifs.read((char*)&itmp, sizeof(int));
        env_ops[i][j]=Eigen::VectorXcd(itmp); 
        for(int k=0; k<env_ops[i][j].size(); k++){
          ifs.read((char*)&c, sizeof(std::complex<double>));
          env_ops[i][j](k)=c;
        }
      }
    }
  }
  bool ProcessTensor_real::is_proper_file_format(std::ifstream &ifs,bool complain){
    char buf[5]; buf[4]='\0';
    ifs.read(buf,4);
    int version;
    ifs.read((char*)&version, sizeof(int));
    
    if(std::string(buf)=="PTOD"){
//      std::cout<<"File format version: "<<version<<std::endl;
      if(version!=4){
        if(complain){
          std::cerr<<"ProcessTensor_real: Expected file format version 4!"<<std::endl;
          exit(1);
        }else{
          return false;
        }
      }
    }else{
      if(complain){
        std::cerr<<"Error reading process tensor (OD): Wrong format!"<<std::endl;
        exit(1);
      }else{
        return false;
      }
    }
    return true;
  }

  bool ProcessTensor_real::is_proper_file_format(const std::string &filename, bool complain){
    std::ifstream ifs(filename.c_str(), std::ios::binary);
    if(!ifs.is_open()){
      if(complain){
        std::cout<<"Error reading ProcessTensor_real: Cannot open file '"<<filename<<"'!"<<std::endl;
        exit(1);
      }else{
        return false;
      }
    }
    return is_proper_file_format(ifs, complain);
  }

  void ProcessTensor_real::read_binary(const std::string &filename){
    std::cout<<"Reading process tensor (OD) '"<<filename<<"'"<<std::endl;
    std::ifstream ifs(filename.c_str(), std::ios::binary);
    if(!ifs.is_open()){
      std::cout<<"Error reading ProcessTensor_real: Cannot open file '"<<filename<<"'!"<<std::endl;  
      exit(1);
    }

    is_proper_file_format(ifs,true);

    dict.read_binary(ifs);
dict.print_beta();std::cout<<std::endl;

    read_env_ops_binary(ifs);

    MPS_real::read_binary(ifs);

    calculate_closures();
//    tgrid.n_tot=a.size();
//    tgrid.n_calc=tgrid.n_tot;
  
    HLU = HermitianLiouvilleBasis(dict.get_N()).get_Matrix();
  }

  void ProcessTensor_real::write_env_ops_binary(std::ostream &ofs)const{
    int itmp=env_ops.size(); ofs.write((char*)&itmp, sizeof(int));
   
    for(size_t i=0; i<env_ops.size(); i++){
      itmp=env_ops[i].size(); ofs.write((char*)&itmp, sizeof(int));
      for(size_t j=0; j<env_ops[i].size(); j++){
        itmp=env_ops[i][j].size(); ofs.write((char*)&itmp, sizeof(int));
        for(int k=0; k<env_ops[i][j].size(); k++){
          std::complex<double> c=env_ops[i][j](k);
          ofs.write((char*)&c, sizeof(std::complex<double>));
        }
      }
    }
  }

  void ProcessTensor_real::write_binary(const std::string &filename)const{
std::cout<<"Writing process tensor to binary '"<<filename<<"'"<<std::endl;
dict.print_beta(); std::cout<<std::endl;
    std::ofstream ofs(filename.c_str(), std::ios::binary);
    char buf[5]="PTOD";
    int version=4;
    ofs.write(buf, 4);
    ofs.write((char*)&version, sizeof(int));
    dict.write_binary(ofs);
    write_env_ops_binary(ofs);
    MPS_real::write_binary(ofs);
  }

  void ProcessTensor_real::set_from_complex(const InfluenceFunctional_OD &IF){
    
    print_timesteps=IF.print_timesteps;
    closure_ops.set_trivial();
    env_ops.clear();
    env_state_filter.reset();
    trafo_chains.clear();

    int N=IF.dict.get_N();
    int NL=N*N;
    HLU = HermitianLiouvilleBasis(N).get_Matrix();
    dict.set_default(N);

    a.resize(IF.a.size());
    Eigen::VectorXcd phases(1); phases(0)=1.;

    for(int n=0; n<(int)a.size(); n++){
if(n>=2)exit(1);
      a[n].resize(NL*NL, IF.a[n].dim_d1, IF.a[n].dim_d2);
      a[n].set_zero();

std::cout<<"TEST1"<<std::endl;
      MPS_Matrix M(NL*NL, IF.a[n].dim_d1, IF.a[n].dim_d2);
      M.set_zero();
      for(int i=0; i<NL; i++){
        for(int j=0; j<NL; j++){
          for(int l=0; l<NL; l++){
            if(dict.beta[i*NL+l]<0)continue;
            for(int d1=0; d1<a[n].dim_d1; d1++){
              for(int d2=0; d2<a[n].dim_d2; d2++){
                M(i*NL+j, d1, d2) += IF.a[n](IF.dict.beta[i*NL+l], d1, d2)*HLU(l,j);
              } 
            }
          }
        } 
      }

std::cout<<"TEST2"<<std::endl;
      MPS_Matrix M2(NL*NL, a[n].dim_d1, a[n].dim_d2);
      M2.set_zero();
      for(int i=0; i<NL; i++){
        for(int j=0; j<NL; j++){
          for(int l=0; l<NL; l++){
            for(int d1=0; d1<a[n].dim_d1; d1++){
              for(int d2=0; d2<a[n].dim_d2; d2++){
                M2(i*NL+j, d1, d2) += HLU.adjoint()(i,l) * M(l*NL+j, d1, d2);
              }
            } 
          }
        }
      }
//std::cout<<"M2:"<<std::endl;for(int d2=0; d2<M2.dim_d2; d2++){
//for(int ij=0; ij<M2.dim_i; ij++){for(int d1=0; d1<M2.dim_d1; d1++){std::cout<<M(ij,d1,d2)<<" ";}std::cout<<std::endl;}std::cout<<std::endl;}


      for(int ij=0; ij<M2.dim_i; ij++){
        for(int d1=0; d1<M2.dim_d1; d1++){
          for(int d2=0; d2<M2.dim_d2; d2++){
            M2(ij, d1, d2)*=phases(d1);
          }
        }
      }
      
std::cout<<"TEST3"<<std::endl;
      Eigen::VectorXcd ph2=Eigen::VectorXcd::Zero(M2.dim_d2);
      double abs2_imag=0;
      for(int d2=0; d2<M2.dim_d2; d2++){
        int ij_max=0, d1_max=0; std::complex<double> val_max=0.;
        for(int ij=0; ij<M2.dim_i; ij++){
          for(int d1=0; d1<M2.dim_d1; d1++){
            if(abs(M2(ij, d1, d2))>abs(val_max)){
              ij_max=ij; d1_max=d1; val_max=M2(ij, d1, d2);
            }
          }
        }
        if(abs(val_max)>1e-32){
          ph2(d2)=val_max/abs(val_max);
        }else{
          ph2(d2)=1.;
        }
      }
std::cout<<"TEST4"<<std::endl;
      for(int d2=0; d2<M2.dim_d2; d2++){
        for(int ij=0; ij<M2.dim_i; ij++){
          for(int d1=0; d1<M2.dim_d1; d1++){
            M2(ij, d1, d2)/=ph2(d2);
            a[n](ij, d1, d2)=M2(ij, d1, d2).real();
            abs2_imag+=M2(ij, d1, d2).imag()*M2(ij, d1, d2).imag();
          }
        }
      }
//std::cout<<"M2:"<<std::endl;for(int d2=0; d2<M2.dim_d2; d2++){
//for(int ij=0; ij<M2.dim_i; ij++){for(int d1=0; d1<M2.dim_d1; d1++){std::cout<<M(ij,d1,d2)<<" ";}std::cout<<std::endl;}std::cout<<std::endl;}

std::cout<<"TEST5"<<std::endl;
      phases=ph2;
      std::cout<<"n="<<n<<" abs2_imag="<<abs2_imag<<std::endl;  
    }//n

    calculate_closures();
  }

//initializers

  void ProcessTensor_real::initialize(int N){
    closure_ops.use=false;
    print_timesteps=false;
    dict.set_default(N);
    HLU = HermitianLiouvilleBasis(N).get_Matrix();
    trafo_chains.clear();
  }

  void ProcessTensor_real::set_trivial(int N, int nmax, double dict_zero, bool verbose){ 
    closure_ops.use=false;
    HLU = HermitianLiouvilleBasis(N).get_Matrix();
    dict.set_default(N);
    a.clear();
    if(nmax>0){
      int NL=dict.get_NL();
      MPS_Matrix_real Mref(NL*NL, 1, 1);
      Mref.set_zero();
      for(int i=0; i<NL; i++){
        Mref(i*NL+i, 0, 0)=1.;
      }
      a.resize(nmax, Mref);
    }

    Eigen::VectorXcd env_id(1); env_id<<1;
    env_ops.clear();
    env_ops.resize(nmax, std::vector<Eigen::VectorXcd>(1, env_id));

    calculate_dict(dict_zero, verbose);
    reduce_to_dict();

    c.clear();
    c.resize(nmax, env_id);
  } 
  
}//namespace
