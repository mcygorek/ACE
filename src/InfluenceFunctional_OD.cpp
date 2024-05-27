#include "InfluenceFunctional_OD.hpp"
#include "InfluenceFunctional_D.hpp"
#include "IF_OD_Abstract.hpp"
#include "ProcessTensor_real.hpp"
#include "SingleBathMode.hpp"
#include "Rep_GIF.hpp"
#include "LiouvilleTools.hpp"
#include "Compress_Trafo_At.hpp"
#include "Sweep_Trafo_Processor.hpp"
#include "ModePropagatorGenerator.hpp"
#include "RankCompressorList.hpp"

namespace ACE{

  const MPS_Matrix & InfluenceFunctional_OD::get_a(int n)const{
    if(n<0||n>=(int)a.size()){
      std::cerr<<"InfluenceFunctional_OD: get_a out of bounds "<<n<<"/"<<a.size()<<std::endl;
      exit(1);
    }
    return a[n];
  }

  const Eigen::VectorXcd & InfluenceFunctional_OD::get_c(int n)const{
    if(n<0||n>=(int)c.size()){
      std::cerr<<"InfluenceFunctional_OD: get_c out of bounds "<<n<<"/"<<c.size()<<std::endl;
      exit(1);
    }
    return c[n];
  }

  const std::vector<Eigen::VectorXcd> & InfluenceFunctional_OD::get_env_ops(int n)const{
    if(n<0||n>=(int)env_ops.size()){
      std::cerr<<"InfluenceFunctional_OD: Trying to access env_ops out of bounds "<<n<<"/"<<env_ops.size()<<std::endl;
      exit(1);
    }
    return env_ops[n];
  }

  void InfluenceFunctional_OD::check_within_limits(int n)const{
    if(get_rank()<n){
        std::cerr<<"InfluenceFunctional_OD::check_within_limits: Not enough steps calculated for influence functional!"<<std::endl;
        std::cerr<<"get_rank(): "<<get_rank()<<" n: "<<n<<std::endl;
        exit(1);
    }
  }
  void InfluenceFunctional_OD::check_env_dims(int n)const{
    if(n<0){ //in this case, check all (via recursion)
      for(int n=0; n<a.size(); n++){
        check_env_dims(n);
      }
      return;
    }

    if(n>=a.size()){
      std::cerr<<"InfluenceFunctional_OD::check_env_dims("<<n<<"): n>=a.size() ("<<n<<" vs. "<<a.size()<<")!"<<std::endl;
      exit(1);
    }

    if(c.size()!=a.size()){
      std::cerr<<"InfluenceFunctional_OD::check_env_dims("<<n<<"): c.size()!=a.size() ("<<c.size()<<" vs. "<<a.size()<<")!"<<std::endl;
      exit(1);
    }
    if(c[n].rows()!=a[n].dim_d2){
      std::cerr<<"InfluenceFunctional_OD::check_env_dims("<<n<<"): c[n].rows()!=a[n].dim_d2 ("<<c[n].rows()<<" vs. "<<a[n].dim_d2<<")!"<<std::endl;
      exit(1);
    }

    if(env_ops.size()<=0){
      return;
    }
    if(env_ops.size()!=a.size()){
      std::cerr<<"InfluenceFunctional_OD::check_env_dims("<<n<<"): env_ops.size()!=a.size() ("<<env_ops.size()<<" vs. "<<a.size()<<")!"<<std::endl;
      exit(1);
    }
    for(int o=0; o<env_ops[n].size(); o++){
      if(env_ops[n][o].rows()!=a[n].dim_d2){
        std::cerr<<"InfluenceFunctional_OD::check_env_dims("<<n<<"): env_ops["<<n<<"]["<<o<<"].rows()!=an["<<n<<"].dim_d2 ("<<env_ops[n][o].rows()<<" vs. "<<a[n].dim_d2<<")!"<<std::endl;
        exit(1);
      }
    }
  }

  void InfluenceFunctional_OD::check_rank_is(int n_max, const std::string &context){
    if(n_max!=get_rank()){
      std::cerr<<context<<": n_max!=get_rank() ("<<n_max<<" vs. "<<get_rank()<<") !"<<std::endl;
      exit(1);
    }
  }
  
  //Define tracking of env_ops with sweeps:
  void InfluenceFunctional_OD::process_low_to_high(int n, const Eigen::MatrixXcd &R){
    if(tgrid.use_rep)cta.set_R_if_correct_n(R,n);
        
    if(n<env_ops.size())for(size_t i=0; i<env_ops[n].size(); i++){
      env_ops[n][i]=R*env_ops[n][i];
    }
    if(c.size()>n && R.cols()==c[n].rows())c[n]=R*c[n];
  }

  void InfluenceFunctional_OD::process_high_to_low(int n, const Eigen::MatrixXcd &L){
    if(tgrid.use_rep)cta.set_L_if_correct_n(L,n);

    //What we need is the inverse of L.
    //Note: L also contains (parts of) the singular values, which 
    //can be obtained from the norms of the columns.
    //L is orthogonal, but not orthonormal. 
    //To get the inverse, we must transpose+conjugate L 
    //and devide by the norms of the columns twice.
    Eigen::VectorXcd vL(L.cols());
    for(int c=0; c<L.cols(); c++)vL(c)=L.col(c).norm();
    Eigen::MatrixXcd Linv=L.adjoint();
    for(int c=0; c<Linv.rows(); c++)Linv.row(c)/=(vL(c)*vL(c));

    if(c.size()>n-1 && Linv.cols()==c[n-1].rows())c[n-1]=Linv*c[n-1];
    if(n-1<env_ops.size())for(size_t i=0; i<env_ops[n-1].size(); i++){
      env_ops[n-1][i]=Linv*env_ops[n-1][i];
    }
  }


  //Genuine IF_OD functions:
  void InfluenceFunctional_OD::calculate_closures(){
    if(a.size()<1){
      std::cerr<<"InfluenceFunctional_OD::calculate_closures: a.size()<1!"<<std::endl;
      exit(1);
    }

    int N=dict.get_N();
//    int NL=dict.get_NL();

    c.resize(a.size());
    c.back().resize(1);
    c.back()(0)=1;

    for(int n=(int)c.size()-2; n>=0; n--){
      c[n]=Eigen::VectorXcd::Zero(a[n].dim_d2);
      for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
          int i_ind=dict(((i*N+i)*N+j)*N+j);
          if(i_ind<0)continue;
          for(int d1=0; d1<a[n+1].dim_d1; d1++){
            for(int d2=0; d2<a[n+1].dim_d2; d2++){
              c[n](d1)+=c[n+1](d2)*a[n+1](i_ind, d1, d2)/((double)N);
            }
          }
        }
      }
    }
  }

  void InfluenceFunctional_OD::calculate_dict(double zero, bool verbose){
    dict.detect(*this, zero);
    if(verbose){
      std::cout<<"InfluenceFunctional_OD: Dictionary(";
      std::cout<<dict.get_reduced_dim()<<"): ";
      dict.print_beta(); 
      std::cout<<std::endl;
    }
  }

  //keep only non-redundent terms with respect to outer indices
  void InfluenceFunctional_OD::reduce_to_dict(){
    if(a.size()<1)return;
    if(a[0].dim_i!=(int)dict.beta.size()){
      std::cerr<<"InfluenceFunctional_OD::reduce_to_dict: a[0].dim_i!=dict.beta.size()!"<<std::endl;
      exit(1);
    }
    for(size_t n=0; n<a.size(); n++){
      dict.reduce_MPS_Matrix(a[n]);
    }
    calculate_closures();
  }

  //build full MPS matrices from reduced (dictionary)
  void InfluenceFunctional_OD::expand_from_dict(){
    if(a.size()<1)return;
    if(a[0].dim_i!=dict.get_reduced_dim()){
      std::cerr<<"InfluenceFunctional_OD::expand_from_dict: a[0].dim_i!=dict.get_reduced_dim()!"<<std::endl;
      exit(1);
    }
    for(size_t n=0; n<a.size(); n++){
      MPS_Matrix m(dict.beta.size(), a[n].dim_d1, a[n].dim_d2);
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


  //expand from diagonal IF
  void InfluenceFunctional_OD::expand_from_diagonal(const InfluenceFunctional_D &other, const HilbertSpaceRotation &hs_rot, double dict_zero){
    if(other.a.size()<1){
      a.clear();
      dict.set_default(2);
      calculate_closures();
      return;
    }
    
    if(hs_rot.used()){
      int N=hs_rot.U.rows();
      dict.set_default(N);
      MPS::copy(other);

std::cout<<"other.lgroups:"; for(int j=0; j<N*N; j++)std::cout<<" "<<other.lgroups[j]; std::cout<<std::endl;
std::cout<<"U:"<<std::endl<<hs_rot.U<<std::endl;

      if(other.lgroups.grp.size()!=N*N){
        std::cerr<<"InfluenceFunctional_OD::expand_from_diagonal:  other.lgroups.grp.size()!=N*N("<<other.lgroups.grp.size()<<" vs. "<<N<<"*"<<N<<")!"<<std::endl;
        exit(1);
      }
      if(other.lgroups.get_Ngrps()!=a[0].dim_i){
        std::cerr<<"InfluenceFunctional_OD::expand_from_diagonal: other.lgroups.get_Ngrps()!=a[0].dim_i ("<<other.lgroups.get_Ngrps()<<" vs. "<<a[0].dim_i<<")!"<<std::endl;
        exit(1);
      }
  
      Eigen::MatrixXcd Map=Eigen::MatrixXcd::Zero(N*N*N*N, a[0].dim_i);
      for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
          for(int k=0; k<N; k++){
            for(int l=0; l<N; l++){
              for(int m=0; m<N; m++){
                for(int o=0; o<N; o++){
                  Map(((i*N+j)*N+k)*N+l, other.lgroups[m*N+o]) +=
                    hs_rot.U(i,m) * hs_rot.U.adjoint()(m,k) *
                    hs_rot.U(l,o) * hs_rot.U.adjoint()(o,j);
                }
              }
            }
          }
        }
      }
       
      for(size_t n=0; n<a.size(); n++){
        MPS_Matrix tmp(N*N*N*N, a[n].dim_d1, a[n].dim_d2);
        tmp.set_zero(); 
        for(int index=0; index<N*N*N*N; index++){
          for(int i=0; i<a[n].dim_i; i++){
            for(int d1=0; d1<a[n].dim_d1; d1++){
              for(int d2=0; d2<a[n].dim_d2; d2++){
                tmp(index, d1, d2) += Map(index, i) * a[n](i, d1, d2);
              }
            }
          }
        }
        a[n].swap(tmp);
      }

      calculate_dict(dict_zero);
      reduce_to_dict();

    }else{
      MPS::copy(other);
      int NL=other.lgroups.sys_dim(); //other.a[0].dim_i;   

#ifdef DEBUG_EXPAND_FROM_DIAGONAL
      std::cout<<"NL: "<<NL<<std::endl;
      for(size_t i=0; i<other.lgroups.grp.size(); i++){
      std::cout<<other.lgroups.grp[i]<<" ";
      }std::cout<<std::endl;
#endif

      dict.beta.clear();
      dict.N=sqrt(NL);
      dict.beta.resize(NL*NL,-1);
      for(int i=0; i<NL; i++){
        dict.beta[i*NL+i]=other.lgroups[i];
      }
      dict.calculate_reduced_dim();

      std::cout<<"expand_from_diagonal: "<<other.lgroups.get_Ngrps()<<" -> "<<dict.get_reduced_dim()<<std::endl;


    }
    calculate_closures();
  }

  void InfluenceFunctional_OD::truncate(int n_tot){
    if(n_tot==(int)a.size()){
      return;
    }else if(n_tot>(int)a.size()||n_tot<1){
      std::cerr<<"InfluenceFunctional_OD::truncate: n_tot>a.size()||n_tot<1!"<<std::endl;
      exit(1);
    }
    MPS_Matrix & last = a[n_tot-1];
 
    MPS_Matrix A(last.dim_i, last.dim_d1, 1);
    A.set_zero();
    for(int i=0; i<last.dim_i; i++){ 
      for(int d1=0; d1<last.dim_d1; d1++){ 
        for(int d2=0; d2<last.dim_d2; d2++){
          A(i,d1,0)+=last(i,d1,d2)*c[n_tot-1](d2);
        }
      }
    } 
    last.swap(A);
    {
      std::vector<MPS_Matrix> a_swp(n_tot);
      for(int n=0; n<n_tot; n++){
        a_swp[n].swap(a[n]);
      }
      a.swap(a_swp);
    }
    
 
    //closures:
    calculate_closures();
/*
    {
      std::vector<Eigen::VectorXcd> c_swp(n_tot);
      for(int n=0; n<n_tot-1; n++){
        c_swp[n](c[n]);
      }
      c.swap(c_swp);
      c.back().resize(1);c.back()(0)=1;
    }
    c.resize(a.size());
    c.back().resize(1);
    c.back()(0)=1;
*/

    //Environment operators:
    if(env_ops.size()>0){
      std::vector<std::vector<Eigen::VectorXcd> > env_ops_swp(n_tot);
      for(int n=0; n<n_tot; n++){
        env_ops_swp[n].swap(env_ops[n]);
      }
      for(int o=0; o<(int)env_ops_swp[n_tot-1].size(); o++){
        env_ops_swp[n_tot-1][o]=Eigen::VectorXcd(1);
        env_ops_swp[n_tot-1][o](0)=(1./0.);
      }
      env_ops.swap(env_ops_swp);
    }
  }

  void InfluenceFunctional_OD::calculate_diagBB(DiagBB &diagBB, RankCompressor &compressor, double dict_zero, int n_extra){

    InfluenceFunctional_D IF_D(tgrid, diagBB, compressor, n_extra);
    expand_from_diagonal(IF_D, diagBB.hs_rot, dict_zero);

    if(n_extra>0)truncate(tgrid.n_tot);
  }


  void InfluenceFunctional_OD::calculate_single_mode(int n_max, double dt, 
                              double t_tot, SingleBathMode &mode, 
                              const Eigen::MatrixXcd &bath_init, int factor){

    int n_tot=t_tot/dt+0.5;
    //ignore n_max for the moment...

    if(n_tot<1){
      std::cerr<<"InfluenceFunctional_OD::calculate_single_mode: n_tot<1!"<<std::endl;
      exit(1);
    }
    if(bath_init.rows()!=mode.get_M()||bath_init.cols()!=mode.get_M()){
      std::cerr<<"InfluenceFunctional_OD::calculate_single_mode: bath_init.size()!=mode.get_M()!"<<std::endl;
      exit(1);
    }


    int N=mode.get_N();
    int NL=N*N;
    int M=mode.get_M();
    int ML=M*M;

    mode.calculateA(dt,factor);

    a.resize(n_tot);
    a[0].resize(NL*NL,1,ML);
    a[0].fill(0.);
    for(int alpha=0; alpha<NL; alpha++){
      for(int beta=0; beta<NL; beta++){
        for(int chi=0; chi<ML; chi++){
          for(int xi0=0; xi0<M; xi0++){
            for(int xi1=0; xi1<M; xi1++){
              a[0](alpha*NL+beta, 0, chi) +=
                mode.A[alpha][beta](chi,xi0*M+xi1) * bath_init(xi0, xi1);
            }
          }
        }
      }
    }
    for(int n=1; n<n_tot; n++){
      a[n].resize(NL*NL, ML, ML);
      for(int alpha=0; alpha<NL; alpha++){
        for(int beta=0; beta<NL; beta++){
          for(int chi=0; chi<ML; chi++){
            for(int chi_=0; chi_<ML; chi_++){
              a[n](alpha*NL+beta, chi_, chi)=mode.A[alpha][beta](chi,chi_);
            }
          }
        }
      }
    }

    MPS_Matrix Mat(NL*NL, a[n_tot-1].dim_d1, 1);
    Mat.fill(0.);
    for(int aa=0; aa<NL*NL; aa++){
      for(int chi=0; chi<a[n_tot-1].dim_d1; chi++){
        for(int xi=0; xi<M; xi++){
          Mat(aa, chi, 0)+=a[n_tot-1](aa, chi, xi*M+xi);
        }
      }
    }
    a[n_tot-1].swap(Mat);
    calculate_closures();
  }

  void InfluenceFunctional_OD::join_diagonal(const InfluenceFunctional_D &other){
    if(other.a.size()<1|| a.size()<1){
      std::cerr<<"InfluenceFunctional_OD::join_diagonal: other.a.size()<1|| a.size()<1!"<<std::endl;
      exit(1);
    }
    int NL=other.a[0].dim_i;
    if(a[0].dim_i!=NL*NL){
      std::cerr<<"InfluenceFunctional_OD::join_diagonal: a[0].dim_i!=NL*NL!"<<std::endl;
      exit(1);
    }

    size_t sz=a.size();
    if(other.a.size()<sz)sz=other.a.size();

    for(size_t n=0; n<sz; n++){
      MPS_Matrix b(NL*NL, a[n].dim_d1*other.a[n].dim_d1,a[n].dim_d2*other.a[n].dim_d2);
      b.fill(0.);
      for(int i=0; i<NL; i++){
        for(int j=0; j<NL; j++){
          for(int d1=0; d1<a[n].dim_d1; d1++){
            for(int d2=0; d2<a[n].dim_d2; d2++){
              for(int od1=0; od1<other.a[n].dim_d1; od1++){
                for(int od2=0; od2<other.a[n].dim_d2; od2++){
  b(i*NL+j, d1*other.a[n].dim_d1 + od1, d2*other.a[n].dim_d2 + od2) =
                 a[n](i*NL+j,d1,d2) * other.a[n](i,od1,od2);
                }
              }
            } 
          }
        }
      }
      a[n].swap(b);
    }
    calculate_closures();
  }

  void InfluenceFunctional_OD::calculate_dt0_ndt0(ModePropagator &mprop, int n_max, double ta, double dt, double dt0, int ndt0){

#ifdef DEBUG
    std::cout<<"DEBUG: InfluenceFunctional_OD: calculate called: n_max="<<n_max<<" ta="<<ta<<" dt="<<dt<<" dt0="<<dt0<<" ndt0="<<ndt0<<std::endl;
#endif
    a.clear();
    if(n_max<1)return;
    a.resize(n_max);
    int N=mprop.get_N_system();
    int M=mprop.get_N_mode();
    int NL=N*N;

    int ML=mprop.rBasis->override_dim( M*M );

    const Eigen::VectorXcd bath_init=
      mprop.rBasis->transform_init(H_Matrix_to_L_Vector(mprop.get_bath_init()));

    //Environment operators
    env_ops.resize(n_max-1); 
    for(int n=0; n<n_max-1; n++){ 
      env_ops[n].resize(mprop.env_ops.size()+1);

      //first one must be identiy:
if(!mprop.rBasis->use()){
      env_ops[n][0]=Eigen::VectorXcd::Zero(ML);
      for(int m=0; m<M; m++)env_ops[n][0](m*M+m)=1;
}else{
      env_ops[n][0]=mprop.rBasis->id_op(M);
}//rBasis

      //others: id*this_op + old_op*id: 
      for(size_t i=0; i<mprop.env_ops.size(); i++){
        if(mprop.env_ops[i].cols()!=M || mprop.env_ops[i].rows()!=M){
          std::cerr<<"mprop.env_ops[i].cols()!=M || mprop.env_ops[i].rows()!=M"<<std::endl;
          exit(1);
        }
if(!mprop.rBasis->use()){
        env_ops[n][i+1]=Eigen::VectorXcd::Zero(ML);
        for(int m1=0; m1<M; m1++){
          for(int m2=0; m2<M; m2++){
            env_ops[n][i+1](m1*M+m2)=mprop.env_ops[i](m2,m1);
          }
        }
}else{
        env_ops[n][i+1]=mprop.rBasis->transform_op(H_Matrix_to_L_Vector(mprop.env_ops[i].transpose()));
}//rBasis
      }
    }

    std::vector<Eigen::VectorXcd> endv(env_ops[0].size(), Eigen::VectorXcd::Zero(1));
    env_ops.push_back(endv);



    //actual MPR:
//Would be better: move treatment of rBasis to FreePropagator
//but for now:
if(!mprop.rBasis->use()){
    for(int n=0; n<n_max; n++){
      a[n].resize(NL*NL, ML, ML); 
      if(n<ndt0){
        mprop.update(ta+n*dt0, dt0); 
      }else{
        mprop.update(ta+ndt0*dt0+(n-ndt0)*dt, dt);
      }
      for(int i=0; i<NL; i++){
        for(int j=0; j<NL; j++){
          Eigen::MatrixXcd &A=mprop.A[i][j];
          for(int d1=0; d1<ML; d1++){
            for(int d2=0; d2<ML; d2++){
              a[n](i*NL+j, d1, d2)=A(d2,d1);
            }
          }
        }
      }
    }
}else{
    Eigen::MatrixXcd A;
    for(int n=0; n<n_max; n++){
      a[n].resize(NL*NL, ML, ML); 
      //only recalculate calculated propagator if previous can't be used
      if(mprop.has_to_recalculate(dt)||n==0){ 
        Eigen::MatrixXcd gen;
        if(n<ndt0){
          gen=mprop.Total_Generator(ta+n*dt0,dt0);
        }else{
          gen=mprop.Total_Generator(ta+ndt0*dt0+(n-ndt0)*dt, dt);
        }
        Eigen::MatrixXcd dgen=Disentangle_Propagator(gen, sqrt(NL));
        Eigen::MatrixXcd B=mprop.rBasis->transform_prop(dgen, sqrt(NL));
        A=B.exp();
      }

/*
      mprop.update(ta+n*dt, dt);
      Eigen::MatrixXcd disM=Disentangle_Propagator(mprop.M, sqrt(NL));
      A=mprop.rBasis->transform_prop(disM, sqrt(NL));
*/

      for(int i=0; i<NL; i++){
        for(int j=0; j<NL; j++){
//          Eigen::MatrixXcd A=mprop.rBasis->transform_prop(mprop.A[i][j].transpose());
          for(int d1=0; d1<ML; d1++){
            for(int d2=0; d2<ML; d2++){
              a[n](i*NL+j, d1, d2)=A(i*ML+d2,j*ML+d1);
            }
          }
        }
      }
    }
}//rBasis

    
    MPS_Matrix ainit(NL*NL, 1, ML);
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

    MPS_Matrix &a_ref=a.back();
    MPS_Matrix afinal(NL*NL, a_ref.dim_d1, 1);
    afinal.fill(0.);
    for(int i=0; i<NL; i++){
      for(int j=0; j<NL; j++){
        for(int d1=0; d1<a_ref.dim_d1; d1++){
if(!mprop.rBasis->use()){
          for(int xi0=0; xi0<M; xi0++){
            afinal(i*NL+j, d1, 0) += a_ref(i*NL+j, d1, xi0*M+xi0);
          }
}else{
          Eigen::VectorXcd id_op=mprop.rBasis->id_op(M);
          for(int d2=0; d2<ML; d2++){
            afinal(i*NL+j, d1, 0) += a_ref(i*NL+j, d1, d2)*id_op(d2);
          }
}//rBasis
        }
      }
    }
    a_ref.swap(afinal);
    
    dict.set_default(N);
    calculate_closures();
  }

  //second order combination of two MPS matrices
  void InfluenceFunctional_OD::single_join_halfdt(
      MPS_Matrix &M, const MPS_Matrix &M1, const MPS_Matrix &M2, 
      const IF_OD_Dictionary &thisdict, const IF_OD_Dictionary &odict, 
      const IF_OD_Dictionary &ndict){

      int NL=thisdict.get_NL();
      std::vector<std::vector<int> > rev=thisdict.get_reverse_beta();
      std::vector<std::vector<int> > rev2=odict.get_reverse_beta();
      std::vector<std::vector<int> > newrev=ndict.get_reverse_beta();


      MPS_Matrix btmp(ndict.get_reduced_dim(), M.dim_d1*M1.dim_d1, M.dim_d2*M1.dim_d2);
      btmp.fill(0.);

      MPS_Matrix b(ndict.get_reduced_dim(), M.dim_d1*M1.dim_d1, M.dim_d2*M2.dim_d2);
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

  //join with other IF_OD under the assumption that the other IF_OD is discretized with time steps dt/2  ->  Symmetric Trotter
  void InfluenceFunctional_OD::join_and_sweep_halfdt(const InfluenceFunctional_OD &other, RankCompressor &compressor, double keep_weight, bool do_sweep){

#ifdef NOSWEEP 
    do_sweep=false;
#endif

    //Implement later: other.a.size()>2*a.size() can be handeled using closures
    if(2*a.size()!=other.a.size()){
      std::cerr<<"InfluenceFunctional_OD::join: 2*a.size()!=other.size()!"<<std::endl; 
      exit(1);
    }
    if(a.size()<1)return;
//    if(a[0].dim_i != other.a[0].dim_i)
    if(dict.beta.size() != other.dict.beta.size()){
      std::cerr<<"InfluenceFunctional_OD::join: a[0].dim_i != other.a[0].dim_i!"<<std::endl; 
      exit(1);
    }
    int NL=sqrt(dict.beta.size());
    if(dict.beta.size()!=NL*NL){
      std::cerr<<"InfluenceFunctional_OD::join: a[0].dim_i != NL*NL!"<<std::endl; 
      exit(1);
    }

    IF_OD_Dictionary newdict(dict); 
    newdict.join(other.dict);

    
    //forward sweep
    for(size_t n=0; n<a.size(); n++){

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
          std::cerr<<"Problem joining environment operators: env_ops["<<n;
          std::cerr<<"].size()!=other.env_ops["<<n<<"].size() (";
          std::cerr<<env_ops[n].size()<<" vs. "<<other.env_ops[n].size()<<")!";
          std::cerr<<std::endl;
          exit(1);
        }
      }else{
      }


      if(tgrid.use_rep && tgrid.rep_unit>0 && tgrid.rep_unit==(int)n-1){
        single_join_halfdt(rep.M, other.a[2*tgrid.rep_unit], other.a[2*tgrid.rep_unit+1], dict, other.dict, newdict);

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

    if(do_sweep&&tgrid.use_rep){
      rep.apply_compress_trafo(cta, true);
    }

    dict=newdict;
    calculate_closures();

  }
  void InfluenceFunctional_OD::join_halfdt(const InfluenceFunctional_OD &other){
    RankCompressor_SVD dummy(0.);
    join_and_sweep_halfdt(other, dummy, 0., false);
  }

  void InfluenceFunctional_OD::normal_form_fwd(){
    Eigen::MatrixXcd lastR=Eigen::MatrixXcd::Identity(1,1);
    for(int n=0; n<a.size()-1; n++){
      a[n].inner_multiply_left(lastR);
      Eigen::MatrixXcd A = a[n].get_Matrix_d1i_d2();
      Eigen::JacobiSVD<Eigen::MatrixXcd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
//std::cout<<"n="<<n<<": A.rows()="<<A.rows()<<" A.cols()="<<A.cols()<<" svd.singularValues().size()="<<svd.singularValues().size()<<" svd.matrixU().rows()="<<svd.matrixU().rows()<<" svd.matrixU().cols()="<<svd.matrixU().cols()<<" svd.matrixV().rows()="<<svd.matrixV().rows()<<" svd.matrixV().cols()="<<svd.matrixV().cols()<<std::endl;
//      std::cout<<n<<": "<<svd.singularValues().transpose()<<std::endl;

      Eigen::VectorXcd diag=svd.singularValues();
      double keep=1.; if(diag.size()>0) keep = svd.singularValues()(0);
      diag/=keep;


      Eigen::MatrixXcd U=svd.matrixU();
      Eigen::MatrixXcd V=svd.matrixV();
      for(int c=0; c<U.cols(); c++){
        std::complex<double> max=0.;
        for(int r=0; r<U.rows(); r++){
          if(abs(U(r,c))>abs(max))max=U(r,c);
        }
        if(abs(max)<1e-12)max=1.;
        max/=abs(max);
        U.col(c)/=max;
        V.col(c)/=max;
      }

      a[n].set_from_Matrix_d1i_d2(U*keep, a[n].dim_i);
      lastR = diag.asDiagonal()*V.adjoint();
      process_low_to_high(n, lastR);
    }
    a.back().inner_multiply_left(lastR);
    calculate_closures();
  }

  void InfluenceFunctional_OD::print_closures(const std::string &fname){
    std::ofstream ofs(fname.c_str());
    for(size_t i=0; i<c.size(); ++i){
      ofs<<c[i].transpose()<<std::endl;
    }
  }

  void InfluenceFunctional_OD::read_env_ops_binary(std::istream &ifs){
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

  void InfluenceFunctional_OD::read_binary(const std::string &filename){
    std::ifstream ifs(filename.c_str(), std::ios::binary);
    if(!ifs.is_open()){
      std::cout<<"Error reading InfluenceFunctional_OD: Cannot open file '"<<filename<<"'!"<<std::endl;  
      exit(1);
    }
    if(ProcessTensor_real::is_proper_file_format(filename,false)){
      ProcessTensor_real PT; 
      PT.read_binary(filename);
      *this=PT.get_IF_OD();
      return;
    }

    char buf[5]; buf[4]='\0';
    ifs.read(buf,4);
    int version;
    ifs.read((char*)&version, sizeof(int));
    
    std::cout<<"Reading process tensor (OD) '"<<filename<<"'"<<std::endl;
    if(std::string(buf)=="PTOD"){
      std::cout<<"File format version: "<<version<<std::endl;
      if(version!=2){
        std::cerr<<"Expected file format version 2!"<<std::endl;
        exit(1);
      }
    }else if(std::string(buf)=="PT__" || std::string(buf)=="PTr_"){ 
      ifs.close();
      std::cout<<"PT file '"<<filename<<"': converting using IF_from_PT"<<std::endl;

//      ProcessTensor PT(filename);
//      InfluenceFunctional_OD IF2=IF_from_PT(PT);
      *this=IF_from_PT(ProcessTensor(filename));

//      InfluenceFunctional_OD IF2=IF_from_PT(ProcessTensor(filename));
/*
      initialize();
      tgrid=IF2.tgrid;
      dict=IF2.dict;
      c=IF2.c;
      env_ops=IF2.env_ops;
      a=IF2.a;
*/
//      calculate_closures();

      return;
    }else{ 
      std::cerr<<"Error reading process tensor (OD): Wrong format!"<<std::endl;
      exit(1);
    }
    dict.read_binary(ifs);
dict.print_beta();std::cout<<std::endl;


    read_env_ops_binary(ifs);

    MPS::read_binary(ifs);

    calculate_closures();

    tgrid.n_tot=a.size();
    tgrid.n_calc=tgrid.n_tot;
  }

  void InfluenceFunctional_OD::write_env_ops_binary(std::ostream &ofs)const{
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

  void InfluenceFunctional_OD::write_binary(const std::string &filename)const{
std::cout<<"Writing process tensor to binary '"<<filename<<"'"<<std::endl;
dict.print_beta(); std::cout<<std::endl;
    std::ofstream ofs(filename.c_str(), std::ios::binary);
    char buf[5]="PTOD";
    int version=2;
    ofs.write(buf, 4);
    ofs.write((char*)&version, sizeof(int));
    dict.write_binary(ofs);
    write_env_ops_binary(ofs);
    MPS::write_binary(ofs);
  }

  void InfluenceFunctional_OD::set_none(int n_max, int N, double dict_zero){
#ifdef OLD_IF_SET_NONE
    int NL=N*N;
    a.resize(n_max);
    for(int i=0; i<n_max; i++){  
      a[i].resize(NL*NL, 1, 1);
      for(int j=0; j<NL; j++){
        a[i](j*NL+j,0,0)=1.;
      }
    }
    calculate_dict(dict_zero);
    reduce_to_dict();
#else
    a.resize(n_max);
    for(int i=0; i<n_max; i++){
      a[i].resize(1, 1, 1);
      a[i](0,0,0)=1.;
    }
    int NL=N*N;
    dict.N=N;
    dict.reduced_dim=1;
    dict.beta.clear();
    dict.beta.resize(NL*NL,-1);
    for(size_t l=0; l<NL; l++){
      dict.beta[l*NL+l]=0;
    }
#endif

    calculate_closures();

    env_ops.clear();
    env_ops.resize(n_max);
    for(size_t i=0; i<env_ops.size(); i++){
      env_ops[i].resize(1,Eigen::MatrixXcd::Identity(1,1));
    }
    rep.set_default(N);
  }

  void InfluenceFunctional_OD::sweep(RankCompressor &compressor, double keep_weight){
    compressor.sweep(*this, keep_weight); 
    calculate_closures();
  }


  void InfluenceFunctional_OD::add_IF_halfdt(InfluenceFunctional_OD &IF2, RankCompressor &compressor){
    
    IF2.check_rank_is(get_rank()*2, "InfluenceFunctional_OD::add_IF");

    if(printdim)std::cout<<"Maxdim: before combination: "<<get_max_dim()<<std::endl;
    if(printdim)std::cout<<"Maxdim: new contribution: "<<IF2.get_max_dim()<<std::endl;

    join_and_sweep_halfdt(IF2, compressor, dict.get_keep_weight());
    if(printdim)std::cout<<"Maxdim: after sweep: "<<get_max_dim()<<std::endl;

    calculate_closures();
  }
   
  void InfluenceFunctional_OD::add_mode(ModePropagatorGenerator &mpg, int k, RankCompressor &compressor, double dict_zero){

    std::cout<<"Add mode "<<k<<"/"<<mpg.get_N_modes()<<std::endl;

    ModePropagatorPtr mpp=mpg.getModePropagator(k);


//    InfluenceFunctional_OD IF2(mpp.ref(), tgrid.n_calc*2, tgrid.ta, tgrid.dt/2., dict_zero);
    InfluenceFunctional_OD IF2(tgrid.construct_half_dt());
    IF2.calculate_dt0_ndt0(*mpp.get(), tgrid.n_calc*2, tgrid.ta, tgrid.dt/2., tgrid.dt0/2., 2);
    IF2.calculate_dict(dict_zero);
    IF2.reduce_to_dict();

    if(tgrid.rep_unit>0){
      rep.expand_ops(*mpp.get());
    }

    add_IF_halfdt(IF2, compressor);
  }

  void InfluenceFunctional_OD::add_modes(ModePropagatorGenerator &mpg, RankCompressor &compressor, double dict_zero){
    for(int k=mpg.first(); k<mpg.get_N_modes(); k=mpg.next(k)){
      add_mode(mpg, k, compressor, dict_zero);
    }
  }

  void InfluenceFunctional_OD::calculate(ModePropagatorGenerator &mpg, RankCompressor &compressor, double dict_zero){

    set_none_from_tgrid(mpg.get_N());    
//    calculate_dict(dict_zero);
//    reduce_to_dict();
    
std::cout<<"TEST: get_rank(): "<<get_rank()<<std::endl;
    add_modes(mpg, compressor, dict_zero);
  }
  
  void InfluenceFunctional_OD::from_Rep_GIF(const Rep_GIF &rep_, int n_max_t){ 
    Rep_GIF tmp_rep=rep_;  // <- avoid problem with aliasing!
    rep=tmp_rep;

    dict=rep.dict;

    a.clear();
    a.resize(n_max_t, rep.M);
    MPS_Matrix tmp(a[0].dim_i, 1, a[0].dim_d2);
    tmp.fill(0.);
    for(int i=0; i<a[0].dim_i; i++){
      for(int d1=0; d1<a[0].dim_d1; d1++){
        for(int d2=0; d2<a[0].dim_d2; d2++){
          tmp(i, 0, d2)+=rep.init(d1)*a[0](i, d1, d2);
        }
      }
    }
    a[0].swap(tmp);

    tmp.resize(a.back().dim_i, a.back().dim_d1,1);
    tmp.fill(0.);
    for(int i=0; i<a.back().dim_i; i++){
      for(int d1=0; d1<a.back().dim_d1; d1++){
        for(int d2=0; d2<a.back().dim_d2; d2++){
          tmp(i, d1, 0)+=a.back()(i, d1, d2) * rep.env_ops[0](d2);
        }
      }
    }
    a.back().swap(tmp);

    c.clear();
    c.resize(n_max_t, rep.env_ops[0]);

    env_ops.clear();
    env_ops.resize(n_max_t, rep.env_ops);
  }


  void InfluenceFunctional_OD::insert_rep(){
    if(tgrid.use_rep){
      if(tgrid.rep_replace){
        std::cout<<"using rep_replace"<<std::endl;
        from_Rep_GIF(rep, tgrid.n_tot);
      }else{
        std::cout<<"using rep but not rep_replace"<<std::endl;
        a.insert(a.begin()+tgrid.rep_unit+1, tgrid.n_rep, rep.M);
        c.insert(c.begin()+tgrid.rep_unit+1, tgrid.n_rep, rep.env_ops[0]);
        if(env_ops[0].size()!=rep.env_ops.size()){
          std::cerr<<"Error: InfluenceFuncitonal_OD::insert_rep: env_ops[0].size()!=rep.env_ops.size()!"<<std::endl;
          exit(1);
        }
//        for(size_t i=0; i<rep.env_ops.size(); i++){
        env_ops.insert(env_ops.begin()+tgrid.rep_unit+1, tgrid.n_rep, rep.env_ops);
//        }
      }
    }
  }
    
  InfluenceFunctional_OD::InfluenceFunctional_OD(const TimeGrid &tgr_, DiagBB &diagBB, RankCompressor &compressor, double dict_zero, int n_extra){
    initialize();
    tgrid=tgr_;
    calculate_diagBB(diagBB, compressor, dict_zero, n_extra);
    cta.n=tgrid.rep_unit;
  }

  InfluenceFunctional_OD::InfluenceFunctional_OD(int n_max, int N){
    tgrid.set_default(n_max);
    initialize();
    set_none(n_max, N);
  }

  InfluenceFunctional_OD::InfluenceFunctional_OD(const TimeGrid &tgr_, int N){
    initialize();
    tgrid=tgr_;
    set_none_from_tgrid(N);

    cta.n=tgrid.rep_unit;
//    cta.use_ortho=compress_trafo_use_ortho;
  }

  InfluenceFunctional_OD::InfluenceFunctional_OD(ModePropagator &mprop, int n_max, double ta, double dt, double dict_zero){
    initialize();
    calculate(mprop, n_max, ta, dt);
    calculate_dict(dict_zero);
    reduce_to_dict();
  }

  InfluenceFunctional_OD::InfluenceFunctional_OD(ModePropagator &mprop, const TimeGrid &tgr_, double dict_zero){
    initialize();
    tgrid=tgr_;
    cta.n=tgrid.rep_unit;

    calculate(mprop, tgrid.n_calc, tgrid.ta, tgrid.dt);
    calculate_dict(dict_zero);
    reduce_to_dict();
  }

  InfluenceFunctional_OD::InfluenceFunctional_OD(ModePropagatorGenerator &mpg, const TimeGrid &tgr_, RankCompressor &compressor, double dict_zero){
    initialize();
    tgrid=tgr_;
    cta.n=tgrid.rep_unit;
    calculate(mpg, compressor, dict_zero);
  }

  InfluenceFunctional_OD::InfluenceFunctional_OD(const std::string &filename){
    initialize();
    read_binary(filename);
  }

  InfluenceFunctional_OD::InfluenceFunctional_OD(){
    initialize();
  }

}//namespace
