#ifndef INFLUENCE_FUNCTIONAL_KRYLOV_DEFINED_H
#define INFLUENCE_FUNCTIONAL_KRYLOV_DEFINED_H

#include "InfluenceFunctional_OD.h"


struct rw_pair{
  Eigen::VectorXcd r;
  std::vector<Eigen::VectorXcd> w;

  void calculate(const std::vector<Eigen::MatrixXcd> &Q, const Eigen::VectorXcd &r_in){
    r=r_in;

    if(r.size()!=Q[0].cols()){
      std::cerr<<"InfluenceFunctional_Krylov::rw_pair: r.size()!=Q[0].cols()!"<<std::endl;
      exit(1);
    }
    w.clear();
    w.resize(Q.size());
    for(int b=0; b<Q.size(); b++){
      w[b]=Q[b]*r;
    }
  }

  rw_pair(const std::vector<Eigen::MatrixXcd> &Q, const Eigen::VectorXcd &r_in) {
    calculate(Q, r_in);    
  }

  rw_pair(const Eigen::VectorXcd &r_in=Eigen::VectorXcd()) : r(r_in){
  }
};

struct rw_pair_list{
  std::vector<rw_pair> list;
   
  size_t size()const{return list.size();}
  rw_pair &operator[](int i){return list[i];}
  void push_back(const rw_pair &rw){list.push_back(rw);}
  void clear(){list.clear();}
  void resize(size_t sz, const rw_pair &rw=rw_pair()){list.resize(sz, rw);}

};

class InfluenceFunctional_Krylov{
public:
  
  std::vector<Eigen::MatrixXcd> Q;
  std::vector<Eigen::VectorXcd> env_ops;
 
  IF_OD_Dictionary dict;
  int ref_n_max, max_k;
  double threshold, dict_zero;
  bool print_timesteps;


  InfluenceFunctional_OD get_IF_OD(int n_max) const{
    if(n_max<1){
      std::cerr<<"InfluenceFunctional_Krylov::get_IF_OD: n_max<1!"<<std::endl;
      exit(1);
    }
    int red_dim=Q.size();
    if(Q.size()<1){
      std::cerr<<"InfluenceFunctional_Krylov::get_IF_OD: Q.size()<1!"<<std::endl;
      exit(1);
    }
    if(dict.get_reduced_dim()!=red_dim){
      std::cerr<<"InfluenceFunctional_Krylov::get_IF_OD: dict.get_reduced_dim()!=Q.size()!"<<std::endl;
      exit(1);
    }

    int rdim=Q[0].rows();
    if(rdim<1){      
      std::cerr<<"InfluenceFunctional_Krylov::get_IF_OD: Q[0].rows()<1!"<<std::endl;
      exit(1);
    }
    if(env_ops.size()<1){
      std::cerr<<"InfluenceFunctional_Krylov::get_IF_OD: env_ops.size()<1!"<<std::endl;
      exit(1);
    }
    if(env_ops[0].size()!=rdim){
      std::cerr<<"InfluenceFunctional_Krylov::get_IF_OD: env_ops[0].size()!=rdim!"<<std::endl;
      exit(1);
    }

    InfluenceFunctional_OD IF;

    IF.dict=dict;
    MPS_Matrix Mref(red_dim, rdim, rdim);
    for(int b=0; b<red_dim; b++){
      for(int m1=0; m1<rdim; m1++){
        for(int m2=0; m2<rdim; m2++){
          Mref(b, m1, m2)=Q[b](m2,m1);
        }
      }
    }

    IF.a.clear();
    IF.a.resize(n_max, Mref);

    //last in chain:
    IF.a.back().resize(Mref.dim_i, Mref.dim_d1, 1); 
    IF.a.back().fill(0.);
    for(int b=0; b<dict.get_reduced_dim(); b++){
      for(int m1=0; m1<rdim; m1++){
        for(int m2=0; m2<rdim; m2++){
          IF.a.back()(b, m1, 0)+=env_ops[0](m2)*Mref(b, m1, m2);
        }
      }
    }
    //first in chain:
    Mref.resize(IF.a[0].dim_i, 1, IF.a[0].dim_d2);
    Mref.fill(0.);
    for(int b=0; b<IF.a[0].dim_i; b++){
      for(int m2=0; m2<IF.a[0].dim_d2; m2++){
        Mref(b, 0, m2)=IF.a[0](b, 0, m2);
      }
    }
    IF.a[0]=Mref;
 
//    IF.calculate_closures();
    IF.c.clear();
    IF.c.resize(n_max+1, env_ops[0]);
    IF.c[n_max]=Eigen::VectorXcd::Identity(1,1);
    //TODO: env_ops!
    return IF;
  }

  operator InfluenceFunctional_OD() const{
    return get_IF_OD(ref_n_max);
  }
 
  void set_none(int n_max, int N){
    ref_n_max=n_max; 
    dict.set_default(N);
//    dict.set_default_diag(N);
    Q.clear();
    Q.resize(dict.get_NL2(),Eigen::MatrixXcd::Zero(1,1));
    for(int i=0; i<N; i++){
      for(int j=0; j<N; j++){
        Q[(i*N+j)*N*N+(i*N+j)]=Eigen::MatrixXcd::Identity(1,1);
      }
    }

    env_ops.clear();
    env_ops.push_back(Eigen::VectorXcd(1));
    env_ops[0](0)=1;
  }


  void add_mode(ModePropagator &mprop, double ta, double dt, const Eigen::MatrixXcd &bath_init){
    double lthresh=log(threshold);

    int N=mprop.get_N_system();
    int NL=N*N;
    if(Q.size()<1){
      std::cerr<<"InfluenceFunctional_Krylov::add_mode: Q.size()<1!"<<std::endl;
      exit(1);
    }
    if(dict.get_NL()!=NL){
      std::cerr<<"InfluenceFunctional_Krylov::add_mode: dict.get_NL()!=NL!"<<std::endl;
      exit(1);
    }
    if(Q.size()!=dict.get_reduced_dim()){
      std::cerr<<"InfluenceFunctional_Krylov::add_mode: Q.size()=dict.get_reduced_dim()!"<<std::endl;
      exit(1);
    }
 
    int ML=Q[0].rows();
    if(ML<1){
      std::cerr<<"InfluenceFunctional_Krylov::add_mode: ML<1!"<<std::endl;
      exit(1);
    }
   
    mprop.update(ta, dt/2);
    int M=mprop.get_N_mode();
    int MLnew=mprop.A[0][0].rows();  
    int MLcombined=MLnew*ML;
    std::vector<Eigen::MatrixXcd> Qbck=Q;
    for(int i=0; i<NL; i++){
      for(int l=0; l<NL; l++){
        Q[i*NL+l]=Eigen::MatrixXcd::Zero(MLcombined, MLcombined);
        for(int j=0; j<NL; j++){
          for(int k=0; k<NL; k++){
            for(int dd1=0; dd1<ML; dd1++){
              for(int dd2=0; dd2<ML; dd2++){
                for(int d1=0; d1<MLnew; d1++){
                  for(int d2=0; d2<MLnew; d2++){
                    for(int d3=0; d3<MLnew; d3++){
                      Q[i*NL+l](dd2*MLnew+d2,dd1*MLnew+d1)+=
                         mprop.A[i][j](d2,d3)*
                         Qbck[j*NL+k](dd2,dd1)*
                         mprop.A[k][l](d3,d1);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

   
    //get rw from matrices:
    Eigen::VectorXcd new_bath_init(MLnew);
    for(int i=0; i<M; i++){
      for(int j=0; j<M; j++){
        new_bath_init(i*M+j)=bath_init(i,j);
      }
    }

    

    Eigen::VectorXcd combined_initial=Eigen::VectorXcd::Zero(MLcombined);
    for(int i=0; i<MLnew; i++)combined_initial(i)=new_bath_init(i);

    int this_max_k=MLcombined;if(max_k>0 && max_k<MLcombined)this_max_k=max_k;
    rw_pair_list rw;
    rw.resize(1);
    rw[0].calculate(Q, combined_initial);
    std::vector<double> weights(1,1.);

    int layer=0;
    bool done=false;
    while(!done){  //layers
      int initial_rw_size=rw.size();
      std::vector<double> new_weights(initial_rw_size, 0.);

      for(int cur_r=0; cur_r<initial_rw_size; cur_r++){
        if(rw.size()>=this_max_k){done=true; break;}
        for(size_t i=0; i<rw[cur_r].w.size(); i++){
          if(rw.size()>=this_max_k){done=true; break;}
          Eigen::VectorXcd w=rw[cur_r].w[i];
        
          //orthogonalize wrt. r:
          double norm=0;
          //loop to reothogonalize multiple times:
          for(int loop=0; loop<2; loop++){ 
            for(size_t j=0; j<rw.size(); j++){
              std::complex<double> overlap=rw[j].r.dot(w);
              w-=overlap*rw[j].r;
              double o2=(overlap.real()*overlap.real()+overlap.imag()*overlap.imag());
              new_weights[j]+=o2*weights[cur_r];
            }
            if(loop==0){ norm=(w.dot(w)).real();}
            w.normalize();
          }
          double nnorm=norm*weights[cur_r];
          if(nnorm<threshold)continue;

          rw.push_back(rw_pair(Q, w));
          new_weights.push_back(nnorm);
        }
      }
      weights=new_weights;
      double totalweight=0;
      for(size_t i=0; i<weights.size(); i++)totalweight+=weights[i];
      for(size_t i=0; i<weights.size(); i++)weights[i]/=totalweight;

      if(rw.size()==initial_rw_size){done=true; break;}
    }


/* 
    std::cout<<"r: "<<std::endl;
    for(size_t i=0; i<rw.size(); i++){
      std::cout<<i<<": "<<rw[i].first.transpose()<<std::endl;
    }
*/
    std::cout<<"r weights: "<<std::endl;
    for(size_t i=0; i<rw.size(); i++){
      std::cout<<i<<": "<<weights[i]<<std::endl;
    }

#ifdef CHECK_ORTHO
    std::cout<<"check orthogonality "<<std::endl;
    for(size_t i=0; i<rw.size(); i++){
      for(size_t j=0; j<rw.size(); j++){
        double norm=abs(rw[i].r.dot(rw[j].r));
        if((i==j && fabs(norm-1.)>1e-12) ){
          std::cout<<"i: "<<i<<" j: "<<j<<": "<<rw[i].r.dot(rw[j].r)-1.<<std::endl;
        } 
        if(i!=j && fabs(norm)>1e-12){
          std::cout<<"i: "<<i<<" j: "<<j<<": "<<rw[i].r.dot(rw[j].r)<<std::endl;
        } 
      }
    }
#endif


    //setup reduced Q:
    for(int b=0; b<dict.get_reduced_dim(); b++){
      Q[b]=Eigen::MatrixXcd::Zero(rw.size(), rw.size());
      for(size_t m1=0; m1<rw.size(); m1++){
        for(size_t m2=0; m2<rw.size(); m2++){
          Q[b](m2, m1)=rw[m2].r.dot(rw[m1].w[b]);
        }
      }
    }

 


    //env_ops:
    Eigen::VectorXcd newtrace=Eigen::VectorXcd::Zero(M*M);
    for(int i=0; i<M; i++)newtrace(i*M+i)=1;
    Eigen::VectorXcd tracev=Vector_otimes(env_ops[0], newtrace);

    env_ops[0]=Eigen::VectorXcd::Zero(rw.size());
//    for(size_t i=0; i<rw.size(); i++)env_ops[0](i)=tracev.dot(rw[i].first);
    for(size_t i=0; i<rw.size(); i++)env_ops[0](i)=rw[i].r.dot(tracev);
  
//    IF_OD_Dictionary dict2;
//    dict2.detect(m,zero);
//    dict.join(dict2);
  }


  void add_modes(ModePropagatorGenerator &mpg, int n_max, double ta, double dt, double thresh, double dictz, double max_k_=0){

    threshold=thresh;
    dict_zero=dictz;
    if(max_k_>0)max_k=max_k_;

    for(int i=0; i<mpg.get_N_modes(); i++){
      std::cout<<"Mode "<<i<<std::endl;
      ModePropagatorPtr mpp=mpg.getModePropagator(i);
      add_mode(mpp.ref(), ta, dt, mpg.get_bath_init(i));
      std::cout<<"Dimension "<<Q[0].rows()<<std::endl;
    }

  }

  InfluenceFunctional_Krylov(int n_max=1, int N=2){
    set_none(n_max, N);
    print_timesteps=false;
    threshold=1e-12;
    dict_zero=1e-20;
    max_k=0;
  }


  InfluenceFunctional_Krylov(ModePropagatorGenerator &mpg, int n_max, double ta, double dt, double thresh=1e-16, double dict_zero_=1e-20, double max_k_=0)
   : threshold(thresh), dict_zero(dict_zero_), max_k(max_k_) {
    if(mpg.get_N_modes()<1){
      std::cerr<<"InfluenceFunctional_Krylov::calculate: mpg.get_N_modes()<1!"<<std::endl;
      exit(1);
    }
    int N=mpg.getModePropagator(0)->get_N_system();
    set_none(n_max,N);
    print_timesteps=false;
    add_modes(mpg, n_max, ta, dt, threshold, dict_zero);
  }
  virtual ~InfluenceFunctional_Krylov(){}
};



#endif
