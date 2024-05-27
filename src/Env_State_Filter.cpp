#include "Env_State_Filter.hpp"
#include "Parameters.hpp"
#include "MPS.hpp"
#include "otimes.hpp"

namespace ACE{

template <typename Scalar_Type>
Eigen::Matrix<Scalar_Type, Eigen::Dynamic, Eigen::Dynamic> order_trafo(
                         Eigen::Matrix<Scalar_Type, Eigen::Dynamic, 1> &q){

  int dim=q.rows();
  Eigen::Matrix<Scalar_Type, Eigen::Dynamic, Eigen::Dynamic> R=
Eigen::Matrix<Scalar_Type, Eigen::Dynamic, Eigen::Dynamic>::Identity(dim, dim);
    
  int i_largest=0; 
  for(int i=1; i<q.size(); i++){
    if(std::abs(q(i))>std::abs(q(i_largest)))i_largest=i;
  }
  R.col(i_largest)=Eigen::VectorXd::Zero(R.rows()); 
  R.col(i_largest)(0)=1.;
  R.col(0)=q;
  R.col(0).normalize();
  for(int i=1; i<R.cols(); i++){
    for(int j=0; j<i; j++){
      Scalar_Type scalar=R.col(j).dot(R.col(i));
      R.col(i)-=scalar*R.col(j);
      R.col(i).normalize();
    }
  }
  return R;
}

template Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> order_trafo(
                   Eigen::Matrix<double, Eigen::Dynamic, 1> &q);

template Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> order_trafo(
                Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> &q);


  void Env_State_Filter::rotate_identity_to_first(std::vector<MPS_Matrix_real> &a, std::vector<std::vector<Eigen::VectorXcd> > &env_ops, int n){

    if(n==a.size()-1)return;
    if(n>a.size()-1){
      std::cerr<<"ProcessTensor_real::rotate_identity_to_first: n>=a.size()-1!"<<std::endl;
      exit(1);
    }

    if(env_ops.size()<a.size()){
      if(a[n].dim_d2!=1){
        std::cerr<<"ProcessTensor_real::rotate_identity_to_first: env_ops.size()<a.size()!"<<std::endl;
        exit(1);
      }else{
        Eigen::VectorXcd template_(1); template_<<1.;
        env_ops=std::vector<std::vector<Eigen::VectorXcd> >(a.size(), 
                            std::vector<Eigen::VectorXcd>(1, template_));
      }
    }

    if(env_ops[n].size()<1){
      if(a[n].dim_d2!=1){
        std::cerr<<"ProcessTensor_real::rotate_identity_to_first: env_ops[n].size()<1!"<<std::endl;
        exit(1);
      }else{
        Eigen::VectorXcd template_(1); template_<<1.;
        env_ops[n].push_back(template_);
      }
    }
    if(env_ops[n][0].rows()!=a[n].dim_d2){
      std::cerr<<"ProcessTensor_real::rotate_identity_to_first: env_ops[n][0].rows()!=a[n].dim_d2!"<<std::endl;
      exit(1);
    }

    if(a[n].dim_d2<2)return;
    if(a[n].dim_d2!=env_ops[n][0].rows()){
      std::cerr<<"ProcessTensor_real::rotate_identity_to_first: a[n].dim_d2!=q.size()!"<<std::endl;
      exit(1);
    }

    Eigen::VectorXd q=env_ops[n][0].real();
    Eigen::MatrixXd R=order_trafo<double>(q);
//print_diff_from_ortho(R, 1e-10);

    {
      MPS_Matrix_real tmp(a[n].dim_i, a[n].dim_d1, a[n].dim_d2);
      tmp.set_zero();
      for(int i=0; i<a[n].dim_i; i++){  
        for(int d1=0; d1<a[n].dim_d1; d1++){
          for(int d2=0; d2<a[n].dim_d2; d2++){
            for(int d=0; d<a[n].dim_d2; d++){
              tmp(i,d1,d2)+=a[n](i,d1,d)*R(d,d2);
            }
          } 
        }
      }
      a[n].swap(tmp);
       
      for(int o=0; o<env_ops[n].size(); o++){
        env_ops[n][o]=R.transpose()*env_ops[n][o];
      }
    }

    {
      MPS_Matrix_real tmp(a[n+1].dim_i, a[n+1].dim_d1, a[n+1].dim_d2);
      tmp.set_zero();
      for(int i=0; i<a[n+1].dim_i; i++){  
        for(int d1=0; d1<a[n+1].dim_d1; d1++){
          for(int d2=0; d2<a[n+1].dim_d2; d2++){
            for(int d=0; d<a[n+1].dim_d1; d++){
              tmp(i,d1,d2)+=R(d,d1)*a[n+1](i,d,d2);
            }
          } 
        }
      }
      a[n+1].swap(tmp);
    }
  }

  void Env_State_Filter::preprocess(std::vector<MPS_Matrix_real> &a1,
                  std::vector<std::vector<Eigen::VectorXcd> > &env_ops1, 
                  int n1,
                  std::vector<MPS_Matrix_real> &a2,
                  std::vector<std::vector<Eigen::VectorXcd> > &env_ops2, 
                  int n2){

    if(mean_field){
      if(!no_rotate_first){
        rotate_identity_to_first(a1, env_ops1, n1);
        rotate_identity_to_first(a2, env_ops2, n2);
      }

      if(n1==0){
        last_selected.clear(); last_selected.push_back(std::make_pair(0,0));
        last_dim1=1; last_dim2=1;
      }else{
        last_selected.swap(selected);
        last_dim1=dim1; last_dim2=dim2;
      }

      if(n1==a1.size()-1){
        selected.clear(); selected.push_back(std::make_pair(0,0));
        dim1=1; dim2=1;
      }else{
        selected.clear();
        dim1=a1[n1].dim_d2; dim2=a2[n2].dim_d2;

        if(nr_MF<1){std::cerr<<"Env_State_Filter: nr_MF<1!"<<std::endl;exit(1);}

/*
        for(int i2=nr_MF; i2<dim2; i2++){
          selected.push_back(std::make_pair(0,i2));
        }
        for(int i1=0; i1<dim1; i1++){
          for(int m=0; m<nr_MF && m<dim2; m++){
           selected.push_back(std::make_pair(i1,m));
          }
        }
*/
        for(int m=0; m<nr_MF && m<dim1; m++){
          for(int i2=0; i2<dim2; i2++){
           selected.push_back(std::make_pair(m,i2));
          }
        }
        for(int i1=nr_MF; i1<dim1; i1++){
          selected.push_back(std::make_pair(i1,0));
        }
      }
    }
  }
  
  void Env_State_Filter::filter(std::vector<MPS_Matrix_real> &a, 
             std::vector<std::vector<Eigen::VectorXcd> > &env_ops, 
             int n){
    
    if(n<0||n>=a.size()||n>=env_ops.size()){
      std::cerr<<"Env_State_Filter::filter: n<0||n>=a.size()||n>=env_ops.size()!"<<std::endl;
      exit(1);
    }

    if(mean_field){
/*
      std::cout<<"n="<<n<<" dim2="<<dim2<<" last_dim2="<<last_dim2;
      std::cout<<" selected.size()="<<selected.size();
      std::cout<<" last_selected.size()="<<last_selected.size();
      std::cout<<" a[n].dim_d1="<<a[n].dim_d1;
      std::cout<<" a[n].dim_d2="<<a[n].dim_d2;
      std::cout<<std::endl;
*/
      if(false){
        last_dim2=dim2=1;
        last_selected.clear(); 
        for(int i=0; i<a[n].dim_d1; i++)last_selected.push_back(std::make_pair(i,0));
        selected.clear(); 
        for(int i=0; i<a[n].dim_d2; i++)selected.push_back(std::make_pair(i,0));
      }

      MPS_Matrix_real tmp(a[n].dim_i, last_selected.size(), selected.size());
      tmp.set_zero();
      for(int i=0; i<tmp.dim_i; i++){
        for(int o1=0; o1<tmp.dim_d1; o1++){
          for(int o2=0; o2<tmp.dim_d2; o2++){
            tmp(i,o1,o2)=a[n](i, 
               last_selected[o1].first * last_dim2 + last_selected[o1].second,
               selected[o2].first * dim2 + selected[o2].second);
          } 
        }
      }
      a[n].swap(tmp);
      for(size_t i=0; i<env_ops[n].size(); i++){
        Eigen::VectorXcd tmp2(selected.size());
        for(int d=0; d<selected.size(); d++){
          tmp2(d)=env_ops[n][i](selected[d].first * dim2 + selected[d].second);
        }
        env_ops[n][i]=tmp2;
      }
    }
  }

  void Env_State_Filter::setup(Parameters &param){
    mean_field=param.get_as_bool("Env_State_Filter_mean_field",false);
    nr_MF=param.get_as_int("Env_State_Filter_mean_field_nr", 1);
    no_rotate_first=param.get_as_bool("Env_State_Filter_no_rotate_first",false);
  }

  Env_State_Filter::Env_State_Filter(Parameters &param){
    setup(param);
  }
  Env_State_Filter::Env_State_Filter(){
    Parameters param;
    setup(param);
  }

}//namespace
