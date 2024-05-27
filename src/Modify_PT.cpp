#include "Modify_PT.hpp"
#include "InfluenceFunctional_OD.hpp"

namespace ACE{
namespace Modify_PT{

// PT with time step dt -> PT with time step (n_coarse*dt)
void coarse_grain(InfluenceFunctional_OD &IF, int n_coarse, double dict_zero){   
  if(n_coarse<1){
    std::cerr<<"Coarse graining PT: n_coarse<1!"<<std::endl;
    exit(1);
  }else if (n_coarse==1){
    return;
  }

  if(IF.a.size()%n_coarse!=0){
    std::cerr<<"Coarse graining PT: IF.a.size()="<<IF.a.size()<<" not divisible by n_coarse="<<n_coarse<<"!"<<std::endl;
    exit(1);
  }

  int intervals=IF.a.size()/n_coarse;

  MPS mps;
  mps.a.resize(intervals);
  int NL=IF.dict.get_NL();
  for(int l=0; l<intervals; l++){
    MPS_Matrix m(NL*NL, IF.a[l*n_coarse].dim_d1, IF.a[l*n_coarse].dim_d2);
    for(int i1=0; i1<NL; i1++){
      for(int i2=0; i2<NL; i2++){
        for(int d1=0; d1<m.dim_d1; d1++){
          for(int d2=0; d2<m.dim_d2; d2++){ 
            m(i1*NL+i2,d1,d2)=IF.a[l*n_coarse](IF.dict.beta[i1*NL+i2],d1,d2);
          }
        }
      }
    }
  
    for(int j=1; j<n_coarse; j++){
      MPS_Matrix m2(NL*NL, m.dim_d1, IF.a[l*n_coarse+j].dim_d2);
      m2.set_zero();
      for(int i1=0; i1<NL; i1++){
        for(int i2=0; i2<NL; i2++){
          for(int i3=0; i3<NL; i3++){
            for(int d1=0; d1<m.dim_d1; d1++){
              for(int d2=0; d2<m2.dim_d2; d2++){
                for(int d=0; d<m.dim_d2; d++){
                  m2(i1*NL+i3,d1,d2) += m(i2*NL+i3,d1,d) * 
                             IF.a[l*n_coarse+j](IF.dict.beta[i1*NL+i2],d,d2);
                }
              }
            } 
          }
        }
      }
      m.swap(m2);
    }
    mps.a[l]=m;
  }
  IF.a.swap(mps.a);

  IF.dict.set_default(IF.dict.get_N());
  IF.calculate_dict(dict_zero);
  IF.reduce_to_dict();
  IF.calculate_closures();

std::cout<<"IF.size() = "<<IF.size()<<std::endl;
std::cout<<"IF.env_ops.size() = "<<IF.env_ops.size()<<std::endl;
  if(IF.env_ops.size()>0){
std::cout<<"IF.env_ops[0].size() = "<<IF.env_ops[0].size()<<std::endl;
    std::vector<std::vector<Eigen::VectorXcd> > env_ops_bck;
    IF.env_ops.swap(env_ops_bck);
    IF.env_ops.clear();
    IF.env_ops.resize(IF.a.size());
    if(env_ops_bck.size()!=IF.env_ops.size()*n_coarse){
      std::cerr<<"Modify_PT: coarse_grain: env_ops_bck.size()!=IF.env_ops.size()*n_coarse!"<<std::endl;
      exit(1);
    }
    for(size_t i=0; i<IF.env_ops.size(); i++){
      IF.env_ops[i]=env_ops_bck[i*n_coarse+(n_coarse-1)];
    }
  }
std::cout<<"IF.env_ops.size() = "<<IF.env_ops.size()<<std::endl;
  IF.rep.set_default(IF.dict.get_N());
}


void apply_system_propagator(InfluenceFunctional_OD &IF, Propagator &prop, double ta, double dt, double dict_zero){
   
  if(IF.dict.get_N()!=prop.get_dim()){
    std::cerr<<"Modify_PT: apply_system_propagator: IF.dict.get_N()!=prop.get_dim()!"<<std::endl;
    exit(1);
  }
  
  int NL=IF.dict.get_NL();
  for(int n=0; n<IF.a.size(); n++){
    double t=ta+n*dt;
    MPS_Matrix &aref=IF.a[n];

    MPS_Matrix m(NL*NL, aref.dim_d1, aref.dim_d2);
    m.set_zero();

    prop.update(t,dt/2.);
    Eigen::MatrixXcd M1=prop.M;
    prop.update(t+dt/2.,dt/2.);
    Eigen::MatrixXcd M2=prop.M;
    
    for(int i1=0; i1<NL; i1++){
      for(int i2=0; i2<NL; i2++){
        for(int i3=0; i3<NL; i3++){
          for(int i4=0; i4<NL; i4++){
            for(int d1=0; d1<aref.dim_d1; d1++){
              for(int d2=0; d2<aref.dim_d2; d2++){
                m(i4*NL+i1, d1, d2)+=M2(i4,i3)*aref(i3*NL+i2, d1, d2)*M1(i2,i1);
              }
            }
          }
        }
      }  
    }
    aref.swap(m);
  }
  IF.dict.set_default(IF.dict.get_N());
  IF.calculate_dict(dict_zero);
  IF.reduce_to_dict();
}

}//namespace
}//namespace


