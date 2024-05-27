#include "Parameters.hpp"
#include "Reader.hpp"
#include "ProcessTensorBuffer.hpp"
#include "DummyException.hpp"

using namespace ACE;

int main(int args, char ** argv){
//  Parameters param(args, argv, true, false);
 try{
  Parameters param(args, argv);

  std::string read_PT=param.get_as_string_check("read_PT");
  std::string write_PT=param.get_as_string_check("write_PT");

//  bool RWA=param.get_as_bool("RWA",false);

  int coarse_grain=param.get_as_size_t("coarse_grain",1);
  int buffer_blocksize=param.get_as_size_t("buffer_blocksize",0);

  TimeGrid tgrid(param);
  if(!param.is_specified("dt")){
    std::cerr<<"Please set parameter 'dt' explicitly!"<<std::endl;
    throw DummyException();
  }
  FreePropagator prop(param);

  ProcessTensorBuffer PTB(read_PT, true);
  if(PTB.get_n_tot()<1){
    std::cerr<<"PTB.get_n_tot()<1!"<<std::endl;
    throw DummyException();
  }
  if(PTB.get_n_tot()%coarse_grain!=0 &&coarse_grain<1){
    std::cerr<<"PTB.get_n_tot()\%coarse_grain!=0"<<std::endl;
    throw DummyException();
  }
  if(PTB.get(0,ForwardPreload).get_N()!=prop.get_dim()){
    std::cerr<<"PT system dim="<<PTB.get(0,ForwardPreload).get_N()<<" != prop.get_dim()="<<prop.get_dim()<<"!"<<std::endl;
    throw DummyException();
  }

  int n_tot_old=PTB.get_n_tot();
  int n_tot_new=PTB.get_n_tot()/coarse_grain;


  ProcessTensorBuffer PTB2;
  PTB2.set_new_file(write_PT,buffer_blocksize);
  PTB2.resize(n_tot_new);

  for(int n2=0; n2<n_tot_new; n2++){
    ProcessTensorElement &e2=PTB2.get(n2,ForwardPreload);
    for(int c=0; c<coarse_grain; c++){
      int n=n2*coarse_grain+c;
      double t=tgrid.get_t(n);
      double dt=tgrid.get_dt(n);

      ProcessTensorElement tmp=PTB.get(n,ForwardPreload);
      tmp.expand_from_dict();
 
      Eigen::MatrixXcd M1,M2;
      if(false){
//        prop.update(0,t+dt/2.);     M1=prop.M;
//        prop.update(t+dt/2.,dt/2.); M2=prop.M;
      }else{
        prop.update(t,dt/2.);       M1=prop.M;
        prop.update(t+dt/2.,dt/2.); M2=prop.M;
      }
      int NL=M1.rows();

      MPS_Matrix m(NL*NL, tmp.M.dim_d1, tmp.M.dim_d2);
      m.set_zero();
      for(int i1=0; i1<NL; i1++){
        for(int i2=0; i2<NL; i2++){
          for(int i3=0; i3<NL; i3++){
            for(int d1=0; d1<tmp.M.dim_d1; d1++){
              for(int d2=0; d2<tmp.M.dim_d2; d2++){
                m(i3*NL+i1, d1, d2)+=tmp.M(i3*NL+i2, d1, d2)*M1(i2,i1);
              }
            }
          }
        }
      }
      tmp.M.swap(m);

      m.set_zero();
      for(int i1=0; i1<NL; i1++){
        for(int i2=0; i2<NL; i2++){
          for(int i3=0; i3<NL; i3++){
            for(int d1=0; d1<tmp.M.dim_d1; d1++){
              for(int d2=0; d2<tmp.M.dim_d2; d2++){
                m(i3*NL+i1, d1, d2)+=M2(i3,i2)*tmp.M(i2*NL+i1, d1, d2);
              }
            }
          }
        }
      }
      tmp.M.swap(m);

      if(c==0){
        e2=tmp;
      }else{
        m.resize(NL*NL, e2.M.dim_d1, tmp.M.dim_d2);
        m.set_zero();
        for(int i1=0; i1<NL; i1++){
          for(int i2=0; i2<NL; i2++){
            for(int i3=0; i3<NL; i3++){
              for(int d1=0; d1<e2.M.dim_d1; d1++){
                for(int d2=0; d2<tmp.M.dim_d1; d2++){
                  for(int d3=0; d3<tmp.M.dim_d2; d3++){
                    m(i3*NL+i1,d1,d3)+=tmp.M(i3*NL+i2,d2,d3)*e2.M(i2*NL+i1,d1,d2);
                  }
                }
              }
            }
          }
        }          
        e2.M.swap(m);
        e2.closure=tmp.closure;
        e2.env_ops=tmp.env_ops;
        e2.forwardNF=tmp.forwardNF;
        e2.backwardNF=tmp.backwardNF;
      }
    }
  } 
 }catch(DummyException &e){
  return 1;
 }
 return 0;
}
