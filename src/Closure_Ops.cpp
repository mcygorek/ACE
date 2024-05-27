#include "Closure_Ops.hpp"
#include "otimes.hpp"
#include "Parameters.hpp"
#include "Reader.hpp"

namespace ACE{

  void Closure_Ops::check_ML(int ML)const{
    for(size_t i=0; i<ops.size(); i++){
      if(ops[i].rows()!=ML){ 
        std::cerr<<"Closure_Ops: ops[i].rows()!=ML (i="<<i<<" ops[i].rows()="<<ops[i].rows()<<" ML="<<ML<<")!"<<std::endl;
        exit(1);
      }
    }
  }

  void Closure_Ops::set_ops_from_matrices(const std::vector<Eigen::MatrixXcd> &mats){
    ops.clear();
    if(mats.size()<1)return;

    int M=mats[0].rows();
    ops.resize(mats.size(),Eigen::VectorXcd::Zero(M*M));

    for(size_t i=0; i<mats.size(); i++){
      if(i>0){
        if(mats[i].rows()!=M){
          std::cerr<<"Closure_Ops: mats[i].rows()!=M!"<<std::endl; 
          exit(1);
        }
      } 
      for(int m1=0; m1<M; m1++){
        for(int m2=0; m2<M; m2++){
          ops[i](m1*M+m2)=mats[i](m2,m1);
        }
      }
    }
  }

  void Closure_Ops::print_info(std::ostream &ofs)const{
    ofs<<"Closure_Ops: ";//<<std::endl;
    if(use)ofs<<"use=true "; else ofs<<"use=false ";
    if(use_env_ops){
      ofs<<"use_env_ops=true env_ops_nr.size()="<<env_ops_nr.size();
    }else{
      ofs<<"use_env_ops=false ops.size()="<<ops.size();
    } 
    ofs<<" H_scale="<<H_scale;
    ofs<<std::endl;
  }

  void Closure_Ops::setup(Parameters &param){
    use=param.get_as_bool("use_env_closure",false);
    //if none are explicitly specified: use get_env_ops() of MPGs. 

    std::vector<std::string> all=param.get_all_strings("add_Closure_Op_explicit");
    if(all.size()>0){
      use=true;
      use_env_ops=false;
      std::vector<Eigen::MatrixXcd> mats(all.size());
      for(size_t i=0; i<all.size(); i++){
        mats[i]=ReadExpression(all[i]);
      } 
      set_ops_from_matrices(mats);
    }else{
      use_env_ops=true;
    }

    Parameters_Entry nrs=param.get("add_Closure_Op");
    env_ops_nr.clear();
    if(nrs.size()>0){
      use=true;
    }
    for(size_t i=0; i<nrs.size(); i++){
      double scale=1.;
      if(nrs[i].size()>1){
        scale=readDouble(nrs[i][1], "add_Closure_Op: scale");
      }

      env_ops_nr.push_back(
        std::make_pair(readSizeT(nrs[i][0], "add_Closure_Op"),scale));
    }
    
    H_scale=param.get_as_double("add_Closure_H_scale",0);
    if(H_scale>0){
      use=true;
    }
  }

  void Closure_Ops::set_trivial(){
    use=false;
    use_env_ops=false;
    H_scale=0;
    ops.clear();
  }
  Closure_Ops::Closure_Ops(){
    Parameters param;
    setup(param);
  }


}//namespace
