#include "CCTEMPO_IF_Nodes.hpp"
#include "CCTEMPO.hpp"
#include "Simulation_CTEMPO.hpp"

namespace ACE{

void CCTEMPO_IF_Nodes::setup(Parameters &param){
  TimeGrid tgrid(param);
  std::vector<int> dims=CCTEMPO::get_dimensions(param);
  int N_Sites=dims.size();
  list.resize(N_Sites);

  // set up b-nodes of CTEMPO networks
  for(int r=0; r<N_Sites; r++){
    int N=dims[r];
    std::string prefix="S"+int_to_string(r)+"_Boson";
    Parameters param2; param2.add_from_prefix(prefix, param);
    if(param2.map.size()<1 || param.get_as_bool("S"+int_to_string(r)+"_use_PT")){ continue; }
   
    Parameters param3=param;
    int n_mem=tgrid.n_mem;
    if(param3.is_specified("S"+int_to_string(r)+"_t_mem") || 
       param3.is_specified("S"+int_to_string(r)+"_n_mem") ){
      int n_mem3=param3.get_as_size_t("S"+int_to_string(r)+"_n_mem",
                 param3.get_as_double("S"+int_to_string(r)+"_t_mem")/tgrid.dt);
      if(n_mem3<tgrid.n_tot){
        n_mem=n_mem3;
        param3.override_param("n_mem", n_mem3);
      }
    }
//    std::cout<<"r="<<r<<" n_mem="<<n_mem<<std::endl;
    list[r].diagBB.setup(param3, prefix);
    list[r].b=Simulation_CTEMPO::initialize_b(n_mem, list[r].diagBB, tgrid.dt);
  }

  //PT-MPO:
  for(int r=0; r<N_Sites; r++){
    std::string prefix="S"+int_to_string(r);
    Parameters paramS; paramS.add_from_prefix(prefix, param);
    //if Sx_Boson is set but not Sx_use_PT, CTEMPO is used instead of PTs
    if(!paramS.get_as_bool("use_PT")){ 
      paramS.erase_with_prefix("Boson");
    } 
    //Need: copy TimeGrid info:
    std::vector<std::string> keys={"ta","te","dt","t_mem","n_mem","threshold"};
    for(size_t i=0; i<keys.size(); i++){
      if(param.is_specified(keys[i])){
        paramS.add_if_not_specified(keys[i],param.get_as_string(keys[i]));
      }
    }
    std::cout<<"Parameters for PT["<<r<<"]:"<<std::endl; paramS.print();
    list[r].PT=std::make_shared<ProcessTensorForwardList>(paramS, dims[r]);
  }

  //Process ...env_same_as_S...
  for(int r=0; r<N_Sites; r++){
    std::string key="S"+int_to_string(r)+"_env_sameas_S";
    if(param.is_specified(key)){
      int r2=param.get_as_size_t(key);
      if(r2<0||r2>=N_Sites){
        std::cerr<<key<<": r2="<<r2<<" out of bounds!"<<std::endl;
        throw DummyException();
      }
      list[r].b=list[r2].b;
    }
  }
}


}//namespace
