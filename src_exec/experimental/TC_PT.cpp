//#include "Trafo_Chain.hpp"
//#include "Parameters.hpp"
//#include "InfluenceFunctional_Repeat.hpp"
//#include "IF_OD_Plan.hpp"
#include "ACE.hpp"

using namespace ACE;

int main(int args, char ** argv){
  Parameters param(args, argv, true);
  std::string print_param=param.get_as_string("print_param");
  if(print_param!="")param.print(print_param);
 
  double dict_zero=param.get_as_double("dict_zero");
  std::vector<std::string> infile=param.get_all_strings("TC");
  if(infile.size()<1){
    std::cerr<<"'TC' needs at least one arguments!"<<std::endl;
    exit(1);
  }

  std::string write_PT=param.get_as_string("write_PT");
  std::string write_repeatPT=param.get_as_string("write_repeatPT");
  if( write_PT=="" && write_repeatPT==""){
    std::cerr<<"write_PT=="" && write_repeatPT==""!"<<std::endl;
    exit(1);
  }
  bool printdim=param.get_as_bool("printdim",true);
  
  Trafo_Chain tc(infile[0]);
  std::cout<<"TC[0]: ";
  tc.print_info();
  if(tc.size()<1){
    std::cerr<<"TC[0].size()<1!"<<std::endl;
    exit(1);
  }
   
  double compress_first=param.get_as_double("compress_first", 0);
  if(compress_first>0){
    std::cout<<"compressing first TC with epsilon="<<compress_first<<std::endl;
    tc.compress(compress_first);
    tc.print_info();
  }
  
  double threshold=param.get_as_double("threshold",0);
  for(size_t i=1; i<infile.size(); i++){
    Trafo_Chain tc2(infile[i]);
    if(compress_first>0)tc2.compress(compress_first);
    std::cout<<"TC["<<i<<"]: ";
    tc2.print_info();

    tc.add_ortho(tc2, threshold);
    std::cout<<"added i="<<i<<std::endl;
    std::cout<<"new overlap: "<<tc.overlap(tc2)<<" of "<<tc.lastdim()<<" or "<<tc2.lastdim()<<std::endl;
    std::cout<<"TC[0]: ";
    tc.print_info();
  }

  std::cout<<"TC preparation completed"<<std::endl;

  TimeGrid tgrid(param);
  IF_OD_Plan plan(param);
  
  if(plan.mpgs.size()<1){
    std::cerr<<"No ModePropagatorGenerator specified!"<<std::endl;
  }


  ModePropagatorGenerator & mpg=*plan.mpgs[0].get();
  if(tc.size()!=mpg.get_N_modes()){
    std::cerr<<"tc.size()!=mpg.get_N_modes() ("<<tc.size()<<" vs. "<<mpg.get_N_modes()<<")!"<<std::endl;
    exit(1);
  }
  std::cout<<"Using ModePropagatorGenerator: '"<<mpg.name()<<"' with N_modes=";
  std::cout<<mpg.get_N_modes()<<std::endl;

  int N=mpg.get_N();
  int NL=N*N;
  

  TimeGrid tgrid_dummy;
  {Parameters param2=param; param2.override_param("te",tgrid.dt*3);tgrid_dummy.setup(param2);}

  //PT has 3 steps: first, 'typical', and last
  ProcessTensor_real PT(N, tgrid_dummy.n_calc, dict_zero);

  TimeGrid tgrid2=tgrid_dummy.construct_half_dt();
  std::shared_ptr<RankCompressor_real> compressor_dummy=std::make_shared<RankCompressor_None_real>();
  for(int k=0; k<mpg.get_N_modes(); k++){
    std::cout<<"Mode "<<k<<"/"<<mpg.get_N_modes()<<std::endl;
    ModePropagatorPtr mpp=mpg.get_ModePropagator(k);

    ProcessTensor_real PT2(*mpp.get(), tgrid2); //, closure_ops);
    PT2.calculate_dict(dict_zero);
    PT2.reduce_to_dict();
    std::cout<<"new mode: dict: ";PT2.dict.print_beta(); std::cout<<std::endl;
    if(printdim)std::cout<<"Maxdim: before combination: "<<PT.get_max_dim()<<std::endl;
    if(printdim)std::cout<<"Maxdim: new contribution: "<<PT2.get_max_dim()<<std::endl;

    PT.join_and_sweep_halfdt(PT2, *compressor_dummy.get(),false);

    if(PT.a[1].dim_d1!=tc.T[k].cols()){
      std::cerr<<"PT.a[1].dim_d1!=tc.T["<<k<<"].cols()!"<<std::endl;
      exit(1);
    }
    if(PT.a[1].dim_d2!=tc.Tinv[k].rows()){
      std::cerr<<"PT.a[1].dim_d2!=tc.Tinv["<<k<<"].rows()!"<<std::endl;
      exit(1);
    }
   
    PT.trafo_d2(tc.Tinv[k],tc.T[k],0);
    PT.trafo_d1(tc.T[k],1);
    PT.trafo_d2(tc.Tinv[k],tc.T[k],1);
    PT.trafo_d1(tc.T[k],2);

    if(printdim)std::cout<<"Maxdim: after sweep: "<<PT.get_max_dim()<<std::endl;
    std::cout<<"dict: ";PT.dict.print_beta(); std::cout<<std::endl;
  }
 
  if(write_repeatPT!=""){
    InfluenceFunctional_OD IF=PT.get_IF_OD();
    InfluenceFunctional_Repeat IFR;
    IFR.dict=IF.dict;
    IFR.a=IF.a;
    IFR.c=IF.c;
    IFR.env_ops=IF.env_ops;
    IFR.write_binary(write_repeatPT);
  }

  if(write_PT!=""){
    PT.make_longer(1,tgrid.n_tot-3);
//    PT.make_longer(1,7);
//    PT.calculate_closures();
//    std::cout<<"printdims: "<<std::endl; PT.print_dims(); 
    PT.write_binary(write_PT);
  }

  return 0;
}
