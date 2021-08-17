#ifndef IF_OD_PLAN_DEFINED_H
#define IF_OD_PLAN_DEFINED_H

#include "Parameters.h"
#include "InitialState.h"
#include "InfluenceFunctional_OD.h"
#include "RankCompressor_Selector.h"
#include "ModePropagatorGenerator.h"
#include "RankCompressor_Selector.h"
#include "Modify_PT.h"

/* 
   Read from Parameter what process tensors are to be calculated.
   Later, the Plan can be executed.
*/

class IF_OD_Plan{
public:
  IF_TimeGrid tgrid;
  bool use_dict;
  double dict_zero;
  int n_coarse;
  int factorization;
  bool compress_trafo_use_ortho;
  int sub_sum_modes;
  Parameters sub_sum_param;

  bool IF_print_timesteps;
  RankCompressor_Ptr compressor;

  int final_sweep_n;
  RankCompressor_Ptr final_sweep_compressor;

  bool use_Gaussian;
  DiagBB Gaussian_DiagBB;

  int dim;
  std::string read_PT, write_PT;
  std::vector<std::string> multi_PT;

  std::vector<Smart_Ptr<ModePropagatorGenerator> > mpgs;
  std::vector<Smart_Ptr<FreePropagator> > PT_apply_prop;

  std::string print_dims_to_file;

  void check_consistency(int check_sys_dim=0)const{
    for(size_t i=0; i<mpgs.size(); i++){
      if(i==0 && check_sys_dim<1){
        check_sys_dim=mpgs[i]->get_N();
      }
      if(mpgs[i]->get_N()!=check_sys_dim){
        std::cerr<<"System part of mode '"<<mpgs[i]->name()<<"' has dimension "<<mpgs[i]->get_N()<<" instead of expected "<<check_sys_dim<<"!"<<std::endl;
        exit(1);
      }
    }

    if(n_coarse>1){
      if(tgrid.n_tot%n_coarse!=0){
        std::cerr<<"IF_OD_Plan: tgrid.n_tot%n_coarse!=0!"<<std::endl;
        exit(1);
      }
    }

    if(use_Gaussian){
      if(read_PT!=""){
        std::cerr<<"'read_PT' cannot be used together with 'use_Gaussian'!"<<std::endl;
        exit(1);
      }
      if(dim!=Gaussian_DiagBB.sys_dim()){
        std::cerr<<"Gaussian_couplings.rows()!=dim!"<<std::endl;
        exit(1);
      }
    }
  }


  void setup(Parameters &param, int check_sys_dim=0){
    //------------------------------------------  
    tgrid.setup(param);
    if(param.get_as_double("print_timegrid_info",true))tgrid.print_info();
    compress_trafo_use_ortho=param.get_as_bool("compress_trafo_use_ortho",false);

    //------------------------------------------  
    compressor=RankCompressor_Selector(param,true);
    use_dict=param.get_as_bool("use_dict",false);
    dict_zero=param.get_as_double("dict_zero", use_dict? 1e-12 : -1.);
    factorization=param.get_as_int("factorization", 0);
   
    sub_sum_modes=param.get_as_int("sub_sum_modes",0);
    sub_sum_param.map.clear();
    sub_sum_param.add_from_prefix("sub_sum",param);

    read_PT=param.get_as_string("read_PT", "");
    write_PT=param.get_as_string("write_PT", "");
    print_dims_to_file=param.get_as_string("print_IF_dims_to_file");
  
    dim=2;
    Eigen::MatrixXcd rho=InitialState(param);
    if(rho.rows()>1){
      dim=rho.rows();
    }

    IF_print_timesteps=param.get_as_bool("IF_print_timesteps",false);
    
    //------------------------------------------  
    n_coarse=param.get_as_size_t("n_coarse",0);

   
    Parameters PT_apply_param; 
    PT_apply_param.add_from_prefix("PT_apply_System_Propagator",param);
    if(!PT_apply_param.is_empty()){
      PT_apply_prop.push_back(new FreePropagator(PT_apply_param));
    }
    if(param.is_specified("PT_apply_System_Hamiltonian")){
      int nrr=param.get_nr_rows("PT_apply_System_Hamiltonian");
      if(PT_apply_prop.size()<1){
        PT_apply_prop.push_back(new FreePropagator());
      }
      for(int i=0; i<nrr; i++){
        Eigen::MatrixXcd H=param.get_as_operator("PT_apply_System_Hamiltonian",
Eigen::MatrixXcd::Zero(1,1), i, 0);
        if(H.rows()>1){
          PT_apply_prop.back()->add_Hamiltonian(H);
        }
      }
    }
    {
    std::vector<std::string> svv=param.get_all_strings("PT_apply_System_Propagator");
     for(size_t i=0; i<svv.size(); i++){
       Parameters param2;
       param2.add_from_file(svv[i]);
       PT_apply_prop.push_back(new FreePropagator(param2));
     }
    }
    //------------------------------------------  

    use_Gaussian=param.get_as_bool("use_Gaussian",false);
    std::string Gaussian_prefix=param.get_as_string("Gaussian_prefix", "Boson");
    if(use_Gaussian){
      Gaussian_DiagBB.setup(param, Gaussian_prefix);

      if(param.is_specified("mem_threshold")){
        double mem_threshold=param.get_as_double("mem_threshold");
        int n_mem_est=Gaussian_DiagBB.estimate_memory_length(tgrid.n_mem, tgrid.dt, mem_threshold, true);
        tgrid.n_mem=n_mem_est;
      }
    }

    //------------------------------------------  

    if(!(use_Gaussian && Gaussian_prefix=="Boson")  && 
       (param.get_as_size_t("Boson_N_modes",0)>0
        || param.get_as_string("Boson_from_table")!="") ){
      mpgs.push_back(new ModePropagatorGenerator_Boson(param));
    }
    if(param.get_as_size_t("MultiSite_N_modes",0)>0
      || param.get_as_string("MultiSite_from_table")!=""){
      mpgs.push_back(new ModePropagatorGenerator_MultiSite(param));
    }
    if(param.get_as_size_t("Potential1D_N_modes",0)>0
      || param.get_as_string("Potential1D_from_table")!=""){
      mpgs.push_back(new ModePropagatorGenerator_Potential1D(param));
    }
    if(param.get_as_size_t("Fermion_N_modes",0)>0
      || param.get_as_string("Fermion_from_table")!=""){
      mpgs.push_back(new ModePropagatorGenerator_Fermion(param));
    }
    if(param.get_as_size_t("RandomSpin_N_modes",0)>0){
      mpgs.push_back(new ModePropagatorGenerator_RandomSpin(param));
    }
    if(param.get_as_size_t("QDPhonon_N_modes",0)>0){
      mpgs.push_back(new ModePropagatorGenerator_QDPhonon(param));
    }
    if(param.is_specified("add_single_mode")){
      std::vector<std::vector<std::string> > lines=param.get("add_single_mode");
      for(size_t r=0; r<lines.size(); r++){
        mpgs.push_back(new ModePropagatorGenerator_SingleMode(lines[r]));
      }
    }
    if(param.is_specified("add_single_mode_from_file")){
      std::vector<std::vector<std::string> > lines=param.get("add_single_mode_from_file");
      for(size_t r=0; r<lines.size(); r++){
        mpgs.push_back(new ModePropagatorGenerator_SingleModeFromFile(lines[r]));
      }
    }
    //------------------------------------------  
    Parameters final_sweep_param; 
    final_sweep_param.add_from_prefix("final_sweep", param);
    if(final_sweep_param.is_empty()){
      final_sweep_n=0;
    }else{
      final_sweep_n=final_sweep_param.get_as_size_t("n",1);
      final_sweep_compressor=RankCompressor_Selector(final_sweep_param, false);
    }

    //------------------------------------------  
    if(param.is_specified("multi_PT")){
      multi_PT=param.get_all_strings("multi_PT");
    }

    //------------------------------------------  
    check_consistency(check_sys_dim);
  }

  //===================================================

  std::vector<Smart_Ptr<InfluenceFunctional_OD> > execute(){
    std::cout<<"Calculating process tensor"<<std::endl;
    check_consistency();
    std::vector<Smart_Ptr<InfluenceFunctional_OD> >IF;
    //------------------------------------------  

    if(read_PT!=""){  
      std::cout<<"Reading process tensor file '"<<read_PT<<"'"<<std::endl;
      IF.push_back(new InfluenceFunctional_OD(read_PT)); 

    }else if(use_Gaussian){
      IF.push_back(new InfluenceFunctional_OD(tgrid, Gaussian_DiagBB, compressor.ref(), dict_zero));

//TODO:    }else if(use_IF_from_SysProp){
//      IF.push_back(IF_from_SysProp(tgrid, prop_IF_from_SysProp));
    }else{
      IF.push_back(new InfluenceFunctional_OD(tgrid, dim));
    }
    //------------------------------------------  

    if(tgrid.use_rep){
      std::cout<<"rep_unit: "<<IF[0]->tgrid.rep_unit<<" n_rep: "<<IF[0]->tgrid.n_rep<<std::endl;
      IF[0]->compress_trafo_use_ortho=compress_trafo_use_ortho;
    }
    IF[0]->print_timesteps=IF_print_timesteps; 
    
    //------------------------------------------ mpgs:
    for(size_t i=0; i<mpgs.size(); i++){
      std::cout<<"Calculating PT for Generator '"<<mpgs[i]->name()<<"'"<<std::endl;
//      IF[0]->add_modes(mpgs[i].ref(), *compressor, dict_zero);

      ModePropagatorGenerator & mpg=mpgs[i].ref();
      for(int k=mpg.first(); k<mpg.get_N_modes(); k=mpg.next(k)){
        std::cout<<"Add mode "<<k<<"/"<<mpg.get_N_modes()<<std::endl;

        ModePropagatorPtr mpp=mpg.getModePropagator(k);
//        InfluenceFunctional_OD IF2(mpp.ref(), tgrid.n_calc*2, tgrid.ta, tgrid.dt/2., dict_zero);
        InfluenceFunctional_OD IF2(tgrid);
        IF2.calculate_dt0_ndt0(mpp.ref(), tgrid.n_calc*2, tgrid.ta, tgrid.dt/2., tgrid.dt0/2., 2);
        IF2.calculate_dict(dict_zero);
        IF2.reduce_to_dict();

     
        //sub_sum:
        int maxdim=IF2.get_max_dim();
//        Parameters param2; 
//        param2.add_to("compress_maxk",Reader::int_to_string(maxdim));
        RankCompressor_Ptr compr2=RankCompressor_Selector(sub_sum_param); 
        if(sub_sum_modes>1){
          IF2.printdim=false;
          for(int count=1; count<sub_sum_modes; count++){ 
//            std::cout<<"Add mode "<<k<<"/"<<mpg.get_N_modes()<<" (inner: "<<count<<"/"<<sub_sum_modes<<", dim: "<<maxdim<<")"<<std::endl;
            k=mpg.next(k); if(!(k<mpg.get_N_modes()))break;
            ModePropagatorPtr mpp2=mpg.getModePropagator(k);
//            InfluenceFunctional_OD IF3(mpp2.ref(), tgrid.n_calc*4, tgrid.ta, tgrid.dt/4., dict_zero);
            InfluenceFunctional_OD IF3(tgrid);
            IF3.calculate_dt0_ndt0(mpp2.ref(), tgrid.n_calc*4, tgrid.ta, tgrid.dt/4., tgrid.dt0/4., 4);
            IF3.calculate_dict(dict_zero, false);
            IF3.reduce_to_dict();

            IF2.add_IF_halfdt(IF3, *compr2);
          }
        }
        IF[0]->add_IF_halfdt(IF2, *compressor);
 
        if(IF[0]->tgrid.rep_unit>0){
          IF[0]->rep.expand_ops(mpp.ref());
        }
      }
    }


    //------------------------------------------  
    //Insert repeated units
    if(tgrid.use_rep){
      if(tgrid.rep_regularize){
        IF[0]->rep.regularize();
      }
#ifdef DEBUG_REP
      std::cout<<"insert: rep.M: "<<IF[0]->rep.M.dim_i<<" "<<IF[0]->rep.M.dim_d1<<" "<<IF[0]->rep.M.dim_d2<<std::endl;
      IF[0]->print_dims(std::cout);
#endif

    //do the actual insertion of the repeated units
      IF[0]->insert_rep();
    }  

    //Possibly multiply part of system propagator and coarse grain time steps
    for(size_t i=0; i<PT_apply_prop.size(); i++){
      Modify_PT::apply_system_propagator(IF[0].ref() , PT_apply_prop[i].ref(), tgrid.ta, tgrid.dt, dict_zero);
    }
    if(n_coarse>1){
      Modify_PT::coarse_grain(IF[0].ref(), n_coarse, dict_zero);
    }

    //------------------------------------------  
    if(IF[0]->size()>0)for(int i=0; i<final_sweep_n; i++){
      std::cout<<"Final sweep: "<<i<<"/"<<final_sweep_n<<":"<<std::endl;
      InfluenceFunctional_OD IF_new(IF[0]->size()*2, IF[0]->get_sys_dim());
      IF[0]->add_IF_halfdt(IF_new, final_sweep_compressor.ref());
    }
    
    //------------------------------------------  
    if(print_dims_to_file!=""){
      std::cout<<"Printing IF dimensions to file '"<<print_dims_to_file<<"'"<<std::endl;
      IF[0]->print_dims(print_dims_to_file);
    }

    //------------------------------------------  
    if(write_PT!=""){
      IF[0]->write_binary(write_PT);
    }

    for(size_t i=0; i<multi_PT.size(); i++){
      std::cout<<"Reading process tensor file ["<<i+1<<"]: '"<<multi_PT[i]<<"'"<<std::endl;
      IF.push_back(new InfluenceFunctional_OD(multi_PT[i]));
    }

    //------------------------------------------  
    return IF;
  }

  IF_OD_Plan(Parameters &param, int check_sys_dim=0){
    setup(param,check_sys_dim);
  }
  IF_OD_Plan(){
    Parameters param;
    setup(param);
  }
};

#endif
