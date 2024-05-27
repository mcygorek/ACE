#include "IF_OD_Plan.hpp"
#include "InitialState.hpp"
#include "InfluenceFunctional_OD.hpp"
#include "ProcessTensor_real.hpp"
#include "RankCompressor_Selector.hpp"
#include "MPG_Selector.hpp"
#include "RankCompressor_Selector.hpp"
#include "Modify_PT.hpp"
#include "ModePropagatorGeneratorList.hpp"

/* 
   Read from Parameter what process tensors are to be calculated.
   Later, the Plan can be executed.
*/
namespace ACE{

  void IF_OD_Plan::check_consistency(int check_sys_dim)const{
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
      if(!!read_PT){
        std::cerr<<"'read_PT' cannot be used together with 'use_Gaussian'!"<<std::endl;
        exit(1);
      }
      if(dim!=Gaussian_DiagBB.sys_dim()){
        std::cerr<<"Gaussian_couplings.rows()!=dim!"<<std::endl;
        exit(1);
      }
    }
  }

  void IF_OD_Plan::setup(Parameters &param, int check_sys_dim){
    //------------------------------------------  
    tgrid.setup(param);
    if(param.get_as_double("print_timegrid_info",true))tgrid.print_info();
    compress_trafo_use_ortho=param.get_as_bool("compress_trafo_use_ortho",false);

    //------------------------------------------  
    compressor=RankCompressor_Selector(param,true);
    use_dict=param.get_as_bool("use_dict",false);
    dict_zero=param.get_as_double("dict_zero", use_dict? 1e-12 : -1.);
    factorization=param.get_as_int("factorization", 0);
 
  
    read_PT.fname=param.get_as_string("read_PT", "");
    read_PT.expand_front=param.get_as_size_t("read_PT_expand_front",0);
    read_PT.expand_back=param.get_as_size_t("read_PT_expand_back",0);

    if(param.is_specified("read_PT_expand")){
      if(param.is_specified("read_PT")){
        std::cerr<<"Please don't specify both 'read_PT' and 'read_PT_expand' at the same time!"<<std::endl;
        exit(1); 
      }
      Parameters_Entry entry=param.get("read_PT_expand");
      if(entry[0].size()<3){
        std::cerr<<"'read_PT_expand' needs 3 parameters: FILENAME EXPAND_FRONT EXPAND_BACK!"<<std::endl;
        exit(1); 
      }
      read_PT=ReadPT_struct(entry[0][0], 
         readDouble(entry[0][1], "read_PT_expand: EXPAND_FRONT"),
         readDouble(entry[0][2], "read_PT_expand: EXPAND_BACK"));
    }

    write_PT=param.get_as_string("write_PT", "");
    print_dims_to_file=param.get_as_string("print_IF_dims_to_file");
  
    dim=2;
    Eigen::MatrixXcd rho=InitialState(param);
    if(rho.rows()>1){
      dim=rho.rows();
    }

    IF_print_timesteps=param.get_as_bool("IF_print_timesteps",false);

    intermediate_sweeps=param.get_as_size_t("intermediate_sweeps",0);  
    {Parameters param2; param2.add_from_prefix("intermediate_sweeps",param);
     intermediate_sweeps_compressor.setup(param2);}

    final_sweeps=param.get_as_size_t("final_sweeps",0);  
    {Parameters param2; param2.add_from_prefix("final_sweeps",param);
     final_sweeps_compressor.setup(param2);}

    //------------------------------------------  
    n_coarse=param.get_as_size_t("n_coarse",0);

   
    Parameters PT_apply_param; 
    PT_apply_param.add_from_prefix("PT_apply_System_Propagator",param);
    if(!PT_apply_param.is_empty()){
      PT_apply_prop.push_back(std::make_shared<FreePropagator>(PT_apply_param));
    }
    if(param.is_specified("PT_apply_System_Hamiltonian")){
      int nrr=param.get_nr_rows("PT_apply_System_Hamiltonian");
      if(PT_apply_prop.size()<1){
        PT_apply_prop.push_back(std::make_shared<FreePropagator>());
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
       PT_apply_prop.push_back(std::make_shared<FreePropagator>(param2));
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
 
      double t_extra=param.get_as_double("t_extra",0);
      n_extra=(t_extra/tgrid.dt+0.5);
    }
    //------------------------------------------  
    use_realPT=param.get_as_bool("use_realPT", false);
    if(use_realPT){
      compressor_real.setup(param);
    }
    closure_ops.setup(param);
    env_state_filter=std::make_shared<Env_State_Filter>(param);

    //------------------------------------------  

    Parameters_Entry param_trafo_chain=param.get("print_trafo_chain");
    for(size_t i=0; i<param_trafo_chain.size(); i++){
      if(param_trafo_chain[i].size()<2){
        std::cerr<<"Usage: print_trafo_chain TIME FILENAME!"<<std::endl;
        exit(1);
      }
      double t=readDouble(param_trafo_chain[i][0],"print_trafo_chain TIME");
      int n=tgrid.get_closest_n(t);
      if(n<1 || n>=tgrid.n_calc-1){
        std::cerr<<"print_trafo_chain: n<1 || n>=tgrid.n_calc-1 ("<<n<<" vs. "<<tgrid.n_calc<<")!"<<std::endl;
        exit(1);
      }
      std::string fname=param_trafo_chain[i][1];
      std::cout<<"Will print trafo_chain at time "<<t<<" (step "<<n<<") to file '"<<fname<<"'."<<std::endl;
      
      trafo_chain_files.push_back(std::make_pair(n, fname));
    }
    //------------------------------------------  
   
    mpgs=MPG_Selector(param);

    //------------------------------------------  
    Parameters final_sweep_param;    //NOTE: this is for the complex IF_OD, not the ProcessTensor_real! (cf. "final_sweeps" with a trailing "s")
    final_sweep_param.add_from_prefix("final_sweep", param);
    if(final_sweep_param.is_empty()){
      final_sweep_n=0;
    }else{
      final_sweep_n=final_sweep_param.get_as_size_t("n",1);
      final_sweep_compressor=RankCompressor_Selector(final_sweep_param, false);
    }

    //------------------------------------------  
    if(param.is_specified("multi_PT")){
      std::vector<std::string> svec=param.get_all_strings("multi_PT");
      for(size_t i=0; i<svec.size(); i++){
        multi_PT.push_back(ReadPT_struct(svec[i]));
      }
    }
    if(param.is_specified("multi_PT_expand")){
      Parameters_Entry entry=param.get("multi_PT_expand");
      for(size_t i=0; i<entry.size(); i++){
        if(entry[i].size()<3){
          std::cerr<<"'multi_PT_expand' needs 3 parameters: FILENAME EXPAND_FRONT EXPAND_BACK!"<<std::endl;
          exit(1); 
        }
        multi_PT.push_back(ReadPT_struct(entry[i][0], 
             readDouble(entry[i][1], "multi_PT_expand: EXPAND_FRONT"),
             readDouble(entry[i][2], "multi_PT_expand: EXPAND_BACK")));
      }
    }

    //------------------------------------------  
    check_consistency(check_sys_dim);
//std::cout<<"TEST:: IF_OD_PLAN setup finished."<<std::endl;
  }

  //===================================================


  std::vector<std::shared_ptr<InfluenceFunctional_OD> > IF_OD_Plan::execute(){
    std::cout<<"Calculating process tensor"<<std::endl;
    check_consistency();
    std::vector<std::shared_ptr<InfluenceFunctional_OD> >IF;
    //------------------------------------------  

    if(!!read_PT && !use_realPT){  
      std::cout<<"Reading process tensor file '"<<read_PT.fname<<"'"<<std::endl;
      IF.push_back(std::make_shared<InfluenceFunctional_OD>(read_PT.fname)); 

      IF[0]->dict.expand(read_PT);

    }else if(use_Gaussian){
      IF.push_back(std::make_shared<InfluenceFunctional_OD>(tgrid, Gaussian_DiagBB, compressor.ref(), dict_zero, n_extra));

    }else{
      IF.push_back(std::make_shared<InfluenceFunctional_OD>(tgrid, dim));
    }
    //------------------------------------------  

    if(tgrid.use_rep){
      std::cout<<"rep_unit: "<<IF[0]->tgrid.rep_unit<<" n_rep: "<<IF[0]->tgrid.n_rep<<std::endl;
      IF[0]->compress_trafo_use_ortho=compress_trafo_use_ortho;
    }
    IF[0]->print_timesteps=IF_print_timesteps; 
    
    //------------------------------------------ mpgs:
    if(use_realPT && (mpgs.size()>0 || !!read_PT) ){ 
 
      ProcessTensor_real PT;
      if(!!read_PT){ 
        PT.read_binary(read_PT.fname);

        PT.dict.expand(read_PT);
      }else{
        PT.initialize(mpgs[0]->get_N());
        PT.set_trivial(mpgs[0]->get_N(), tgrid.n_calc, dict_zero);
        PT.closure_ops=closure_ops;
      }
      PT.print_timesteps=IF_print_timesteps;
      PT.env_state_filter=env_state_filter;

//Pass Trafo_Chain info to PT
      for(size_t i=0; i<trafo_chain_files.size(); i++){
        PT.trafo_chains.push_back(std::make_pair(Trafo_Chain(), trafo_chain_files[i].first));
      }


//MPG loop
      bool printdim=true;
      TimeGrid tgrid2=tgrid.construct_half_dt();
      for(size_t i=0; i<mpgs.size(); i++){
        std::cout<<"Calculating PT for Generator '"<<mpgs[i]->name()<<"' (realPT)"<<std::endl;
 
        ModePropagatorGenerator &mpg=*mpgs[i].get();
        for(int k=mpg.first(); k<mpg.get_N_modes(); k=mpg.next(k)){
          std::cout<<"Mode "<<k<<"/"<<mpg.get_N_modes()<<std::endl;

          ModePropagatorPtr mpp=mpg.getModePropagator(k);

          ProcessTensor_real PT2(*mpp.get(), tgrid2, closure_ops);
          PT2.calculate_dict(dict_zero);
          PT2.reduce_to_dict();
          std::cout<<"new mode: dict: ";PT2.dict.print_beta(); std::cout<<std::endl;
          if(printdim)std::cout<<"Maxdim: before combination: "<<PT.get_max_dim()<<std::endl;
          if(printdim)std::cout<<"Maxdim: new contribution: "<<PT2.get_max_dim()<<std::endl;

          PT.join_and_sweep_halfdt(PT2, compressor_real);
          if(printdim)std::cout<<"Maxdim: after sweep: "<<PT.get_max_dim()<<std::endl;

          std::cout<<"dict: ";PT.dict.print_beta(); std::cout<<std::endl;
           
          for(int i=0; i<intermediate_sweeps; i++){
            std::cout<<"Performing intermediate sweep "<<i<<"/"<<intermediate_sweeps<<std::endl;
            if(intermediate_sweeps_compressor.has_effect()){
              if(i==intermediate_sweeps-1){
                PT.sweep_low_high_low(intermediate_sweeps_compressor);
              }else{
                PT.sweep_low_high_low(compressor_real);
              }
            }else{
              PT.sweep_low_high_low(compressor_real);
            }
            
            for(size_t i=0; i<PT.trafo_chains.size(); i++){
              Trafo_Chain & tc = PT.trafo_chains[i].first;
              int sz=tc.T.size();
              if(tc.T.size()<2){ std::cerr<<"intermediate sweep: tc.T.size()<2!"<<std::endl; exit(1);}
              if(tc.T.size()!=tc.Tinv.size()){ std::cerr<<"intermediate sweep: tc.T.size()!=tc.Tinv.size()!"<<std::endl; exit(1);}
              if(tc.T[sz-1].cols()!=tc.T[sz-2].rows()){std::cerr<<"intermediate sweep: tc.T[sz-1].cols()!=tc.T[sz-2].rows()!"<<std::endl; exit(1);}
              if(tc.Tinv[sz-1].rows()!=tc.Tinv[sz-2].cols()){std::cerr<<"intermediate sweep: tc.Tinv[sz-1].rows()!=tc.Tinv[sz-2].cols()!"<<std::endl; exit(1);}

              tc.T[sz-2]=tc.T[sz-1]*tc.T[sz-2];
              tc.T.pop_back();
              
              tc.Tinv[sz-2]=tc.Tinv[sz-2]*tc.Tinv[sz-1];
              tc.Tinv.pop_back();
            }
          }
        }
      }

      for(int i=0; i<final_sweeps; i++){
        std::cout<<"Performing final sweep "<<i<<"/"<<final_sweeps<<std::endl;
        if(final_sweeps_compressor.has_effect() && i==final_sweeps-1){
          PT.sweep_low_high_low(final_sweeps_compressor);
        }else{
          PT.sweep_low_high_low(compressor_real);
        }
        
        for(size_t i=0; i<PT.trafo_chains.size(); i++){
          Trafo_Chain & tc = PT.trafo_chains[i].first;
          int sz=tc.T.size();
          if(tc.T.size()<2){ std::cerr<<"final sweep: tc.T.size()<2!"<<std::endl; exit(1);}
          if(tc.T.size()!=tc.Tinv.size()){ std::cerr<<"final sweep: tc.T.size()!=tc.Tinv.size()!"<<std::endl; exit(1);}
          if(tc.T[sz-1].cols()!=tc.T[sz-2].rows()){std::cerr<<"final sweep: tc.T[sz-1].cols()!=tc.T[sz-2].rows()!"<<std::endl; exit(1);}
          if(tc.Tinv[sz-1].rows()!=tc.Tinv[sz-2].cols()){std::cerr<<"final sweep: tc.Tinv[sz-1].rows()!=tc.Tinv[sz-2].cols()!"<<std::endl; exit(1);}

          tc.T[sz-2]=tc.T[sz-1]*tc.T[sz-2];
          tc.T.pop_back();
          
          tc.Tinv[sz-2]=tc.Tinv[sz-2]*tc.Tinv[sz-1];
          tc.Tinv.pop_back();
        }
      }


//Print Trafo_Chain
      for(size_t i=0; i<PT.trafo_chains.size(); i++){
        PT.trafo_chains[i].first.print_info();
//        PT.trafo_chains[i].first.compress(1e-20);
//        PT.trafo_chains[i].first.print_info();
        PT.trafo_chains[i].first.write(trafo_chain_files[i].second);
      }

#ifdef FINAL_ROTATE_IDENTITY_FIRST
      if(PT.env_state_filter.get_count()!=0 && PT.env_state_filter->mean_field){
        for(size_t n=0; n<PT.a.size(); n++){
          std::cout<<"called: apply: n="<<n<<std::endl;
          PT.env_state_filter->rotate_identity_to_first(PT.a, PT.env_ops, n);
        }
      }
#endif

      if(write_PT!=""){
        PT.write_binary(write_PT);
      }
      IF[0]=std::make_shared<InfluenceFunctional_OD>( PT.get_IF_OD() );
    }else{
      for(size_t i=0; i<mpgs.size(); i++){
        std::cout<<"Calculating PT for Generator '"<<mpgs[i]->name()<<"'"<<std::endl;
        IF[0]->add_modes(*mpgs[i].get(), *compressor, dict_zero);
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
      Modify_PT::apply_system_propagator(*IF[0].get() , *PT_apply_prop[i].get(), tgrid.ta, tgrid.dt, dict_zero);
    }
    if(n_coarse>1){
      Modify_PT::coarse_grain(*IF[0].get(), n_coarse, dict_zero);
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
    if(write_PT!="" && !use_realPT){
      IF[0]->write_binary(write_PT);
    }

    for(size_t i=0; i<multi_PT.size(); i++){
      std::cout<<"Reading process tensor file ["<<i+1<<"]: '"<<multi_PT[i].fname<<"'"<<std::endl;
      IF.push_back(std::make_shared<InfluenceFunctional_OD>(multi_PT[i].fname));
      IF.back()->dict.expand(multi_PT[i]);
    }

    //------------------------------------------  
    return IF;
  }

  IF_OD_Plan::IF_OD_Plan(Parameters &param, int check_sys_dim){
    setup(param,check_sys_dim);
  }
  IF_OD_Plan::IF_OD_Plan(){
    Parameters param;
    setup(param);
  }

}//namespace
