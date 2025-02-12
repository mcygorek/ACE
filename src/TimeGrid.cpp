#include "TimeGrid.hpp"
#include "Parameters.hpp"

namespace ACE{

  bool TimeGrid::n_in_range(int n)const{
    if(n<0)return false;
    if( !rep_infinite && n>=n_tot)return false;
    return true;
  }
  void TimeGrid::complain_if_n_out_of_range(int n)const{
    if(!n_in_range(n)){
      std::cerr<<"n out of range: "<<n<<"/"<<n_tot<<std::endl;
      exit(1);
    }
  }
  double TimeGrid::get_t(int n)const{  
    if(n<=0)return ta;
    if(n<=ndt0)return ta+n*dt0;
    else return ta+ndt0*dt0+(n-ndt0)*dt;  
  }
  double TimeGrid::get_dt(int n)const{  
    if(n<ndt0)return dt0;
    else return dt;  
  }
  std::vector<double> TimeGrid::get_all()const{
    std::vector<double> times(n_tot+1);
    for(int n=0; n<=n_tot; n++){
      times[n]=get_t(n);
    }
    return times;
  }

  int TimeGrid::get_closest_n(double t)const{
    double te=get_t_tot();
    double eps=dt*1e-6;
    if(t<ta-eps)return -1;
    if(t>=te+eps)return -1;
    double t1=t-ta;
    if(t1<dt0/2.)return 0;
    if(t1<(ndt0+0.5)*dt0)return t1/dt0;
    return (t1-ndt0*dt0)/dt+ndt0+0.5;
  }
  std::vector<int> TimeGrid::get_interval_set(double interval_width){
    std::vector<int> list(1,0); 
    for(double t=ta; t<=get_t_tot(); t+=interval_width){
      int n_next=get_closest_n(t);
      if(n_next<0)break;
      if(n_next==list.back())continue;
      list.push_back(n_next);
    }
    return list;
  }

  void TimeGrid::print_info()const{
    std::cout<<"timegrid: ta: "<<ta<<" dt: "<<dt<<" te: "<<get_t(n_tot)<<std::endl;
    std::cout<<"dt0: "<<dt0<<" ndt0: "<<ndt0<<std::endl;
    std::cout<<"n_tot: "<<n_tot<<" n_calc: "<<n_calc<<" n_mem: "<<n_mem<<" n_rep: "<<n_rep<<std::endl;
    std::cout<<"use_rep: "<<(use_rep?"true":"false")<<" ";
    std::cout<<" rep_replace: "<<(rep_replace?"true":"false")<<" ";
    std::cout<<" rep_regularize: "<<(rep_regularize?"true":"false")<<" ";
    std::cout<<" rep_unit: "<<rep_unit<<std::endl;
  }

  void TimeGrid::setup(Parameters &param){
    IF_print_timesteps=param.get_as_bool("IF_print_timesteps",false);

    ta=param.get_as_double("ta", 0);
    dt=param.get_as_double("dt", 1e-2);
    dt0=param.get_as_double("dt0", dt);
    ndt0=param.get_as_int("ndt0",1);
    if(dt0<1e-16){
      std::cerr<<"dt0 has to be larger than 0!"<<std::endl;
      exit(1);
    }
    double te=param.get_as_double("te", 10);
//    int n_max_t=(te-ta)/dt+0.5;
//    int n_max_t=(te-ta-dt0+1e-12*dt)/dt+1;
//    int n_max_t=(te-ta-dt0+1e-6*dt)/dt+1;
    int n_max_t=round((te-ta-ndt0*dt0)/dt)+ndt0;
    if(te-ta<=ndt0*dt0)n_max_t=round((te-ta)/dt0);

    if(param.get_as_string("t_mem")=="auto" ||  param.get_as_string("n_mem")=="auto" || param.get_as_int("n_mem")==-2){
      n_mem=-2;
    }else{
      double t_mem=param.get_as_double("t_mem", te-ta);
      n_mem=param.get_as_double("n_mem", t_mem/dt+0.5);
      if(n_mem>n_max_t)n_mem=n_max_t;
    }
    
    n_tot=n_max_t;
    n_calc=n_tot;

    rep_unit=-1;
    use_rep=false;
    n_rep=0;

    if(param.is_specified("t_rep")){
      double t_rep=param.get_as_double("t_rep"); //calculate until t_rep; construct repetitive unit and insert until n_tot
      n_rep=n_tot-( t_rep/dt+0.5);
      if(n_rep>0){  
        use_rep=true;
        n_calc=n_tot-n_rep;
        rep_unit=n_calc/2;
      }
    }

    rep_replace=param.get_as_bool("rep_replace", false);
    rep_regularize=param.get_as_bool("rep_regularize", false);

  }
  void TimeGrid::half_dt(int fac){
    dt/=fac;
    dt0/=fac;
    ndt0*=fac;
    n_tot*=fac;
    if(n_mem>-2)n_mem*=fac;
    n_calc*=fac;
    n_rep*=fac;
    if(rep_unit>=0)rep_unit*=fac;
  }
  TimeGrid TimeGrid::construct_half_dt(int fac)const{
    TimeGrid tgr(*this);
    tgr.half_dt(fac);
    return tgr;
  }
  TimeGrid TimeGrid::get_no_rep()const{
    TimeGrid tgr(*this);
    if(use_rep || rep_unit>0){
      tgr.rep_unit=-1;
      tgr.use_rep=false;
      tgr.n_calc=n_rep;
      tgr.n_mem=n_rep;
      tgr.n_tot=n_rep;
      tgr.n_rep=0;
    }
    return tgr;
  } 
  void TimeGrid::set_default(int nmax_, double dt_, double ta_){
    Parameters param;
    param.add_to("ta",ta_);
    param.add_to("dt",dt_);
    param.add_to("te",nmax_*dt_);
    setup(param);
  }
  void TimeGrid::setup_coarse(Parameters &param){
    double old_dt=param.get_as_double("dt",1e-2);
    int n_coarse=param.get_as_size_t("n_coarse",0);
    if(n_coarse<2){
      setup(param);
      return;
    }
    Parameters param2=param;
    param2.override_param("dt",old_dt*n_coarse);
    param2.override_param("dt0",param.get_as_double("dt0",old_dt)+(n_coarse-1)*old_dt);
    if(param.is_specified("n_mem")){
      param2.override_param("n_mem",param.get_as_double("n_mem")/n_coarse);
    }
    setup(param2);
  }

}//namespace
