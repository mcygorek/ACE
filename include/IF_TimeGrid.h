#ifndef IF_TIMEGRID_DEFINED_H
#define IF_TIMEGRID_DEFINED_H
#include "Parameters.h"

/* Structure to store the parametrization of the time grid as well
   as different parameters determining how exactly the influence functionals
   are to be calculated */

struct IF_TimeGrid{
  double ta;  //initial time
  double dt;  //time step width
  
  int n_tot;  //total length of MPO
  int n_calc; //legth of MPO to be calculated (n_tot may be padded by some repetitive element)
  int n_mem;  //memory truncation (used to speed up calculation)
 
  bool use_rep; //turns on the use of prepetitive units
  bool rep_replace; //replace all chain elements by repetitive unit
  bool rep_regularize; //try to make repetitive element more "physical"
  bool rep_infinite; //repeat units indefinitely
  int n_rep;  //how many repetetive units are to be included
  int rep_unit; //at which point int the MPO chain is the repetitive unit(s) to be inserted

  bool IF_print_timesteps; //print time step during calculation of IF



  bool n_in_range(int n)const{
    if(n<0)return false;
    if( !rep_infinite && n>=n_tot)return false;
    return true;
  }
  void complain_if_n_out_of_range(int n)const{
    if(!n_in_range(n)){
      std::cerr<<"n out of range: "<<n<<"/"<<n_tot<<std::endl;
      exit(1);
    }
  }
  double get_t(int n)const{  return ta+n*dt;  }
  double get_t_tot()const{ return get_t(n_tot); }

  void print_info()const{
    std::cout<<"timegrid: ta: "<<ta<<" dt: "<<dt<<" te: "<<get_t(n_tot)<<std::endl;
    std::cout<<"n_tot: "<<n_tot<<" n_calc: "<<n_calc<<" n_mem: "<<n_mem<<" n_rep: "<<n_rep<<std::endl;
    std::cout<<"use_rep: "<<(use_rep?"true":"false")<<" ";
    std::cout<<" rep_replace: "<<(rep_replace?"true":"false")<<" ";
    std::cout<<" rep_regularize: "<<(rep_regularize?"true":"false")<<" ";
    std::cout<<" rep_unit: "<<rep_unit<<std::endl;
  }
  void setup(Parameters &param){
    IF_print_timesteps=param.get_as_bool("IF_print_timesteps",false);

    ta=param.get_as_double("ta", 0);
    dt=param.get_as_double("dt", 1e-2);
    double te=param.get_as_double("te", 10);
    int n_max_t=(te-ta)/dt+0.5;

    double t_mem=param.get_as_double("t_mem", te-ta);
    n_mem=param.get_as_double("n_mem", t_mem/dt+0.5);
    if(n_mem>n_max_t)n_mem=n_max_t;
    
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
  void set_default(int nmax=0){
    Parameters param;
    param.add_to("ta",0);
    param.add_to("dt",1);
    param.add_to("te",nmax);
    setup(param);
  }
  IF_TimeGrid(Parameters &param){
    setup(param);
  }
  IF_TimeGrid(int nmax=0){
    set_default(nmax);
  }
};


#endif
