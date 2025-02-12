#pragma once
#ifndef IF_TIMEGRID_DEFINED_H
#define IF_TIMEGRID_DEFINED_H
//#include "Parameters.hpp"

#include <vector>

/* Structure to store the parametrization of the time grid as well
   as different parameters determining how exactly the influence functionals
   are to be calculated */

namespace ACE{
class Parameters;

struct TimeGrid{
  double ta;   //initial time
  double dt;   //time step width
  double dt0;  //first time step width
  int ndt0;    //number of times first time step is applied (required for symmetric trotter formulas)
  
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

  bool n_in_range(int n)const;
  
  void complain_if_n_out_of_range(int n)const;

  double get_t(int n)const;
  
  double get_dt(int n)const;
  
  inline double get_t_tot()const{ return get_t(n_tot); }

  std::vector<double> get_all()const;

  int get_closest_n(double t)const;
  
  std::vector<int> get_interval_set(double interval_width);

  void print_info()const;
  
  void setup(Parameters &param);

  void half_dt(int fac=2);
  
  TimeGrid construct_half_dt(int fac=2)const;
  
  TimeGrid get_no_rep()const;
   
  void set_default(int nmax=0, double dt=0.01, double ta=0.);
  
  void setup_coarse(Parameters &param);
  
  inline TimeGrid(Parameters &param){
    setup(param);
  }
  inline TimeGrid(int nmax=0, double dt=0.01, double ta=0.){
    set_default(nmax, dt, ta);
  }
};

}//namespace
#endif
