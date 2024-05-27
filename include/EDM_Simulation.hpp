#ifndef ACE_EDM_SIMULATION_DEFINED_H
#define ACE_EDM_SIMULATION_DEFINED_H

#include "EDM_Initial_State.hpp"
#include "EDM_Filter.hpp"
#include "EDM_OutputPrinter.hpp"
#include "EDM_OneBodyTerm_FreePropagator.hpp"
#include "EDM_TwoBodyTerm.hpp"
//#include "EDM_OneBodyObservable.hpp"
#include "Eigen_fwd.hpp"
#include "TimeGrid.hpp"
#include "Parameters.hpp"
#include <memory>

namespace ACE{

class EDM_Simulation{
public:
  std::vector<std::pair<std::shared_ptr<EDM_OneBodyTerm>, int> >  oneBodyTerms;
  std::vector<std::pair<std::shared_ptr<EDM_TwoBodyTerm>, std::pair<int,int> > >  twoBodyTerms;
 
  std::vector<std::pair<Eigen::VectorXcd, int> > oneBodyObservables;
  EDM_OutputPrinter printer;
  EDM_OutputPrinter debug_printer;
  bool print_timesteps;

  bool use_symmetric_Trotter;
  bool propagate_alternate;

  //function
  void extract_Observables(double t, const EDM_State &state);

  void apply_OneBodyTerms(int n, const TimeGrid &tgrid, EDM_State &state, const EDM_Filter &filter, bool reverse_order=false);
  void apply_TwoBodyTerms(int n, const TimeGrid &tgrid, EDM_State &state, const EDM_Filter &filter, bool reverse_order=false);
  void one_step(int n, const TimeGrid &tgrid, EDM_State &state, const EDM_Filter &filter);
  void run(const TimeGrid &tgrid, const EDM_State &initial, const EDM_Filter &filter);

  void print_debug(double t, const EDM_State &state);

  //initialization
  void clear();
  void setup(Parameters &param);

/*  inline void add_OneBodyTerm(std::shared_ptr<EDM_OneBodyTerm> term, int site, const EDM_Index &Ldim){
    Ldim.check_in_range(site);
    oneBodyTerms.push_back(std::make_pair(term, site));
  }*/

  EDM_Simulation(){}
  EDM_Simulation(Parameters &param){
    setup(param);
  }
};
}
#endif
