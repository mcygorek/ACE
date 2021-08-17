#include "ACE.h"

int main(int args, char** argv){

  Parameters param(args, argv, true);
  if(param.is_specified("print_param")){
    param.print(param.get_as_string("print_param"));
  }
  std::string outfile=param.get_as_string("outfile", "ACE.out");

  IF_TimeGrid tgrid(param);

  Eigen::MatrixXcd initial_rho=InitialState(param);
  FreePropagator fprop(param);

  //Check if system dimensions are compatible
  fprop.set_dim(initial_rho, "Mismatch in dimensions between initial system state and system propagator!\nPlease specify 'initial' and check whether the dimensions agree with 'add_Hamiltonian', 'add_Lindblad', ...");
  fprop.set_dim(Output_Ops(param).get_dim(), "Mismatch in dimensions between observables system propagator!\nPlease specify 'add_Output' and check whether the dimensions agree with 'add_Hamiltonian', 'add_Lindblad', ...");

  
  std::vector<Smart_Ptr<InfluenceFunctional_OD> >IF=IF_from_Parameters(param, fprop.get_dim());

  if(param.get_as_string("print_dims")!=""){
    std::string str=param.get_as_string("print_dims");
    std::ofstream ofs_dims(str.c_str());
    IF[0]->print_dims(ofs_dims);
  }


  std::cout<<"Calculating"<<std::endl;

  //Be aware that IF may be generated with coarse graining. Then, run dynamics on coarser time grid:
  IF_TimeGrid tgrid_coarse;  
  tgrid_coarse.setup_coarse(param);

  Simulation_OD sim;
  sim.setup_output(param);
  sim.use_symmetric_Trotter=param.get_as_bool("use_symmetric_Trotter",false);
  if(sim.use_symmetric_Trotter){
    std::cout<<"Using symmetric Trotter decomposition"<<std::endl;
  }
  
  sim.print_timestep=param.get_as_bool("print_timestep", param.get_as_bool("print_timesteps", false));

  sim.run(fprop, IF, tgrid_coarse, initial_rho);

  sim.print_results(outfile);


//  if(param.is_specified("print_pulse_FT"))Pulse_Printer(param, fprop);

  return 0;
}

