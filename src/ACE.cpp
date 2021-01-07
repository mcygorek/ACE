#include "Simulation_OD.h"

int main(int args, char** argv){

  Parameters param(args, argv, true);
  if(param.is_specified("print_param")){
    param.print(param.get_as_string("print_param"));
  }
  std::string outfile=param.get_as_string("outfile", "ACE.out");

  double ta=param.get_as_double("ta", 0);
  double dt=param.get_as_double("dt", 1e-2);
  double te=param.get_as_double("te", 10);


  FreePropagator fprop(param);

  std::vector<Smart_Ptr<InfluenceFunctional_OD> >IF=IF_from_Parameters(param);


  Eigen::MatrixXcd initial_rho=InitialState(param);

  Simulation_OD sim;
  sim.setup_output(param);

  cavityfy(fprop, IF, initial_rho, sim.output_Op, param);

  if(param.get_as_string("print_dims")!=""){
    std::string str=param.get_as_string("print_dims");
    std::ofstream ofs_dims(str.c_str());
    IF[0]->print_dims(ofs_dims);
  }

  std::cout<<"Calculating"<<std::endl;
  sim.print_timestep=param.get_as_bool("print_timestep", param.get_as_bool("print_timesteps", false));
  sim.run(fprop, IF, ta, dt, te, initial_rho);
//  sim.run_nobath(fprop, ta, dt, te, initial_rho);
  sim.print_results(outfile);


  if(param.is_specified("print_pulse_FT"))Pulse_Printer(param, fprop);

  return 0;
}

