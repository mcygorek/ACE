#include "Simulation.h"


int main(int args, char** argv){

  Parameters param(args, argv, true);
  if(param.is_specified("print_param")){
    param.print(param.get_as_string("print_param"));
  }
  std::string outfile=param.get_as_string("outfile", "results.dat");


  double ta=param.get_as_double("ta", 0);
  double dt=param.get_as_double("dt", 0.5);
  double te=param.get_as_double("te", 30);
  int n_max=param.get_as_double("n_max", 5);

  bool silent=param.get_as_bool("silent",false);
  bool use_bath=param.get_as_bool("use_bath", true);
  bool noSubPS= !param.get_as_bool("subtract_polaron_shift", true);
  double temperature=param.get_as_double("temperature", 10.);

  Operators2x2  op;

//  InitialState initial_rho(param);
  Eigen::MatrixXcd initial_rho=InitialState(param);


  ///Bath-free propagator:
  FreePropagator fprop(param);


  ///Bath
  RealFunctionPtr SD=new RealFunction_Zero_Class();
  if(use_bath)SD=new SpectralDensity_QD();

  InfluenceFunctional IF(n_max, dt, op.ketbra(1,1), SD, temperature, noSubPS);

  if(param.is_specified("print_IF"))IF.print(param.get_as_string("print_IF"));


  Simulation sim;
  sim.add_output_Op( op.ketbra(1,1) );
  sim.add_output_Op( op.ketbra(0,0) );
  sim.add_output_Op( op.ketbra(0,1) );

  if(use_bath){
    sim.run(fprop, IF, ta, dt, te, initial_rho, silent);
  }else{
    sim.run_nobath(fprop, ta, dt, te, initial_rho, silent);
  }
  sim.print_results(outfile);


  return 0;
}

