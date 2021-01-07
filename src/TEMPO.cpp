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
  double temperature=param.get_as_double("temperature", 10.);
   

  Eigen::MatrixXcd initial_rho=InitialState(param);

  ///Bath-free propagator:
  FreePropagator fprop(param);

  ///Bath
  RealFunctionPtr SD=new RealFunction_Zero_Class();
  if(use_bath)SD=new SpectralDensity_QD();

  Operators2x2 op;
  InfluenceFunctional_Vector IF(n_max, dt, op.ketbra(1,1), SD, temperature);

/*
  if(initial_rho.rows()!=2){
    int MTLS=param.get_as_size_t("expand_MTLS");
    if(MTLS>1){
      int dim=pow(2, MTLS);
      if(initial_rho.rows()!=dim){
        std::cerr<<"Error using 'expand_MTLS': Mismatch in dimensions between InfluenceFunctional and initial state!"<<std::endl; 
        exit(1);
      }
//TODO: expand!
    }else{ 
      std::cerr<<"Error: Mismatch in dimensions between InfluenceFunctional and initial state!"<<std::endl; 
      exit(1);
    }
  }
*/


  Simulation_TEMPO sim;
  sim.setup_output(param);
  sim.compressor=RankCompressor_Selector(param);


  if(use_bath==false){
    sim.run_nobath(fprop, ta, dt, te, initial_rho, silent);
  }else{
    sim.run(fprop, IF, ta, dt, te, initial_rho, silent);
  }

  sim.print_results(outfile);


  return 0;
}

