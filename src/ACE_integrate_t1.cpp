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
  double t1_max=param.get_as_double("t1_max", te/2.);
  if(t1_max>=te){ 
    std::cerr<<"t1_max must be less than te!"<<std::endl;
    exit(1);
  }
  int n_max_t1=(t1_max-ta)/dt;
  int n_max=(te-ta)/dt;
  int n_max_tau=n_max-n_max_t1;



  std::vector<Smart_Ptr<InfluenceFunctional_OD> >IF=IF_from_Parameters(param);



  if(param.is_specified("outfile_noop")){
    std::cout<<"Performing reference calculation..."<<std::endl;

    Eigen::MatrixXcd initial_rho=InitialState(param);
    FreePropagator fprop(param);

    Simulation_OD sim;
    sim.setup_output(param);
    cavityfy(fprop, IF, initial_rho, sim.output_Op, param);

    sim.print_timestep=param.get_as_bool("print_timestep", false);
    sim.run(fprop, IF, ta, dt, te, initial_rho);
    sim.print_results(param.get_as_string("outfile_noop"));
  }

  

  std::vector<std::vector<std::string> > loop_ops=param.get("loop_Operator");
  if(loop_ops.size()<1){
    std::cerr<<"For ACE_integrate_t1: Need to specify parameter 'loop_Operator'!"<<std::endl;
    exit(1);
  }

  std::vector<std::pair<double, std::string> > print_intermediate;
  { 
    std::vector<std::vector<std::string> > svv=param.get("print_fixed_t1");
    for(size_t i=0; i<svv.size(); i++){
      if(svv[i].size()<2){
        std::cerr<<"Usage: -print_fixed_t1 TIME FILENAME"<<std::endl;
        exit(1);
      }
      double t=Reader::readDouble(svv[i][0],"-print_fixed_t1 TIME FILENAME");
      print_intermediate.push_back(std::make_pair(t, svv[i][1]));
    }
  }

  Simulation_Results results;
  for(int nt1=0; nt1<=n_max_t1; nt1++){
    double t1=ta+nt1*dt;
    std::cout<<"t1: "<<t1<<" te: "<<te<<std::endl;

    Parameters param2(param);

    for(size_t i=0; i<loop_ops.size(); i++){ 
      std::stringstream ss; ss<<t1<<" "<<loop_ops[i][0];
      param2.add_to("apply_Operator", ss.str());   
    }

    Eigen::MatrixXcd initial_rho=InitialState(param2);
    FreePropagator fprop(param2);

    Simulation_OD sim;
    sim.setup_output(param2);
    cavityfy(fprop, IF, initial_rho, sim.output_Op, param2);

    sim.print_timestep=param2.get_as_bool("print_timestep", false);
    sim.run(fprop, IF, ta, dt, te, initial_rho);
  
    if(nt1==0){
      results.list.resize(n_max_tau);
      for(int i=0; i<n_max_tau; i++){
        results.list[i].first=i*dt;
        results.list[i].second.resize(sim.results.list[i].second.size());
        for(size_t j=0; j<sim.results.list[i].second.size(); j++){
          results.list[i].second[j]=sim.results.list[i].second[j]*dt;
        }
      }
    }else{
      for(int i=0; i<n_max_tau; i++){
        for(size_t j=0; j<sim.results.list[i].second.size(); j++){
          results.list[i].second[j]+=sim.results.list[i+nt1].second[j]*dt;
        }
      }
    }
 
    for(size_t i=0; i<print_intermediate.size(); i++){
      if(abs(t1-print_intermediate[i].first)<1e-12){
        param2.print("TEST.TXT");
        sim.print_results(print_intermediate[i].second);
      }
    }
  }
 
  results.print(outfile);


  return 0;
}

