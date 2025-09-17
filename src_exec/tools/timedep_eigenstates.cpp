#include "FreePropagator.hpp"
#include "Parameters.hpp"
#include "TimeGrid.hpp"
#include "Simulation_Results.hpp"
#include <iostream>

using namespace ACE;

int main(int args, char** argv){

 try{
  Parameters param(args, argv, true);
  TimeGrid tgrid(param);
  FreePropagator fprop(param);

  std::string outfile=param.get_as_string_check("outfile");
  // add .ds to outfile
  outfile += ".ds";
  std::cout << "Writing to " << outfile << std::endl;

  Simulation_Results Eigenv_res;
  Eigenv_res.resize(tgrid.n_tot+1);

  bool print_eigenstates=param.get_as_bool("print_eigenstates",false);
  bool use_nonhermitian=param.get_as_bool("use_nonhermitian",false);

  // fill up results object for each time step
  for(int n=0; n<tgrid.n_tot+1; n++){
    double t=tgrid.get_t(n);

    Eigen::MatrixXcd H=fprop.get_Htot(t);
    if(use_nonhermitian){
      H=fprop.get_H_nonhermitian(t);
    }
    int N=H.rows();
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver(H);


    // Sort:
    std::vector<std::pair<int,std::complex<double> > > srt(N);
    for(int i=0; i<N; i++){  
      srt[i].first=i;
      srt[i].second=solver.eigenvalues()(i);
    }
    std::sort(srt.begin(), srt.end(), 
      [](const std::pair<int,std::complex<double> > & a,
         const std::pair<int,std::complex<double> > & b){
           return  a.second.real() < b.second.real();   });


    // Collect for writing:
    Eigenv_res[n].first = t;
    // pushes all elements of eigenvector i, ...
    for(int i=0; i<N; i++){
      Eigenv_res[n].second.push_back(solver.eigenvalues()(srt[i].first));
    } 
    // ... then all elements of eigenvector i+1, etc.
    if(print_eigenstates){
      for(int i=0; i<N; i++){
        for (int j=0; j<N; j++){
          Eigenv_res[n].second.push_back(solver.eigenvectors()(j,srt[i].first));
        }
      }
    }
  }

  Eigenv_res.print(outfile);
 }catch(std::exception &e){
 }
  return 0;
}
