#include "ACE.hpp"

using namespace ACE;

int main(int args, char** argv){

  Parameters param(args, argv, true);

  std::string outfile=param.get_as_string_check("outfile");
  std::string Gaussian_prefix=param.get_as_string("Gaussian_prefix", "Boson");
  double mem_threshold=param.get_as_double("mem_threshold",-1);
  
  DiagBB diagBB(param, Gaussian_prefix);
 
  TimeGrid tgrid(param);



#ifdef K_EXISTING_FUNCTIONS

  int n_mem=diagBB.estimate_memory_length(tgrid.n_calc, tgrid.dt, mem_threshold, true);
  std::cout<<"n_mem: "<<n_mem<<std::endl;
 
  diagBB.print_K(outfile, n_mem, tgrid.dt);

#else 

/*
  std::ofstream ofs(outfile.c_str());
  double ref=0.;
  int below_thr=-1;
  for(int n=0; n<tgrid.n_calc; n++){
    double t=n*tgrid.dt;
    std::complex<double> K=diagBB.calculate_K(n, tgrid.dt);
    if(n==0 || n==1){ 
      if(abs(K)>ref)ref=abs(K);
    }else if(mem_threshold>0.){
      if(abs(K)>mem_threshold*ref){
        below_thr=-1.;
      }else{
        if(below_thr<0){
          below_thr=n;
        }
      }
    }
    ofs<<t<<" "<<K.real()<<" "<<K.imag()<<std::endl;
  } 

  if(mem_threshold>0. && below_thr>=0){
    std::cout<<"below threshold from time: "<<below_thr*tgrid.dt<<std::endl;
  }
*/

// Find memory time: relative to maximal (max) value around K(t=0).
// K can be oscillatory. So, after finding first |K| < threshold * max, go at least 2 twice as long to see if values larger than that are found


  std::ofstream ofs(outfile.c_str());
  double max=0.;
  int below_thr=-1;
  int first_below_thr=-1;
  int n_break=-1;

  std::cout<<"Analyzing memory kernel"<<std::endl;
  for(int n=0; n<tgrid.n_calc; n++){
    double t=n*tgrid.dt;
    std::complex<double> K=diagBB.calculate_K(n, tgrid.dt); ///(tgrid.dt*tgrid.dt);

    ofs<<t<<" "<<K.real()<<" "<<K.imag()<<std::endl;

    if(mem_threshold>0.){
      if(n<=10){ 
        if(abs(K)>max)max=abs(K);
      }else{
        if(abs(K)>mem_threshold*max){
          below_thr=-1.;
        }else{
          if(below_thr<0){
            below_thr=n;
            if(first_below_thr<0)first_below_thr=n;
          }
        }
        //test break criterion:
        if(below_thr>0 && n>below_thr+2*first_below_thr>0){
          n_break=n;
          break;
        }
      }
    }
  } 

  if(mem_threshold>0){
    if(n_break>=0){
      std::cout<<"First time below threshold: "<<first_below_thr*tgrid.dt<<std::endl;
      std::cout<<"Estimated memory time: "<<n_break*tgrid.dt<<std::endl;
    }else{
      std::cout<<"No suitable cutoff found within n_max="<<tgrid.n_calc<<std::endl;
    }
  }
 
#endif


  return 0;
}
