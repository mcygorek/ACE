#include "Parameters.hpp"
#include "slowFT.hpp"
#include "ReadTable.hpp"
#include "Constants.hpp"
#include "Simulation_Results.hpp"
#include "discreteFT.hpp"
#include <iostream>
#include <fstream>

using namespace ACE;

int main(int args, char ** argv){
  Parameters param(args, argv);

  std::string outfile=param.get_as_string_check("outfile");
  std::string infile=param.get_as_string_check("infile");
  int col=param.get_as_size_t("col",2);  //number of column starts with 1 for time, then come real and imaginary parts as pairs
  if(col%2!=0 || col<2){
    std::cerr<<"col="<<col<<" must be an even number >=2!"<<std::endl;
    exit(1);
  }

  int sig=param.get_as_int("sign",1);
  int integrate_mode=param.get_as_int("integrate_mode",1);
  double t_zero=param.get_as_double("t_zero",0.);

  if(param.get_as_bool("use_FFT",true)==true){
    Simulation_Results in(infile);
    if(in.size()<4){
      std::cerr<<"'"<<infile<<"' contains less than 4 usable lines!"<<std::endl;
      exit(1);
    } 
    if(param.is_specified("t_start")){
      double t_start=param.get_as_double_check("t_start");
      int first_start=in.size();
      for(int i=0; i<in.size(); i++){
        if(in[i].first>=t_start){
          first_start=i;
          break;
        }
      }  
      if(first_start>=(int)in.size()-4){
        std::cerr<<"t_start="<<t_start<<" exceeds domain of file '"<<infile<<"'!"<<std::endl;
        exit(1);
      }
      in.list.erase(in.list.begin(), in.list.begin()+first_start);
    }

    if(param.get_as_bool("subtract_final")){
      for(size_t i=0; i<in.size(); i++){
        if(in[i].second.size()<=(col-2)/2 || in.back().second.size()<=(col-2)/2){
          std::cerr<<"Not enough columns in file '"<<infile<<"': ";
          std::cerr<<"in["<<i<<"].second.size()="<<in[i].second.size()<<" vs. "<<col+1<<std::endl;
          exit(1);
        }
        in[i].second[(col-2)/2]-=in.back().second[(col-2)/2];
      }
    }    

    int N=param.get_as_int("Ndiscr",0);
    if(N>1){
      if(in.size()>N){ //truncate input
        in.resize(N);
      }else if(in.size()<N){//zero-pad (in addition to expanding to powers of 2)
        in.list.resize(N,std::pair<double, std::vector<std::complex<double> > >
(0.,std::vector<std::complex<double> >(col+1, 0.)));   
      }
    }

    //shift time zero:
    double ta=in[0].first;
    double dt=in[1].first-in[0].first;
    for(int i=0; i<(int)in.size(); i++){
       in[i].first=ta-t_zero+i*dt;
    }
 
    Simulation_Results out=resultsFFT(in, (col-2)/2, sig, integrate_mode);
    out.print(outfile);
    return 0;
  }


  double wa=param.get_as_double("wa",param.get_as_double("Ea")/hbar_in_meV_ps); 
  double wb=param.get_as_double("wb",param.get_as_double("Eb")/hbar_in_meV_ps); 
  if(param.is_specified("we")){
    std::cout<<"Use 'wb' instead of 'we'!"<<std::endl;
    exit(1);
  }
  int Ndiscr=param.get_as_size_t("Ndiscr",1000);
//  double ta=param.get_as_double("ta");
//  double dt=param.get_as_double_check("dt");
  int Nsubdiv=param.get_as_size_t("Nsubdiv",1);
  double dw=(wb-wa)/Ndiscr;

  if(!(wb>wa)){
    std::cerr<<"wb not greater than wa!"<<std::endl;
    exit(1);
  }


  ReadTable tab(infile, 0, col/2, col/2+1);
  int N=tab.table.size();
  if(N<3){
    std::cerr<<"Simpson integration needs at least 3 data point!"<<std::endl;
    exit(1);
  }
  double ta=tab.table[0][0]-t_zero;
  double dt=tab.table[1][0]-tab.table[0][0];
  double te=tab.table.back()[0]-t_zero;

  std::vector<std::complex<double> > in(tab.table.size(), 0.);
  for(size_t i=0; i<tab.table.size(); i++){
    in[i]=std::complex<double>( tab.table[i][1], tab.table[i][2] );
  }


  if(param.is_specified("t_start")){
    double t_start=param.get_as_double_check("t_start")-t_zero;
    int first_start=round((t_start-ta)/dt);
    if(first_start>=(int)in.size()-4){
      std::cerr<<"t_start="<<t_start<<" exceeds domain of file '"<<infile<<"'!"<<std::endl;
      exit(1);
    }
    in.erase(in.begin(), in.begin()+first_start);
    ta=t_start;
  }


  std::vector<std::complex<double> > result=slowFT(in, wa, wb, Ndiscr, ta, dt, sig, Nsubdiv);

  std::cout<<"wa="<<wa<<" wb="<<wb<<" dw="<<dw<<" Ndiscr="<<Ndiscr<<std::endl;
  std::cout<<"result.size()="<<result.size()<<std::endl;

  std::ofstream ofs(outfile.c_str());
  for(size_t i=0; i<result.size(); i++){
    ofs<<wa+i*dw<<" "<<result[i].real()<<" "<<result[i].imag()<<std::endl;
  }  
}
