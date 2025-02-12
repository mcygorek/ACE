#include <iostream>
#include <fstream>
#include "Parameters.hpp"
#include "TimeGrid.hpp"
#include <map>
#include <random>

using namespace ACE;

int main(int args, char **argv){
  Parameters param(args, argv);

  TimeGrid tgrid(param);
  double dist_mean=param.get_as_double("dist_mean",0.);
  double dist_std=param.get_as_double("dist_std",1.);
  double histowidth=param.get_as_double("histo_width",0.1);
  std::string print_histo=param.get_as_string("print_histo","");
  std::string outfile=param.get_as_string("outfile","noise.pulse");

  double revert=param.get_as_double("revert",0.);
  double sigma=param.get_as_double("sigma",1.);

/*
  enum MODE {RW, OU} mode=RW; //RW=random walk, OU=Ornstein-Uhlenbeck
  if(param.get_as_bool("use_OU")){
    mode=OU;
  }*/


  std::mt19937 mt;
  int seed=param.get_as_int("seed",0);
  if(seed>0){
    std::seed_seq seq{seed};
    mt.seed(seq);
  }else{
    std::random_device rd;
    std::seed_seq seq{rd()};
    mt.seed(seq);
  }
  
  std::normal_distribution<double> dist(dist_mean, dist_std);

  if(param.get_as_bool("only_print_single")){
    std::cout<<(double)dist(mt)<<std::endl;
    return 0;
  }

  std::map<int, int> hist;
  std::ofstream ofs(outfile.c_str());
  double val=0;
  if(param.get_as_bool("start_from_distr")){
    val=dist(mt)*sigma;
  }
  for(int n=0; n<=tgrid.n_tot; n++){
    double t=tgrid.get_t(n);
    double dt=tgrid.get_dt(n);
    double rand=dist(mt);
     
    val+=-revert*val*dt;
    val+=rand*sigma*sqrt(2.*dt*revert);
    
    ofs<<t<<" "<<val<<" 0"<<std::endl;
    if(print_histo!=""){
      ++hist[std::round(val/histowidth)];
    }
  }

  if(print_histo!=""){
    std::ofstream ofs(print_histo.c_str());
    for(const std::pair<int,int> & h : hist){
//      ofs<<h.first*histowidth<<" "<<h.second*histowidth<<std::endl;
      ofs<<h.first*histowidth<<" "<<h.second/histowidth/(tgrid.n_tot+1.)<<std::endl;
    }
  }

  return 0;
}
