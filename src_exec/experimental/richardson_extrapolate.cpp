#include "ACE.hpp"
#include "ReadTable.hpp"

/*
Richardson extrapolation: 
Combine two ACE calculations with different time step widths dt.
If ACE has an error ~dt^2, Richardson extrapolation will lead to an error ~dt^3

Usage: e.g.:
> richardson_extrapolate -infile dt5e-2.out dt1e-1.out  -outfile richardson.out

*/

using namespace ACE;

class ACE_Outfile : public Printable{
public:
  int nmax;
  int entries;  // number of complex entries: (columns-1)/2  
  double dt;
  double ta;
  std::vector<std::vector<std::complex<double> > > table;

  virtual std::ostream & print(std::ostream & os=std::cout)const{
    os<<"nmax: "<<nmax<<" ta: "<<ta<<" dt: "<<dt<<" entries: "<<entries;
    return os;
  }

  void read(const std::string &fname){
    int dcols=ReadTable::colums_in_first_row(fname);
    if(dcols%2!=1){
      std::cerr<<"ACE output file '"<<fname<<"' should have an odd number of columns!"<<std::endl;
      exit(1);
    }
    if(dcols<3){
      std::cerr<<"ACE output file '"<<fname<<"' should have at least 3 columns!"<<std::endl;
      exit(1);
    }
    entries=(dcols-1)/2.;
    ReadTable RT(fname); 
    if(RT.size()<2){
      std::cerr<<"ACE output file '"<<fname<<"' should have at least 2 lines!"<<std::endl;
      exit(1);
    }
    nmax=RT.size();
    ta=RT.table[0][0];
    dt=RT.table[1][0]-ta;

    table.clear();
    table.resize(RT.size(), std::vector<std::complex<double> >(entries, 0.));
    for(size_t n=0; n<RT.size(); n++){
      for(size_t e=0; e<table[n].size(); e++){
        table[n][e]=std::complex<double>(RT.table[n][2*e+1], RT.table[n][2*e+2]);
      }
    }
  }
  ACE_Outfile(){}
  ACE_Outfile(const std::string &fname){
    read(fname);
  }
};

std::vector<std::vector<std::complex<double> > > expand_interpolate(
     const std::vector<std::vector<std::complex<double> > > & tab_in, 
                                              int samples, int order){
  std::vector<std::vector<std::complex<double> > > tab_expand;

  if(tab_in.size()<order+1){
    std::cerr<<"Error: expand_interpolate: tab_in.size()<order+1!"<<std::endl;
    exit(1);
  }

  int n_full_intervals=(tab_in.size()-1)/order;
  int n_res=tab_in.size()-(n_full_intervals*order+1);

std::cout<<"tab_in.size(): "<<tab_in.size()<<std::endl;
std::cout<<"n_full_intervals: "<<n_full_intervals<<std::endl;
std::cout<<"n_res: "<<n_res<<std::endl;
  
  for(int n=0; n<n_full_intervals; n++){
    for(int s=0; s<samples*order; s++){
      std::vector<std::complex<double> > tmp(tab_in[0].size(), 0.);
      for(size_t e=0; e<tmp.size(); e++){
        for(int j=0; j<order+1; j++){
          double L=1;
          for(int k=0; k<order+1; k++){
            if(k==j)continue;
            L*=(((double)s/(double)samples)-(double)k)/((double)j-(double)k);
          }
          tmp[e]+=tab_in[n*order+j][e]*L;
        }
      }
      tab_expand.push_back(tmp);
    }
  }
  if(n_res!=0){
    for(int s=samples*(order-n_res); s<samples*order; s++){
      std::vector<std::complex<double> > tmp(tab_in[0].size(), 0.);
      for(size_t e=0; e<tmp.size(); e++){
        for(int j=0; j<order+1; j++){
          double L=1;
          for(int k=0; k<order+1; k++){
            if(k==j)continue;
            L*=(((double)s/(double)samples)-(double)k)/((double)j-(double)k);
          }
//std::cout<<"TEST: s: "<<s<<" j: "<<j<<" arg: "<<(n_full_intervals-1)*order+n_res+j<<" tab_in.size(): "<<tab_in.size()<<std::endl;
          tmp[e]+=tab_in[(n_full_intervals-1)*order+n_res+j][e]*L;
        }
      }
      tab_expand.push_back(tmp);
    }
  }
  tab_expand.push_back(tab_in.back());


  return tab_expand;
}

int main(int args, char **argv){
  Parameters param(args, argv, true);

  param.complain_if_not_specified("infile");
  std::string outfile=param.get_as_string_check("outfile");
  int richardson_n=param.get_as_int("richardson_n",2);

  std::vector<std::string> infiles=param.get_all_strings("infile");
  if(infiles.size()!=2){
    std::cerr<<"Richardson extrapolation implemented only for two input files!"<<std::endl;
  }

  std::vector<ACE_Outfile> aof;
  for(size_t i=0; i<infiles.size(); i++){
    aof.push_back(ACE_Outfile(infiles[i]));
  }
  for(size_t i=0; i<aof.size(); i++){
    std::cout<<aof[i]<<std::endl;
  }

  if(aof[0].entries!=aof[1].entries){
    std::cerr<<"aof[0].entries!=aof[1].entries!"<<std::endl;
    exit(1);
  }
 
  int richardson_t=(aof[1].dt/aof[0].dt+0.5);
  std::cout<<"richardson parameters: n="<<richardson_n<<" t="<<richardson_t<<std::endl;
  if(richardson_t<2){
    std::cerr<<"richardson_t has to be at least 2!"<<std::endl;
    exit(1);
  }

  if(fabs(richardson_t * aof[1].dt- aof[0].dt)<1e-16){
    std::cerr<<"richardson_t * dt1 != dt2 !"<<std::endl;
    exit(1);
  }
  if(aof[0].nmax-1 != richardson_t*(aof[1].nmax-1)){
    std::cerr<<"aof[0].nmax-1 != richardson_t*(aof[1].nmax-1) !"<<std::endl;
    exit(1);
  }
   
  double rtn=pow(richardson_t, richardson_n);

  
  int interpolate_order=param.get_as_size_t("interpolate_order",0);

  if(interpolate_order<1){
    std::ofstream ofs(outfile.c_str());
    for(int i=0; i<aof[1].nmax; i++){
      double time=aof[1].ta+i*aof[1].dt;
      ofs<<time;
      for(int e=0; e<aof[1].entries; e++){
        std::complex<double> A0=aof[0].table[i*richardson_t][e];
        std::complex<double> A1=aof[1].table[i][e];
        std::complex<double> value=(rtn*A0 -1.*A1)/(rtn-1.);
        ofs<<" "<<value.real()<<" "<<value.imag();
      }
      ofs<<std::endl;
    }
  }else{

    std::vector<std::vector<std::complex<double> > > tab1=expand_interpolate(
                               aof[1].table, richardson_t, interpolate_order);

    

    std::string interpolate_file=param.get_as_string("print_interpolated");
    if(interpolate_file!=""){
      std::ofstream ofs(interpolate_file.c_str());
      for(size_t i=0; i<tab1.size(); i++){
        double time=aof[1].ta+i*aof[1].dt/(double)richardson_t;
        ofs<<time;
        for(size_t e=0; e<tab1[i].size(); e++){
          ofs<<" "<<tab1[i][e].real()<<" "<<tab1[i][e].imag();
        }
        ofs<<std::endl;
      }
    }


    std::ofstream ofs(outfile.c_str());
    for(int i=0; i<aof[0].nmax; i++){
      double time=aof[0].ta+i*aof[0].dt;
      ofs<<time;
      for(int e=0; e<aof[0].entries; e++){
        std::complex<double> A0=aof[0].table[i][e];
        std::complex<double> A1=tab1[i][e];
        std::complex<double> value=(rtn*A0 -1.*A1)/(rtn-1.);
        ofs<<" "<<value.real()<<" "<<value.imag();
      }
      ofs<<std::endl;
    }
  }
 
  return 0;
}
