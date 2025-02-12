#include "MPG_Discretization.hpp"
#include "Parameters.hpp"
#include "EnergyRange.hpp"
#include "SpectralDensity_Selector.hpp"
#include "Reader.hpp"
#include "ReadTable.hpp"
#include "DummyException.hpp"

namespace ACE{

  void MPG_Discretization::check_bounds(const std::string fct, int k)const{
    if( (int)E.size() != N ){
      std::cerr<<"Error: MPG_Discretization::check_bounds: inconsistent sizes: E.size() != N!"<<std::endl;
      throw DummyException();
    }
    if( (int)dE.size() != N ){
      std::cerr<<"Error: MPG_Discretization::check_bounds: inconsistent sizes: dE.size() != N!"<<std::endl;
      throw DummyException();
    }
    if(k<0||k>=N){
      std::cerr<<"Error: MPG_Discretization::"<<fct<<" out of bounds: "<<k<<"/N"<<std::endl; 
      throw DummyException();
    }
  }

  double MPG_Discretization::get_E(int k)const{
    check_bounds("get_E",k);
    return E[k];
  }
  double MPG_Discretization::get_dE(int k)const{
    check_bounds("get_dE",k);
    return dE[k];
  }
  double MPG_Discretization::get_omega(int k)const{
    check_bounds("get_omega",k);
    return E[k]/hbar_in_meV_ps;
  }
  double MPG_Discretization::get_domega(int k)const{
    check_bounds("get_domega",k);
    return dE[k]/hbar_in_meV_ps;
  }


  std::string MPG_Discretization::add_name(const std::string &mpgname, const std::string &str)const{
    if(mpgname=="")return str;
    return std::string(mpgname+"_"+str);
  }  

  void MPG_Discretization::setup(Parameters &param, const std::string &mpgname){
    E.clear(); dE.clear();
    N=param.get_as_size_t(add_name(mpgname,"N_modes"));
    if(N<1)return;


    EnergyRange range(param, mpgname);
    //Discretize so that E_max and E_min are sample points,
    //as opposed to sampling mid-points of intervals: 
    //('false' is more consistent with continuum limit)
    bool sample_range_ends=param.get_as_bool(add_name(mpgname,"sample_range_ends"),false); 
    if(N==1){sample_range_ends=false;}//avoid division by zero in this case

    double interval=range.E_range()/N;
    if(sample_range_ends){ 
      interval=range.E_range()/(N-1);
    }
     
    dE.clear();
    dE.resize(N,interval);
    E.resize(N);
    for(int k=0; k<N; k++){
      E[k]=range.E_min()+(sample_range_ends? 0 : 0.5*interval)+k*interval;
    }

    //use twice as many modes: for every energy E, add an energy -E
    if(param.get_as_bool(add_name(mpgname,"E_duplicate_negative"),false)){
      for(int k=0; k<N; k++){
        E.push_back(-E[k]);
        dE.push_back(dE[k]);
      }
      N*=2;
    }
  }




//MPG_Discretization for a one-parameter interaction
  void MPG_Discretization_E_g::check_bounds(const std::string fct, int k)const{
    MPG_Discretization::check_bounds(fct, k);
    if( (int)g.size() != N ){
      std::cerr<<"Error: MPG_Discretization_E_g::check_bounds: inconsistent sizes: g.size() != N!"<<std::endl;
      throw DummyException();
    }
  }
  void MPG_Discretization_E_g::zero_pad(int N_new){
    N = N_new;
    E.resize(N_new, 0.);
    dE.resize(N_new, 0.);
    g.resize(N_new, 0.);
  }

  double MPG_Discretization_E_g::get_g(int k)const{
    check_bounds("get_g",k);
    return g[k];
  }

  double MPG_Discretization_E_g::rate_from_g(double g, double DOS){
    return 2.*M_PI*hbar_in_meV_ps*g*g*DOS;
  }
  double MPG_Discretization_E_g::g_from_rate(double rate, double DOS){
    return sqrt(rate/(2.*M_PI*hbar_in_meV_ps*DOS));
  }
  bool MPG_Discretization_E_g::compare_abs(const std::pair<double,int> &p1,
                          const std::pair<double,int> &p2){
    return fabs(p1.first)>fabs(p2.first);
  }
  bool MPG_Discretization_E_g::compare_abs_less(const std::pair<double,int> &p1,
                          const std::pair<double,int> &p2){
    return fabs(p1.first)<fabs(p2.first);
  }
  bool MPG_Discretization_E_g::compare_greater(const std::pair<double,int> &p1,
                          const std::pair<double,int> &p2){
    return (p1.first)>(p2.first);
  }
  bool MPG_Discretization_E_g::compare_less(const std::pair<double,int> &p1,
                          const std::pair<double,int> &p2){
    return (p1.first)<(p2.first);
  }

  void MPG_Discretization_E_g::print(std::ostream &os)const{
    for(int i=0; i<N; i++){
      os<<get_E(i)<<" "<<get_g(i)<<std::endl;
    }
  } 
  void MPG_Discretization_E_g::print(const std::string &fname)const{
    std::ofstream ofs(fname.c_str());
    print(ofs);
  }

  
  void MPG_Discretization_E_g::setup(Parameters &param, const std::string &mpgname){
    MPG_Discretization::setup(param, mpgname);

    std::string sort_default="none"; //"E_abs";
 

    //construct sets {E,g} by interpolating spectral density J(omega) from file
//    if(param.is_specified(add_name(mpgname,"J_from_file")) ||
//       param.is_specified(add_name(mpgname,"J_type")) ||
//       param.is_specified(add_name(mpgname,"rate"))){
      if(SpectralDensity_Selector::is_specified(param, add_name(mpgname,"J"))){

 
//      SpectralDensity_Selector SD_select(param, add_name(mpgname,"J"));
//      RealFunctionPtr J=SD_select;
      RealFunctionPtr J=SpectralDensity_Selector(param, add_name(mpgname,"J"));

      g.clear(); g.resize(N);
      for(int k=0; k<N; k++){
        double w=get_E(k)/hbar_in_meV_ps;
        double dw=get_dE(k)/hbar_in_meV_ps;
        g[k]=sqrt( J->f(w) * dw );
      }
    //read sets {E,g} from a file
    }else if(param.is_specified(add_name(mpgname,"E_g_from_file"))||
             param.is_specified(add_name(mpgname,"omega_g_from_file")) ){

      sort_default="none";     
      E.clear(); g.clear();
      
      double hbar_scale=hbar_in_meV_ps;
      std::string ftparam=add_name(mpgname,"omega_g_from_file");
      if(param.is_specified(add_name(mpgname,"E_g_from_file"))){
        ftparam=add_name(mpgname,"E_g_from_file");
        hbar_scale=1.;
      }
      std::vector<std::vector<std::string> > entries=param.get(ftparam);
      for(size_t r=0; r<entries.size(); r++){
        if(entries[r].size()<1){
          std::cerr<<"Error reading parameter '"<<ftparam<<"': No argument!"<<std::endl;
          throw DummyException();
        }
        std::string fname=entries[r][0];
        //which column of the file will be the 'E' and which will be the 'g'?
        int colE=0, colg=1;
        if(entries[r].size()>1){
          colE=readSizeT(entries[r][1],std::string(ftparam+": column_E"));
        }
        if(entries[r].size()>2){
          colg=readSizeT(entries[r][2],std::string(ftparam+": column_g"));
        }

        ReadTable tab(entries[r][0], colE, colg);
        for(size_t tr=0; tr<tab.size(); tr++){
          E.push_back(tab[tr][colE]*hbar_scale);
          g.push_back(tab[tr][colg]);
std::cout<<"E="<<E.back()<<" g="<<g.back()<<std::endl;
        }
      }
      N=E.size();
      dE.clear(); dE.resize(E.size());


    }else{
      //flat spectral density 
      double g_=param.get_as_double(add_name(mpgname,"g"), 0.1);
      g.clear(); g.resize(N, g_);
 
    }
    
    //---------------------------------------------------
    {
    double shift_E=0.;
    int nrr=param.get_nr_rows(add_name(mpgname,"shift_E"));
    for(int i=0; i<nrr; i++){
      shift_E+=param.get_as_double(add_name(mpgname,"shift_E"),0,i,0);
    }
    nrr=param.get_nr_rows(add_name(mpgname,"shift_omega"));
    for(int i=0; i<nrr; i++){
      shift_E+=hbar_in_meV_ps*param.get_as_double(add_name(mpgname,"shift_omega"),0,i,0);
    }
    for(size_t k=0; k<E.size(); k++){
      E[k]+=shift_E;
    }
    }
    //---------------------------------------------------

    if(param.is_specified(add_name(mpgname,"g_scale"))){
      double g_scale=param.get_as_double(add_name(mpgname,"g_scale"));
      for(size_t k=0; k<g.size(); k++){
        g[k]*=g_scale;
      }
    }

    //---------------------------------------------------

    std::string printEg=param.get_as_string(add_name(mpgname,"print_E_g"));
    if(printEg!="")print(printEg);
    if(param.get_as_bool(add_name(mpgname,"stop_after_print_E_g"),false)){
      std::cout<<"File '"<<printEg<<"' written."<<std::endl;
      std::cout<<"Stopping due to option '"<<add_name(mpgname,"stop_after_print_E_g")<<"'."<<std::endl;
      throw DummyException();
    }
    
    /* sorting: modes: 
       "none":  no sorting
       "E_abs":  decreasing modulus of energies |E[k]|
       "g_abs":  increasing modulus of couplings |g[k]|
       "rabi":   In the problem of Rabi rotations in a driven, off-resonant
                 two-level system, the amplitude of oscillations goes with
                 sin^2(theta/2) = (1-cos(theta))/2 where
                 cos(theta)= delta/sqrt(delta^2 + (hbar g)^2) 
                 -> sort according to expected amplitudes.
                 It turns out this is monotonically related to the angle
                 tan(theta)=hbar|g[k]|/|E[k]|
       "random": randomly shuffle list
    */
    std::string sort_mode=param.get_as_string(add_name(mpgname,"sort"),sort_default);
    std::vector<std::pair<double,int> > sorter(N);
    bool do_sort=false;

    if(sort_mode=="reverse"){
      do_sort=true;
      for(size_t i=0; i<sorter.size(); i++){
        sorter[i]=std::make_pair(E[i], i);
      }
      std::reverse(sorter.begin(), sorter.end());

    }else if(sort_mode=="E_abs"||sort_mode=="E_abs_decreasing"){
      do_sort=true;
      for(size_t i=0; i<sorter.size(); i++){
        sorter[i]=std::make_pair(E[i], i);
      }
      sort(sorter.begin(), sorter.end(), compare_abs);

    }else if(sort_mode=="E_abs_reverse"||sort_mode=="E_abs_increasing"){
      do_sort=true;
      for(size_t i=0; i<sorter.size(); i++){
        sorter[i]=std::make_pair(E[i], i);
      }
      sort(sorter.begin(), sorter.end(), compare_abs_less);

    }else if(sort_mode=="E_half_increasing"){
      do_sort=true;
      for(size_t i=0; i<sorter.size(); i++){
        sorter[i]=std::make_pair(E[i], i);
      }
      sort(sorter.begin(), sorter.end(), compare_less);
      int half=sorter.size()/2;
      std::reverse(sorter.begin()+half, sorter.end());

    }else if(sort_mode=="E_half_decreasing"){
      do_sort=true;
      for(size_t i=0; i<sorter.size(); i++){
        sorter[i]=std::make_pair(E[i], i);
      }
      sort(sorter.begin(), sorter.end(), compare_greater);
      int half=sorter.size()/2;
      std::reverse(sorter.begin()+half, sorter.end());

    }else if(sort_mode=="g_abs"){
      do_sort=true;
      for(size_t i=0; i<sorter.size(); i++){
        sorter[i]=std::make_pair(g[i], i);
      }
      sort(sorter.begin(), sorter.end(), compare_abs_less);

    }else if(sort_mode=="rabi"||sort_mode=="Rabi"){
      do_sort=true;
      for(size_t i=0; i<sorter.size(); i++){
/*
        double hg=hbar_in_meV_ps*g[i];
        double denom=E[i]*E[i]+hg*hg; 
        double ampl=1./2.;
        if(denom>1e-6)ampl=(1.-fabs(E[i])/sqrt(denom))/2.;
        sorter[i]=std::make_pair( denom, i);
*/
        if(fabs(E[i])<1e-6){
          sorter[i]=std::make_pair(1e12, i);
        }else{
          sorter[i]=std::make_pair(fabs(g[i]/E[i]), i);
        }
      }
      sort(sorter.begin(), sorter.end(), compare_abs_less);

    }else if(sort_mode=="random"){
      do_sort=true;
      for(size_t i=0; i<sorter.size(); i++){
        sorter[i]=std::make_pair(E[i], i);
      }
      std::random_shuffle(sorter.begin(), sorter.end());

    }else if(sort_mode=="none"||sort_mode==""){

    }else{
      std::cerr<<"Cannot understand sort mode '"<<sort_mode<<"'. Please try: 'none', 'E_abs', 'Rabi'"<<std::endl;
      throw DummyException();
    }


    if(do_sort){
      std::cout<<mpgname<<": sort mode: "<<sort_mode<<std::endl;

#ifdef DEBUG_MODE_SORTER
      std::cout<<"(E, g) before sorting:"<<std::endl;
      for(size_t i=0; i<sorter.size(); i++){
        std::cout<<"("<<E[i]<<", "<<g[i]<<"), ";
      }
      std::cout<<std::endl;

      std::cout<<"Sorter:"<<std::endl;
      for(size_t i=0; i<sorter.size(); i++){
        std::cout<<"("<<sorter[i].first<<", "<<sorter[i].second<<"), ";
      }
      std::cout<<std::endl;
#endif

      std::vector<double> E2(N);
      for(size_t i=0; i<sorter.size(); i++)E2[i]=E[sorter[i].second];
      E.swap(E2);

      std::vector<double> g2(N);
      for(size_t i=0; i<sorter.size(); i++)g2[i]=g[sorter[i].second];
      g.swap(g2);

#ifdef DEBUG_MODE_SORTER
      std::cout<<"(E, g) after sorting:"<<std::endl;
      for(size_t i=0; i<sorter.size(); i++){
        std::cout<<"("<<E[i]<<", "<<g[i]<<"), ";
      }
      std::cout<<std::endl;
#endif
    }
    
  }


}//namespace
