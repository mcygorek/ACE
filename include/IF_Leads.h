#ifndef IF_LEADS_DEFINED_H
#define IF_LEADS_DEFINED_H

#include <algorithm>
#include "InfluenceFunctional_OD.h"
#include "SingleBathMode.h"

class IF_Leads: public InfluenceFunctional_OD{
public:
  //Energy [meV] and coupling [ps^{-1}] of k-th mode  
  std::vector<std::pair<double, double> > E_g;  
  double EFermi, temperature;


  int get_Nmodes()const{
    return E_g.size();
  }
  static double rate_from_g(double g, int N, double Ediff){
    return 2.*M_PI*Constants::hbar_in_meV_ps*g*g*(N-1)/Ediff;
  }
  static double g_from_rate(double rate, int N, double Ediff){
    return sqrt(rate/(2.*M_PI*Constants::hbar_in_meV_ps*(N-1)/Ediff));
  }

  static bool compare_abs(const std::pair<double,double> &p1,
                          const std::pair<double,double> &p2){
    return fabs(p1.first)>fabs(p2.first);
  } 
  void setup(double EF_, double T, int Nmod, double Emin, double Emax, double g_){
    EFermi=EF_;
    temperature=T;

    if(Nmod<0){ E_g.clear(); return;}
    if(Nmod==1){
      E_g.resize(1); 
      E_g[0].first=(Emax+Emin)/2.;
      E_g[0].second=g_;
      return;
    }
    
    E_g.resize(Nmod);
    double dE=(Emax-Emin)/(Nmod-1);
    for(int i=0; i<Nmod; i++){
      E_g[i].first=Emin+i*dE;
      E_g[i].second=g_;
    }
#ifdef IF_LEADS_NOSORT 
    sort(E_g.begin(), E_g.end(), compare_abs);
#endif
  }
  void calculate(int n_max, double dt, RankCompressor &compressor, double dict_zero=1e-12){
    bool printdim=false; 
#ifdef PRINT_MAX_DIM 
    printdim=true;
#endif

//    std::cout<<"NONE:"<<std::endl;
    InfluenceFunctional_OD::set_none(n_max, 2);    
    calculate_dict(dict_zero);
    reduce_to_dict();
//dict.set_default(2);
   
    Operators2x2 op;
    Eigen::MatrixXcd HB_diag=otimes(op.id(), op.ketbra(1,1));
    Eigen::MatrixXcd HB_base=otimes(op.ketbra(0,1), op.ketbra(1,0)) +
                             otimes(op.ketbra(1,0), op.ketbra(0,1));

    for(int i=0; i<get_Nmodes(); i++){
      std::cout<<"Lead mode "<<i<<"/"<<get_Nmodes()<<std::endl;
      if(printdim)std::cout<<"Maxdim: before combination: "<<get_max_dim()<<std::endl;

      Eigen::MatrixXcd HB = E_g[i].first*HB_diag
           +Constants::hbar_in_meV_ps*E_g[i].second*HB_base;
      Eigen::MatrixXcd bath_init=op.ketbra(1,1);


      SingleBathMode mode(2, 2, HB);
      InfluenceFunctional_OD IF2;
//      IF2.calculate_single_mode(n_max, dt, n_max*dt, mode, bath_init, factorization);
      IF2.calculate_single_mode(n_max*2, dt/2., n_max*dt, mode, bath_init, factorization);
      if(printdim)std::cout<<"Maxdim: new contribution: "<<IF2.get_max_dim()<<std::endl;
      IF2.calculate_dict(dict_zero);
      IF2.reduce_to_dict();
//      join_and_sweep(IF2, compressor, dict.get_keep_weight());
      join_and_sweep_halfdt(IF2, compressor, dict.get_keep_weight());



/*
      ModePropagator mprop(2,2,HB);
      join_and_sweep(mprop, 0, dt, bath_init, compressor, dict.get_keep_weight());
*/


      if(printdim)std::cout<<"Maxdim: after sweep: "<<get_max_dim()<<std::endl;
    }

    calculate_closures();
  }
  

  IF_Leads(double EF_, double T, int Nmod, double Emin, double Emax, double g_,
           int n_max, double dt, RankCompressor &compressor, 
           double dict_zero=1e-12, int factor=0){
    setup(EF_, T, Nmod, Emin, Emax, g_);
    factorization=factor;
    calculate(n_max, dt, compressor, dict_zero);
  }
  IF_Leads(){
  }

};

#endif
