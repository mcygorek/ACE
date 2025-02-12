#include "discreteFT.hpp"
#include "Simulation_Results.hpp"
#include <cmath>
#include <iostream>

namespace ACE{

void discreteFT(std::complex<double> *in, std::complex<double> *out, int N, int stride, int sign){
  if(N<2){ 
    out[0]=in[0];
  }else{
    discreteFT(in, out, N/2, 2*stride, sign);
    discreteFT(in+stride, out+N/2, N/2, 2*stride, sign);
    for(int k=0; k<N/2; k++){
      std::complex<double> p=out[k];
      std::complex<double> q=exp(std::complex<double>(0,sign*2.*M_PI*k/N))*out[k+N/2];
      out[k]=p+q;
      out[k+N/2]=p-q;
    }
  }
}

std::vector<std::complex<double> > discreteFT(
     const std::vector<std::complex<double> > & in, int sign){
  
  int M=ceil( log((double)in.size())/log(2.) ); 
  int N=pow(2, M);

//  std::cerr<<"Zero padding: "<<in.size()<<" -> "<<N<<std::endl;

  std::complex<double> *in_array=new std::complex<double>[N];

  for(size_t i=0; i<in.size(); i++){
    in_array[i]=in[i];
  }
  for(size_t i=in.size(); i<N; i++){
    in_array[i]=0.;
  }

  std::vector<std::complex<double> > out(N); 
  std::complex<double> *out_array=&out[0];
  discreteFT(in_array, out_array, N, 1, sign);

  delete[] in_array;
 
  return out; 
}

//Corrections for integrals  (cf. "Numerical Recipes" Chap. 13.9)
double FFT_trapezoidal_correction_W(double t){
  if(abs(t)>1e-6){
    return 2.*(1.-cos(t))/(t*t);
  }else{
    double t2=t*t; double t4=t2*t2; double t6=t4*t2;
    return 1.-t2/12.+t4/360.-t6/20160.;   
  }
}
std::complex<double> FFT_trapezoidal_correction_a0(double t) { 
  if(abs(t)>1e-6){ 
    double t2=t*t;
    return std::complex<double>(-(1.-cos(t))/t2, (t-sin(t))/t2);
  }else{
    double t2=t*t; double t4=t2*t2; double t6=t4*t2;
    return std::complex<double>(
            -1./2.+t2/24.-t4/720.+t6/40320.,
            t*(1./6.-t2/120.+t4/5040.-t6/362880.));
  }
}

Simulation_Results resultsFFT(const Simulation_Results & in, int col, int sign, int integrate_mode){

  if(in.size()<4){
    std::cerr<<"resultsFFT: in.size()="<<in.size()<<"<4!"<<std::endl;
    exit(1);
  }
  double t_min=in[0].first;
  double dt=in[1].first-in[0].first;

  double omega_range=(2.*M_PI)/dt;
  int N_old=in.size();
  int N_new=pow(2, ceil(log((double)N_old)/log(2.)) );
  double domega=omega_range/N_new;

  std::cout<<"resultsFFT: in.size()="<<(int)in.size();
  std::cout<<" N_new="<<N_new;
  std::cout<<" dt="<<dt;
  std::cout<<" domega="<<domega;
  std::cout<<" integrate_mode="<<integrate_mode;
//  std::cout<<" orig_omega_range="<<omega_range;
//  std::cout<<" my_omega_range="<<my_omega_range;
  std::cout<<std::endl;


  std::complex<double> *in_array=new std::complex<double>[N_new];
  for(size_t i=0; i<N_new; i++){
    if(i<in.size()){
      if(in[i].second.size()<col){
        std::cerr<<"resultsFFT: in["<<i<<"].second.size()="<<(int)in[i].second.size()<<"<col="<<col<<"!"<<std::endl;
        exit(1);
      }
      in_array[i]=in[i].second[col];
    }else{
      in_array[i]=0.;
    }
  }

  std::complex<double> *out_array=new std::complex<double>[N_new];
  discreteFT(in_array, out_array, N_new, 1, sign);

  std::complex<double> low_end[1]; low_end[0]=in_array[0];
  std::complex<double> high_end[1]; high_end[0]=in_array[N_old-1];
  delete[] in_array;

  Simulation_Results out(N_new);

  for(int i=0; i<N_new; i++){ 
    std::complex<double> c=out_array[i];
    double w = i*domega ;   //reorder: first negative, then positive freqs.

    int j=i+N_new/2;   
    if(i>=N_new/2){
      w=(i-N_new)*domega;
      j=i-N_new/2;
    }

    if(integrate_mode==1){

//      std::function<double(double)> W { [](double t){ 
//        if(abs(t)>1e-9){return 2.*(1.-cos(t))/(t*t);}
//        else{double t2=t*t; double t4=t2*t2; double t6=t4*t2;
//          return 1.-t2/12.+t4/360.-t6/20160.;   }
//      }}; 
//      std::function<std::complex<double>(double)> a0 { [](double t){ 
//        if(abs(t)>1e-9){ double t2=t*t;
//          return std::complex<double>(-(1.-cos(t))/t2, (t-sin(t))/t2);}
//        else{double t2=t*t; double t4=t2*t2; double t6=t4*t2;
//          return std::complex<double>(
//            -1./2.+t2/24.-t4/720.+t6/40320.,
//            1./6.-t2/120.+t4/5040.-t6/362880.);}
//      }}; 

    
      c*=FFT_trapezoidal_correction_W(w*dt);
      c+=FFT_trapezoidal_correction_a0(w*dt)*low_end[0];
      c+=exp(std::complex<double>(0,(sign*w*dt*N_old)))*
           (std::conj(FFT_trapezoidal_correction_a0(w*dt))*high_end[0]);

    }

    c*=dt*exp(std::complex<double>(0,(sign*w*t_min)));
    out[j].first=w;
    out[j].second=std::vector<std::complex<double> >(1, c);
  }

  delete[] out_array;

  return out;
}

}//namespace
