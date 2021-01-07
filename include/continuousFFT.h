#include <fftw3.h> 
#include <complex>
#include <vector>
#include <cstdlib>
#include <iostream>

//Note: Convention: forward: multiplication with e^{-i omega t}
void DFT(const std::vector<std::complex<double> > &vin, 
               std::vector<std::complex<double> > &vout,
               bool forward=true, int zeropad=0){

  int N=vin.size()+zeropad;
  fftw_complex *in =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  fftw_complex *out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

  for(size_t i=0; i<vin.size(); i++){
//    std::complex<double> expfac=1.;exp(std::complex<double>(0,-i*wadt));
    std::complex<double> c=vin[i]; //expfac*vin[i];
    in[i][0]=c.real();
    in[i][1]=c.imag();
  }
  for(int z=0; z<zeropad; z++){
    in[vin.size()+z][0]=0.;
    in[vin.size()+z][1]=0.;
  }

  int sign=FFTW_FORWARD;
  if(!forward)sign=FFTW_BACKWARD;

  fftw_plan p = fftw_plan_dft_1d(N, in, out, sign, FFTW_ESTIMATE);
  fftw_execute(p); 
  fftw_destroy_plan(p);
  fftw_free(in); 
  
  vout.clear();
  vout.resize(N);
  for(size_t i=0; i<vout.size(); i++){
    vout[i]=std::complex<double>(out[i][0], out[i][1]);
  }
  fftw_free(out);
}

double continuousFFT_dw(double dt, int N){
  return (2.*M_PI)/(N*dt);  //w= n*dw; n=0,1,..., N/2-1
}
//wa, so that 0 of frequency is in the center
double continuousFFT_center_wa(double dt, int N){
  double dw=continuousFFT_dw(dt,N);
  double wa=-dw*N/2;
  if(N%2==1)wa=-dw*(N-1)/2;
  return wa;
}



/** 
vin: input data
vout: output data
ta: time corresponding to first input data point
dt: time intervals
N: number of discretization points for omega 
(larger N gives finer grid; vin.size() determines maximal omega)
wa: omega corresponding to first output data point
dw: get size of the omega intervals
sign: forward corresponds to integral with e^{-i omega t}
*/

void continuousFFT(const std::vector<std::complex<double> > &vin,
               std::vector<std::complex<double> > &vout,
               double ta, double dt, int N, double wa, double *dw_=NULL, 
               bool forward=true){

  if(vin.size()<2){
    std::cerr<<"continuousFFT with input data of size "<<vin.size()<<"!"<<std::endl;
    exit(1);
  }
  if(N<vin.size()){
    std::cerr<<"continuousFFT with N < input data size !"<<std::endl;
    exit(1);
  }

  int M=vin.size()-1;
  double dw=continuousFFT_dw(dt, N);
  if(dw_!=NULL)*dw_=dw;


  //shift input and perform DFT

  fftw_complex *in =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  fftw_complex *out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

  for(size_t i=0; i<vin.size(); i++){
    std::complex<double> c=exp(std::complex<double>(0., -wa*(ta+i*dt) ))*vin[i];
    in[i][0]=c.real();
    in[i][1]=c.imag();
  }
  for(int z=0; z<N-vin.size(); z++){
    in[vin.size()+z][0]=0.;
    in[vin.size()+z][1]=0.;
  }
  int sign=FFTW_FORWARD;
  if(!forward)sign=FFTW_BACKWARD;

  fftw_plan p = fftw_plan_dft_1d(N, in, out, sign, FFTW_ESTIMATE);
  fftw_execute(p); 
  fftw_destroy_plan(p);
  fftw_free(in); 
  
  vout.clear();
  vout.resize(N);
  for(size_t i=0; i<vout.size(); i++){
    vout[i]=std::complex<double>(out[i][0], out[i][1]);
  }
  fftw_free(out);


/*
  for(int i=0; i<vin.size(); i++){
    vin[i]*=exp(std::complex<double>(0., -wa*(ta+i*dt) ));
  }
  DFT(vin, vout, forward, N-vin.size()); //, wa*dt);
*/


  std::complex<double> dexp=exp(std::complex<double>(0., -dw*ta));
  std::complex<double> expfac=1.;
  std::complex<double> dexpfin=exp(std::complex<double>(0., -dw*M*dt));
  std::complex<double> expfin=1;
  for(size_t i=0; i<vout.size(); i++){
    //setting up correction factors:
    double theta=i*dw*dt;
    double W;
    std::complex<double> a0;
    if(theta>1e-4){
      W=2.*(1.-cos(theta))/(theta*theta);
      a0=std::complex<double>( (-1.+cos(theta))/(theta*theta), 
                            -  (theta-sin(theta))/(theta*theta));
    }else{
      double t2=theta*theta;
      W=1.+t2*(-1./12.+t2*(1./360.-t2/20160.));
      a0=std::complex<double>( 
          -1./2.+t2*(1./24.+t2*(-1./720.+t2/40320)),
        -  theta*(1./6.+t2*(-1./120.+t2*(1./5040.-t2/362880.))));
    }
    //a0=0; expfac=1.; W=1.;
    
//    expfac=exp(std::complex<double>(0., -i*dw*ta));
    //applying corrections
    vout[i]=dt*expfac*(W*vout[i] + a0*vin[0] +expfin*std::conj(a0)*vin.back());

    expfac*=dexp;
    expfin*=dexpfin;
  }

}

void continuousFFT(const std::vector<std::pair<double, std::vector<std::complex<double> > > > &vin,
               int col,
               std::vector<std::complex<double> > &vout,
               double ta, double dt, int N, double wa, double *dw_=NULL, 
               bool forward=true){

  if(vin.size()<1){
    std::cerr<<"continuousFFT vin.size()<1!"<<std::endl;
    exit(1);
  }
  std::vector<std::complex<double> > v(vin.size());
  for(size_t i=0; i<vin.size(); i++){
    if(vin[i].second.size()<=col){
      std::cerr<<"continuousFFT vin["<<i<<"].second.size() <= col!"<<std::endl;
      exit(1);
    }
    v[i]=vin[i].second[col];
  }
  continuousFFT(v, vout, ta, dt, N, wa, dw_, forward);
}
