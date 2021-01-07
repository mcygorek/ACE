#ifndef FREE_PROPAGATOR_DEFINED_H
#define FREE_PROPAGATOR_DEFINED_H

#include "Constants.h"
#include "Function.h"
#include "Pulse.h"
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>
#include <vector>
#include "Propagator.h"
#include "Parameters.h"
#include "Operators.h"
#include "ReadOperator.h"
#include "MultitimeOp.h"

class FreePropagator: public Propagator{
public:
  /// Propagator in Liouville space. To be updated in each time step:
  //Eigen::MatrixXcd M;  <- from base class Propagator 

  //These are the necessary ingredients to compose M:

  ///time-independent part of Hamiltonian ( NxN matrix )
  Eigen::MatrixXcd const_H;

  ///time-dependent parts of Hamiltonian
  std::vector<std::pair<ComplexFunctionPtr, Eigen::MatrixXcd> > timedep_H;

  /** N^2 x N^2 matrix describing the mapping:  
      d/dt rho = -i/hbar [H,rho] + nnonH {rho}  (as matrix-vector mult.) */
  Eigen::MatrixXcd nonH;


  ///for fast changing time-dependent Hamiltonian
  int Nintermediate;
  
  bool precalculated;


  ///Operators for multitime correlation functions
  std::vector<MultitimeOp> multitime_op;


  //Functions:
  ///get Dimension of Hamiltonian
  virtual int get_dim()const{return const_H.rows();}
  virtual bool is_time_independent() const{return timedep_H.size()<1;}
  
  ///sanity checks
  virtual void check_dimensions()const{
    if(const_H.rows()!=const_H.cols()){
      std::cerr<<"Hamiltonian::check_dimensions: const_H.rows()!=const_H.cols()!"<<std::endl;
      exit(1);
    }
    for(size_t i=0; i<timedep_H.size(); i++){
      if(timedep_H[i].second.rows()!=const_H.cols()){
        std::cerr<<"Hamiltonian::check_dimensions: timedep_H[i].second.rows()!=const_H.cols() for i="<<i<<"!"<<std::endl;
        exit(1);
      }
      if(timedep_H[i].second.rows()!=timedep_H[i].second.cols()){
        std::cerr<<"Hamiltonian::check_dimensions: timedep_H[i].second.rows()!=timedep_H[i].second.cols() for i="<<i<<"!"<<std::endl;
        exit(1);
      }
    }
  }

  virtual Eigen::MatrixXcd get_Htot(double t)const{
    Eigen::MatrixXcd matH=const_H;
    for(size_t i=0; i<timedep_H.size(); i++){
      for(int j=0; j<matH.rows(); j++){
        matH(j,j)+=(timedep_H[i].first->f(t)).real() * timedep_H[i].second(j,j);
        for(int k=0; k<j; k++){
          matH(j,k)+=timedep_H[i].first->f(t) * timedep_H[i].second(j,k);
          matH(k,j)+=std::conj(timedep_H[i].first->f(t)) * timedep_H[i].second(k,j);
        }
      }
    }
    return matH;
  }

  /* returns a matrix in Liouville space whose ".exp()" gives the time evolution operator \mathcal{M} in Liouville space */
  static Eigen::MatrixXcd Hamiltonian_to_Liouville_Generator(const Eigen::MatrixXcd &H, double dt){

    int dim=H.rows();
    std::complex<double> gamma(0.,-dt/Constants::hbar_in_meV_ps);

    Eigen::MatrixXcd gen=Eigen::MatrixXcd::Zero(dim*dim, dim*dim);

    for(size_t i0=0; i0<dim; i0++){
      for(size_t i1=0; i1<dim; i1++){
        for(size_t i2=0; i2<dim; i2++){
            gen(i0*dim+i1, i2*dim+i1) += gamma*H(i0,i2);
        }
      }
    }
    for(size_t i0=0; i0<dim; i0++){
      for(size_t i1=0; i1<dim; i1++){
        for(size_t i3=0; i3<dim; i3++){
          gen(i0*dim+i1, i0*dim+i3) += -gamma*H(i3,i1);
        }
      }
    } 
    return gen;
  }
  void update_single(double t, double dt){
    size_t dim=get_dim();
    Eigen::MatrixXcd gen=Eigen::MatrixXcd::Zero(dim*dim, dim*dim);
    if(nonH.rows()>1)gen=nonH*dt;

    gen+=Hamiltonian_to_Liouville_Generator(get_Htot(t),  dt);
    
    M=gen.exp();

  }
  virtual void update(double t, double dt){
    if(precalculated && timedep_H.size() <1){
       //skip calculation
    }else{
      if(Nintermediate<1){ 
        update_single(t, dt);
      }
      double dt2=dt/(Nintermediate+1);
      int NL=get_dim()*get_dim();
      Eigen::MatrixXcd N=Eigen::MatrixXcd::Identity(NL,NL);
      for(size_t i=0; i<Nintermediate+1; i++){
        update_single(t+i*dt2, dt2);
        N*=M;
      }
      M=N;
      precalculated=true;
    }
 
    //apply multitime operators
    for(size_t o=0; o<multitime_op.size(); o++){
      if(multitime_op[o].is_now(t,dt)){
        multitime_op[o].apply(M);
        precalculated=false;
      }
    }
  }

  void set_Nintermediate(int N){ Nintermediate=N; }
  void set_dim(int dim){
    if(const_H.rows()<1){
      const_H=Eigen::MatrixXcd::Zero(dim,dim);
      return;
    }
    if(const_H.rows()!=dim || const_H.cols()!=dim){
      std::cerr<<"FreePropagator::set_dim: Mismatch with dimension of Hamiltonian!"<<std::endl;
      exit(1);
    }
  }
  void set_Hamiltonian(const Eigen::MatrixXcd & H){
    const_H=H;
    precalculated=false;
  }
  void add_Hamiltonian(const Eigen::MatrixXcd & H){
    if(const_H.rows()<1){
      const_H=H;
    }else{
      if(H.rows()!=const_H.rows() || H.cols()!=const_H.rows()){
        std::cerr<<"FreePropagator::add_Hamiltonian: Mismatch in dimensions!"<<std::endl;
        exit(1);
      }
      const_H+=H;
    } 
    precalculated=false;
  }
  void add_Hamiltonian(ComplexFunctionPtr &f, const Eigen::MatrixXcd &A){
    if(const_H.rows()<1)set_dim(A.rows());
    if(A.rows()!=const_H.rows()||A.cols()!=const_H.rows()){
        std::cerr<<"FreePropagator::add_Hamiltonian: Mismatch in dimensions!"<<std::endl;
        exit(1);
    }
    timedep_H.push_back(std::make_pair(f, A));
  }

  void add_Lindblad(double gamma, const Eigen::MatrixXcd & L){
    if(nonH.rows()==0){
      nonH=Eigen::MatrixXcd::Zero(L.rows()*L.rows(), L.cols()*L.cols());
    }
    size_t dim=get_dim();
    if(L.rows()!=L.cols()){
      std::cerr<<"FreePropagator_nonHamiltonian::add_Lindblad: L.rows()!=L.cols()"<<std::endl;
      exit(1);
    }
    if(L.rows()!=dim){
      std::cerr<<"FreePropagator_nonHamiltonian::add_Lindblad: L.rows()!=H.get_dim()"<<std::endl;
      exit(1);
    }

    double g=gamma;
    Eigen::MatrixXcd Ldag=L.adjoint();
    Eigen::MatrixXcd L2=Ldag*L;

    for(size_t i0=0; i0<dim; i0++){
      for(size_t i1=0; i1<dim; i1++){
        for(size_t i2=0; i2<dim; i2++){
          for(size_t i3=0; i3<dim; i3++){
            nonH(i0*dim+i1, i2*dim+i3) += g*L(i0,i2)*Ldag(i3,i1);
          }
        }
      }
    }
    for(size_t i0=0; i0<dim; i0++){
      for(size_t i1=0; i1<dim; i1++){
        for(size_t i2=0; i2<dim; i2++){
          nonH(i0*dim+i1, i2*dim+i1) += -0.5*g*L2(i0,i2);
        }
      }
    }
    for(size_t i0=0; i0<dim; i0++){
      for(size_t i1=0; i1<dim; i1++){
        for(size_t i3=0; i3<dim; i3++){
          nonH(i0*dim+i1, i0*dim+i3) += -0.5*g*L2(i3,i1);
        }
      }
    }
    precalculated=false;
  }
  void add_diagonal_loss(double gamma, int index){
    size_t dim=get_dim();
    if(index<0||index>=dim){
      std::cerr<<"FreePropagator_nonHamiltonian::add_diagonal_loss: index<0||index>=dim"<<std::endl;
      exit(1);
    }
    for(int k=0; k<dim; k++){
      nonH( index*dim+k , index*dim+k) -= gamma/2.;
    }
    for(int k=0; k<dim; k++){
      nonH( k*dim+index , k*dim+index) -= gamma/2.;
    }
    precalculated=false;
  }
  void add_diagonal_loss(const std::pair<double, int> & dloss){
    add_diagonal_loss(dloss.first, dloss.second);
  }
  void add_MultitimeOp(double t, Eigen::MatrixXcd m1, Eigen::MatrixXcd m2){
    multitime_op.push_back(MultitimeOp(t, m1, m2));
  }
  void printPulses(const std::string &filename, double ta, double dt, double te){
    if(timedep_H.size()<1)return;
    std::ofstream ofs(filename.c_str());
    int Nt=(int) ((double)(te-ta)/dt);
std::cout<<"Nt: "<<Nt<<" ta: "<<ta<<" te: "<<te<<" dt: "<<dt<<std::endl;
    for(size_t ti=0; ti<=Nt; ti++){
      double t=ta+ti*dt;
      ofs<<t;
      for(size_t i=0; i<timedep_H.size(); i++){
        std::complex<double> c=timedep_H[i].first->f(t);
        ofs<<" "<<c.real()<<" "<<c.imag();
      } 
      ofs<<std::endl;
    }
  }
  void initialize(){
    Nintermediate=0;
    precalculated=false;
  }
  void setup(Parameters &param){
    //Two-level systems
    bool was_set=false;
    Operators2x2 op2;
    const_H=Eigen::MatrixXcd::Zero(2,2);
    Nintermediate=param.get_as_size_t("Nintermediate", 0);

    if(param.is_specified("TLS_EX")){
      const_H+=param.get_as_double("TLS_EX",0.)*op2.ketbra(1,1);
      was_set=true;
    }
    if(param.is_specified("TLS_const_driving")){
      int nr=param.get_nr_rows("TLS_const_driving");
      for(int r=0; r<nr; r++){
        std::vector<double> dvec=param.get_row_doubles("TLS_const_driving", r, 1);
        if(dvec.size()<2){
          const_H+=dvec[0]*Constants::hbar_in_meV_ps*0.5*op2.sigma_x(); 
        }else{
          ComplexFunctionPtr pulse=new Pulse_Constant(dvec[0], dvec[1]);
          add_Hamiltonian(pulse, Constants::hbar_in_meV_ps*0.5*op2.sigma_x());
        }
      }
      was_set=true;
    }
    if(param.is_specified("TLS_B_eff")){
      std::vector<double> dvec=param.get_row_doubles("TLS_B_eff", 0, 3);
      const_H+= dvec[0]*0.5*op2.sigma_x()
              + dvec[1]*0.5*op2.sigma_y()
              + dvec[2]*0.5*op2.sigma_z();
      was_set=true;
    }
    if(param.is_specified("TLS_gamma_rad")){
      add_Lindblad(param.get_as_double("TLS_gamma_rad"), op2.ketbra(0,1));
      was_set=true;
    }
  
    //Pulses:
    std::vector<std::vector<std::string> >  svv;
    if(param.is_specified("TLS_pulse_gauss")){
      svv=param.get("TLS_pulse_gauss");
      for(size_t r=0; r<svv.size(); r++){
        ComplexFunctionPtr pulse=new Pulse_Gauss(svv[r]);
        add_Hamiltonian(pulse, 0.5*op2.sigma_x());
      }
      was_set=true;
    }
    if(param.is_specified("TLS_pulse_gauss_FWHM_E")){
      svv=param.get("TLS_pulse_gauss_FWHM_E");
      for(size_t r=0; r<svv.size(); r++){
        ComplexFunctionPtr pulse=new Pulse_Gauss(svv[r]);
        ((Pulse_Gauss &)pulse.ref()).FWHM=(Constants::hbar_in_meV_ps*8.*log(2.))/((Pulse_Gauss &)pulse.ref()).FWHM;
        add_Hamiltonian(pulse, 0.5*op2.sigma_x());
      }
      was_set=true;
    }
    if(param.is_specified("TLS_pulse_rect")){
      svv=param.get("TLS_pulse_rect");
      for(size_t r=0; r<svv.size(); r++){
        ComplexFunctionPtr pulse=new Pulse_Rect(svv[r]);
        add_Hamiltonian(pulse, 0.5*op2.sigma_x());
      }
      was_set=true;
    }

    //Multi-two-level-systems
    int MTLS=param.get_as_size_t("MTLS",0);
    int MDIM=loop_pow(2, MTLS);
    if(param.is_specified("MTLS")){
      if(was_set){
        std::cerr<<"Please check parameters: cannot use 'TLS_*' and 'MTLS_*' at the same time!"<<std::endl;
        exit(1); 
      }
      const_H=Eigen::MatrixXcd::Zero(MDIM,MDIM);
      was_set=true;
    }
    if(param.is_specified("MTLS_EX")){
      if(MTLS<1){std::cerr<<"Please specify MTLS!"<<std::endl;exit(1);}
      Parameters_Entry pe=param.get("MTLS_EX");
      for(size_t i=0; i<pe.size(); i++){
        if(pe[i].size()<2){
          std::cerr<<"Usage: MTLS_EX index(<MTLS) value!"<<std::endl;exit(1);
        }
        size_t index=Reader::readSizeT(pe[i][0],"MTLS_EX index", MTLS);
        double val=Reader::readDouble(pe[i][1],"MTLS_EX value");
  
        const_H+=ExpandSingleOp(loop_pow(2, index), val*op2.ketbra(1,1),
                                loop_pow(2, MTLS-1-index));
        was_set=true;
      }
    }
    if(param.is_specified("MTLS_const_driving")){
      if(MTLS<1){std::cerr<<"Please specify MTLS!"<<std::endl;exit(1);}
      Parameters_Entry pe=param.get("MTLS_const_driving");
      for(size_t i=0; i<pe.size(); i++){
        if(pe[i].size()<2){
          std::cerr<<"Usage: MTLS_const_driving index(<MTLS) value!"<<std::endl;exit(1);
        }
        double val=Reader::readDouble(pe[i][1],"MTLS_const_driving value");
        Eigen::MatrixXcd mat=val*Constants::hbar_in_meV_ps*0.5*op2.sigma_x();

        if(pe[i][0]=="*"){
          for(size_t i=0; i<MTLS; i++){
            const_H+=ExpandSingleOp(loop_pow(2, i), mat, loop_pow(2, MTLS-1-i));
            was_set=true;
          }
        }else{
          size_t index=Reader::readSizeT(pe[i][0],"MTLS_const_driving index", MTLS);
          const_H+=ExpandSingleOp(loop_pow(2, index), mat, loop_pow(2, MTLS-1-index));
          was_set=true;
        }
      }
    }

    if(param.is_specified("add_Hamiltonian")){
      for(int i=0; i<param.get_nr_rows("add_Hamiltonian"); i++){
        std::string str=param.get_as_single_string("add_Hamiltonian", i);
std::cout<<"add_Hamiltonian: "<<str<<std::endl;
        Eigen::MatrixXcd M=ReadExpression(str);
        if(was_set && const_H.rows()!=M.rows()){
          std::cerr<<"Error in 'add_Hamiltonian': wrong dimensions!"<<std::endl;
          exit(1);
        }
        if(was_set)const_H+=M; 
        else const_H=M;
        was_set=true;
      }

      std::cout<<"Hamiltonian: "<<std::endl<<const_H<<std::endl;
    }

    if(param.is_specified("add_Lindblad")){
      for(int i=0; i<param.get_nr_rows("add_Lindblad"); i++){
        std::string str=param.get_as_single_string("add_Lindblad", i);
std::cout<<"add_Lindblad: "<<str<<std::endl;
        std::stringstream ss(str);
        double rate; ss>>rate;
        if(ss.fail()){
          std::cerr<<"Error: add_Lindblad needs a rate argument!"<<std::endl;
          exit(1);
        }
        Eigen::MatrixXcd M=ReadExpression(ss.str().substr(ss.tellg()));
        add_Lindblad(rate, M);
      }
    }
    if(param.is_specified("apply_Operator")){
      for(int i=0; i<param.get_nr_rows("apply_Operator"); i++){
        std::string str=param.get_as_single_string("apply_Operator", i);
std::cout<<"apply_Operator: "<<str<<std::endl;
        std::stringstream ss(str);
        double time; ss>>time;
        if(ss.fail()){
          std::cerr<<"Error: apply_Operator needs a time argument!"<<std::endl;
          exit(1);
        }
std::cout<<"time: "<<time<<std::endl;
std::cout<<"'"<<ss.str()<<std::endl;
        Eigen::MatrixXcd M=ReadExpression(ss.str().substr(ss.tellg()));
std::cout<<"M: "<<std::endl<<M<<std::endl;

        multitime_op.push_back(MultitimeOp(time, M, M));
      }
    }
    if(param.is_specified("apply_Operator_left")){
      for(int i=0; i<param.get_nr_rows("apply_Operator_left"); i++){
        std::string str=param.get_as_single_string("apply_Operator_left", i);
std::cout<<"apply_Operator_left: "<<str<<std::endl;
        std::stringstream ss(str);
        double time; ss>>time;
        if(ss.fail()){
          std::cerr<<"Error: apply_Operator needs a time argument!"<<std::endl;
          exit(1);
        }
std::cout<<"time: "<<time<<std::endl;
std::cout<<"'"<<ss.str()<<std::endl;
        Eigen::MatrixXcd M=ReadExpression(ss.str().substr(ss.tellg()));
std::cout<<"M: "<<std::endl<<M<<std::endl;

        multitime_op.push_back(MultitimeOp(time, M, Eigen::MatrixXcd::Identity(M.rows(), M.cols())));
      }
    }
    //pulses with general operators
    if(param.is_specified("add_Pulse")){
      std::vector<std::vector<std::string> >  svv=param.get("add_Pulse");
      for(size_t r=0; r<svv.size(); r++){
        if(svv[r].size()<1)continue;
        std::string type=svv[r][0];

        if(type=="Gauss"){ //center, FWHM, area, detuning
          std::string errmsg="Usage: add_Pulse Gauss CENTER FWHM AREA DETUNING OPERATOR!";
          if(svv[r].size()<6){
            std::cerr<<"Not enough parameters: "<<errmsg<<std::endl;
            exit(1);
          }
          ComplexFunctionPtr pulse=new Pulse_Gauss(
            Reader::readDouble(svv[r][1],errmsg+"\nPulse: Gauss: center(time)"),
            Reader::readDouble(svv[r][2],errmsg+"\nPulse: Gauss: FWHM(time)"),
            Reader::readDouble(svv[r][3],errmsg+"\nPulse: Gauss: area(pi)"),
            Reader::readDouble(svv[r][4],errmsg+"\nPulse: Gauss: detuning(meV)"));
   
          Eigen::MatrixXcd M=ReadExpression(svv[r][5]);
          add_Hamiltonian(pulse, 0.5*M);
        }else{
          std::cerr<<"Cannot recognize pulse type '"<<type<<"'!"<<std::endl;
          exit(1);
        }
      }
    } 



  
    if(param.get_as_string("print_pulses")!=""){
      printPulses(param.get_as_string("print_pulses"), 
                  param.get_as_double("ta",0),
                  param.get_as_double("dt",0.1),
                  param.get_as_double("te",20) );
    }
  }
  FreePropagator(){
    initialize();
  }
  FreePropagator(Parameters &param){
    initialize();
    setup(param);
  }
};

#endif
