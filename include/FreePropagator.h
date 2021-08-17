#ifndef FREE_PROPAGATOR_DEFINED_H
#define FREE_PROPAGATOR_DEFINED_H

#include "Constants.h"
#include "Pulse.h"
#include "TimedepMatrix.h"
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>
#include <vector>
#include "Propagator.h"
#include "Parameters.h"
#include "Operators.h"
#include "ReadOperator.h"
#include "HilbertSpaceRotation.h"
#include "MultitimeOp.h"
#include "InitialState.h"

class FreePropagator: public Propagator{
public:
  /// Propagator in Liouville space. To be updated in each time step:
  //Eigen::MatrixXcd M;  <- from base class Propagator 


  //These are the necessary ingredients to compose M:

  ///time-independent part of Hamiltonian ( NxN matrix )
  Eigen::MatrixXcd const_H;

  ///time-dependent parts of Hamiltonian
  std::vector<TimedepMatrixPtr> timedep_H;

  /** N^2 x N^2 matrix describing the mapping:  
      d/dt rho = -i/hbar [H,rho] + nnonH {rho}  (as matrix-vector mult.) */
  Eigen::MatrixXcd nonH;

  ///apply Hilbert space rotation _after_ calculation:
  HilbertSpaceRotation hs_rot;

  ///for fast changing time-dependent Hamiltonian
  int Nintermediate;
  
  //Used for consistency checks , cf. set_dim;
  bool dim_fixed;

  //For time-independent dynamics: no need to recalculate propagator in every time step. Then, required: info whether or not the propagator was calculated before:
  bool precalculated;
  double precalculated_dt;
  bool never_update;


  ///Operators for multitime correlation functions
  std::vector<MultitimeOp> multitime_op;


  //Functions:

  virtual bool is_time_independent() const{ 
    return timedep_H.size()<1;
  }
  void set_Nintermediate(int N){ 
    Nintermediate=N; 
  }

  //get Dimension of Hamiltonian
  virtual int get_dim()const{
    return const_H.rows();
  }
  //Either set dimension with a corresponding (zero) Hamiltonian matrix, or make sure adding a new term is consistent with existing dimension
  int set_dim(int dim, const std::string &error_comment=""){
    if(!dim_fixed){
      const_H=Eigen::MatrixXcd::Zero(dim,dim);
      dim_fixed=true;
      return dim;
    }else{
      if(const_H.rows()!=dim || const_H.cols()!=dim){
        std::cerr<<"FreePropagator::set_dim: Mismatch with dimension of Hamiltonian!"<<std::endl;
        if(error_comment!=""){
          std::cerr<<"Comment: "<<error_comment<<std::endl;
        }
        exit(1);
      }
    
      return const_H.rows();
    }
  } 
  int set_dim(const Eigen::MatrixXcd M, const std::string &error_comment=""){
    set_dim(M.cols(),error_comment);
    return set_dim(M.rows(), error_comment);
  }
  

  ///sanity checks
  virtual void check_dimensions()const{
    if(const_H.rows()!=const_H.cols()){
      std::cerr<<"Hamiltonian::check_dimensions: const_H.rows()!=const_H.cols()!"<<std::endl;
      exit(1);
    }
    for(size_t i=0; i<timedep_H.size(); i++){
      if(timedep_H[i]->get_dim()!=const_H.cols()){
        std::cerr<<"Hamiltonian::check_dimensions: timedep_H[i].get_dim()!=const_H.cols() for i="<<i<<"!"<<std::endl;
        exit(1);
      }
    } 

    for(size_t i=0; i<multitime_op.size(); i++){
      if(multitime_op[i].op.rows()!=const_H.rows()*const_H.rows()){
        std::cerr<<"Hamiltonian::check_dimensions: multitime_op[i].op.rows()!=const_H.rows()*const_H.rows() for i="<<i<<"!"<<std::endl;
        exit(1);
      }
    }
  }

  virtual Eigen::MatrixXcd get_Htot(double t)const{
    //constant part:
    Eigen::MatrixXcd matH=const_H;
    //add time-dependent parts:
    for(size_t i=0; i<timedep_H.size(); i++){
      matH += timedep_H[i]->f(t);
    }
    return matH;
  }

  /* returns a matrix in Liouville space whose ".exp()" gives the time evolution operator \mathcal{M} in Liouville space */
  static Eigen::MatrixXcd Hamiltonian_to_Liouville_Generator(const Eigen::MatrixXcd &H, double dt){

    int dim=H.rows();
    std::complex<double> gamma(0.,-dt/Constants::hbar_in_meV_ps);

    Eigen::MatrixXcd gen=Eigen::MatrixXcd::Zero(dim*dim, dim*dim);

    for(int i0=0; i0<dim; i0++){
      for(int i1=0; i1<dim; i1++){
        for(int i2=0; i2<dim; i2++){
            gen(i0*dim+i1, i2*dim+i1) += gamma*H(i0,i2);
        }
      }
    }
    for(int i0=0; i0<dim; i0++){
      for(int i1=0; i1<dim; i1++){
        for(int i3=0; i3<dim; i3++){
          gen(i0*dim+i1, i0*dim+i3) += -gamma*H(i3,i1);
        }
      }
    } 
    return gen;
  }
  bool has_to_recalculate(double dt)const{
    if(precalculated && fabs(precalculated_dt-dt)<1e-16 && timedep_H.size() <1 || never_update){
      return false;
    }else{
      return true;
    }
  }
  virtual Eigen::MatrixXcd Total_Generator(double t, double dt){
    size_t dim=get_dim();
    Eigen::MatrixXcd gen=Eigen::MatrixXcd::Zero(dim*dim, dim*dim);
    if(nonH.rows()>1)gen=nonH*dt;
    gen+=Hamiltonian_to_Liouville_Generator(get_Htot(t),  dt);
    return gen; 
  }

  virtual void update_single(double t, double dt){
    Eigen::MatrixXcd gen=Total_Generator(t,dt);
    M=gen.exp();
    M=hs_rot.apply_Liouville(M);
  }

  virtual void update(double t, double dt){
#ifdef DEBUG
std::cout<<"DEBUG: FreePropagator: update called: t="<<t<<std::endl;
//if(precalculated)std::cout<<"DEBUG: precalculated=true";
//else std::cout<<"DEBUG: precalculated=false";
//std::cout<<" precalculated_dt-dt="<<precalculated_dt-dt<<" timedep_H.size()"<<std::endl;
#endif

    if(has_to_recalculate(dt)){
      if(Nintermediate<1){ 
        update_single(t, dt);
      }else{
        double dt2=dt/(Nintermediate+1);
        int NL=get_dim()*get_dim();
        Eigen::MatrixXcd N=Eigen::MatrixXcd::Identity(NL,NL);
        for(int i=0; i<Nintermediate+1; i++){
          update_single(t+i*dt2, dt2);
          N*=M;
        }
        M=N;
      }
      precalculated=true;
      precalculated_dt=dt;
    }else{
#ifdef DEBUG
std::cout<<"DEBUG: skipping calculation (precalculated)"<<t<<std::endl;
#endif
       //skip calculation
    }
 
    //apply multitime operators
    for(size_t o=0; o<multitime_op.size(); o++){
      if(multitime_op[o].is_now(t,dt)){
        multitime_op[o].apply(M);
        precalculated=false;
      }
    }
  }


  void set_Hamiltonian(const Eigen::MatrixXcd & H){
    const_H=H;
    precalculated=false;
    dim_fixed=true;
    check_dimensions();
  }
  void add_Hamiltonian(const Eigen::MatrixXcd & H){
    if(H.rows()!=H.cols()){
      std::cerr<<"FreePropagator::add_Hamiltonian: H not square!"<<std::endl;
      exit(1);
    }
    set_dim(H);
    const_H+=H;
    precalculated=false;
  }


  void add_Hamiltonian(ComplexFunctionPtr &f, const Eigen::MatrixXcd &A){
    if(A.rows()!=A.cols()){
        std::cerr<<"FreePropagator::add_Hamiltonian (time dependent): Operator not square!"<<std::endl;
        exit(1);
    }
    set_dim(A,"FreePropagator::add_Hamiltonian (time dependent): Mismatch in dimensions!");
    timedep_H.push_back(new TimedepMatrix_Const_Shape_Hermitian(f, A));
  }
  
  void add_TimedepMatrix(TimedepMatrixPtr & mp){
    set_dim(mp->get_dim(), "FreePropagator::add_Timedep: mp");
    timedep_H.push_back(mp);
  }

  void add_Lindblad(double gamma, const Eigen::MatrixXcd & L){
    if(nonH.rows()==0){
      nonH=Eigen::MatrixXcd::Zero(L.rows()*L.rows(), L.cols()*L.cols());
    }
    if(L.rows()!=L.cols()){
      std::cerr<<"FreePropagator::add_Lindblad: L.rows()!=L.cols()"<<std::endl;
      exit(1);
    }
   
    size_t dim=set_dim(L, "FreePropagator::add_Lindblad: L.rows()!=H.get_dim()");
    

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
    int dim=get_dim();
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
    set_dim(m1); set_dim(m2); 
    multitime_op.push_back(MultitimeOp(t, m1, m2));
  }

/*
  void printPulses(const std::string &filename, double ta, double dt, double te){
    if(timedep_H.size()<1)return;
    std::ofstream ofs(filename.c_str());
    int Nt=(int) ((double)(te-ta)/dt);
std::cout<<"Nt: "<<Nt<<" ta: "<<ta<<" te: "<<te<<" dt: "<<dt<<std::endl;
    for(int ti=0; ti<=Nt; ti++){
      double t=ta+ti*dt;
      ofs<<t;
      for(size_t i=0; i<timedep_H.size(); i++){
        std::complex<double> c=timedep_H[i].first->f(t);
        ofs<<" "<<c.real()<<" "<<c.imag();
      } 
      ofs<<std::endl;
    }
  }
*/
  void initialize(){
    dim_fixed=false;
    Nintermediate=0;
    precalculated=false;
    never_update=false;
  }
  void setup(Parameters &param){
    initialize();
    Nintermediate=param.get_as_size_t("Nintermediate", 0);

    //Special case: magnetic field of spin 1/2
    if(param.is_specified("TLS_B_eff")){
      Operators2x2 op2;
      set_dim(2);
      std::vector<double> dvec=param.get_row_doubles("TLS_B_eff", 0, 3);
      const_H+= dvec[0]*0.5*op2.sigma_x()
              + dvec[1]*0.5*op2.sigma_y()
              + dvec[2]*0.5*op2.sigma_z();
    }
  
    //Pulses:
    std::vector<std::vector<std::string> >  svv;


    if(param.is_specified("add_Hamiltonian")){
      for(int i=0; i<param.get_nr_rows("add_Hamiltonian"); i++){
        std::string str=param.get_as_single_string("add_Hamiltonian", i);
std::cout<<"add_Hamiltonian: "<<str<<std::endl;
        add_Hamiltonian(ReadExpression(str));
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

        add_MultitimeOp(time, M, M);
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

        add_MultitimeOp(time, M, Eigen::MatrixXcd::Identity(M.rows(), M.cols()));
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
          if(svv[r].size()<4){
            std::cerr<<"Not enough parameters: "<<errmsg<<std::endl;
            exit(1);
          }
          //defaults:
          if(svv[r].size()<5)svv[r].push_back("0.");
          if(svv[r].size()<6)svv[r].push_back("{|1><0|_2}");
       
          ComplexFunctionPtr pulse=new Pulse_Gauss(
            Reader::readDouble(svv[r][1],errmsg+"\nPulse: Gauss: center(time)"),
            Reader::readDouble(svv[r][2],errmsg+"\nPulse: Gauss: FWHM(time)"),
            Reader::readDouble(svv[r][3],errmsg+"\nPulse: Gauss: area(pi)"),
            Reader::readDouble(svv[r][4],errmsg+"\nPulse: Gauss: detuning(meV)"));
   
          Eigen::MatrixXcd M=ReadExpression(svv[r][5]);
          set_dim(M.rows(), "at 'add_Pulse'");
          add_Hamiltonian(pulse, (Constants::hbar_in_meV_ps/2.)*M);
        }else{
          std::cerr<<"Cannot recognize pulse type '"<<type<<"'!"<<std::endl;
          exit(1);
        }
      }
    } 

/*
    if(param.get_as_string("print_pulses")!=""){
      printPulses(param.get_as_string("print_pulses"), 
                  param.get_as_double("ta",0),
                  param.get_as_double("dt",0.1),
                  param.get_as_double("te",20) );
    }
*/

/*
    if(!const_H_was_set){
      Eigen::MatrixXcd rhoinit=InitialState(param);
      set_dim(rhoinit.rows());
    }
*/
  }

  FreePropagator(){
    initialize();
  }
  FreePropagator(int dim){
    initialize();
    set_dim(dim);
  }
  FreePropagator(Parameters &param){
    setup(param);
  }
};

#endif
