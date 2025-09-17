#ifndef FREE_PROPAGATOR_DEFINED_H
#define FREE_PROPAGATOR_DEFINED_H

#include <vector>
#include <iosfwd>
#include "Eigen_fwd.hpp"
#include "Propagator.hpp"
#include "Parameters.hpp"
#include "TimedepMatrix.hpp"
#include "HilbertSpaceRotation.hpp"
#include "MultitimeOp.hpp"

namespace ACE{
class Parameters;

class FreePropagator: public Propagator{
public:
  /// Propagator in Liouville space. To be updated in each time step:
  //Eigen::MatrixXcd M;  <- from base class Propagator 

  //These are the necessary ingredients to compose M:

  ///time-independent part of Hamiltonian ( NxN matrix )
  Eigen::MatrixXcd const_H;

  ///time-dependent parts of Hamiltonian
  std::vector<TimedepMatrixPtr> timedep_H;
  std::vector<TimedepMatrixPtr> timedep_H_forward;
  std::vector<TimedepMatrixPtr> timedep_H_backward;

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
  virtual int get_dim()const;
  
  //Either set dimension with a corresponding (zero) Hamiltonian matrix, or make sure adding a new term is consistent with existing dimension
  int set_dim(int dim, const std::string &error_comment="");
  int set_dim(const Eigen::MatrixXcd M, const std::string &error_comment="");

  ///sanity checks
  virtual void check_dimensions()const;

  virtual Eigen::MatrixXcd get_Htot(double t)const;
  virtual Eigen::MatrixXcd get_H_nonhermitian(double t, double dt=1e-5);
  virtual Eigen::MatrixXcd get_H_forward(double t)const;
  virtual Eigen::MatrixXcd get_H_backward(double t)const;

  /* returns a matrix in Liouville space whose ".exp()" gives the time evolution operator \mathcal{M} in Liouville space */
  static Eigen::MatrixXcd Hamiltonian_to_Liouville_Generator(const Eigen::MatrixXcd &H, double dt);
  static Eigen::MatrixXcd forward_Hamiltonian_to_Liouville_Generator(const Eigen::MatrixXcd &H, double dt);
  static Eigen::MatrixXcd backward_Hamiltonian_to_Liouville_Generator(const Eigen::MatrixXcd &H, double dt);
  
  bool has_to_recalculate(double dt)const;
  
  virtual Eigen::MatrixXcd Total_Generator(double t, double dt);

  virtual void update_single(double t, double dt);

  virtual void update(double t, double dt);

  void set_Hamiltonian(const Eigen::MatrixXcd & H);
  
  void add_Hamiltonian(const Eigen::MatrixXcd & H);

  void add_Pulse(ComplexFunctionPtr &f, const Eigen::MatrixXcd &A);
  void add_Pulse(const std::pair<std::vector<double>,std::vector<std::complex<double> > > & shape, const Eigen::MatrixXcd &A);

  void add_forward_Pulse(ComplexFunctionPtr &f, const Eigen::MatrixXcd &A);
  void add_forward_Pulse(const std::pair<std::vector<double>,std::vector<std::complex<double> > > & shape, const Eigen::MatrixXcd &A);

  void add_backward_Pulse(ComplexFunctionPtr &f, const Eigen::MatrixXcd &A);
  void add_backward_Pulse(const std::pair<std::vector<double>,std::vector<std::complex<double> > > & shape, const Eigen::MatrixXcd &A);
  
  void add_TimedepMatrix(TimedepMatrixPtr & mp);

  void add_Lindblad(double gamma, const Eigen::MatrixXcd & L, const Eigen::MatrixXcd & Ldag=Eigen::MatrixXcd());
  
  void add_diagonal_loss(double gamma, int index);
  
  inline void add_diagonal_loss(const std::pair<double, int> & dloss){
    add_diagonal_loss(dloss.first, dloss.second);
  }

  void add_MultitimeOp(double t, Eigen::MatrixXcd m1, Eigen::MatrixXcd m2, bool apply_before=false);

  void print_pulses(Parameters &param);

  inline void initialize(){
    dim_fixed=false;
    Nintermediate=0;
    precalculated=false;
    never_update=false;
    nonH=Eigen::MatrixXcd();
    const_H=Eigen::MatrixXcd();
    timedep_H=std::vector<TimedepMatrixPtr>();
    timedep_H_forward=std::vector<TimedepMatrixPtr>();
    timedep_H_backward=std::vector<TimedepMatrixPtr>();
  }
  void setup(Parameters &param, int setdim=-1);

  inline FreePropagator(){
    initialize();
  }
  inline FreePropagator(int dim){
    initialize();
    set_dim(dim);
  }
  inline FreePropagator(Parameters &param, int setdim=-1){
    setup(param, setdim);
  }
  inline FreePropagator(const std::string &filename, int setdim=-1){
    Parameters param(filename);
    setup(param, setdim);
  }
};

}//namespace
#endif
