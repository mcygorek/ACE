#ifndef PROCESS_TENSOR_REAL_DEFINED_H
#define PROCESS_TENSOR_REAL_DEFINED_H

#include "IF_OD_Dictionary.hpp"
#include <vector>
#include <memory>
#include <Eigen/Dense>
#include "InfluenceFunctional_OD.hpp"
#include "HermitianLiouvilleBasis.hpp"
#include "Closure_Ops.hpp"
#include "Env_State_Filter.hpp"
#include "Trafo_Chain.hpp"
#include "LiouvilleTools.hpp"

namespace ACE{

class ProcessTensor_real: public MPS_real, public Sweep_Trafo_Processor_real{
public: 
  IF_OD_Dictionary dict;
  Eigen::MatrixXcd HLU;

  bool print_timesteps;
  Closure_Ops closure_ops;

  //closures:
  std::vector<Eigen::VectorXcd> c;
  //Environment operators:
  std::vector<std::vector<Eigen::VectorXcd> > env_ops;

  //Manual selection of environment states:
  std::shared_ptr<Env_State_Filter> env_state_filter;

  //Trafo_Chain:
  std::vector<std::pair<Trafo_Chain, int> > trafo_chains;

  inline virtual int get_sys_dim()const{
    return dict.get_N();
  }

  MPS_Matrix get_a_phys(int n)const;

  void trafo_d1(const Eigen::MatrixXd &T, int n);
  
  void trafo_d2(const Eigen::MatrixXd &T, const Eigen::MatrixXd &Tback, int n);

  static std::vector<std::vector<Eigen::MatrixXd> >  propA_to_real(
        const std::vector<std::vector<Eigen::MatrixXcd> > & propA, 
        const Eigen::MatrixXcd &HLBasisN, const Eigen::MatrixXcd &HLBasisM,
        double *res=NULL);
 
  void calculate_closures();

  void calculate_dict(double zero=1e-12, bool verbose=true);
  
  //keep only non-redundent terms with respect to outer indices
  void reduce_to_dict();

  //build full MPS matrices from reduced (dictionary)
  void expand_from_dict();

  void make_longer(int n_template, int how_many);

  void calculate_dt0_ndt0(ModePropagator &mprop, int n_max, double ta, double dt, double dt0, int ndt0);

  void calculate(ModePropagator &mprop, int n_max, double ta, double dt, 
                            const Closure_Ops &closure_ops_=Closure_Ops());
  
  void calculate(ModePropagator &mprop, const TimeGrid &tgrid, 
                            const Closure_Ops &closure_ops_=Closure_Ops());

  //second order combination of two MPS matrices
  void single_join_halfdt(
      MPS_Matrix_real &M, const MPS_Matrix_real &M1, const MPS_Matrix_real &M2, 
      const IF_OD_Dictionary &thisdict, const IF_OD_Dictionary &odict, 
      const IF_OD_Dictionary &ndict);

  void join_and_sweep_halfdt(ProcessTensor_real &other, 
                             RankCompressor_real &compressor,
                             bool do_sweep=true, double keep_weight=0.,
                             bool do_join=true);

  void sweep_low_high_low(RankCompressor_real &compressor, double keep_weight=0.);

  //Define tracking of env_ops with sweeps:
  virtual void process_low_to_high(int n, const Eigen::MatrixXd &R);
  virtual void process_high_to_low(int n, const Eigen::MatrixXd &L);

  void sweep(RankCompressor_real &compressor, bool do_sweep=true, double keep_weight=0.);
  
  InfluenceFunctional_OD get_IF_OD()const;

//read - write 
  void read_env_ops_binary(std::istream &ifs);
  virtual void read_binary(const std::string &filename);
  static bool is_proper_file_format(std::ifstream &ifs, bool complain=false);
  static bool is_proper_file_format(const std::string &fname, bool complain=false);
  
  void write_env_ops_binary(std::ostream &ofs)const;
  virtual void write_binary(const std::string &filename)const;

  void set_from_complex(const InfluenceFunctional_OD &IF);

//initializers
  void initialize(int N=2);
  
  void set_trivial(int N, int nmax, double dict_zero=0., bool verbose=false); 
  
  inline ProcessTensor_real(ModePropagator &mprop, const TimeGrid &tgrid, 
                            const Closure_Ops &closure_ops_=Closure_Ops()){
    initialize();
    calculate(mprop, tgrid, closure_ops_);
  }
  inline ProcessTensor_real(int N=2, int nmax=0, double dict_zero=0.){
    initialize(N);
    set_trivial(N, nmax, dict_zero);
  }
  inline virtual ~ProcessTensor_real(){}
};

}//namespace
#endif
