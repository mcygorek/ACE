#ifndef INFLUENCE_FUNCTIONAL_OD_DEFINED_H
#define INFLUENCE_FUNCTIONAL_OD_DEFINED_H

#include "InfluenceFunctional_D.hpp"
#include "IF_OD_Abstract.hpp"
#include "LiouvilleTools.hpp"
#include "Rep_GIF.hpp"
#include "Compress_Trafo_At.hpp"
#include "Sweep_Trafo_Processor.hpp"
#include "ModePropagatorGenerator.hpp"
#include "RankCompressorList.hpp"
#include "IF_from_PT.hpp"

namespace ACE{
class SingleBathMode;
class Rep_GIF;

class InfluenceFunctional_OD: public IF_OD_Abstract, public MPS, public Sweep_Trafo_Processor{
public:

  //closures:
  std::vector<Eigen::VectorXcd> c;
  //Environment operators:
  std::vector<std::vector<Eigen::VectorXcd> > env_ops;

  
  int factorization;  //0: full exp; 1: e^{diag}(1-i/hbar Offdiag)
  bool print_timesteps;
  bool printdim;
 
  //for calculation of repetitive unit
  Rep_GIF rep; 
  bool compress_trafo_use_ortho;
  Compress_Trafo_At cta;
  

  //Implementation of IF_OD_Abstract
  virtual const MPS_Matrix & get_a(int n)const;

  virtual const Eigen::VectorXcd & get_c(int n)const;
  
  virtual const std::vector<Eigen::VectorXcd> &get_env_ops(int n)const;

  virtual void check_within_limits(int n)const;

  void check_env_dims(int n=-1)const;
  
  void check_rank_is(int n_max, const std::string &context);
  
  //Define tracking of env_ops with sweeps:
  virtual void process_low_to_high(int n, const Eigen::MatrixXcd &R);

  virtual void process_high_to_low(int n, const Eigen::MatrixXcd &L);


  //Genuine IF_OD functions:
  void calculate_closures();

  void calculate_dict(double zero=1e-12, bool verbose=true);
  
  //keep only non-redundent terms with respect to outer indices
  void reduce_to_dict();

  //build full MPS matrices from reduced (dictionary)
  void expand_from_dict();

  //expand from diagonal IF
  void expand_from_diagonal(const InfluenceFunctional_D &other, const HilbertSpaceRotation &hs_rot, double dict_zero=1e-20);

  void truncate(int n_tot);

  void calculate_diagBB(DiagBB &diagBB, RankCompressor &compressor, double dict_zero=1e-20, int n_extra=0);

  virtual void calculate_single_mode(int n_max, double dt, double t_tot, 
      SingleBathMode &mode, const Eigen::MatrixXcd &bath_init, int factor=0);

  void join_diagonal(const InfluenceFunctional_D &other);

  void calculate_dt0_ndt0(ModePropagator &mprop, int n_max, double ta, double dt, double dt0, int ndt0);

  inline void calculate(ModePropagator &mprop, int n_max, double ta, double dt){
    calculate_dt0_ndt0(mprop, n_max, ta, dt, dt, 0);
  }

  //second order combination of two MPS matrices
  void single_join_halfdt(
      MPS_Matrix &M, const MPS_Matrix &M1, const MPS_Matrix &M2, 
      const IF_OD_Dictionary &thisdict, const IF_OD_Dictionary &odict, 
      const IF_OD_Dictionary &ndict);

  //join with other IF_OD under the assumption that the other IF_OD is discretized with time steps dt/2  ->  Symmetric Trotter
  void join_and_sweep_halfdt(const InfluenceFunctional_OD &other, RankCompressor &compressor, double keep_weight=0., bool do_sweep=true);

  void join_halfdt(const InfluenceFunctional_OD &other);

  void normal_form_fwd();

  void print_closures(const std::string &fname);

  void read_env_ops_binary(std::istream &ifs);
  
  virtual void read_binary(const std::string &filename);
  
  void write_env_ops_binary(std::ostream &ofs)const;
  
  virtual void write_binary(const std::string &filename)const;

  void set_none(int n_max, int N, double dict_zero=1e-20);
  
  inline void set_none_from_tgrid(int N, double dict_zero=1e-20){
    set_none(tgrid.n_calc, N, dict_zero);
  }


  virtual void sweep(RankCompressor &compressor, double keep_weight);

  void add_IF_halfdt(InfluenceFunctional_OD &IF2, RankCompressor &compressor);
    
  void add_mode(ModePropagatorGenerator &mpg, int k, RankCompressor &compressor, double dict_zero=1e-12);

  void add_modes(ModePropagatorGenerator &mpg, RankCompressor &compressor, double dict_zero=1e-12);

  void calculate(ModePropagatorGenerator &mpg, RankCompressor &compressor, double dict_zero=1e-12);

  inline void initialize(){
    factorization=0;
    print_timesteps=false;
    compress_trafo_use_ortho=false;
    printdim=true;
  }

  void from_Rep_GIF(const Rep_GIF &rep_, int n_max_t);

  void insert_rep();
    
  InfluenceFunctional_OD(const TimeGrid &tgr_, DiagBB &diagBB, RankCompressor &compressor, double dict_zero=1e-20, int n_extra=0);
  
  InfluenceFunctional_OD(int n_max, int N);
  
  InfluenceFunctional_OD(const TimeGrid &tgr_, int N=2);
  
  InfluenceFunctional_OD(ModePropagator &mprop, int n_max, double ta, double dt, double dict_zero=1e-12);

  InfluenceFunctional_OD(ModePropagator &mprop, const TimeGrid &tgr_, double dict_zero=1e-12);

  InfluenceFunctional_OD(ModePropagatorGenerator &mpg, const TimeGrid &tgr_, RankCompressor &compressor, double dict_zero=1e-12);
  
  InfluenceFunctional_OD(const std::string &filename);

  InfluenceFunctional_OD();
  
  inline virtual ~InfluenceFunctional_OD(){}
};

}//namespace
#endif
