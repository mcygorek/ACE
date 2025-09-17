#ifndef ACE_TRUNCATION_LAYOUT_DEFINED_H
#define ACE_TRUNCATION_LAYOUT_DEFINED_H

#include "TruncatedSVD.hpp"
#include "Parameters.hpp"

namespace ACE{

/* New scheme:

While a TruncatedSVD stores essentially just the threshold and maximal number of SVDs (maxk) and does the corresponding truncation, the TruncationLayout provides such TruncatedSVDs.

The concept is:

- there is a base threshold and base maximal number 
  (Parameter names "threshold" and "maxk")

- for a forward sweep, we use the effective threshold:
  thr = threshold * r_forward
  (Parameter "forward_threshold_ratio")

- for a backward sweep, we use the effective threshold:
  thr = threshold * r_backward
  (Parameter "backward_threshold_ratio")

- when using the preselection scheme (see [Phys. Rev. X 14, 011010 (2024)]),e.g.
  thr = threshold * r_{for/back}ward * r_select
  (Parameter "select_threshold_ratio")

- sometimes, we want to ramp up the threshold, i.e., use a smaller threshold for
  combining the first lines and have the last compression be the base threshold
  (to reduce error accumulation and hence also inner bond dimensions). We use:
  thr = threshold * r_{for/back}ward * r_select * r_line,
  where r_line exponentially interpolates:
  r_line(current_line, max_lines) 
     = exp(-ln(range_factor)*(max_lines-1-current_line)/(max_lines-1)
  e.g: r_line(0, 10) = 1/range_factor
       r_line(9, 10) = 1
  (Parameter "threshold_range_factor")

- setting "intermediate_sweep_n" to values > 0 triggers additional line sweeps
  after a combination step. Convention is to use the same TruncatedSVD as
  in the previous sweep.

- setting "final_sweep_n" to values > 0 triggers multiple sweeps after
  the final PT-MPO has been calculated to "clean up" additional degrees of
  freedom. The value can be integer or half of an integer (e.g., 2.5). 
  The latter is used when we want to achieve a certain canonical form, 
  e.g., ending specifically with either a forward sweep or a backward sweep.
  Final sweeps are done with the base threshold and maxk except for the very
  last (single direction) sweep, where we use the values set by 
  "final_threshold" and "final_maxk"

*/

template <typename T> class TruncationLayout_T{
public:

  double base_threshold;
  int base_maxk;
  int base_mink;
  double keep;  //scaling factor: keep <= 0: keep largest singular value.
  double base_Tikhonov;
  
  double forward_threshold_ratio;
  double backward_threshold_ratio;
  double select_threshold_ratio;
  double threshold_range_factor;

  bool use_QR;
 
  int intermediate_sweep_n;
  int final_sweep_n;
  bool final_sweep_half;
  double final_sweep_threshold;
  int final_sweep_maxk;

  //setup at initialization

  inline double get_line_factor(double current_line, int max_lines)const{
    bool debug=false;
if(debug){std::cout<<"DEBUG: current_line="<<current_line<<" max_lines="<<max_lines<<std::endl;}

    if(fabs(threshold_range_factor-1.)<1e-6)return 1.;  // no range

    double lin_fac=1;
    if(max_lines>1){
      lin_fac=1.-current_line/(max_lines-1.);
    }
    if(current_line<0 || max_lines<0){ //negative arguments -> get base 
      return 1.;
    }

    if(abs(lin_fac-1.)<1e-6){
if(debug){std::cout<<1./threshold_range_factor<<std::endl;}
      return 1./threshold_range_factor;

    }else if(lin_fac<0. || lin_fac>1.){ //stop weird things from happening when out of range
if(debug){std::cout<<1.<<std::endl;}
      return 1.;

    }else{
if(debug){std::cout<<exp(-log(threshold_range_factor)*lin_fac)<<std::endl;}
//      return exp(-log(threshold_range_factor)*lin_fac);
      return pow(1./threshold_range_factor, lin_fac);
    }
  }
  
  inline TruncatedSVD_T<T> get_base()const{
    return TruncatedSVD_T<T>(base_threshold, base_maxk, base_mink, keep, false, base_Tikhonov);
  }
  inline TruncatedSVD_T<T> get_QR()const{
    return TruncatedSVD_T<T>(0, 0, 0, keep, true);
  }
  inline TruncatedSVD_T<T> get_base_line(double current_line, int max_lines)const{
    double thr=base_threshold*get_line_factor(current_line,max_lines);
    return TruncatedSVD_T<T>(thr, base_maxk, base_mink, keep, false, base_Tikhonov);
  }
  inline TruncatedSVD_T<T> get_forward(double current_line, int max_lines, bool select=false )const{
    double thr=base_threshold*forward_threshold_ratio;
    thr*=get_line_factor(current_line,max_lines);
    if(select){ thr*=select_threshold_ratio; }
    return TruncatedSVD_T<T>(thr, base_maxk, base_mink, keep, false, base_Tikhonov);
  }
  inline TruncatedSVD_T<T> get_backward(double current_line, int max_lines, bool select=false )const{
    double thr=base_threshold*backward_threshold_ratio;
    thr*=get_line_factor(current_line,max_lines);
    if(select){ thr*=select_threshold_ratio; }
    return TruncatedSVD_T<T>(thr, base_maxk, base_mink, keep, false, base_Tikhonov);
  }
  inline TruncatedSVD_T<T> get_final_sweep()const{
    return TruncatedSVD_T<T>(final_sweep_threshold, final_sweep_maxk, base_mink, keep, false, base_Tikhonov);
  }
  inline int get_intermediate_sweep_n()const{
    return intermediate_sweep_n;
  }
  inline int get_final_sweep_n()const{
    return final_sweep_n;
  }
  inline bool get_final_sweep_half()const{
    return final_sweep_half;
  }

  virtual void print_info(std::ostream &ofs=std::cout)const;

  void setup(Parameters &param);
  
  TruncationLayout_T(Parameters &param){
    setup(param);
  }
  TruncationLayout_T(){
    Parameters param;
    setup(param);
  }
};

typedef TruncationLayout_T<std::complex<double> > TruncationLayout;
typedef TruncationLayout_T<double> TruncationLayout_real;

}//namespace
#endif
