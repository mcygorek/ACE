#ifndef ACE_PROCESS_TENSOR_STREAM_DEFINED_H
#define ACE_PROCESS_TENSOR_STREAM_DEFINED_H

#include "ProcessTensorStream_ro.hpp"
#include "ProcessTensorStream_wo.hpp"
#include "ModePropagatorGenerator.hpp"
#include "DiagBB.hpp"
#include "TimeGrid.hpp"

namespace ACE{
 namespace ProcessTensorStream{
  //contains information about how to shift and extend PTs (in time) before combination
  struct ShiftExtend{   //For combination of PT with different lengths:
    int shift_second;   //Start of second PT is shifted wrt. first PT
    int truncate_at;    //Overall length (first PT will be expanded by 1s). 
                        //If < 0, length corresponds to shifted second PT 
    ShiftExtend() : shift_second(0), truncate_at(-1){
    }
  };
  

  int get_length(const std::string & file);
 
  void check_consistency(const std::string & file);

  void set_trivial(const std::string & file_out, int n_max, int sysdim);

  void copy(const std::string & file_out, const std::string & file_in);

  void calculate_closures(const std::string & file_out, const std::string & file_in);

  void sweep_forward(const std::string & file_out, 
                     const std::string & file_in, 
                     const TruncatedSVD &trunc, int verbosity,
                     int range_start=0, int range_end=-1);

  void sweep_backward(const std::string & file_out, 
                     const std::string & file_in, 
                     const TruncatedSVD &trunc, int verbosity,
                     int range_start=0, int range_end=-1);

  void join_select_and_sweep_backward(const std::string & file_out,
               const std::string & file_in1, const std::string & file_in2,
               const TruncatedSVD &trunc, int verbosity, 
               ShiftExtend shift_extend=ShiftExtend());

  void join_select_and_sweep_forward(const std::string & file_out,
               const std::string & file_in1, const std::string & file_in2,
               const TruncatedSVD &trunc, int verbosity, 
               ShiftExtend shift_extend=ShiftExtend());


  void set_from_DiagBB_single_line(const std::string & file_out,
                    DiagBB &diagBB, double dt, int n);

  void set_from_DiagBB(const std::string & file_out,
                    DiagBB &diagBB, const TimeGrid &tgrid, 
                    const TruncatedSVD & trunc, int intermediate_sweep_n,
                    double dict_zero, int verbosity);

  void set_from_DiagBB_reverse(const std::string & file_out,
                    DiagBB &diagBB, const TimeGrid &tgrid, 
                    const TruncatedSVD & trunc, int intermediate_sweep_n,
                    double dict_zero, int verbosity);

  void set_from_DiagBB_log(const std::string & file_out,
                    DiagBB &diagBB, const TimeGrid &tgrid, 
                    TruncatedSVD trunc, int intermediate_sweep_n,
                    double dict_zero, int verbosity);

  void set_from_DiagBB_log_reverse(const std::string & file_out,
                    DiagBB &diagBB, const TimeGrid &tgrid, 
                    const TruncatedSVD & trunc, int intermediate_sweep_n,
                    double dict_zero, int verbosity);


  void set_from_ModePropagator(const std::string & file_out, 
         ModePropagator &mprop, const TimeGrid &tgrid, double dict_zero);

  void add_modes(const std::string & file_out, 
          ModePropagatorGenerator &mpg, const TimeGrid &tgrid, 
          const TruncatedSVD & trunc, int intermediate_sweep_n,
          double dict_zero, int verbosity);

  void add_modes_reverse(const std::string & file_out, 
          ModePropagatorGenerator &mpg, const TimeGrid &tgrid, 
          const TruncatedSVD & trunc, int intermediate_sweep_n,
          double dict_zero, int verbosity);


   void add_modes_tree_get(int level, int first_elem, 
          const std::string & file_out,
          ModePropagatorGenerator &mpg, const TimeGrid &tgrid,
          const TruncatedSVD & trunc, int intermediate_sweep_n,
          double dict_zero, int verbosity);

  void add_modes_tree(const std::string & file_out, 
          ModePropagatorGenerator &mpg, const TimeGrid &tgrid, 
          const TruncatedSVD & trunc, int intermediate_sweep_n,
          double dict_zero, int verbosity);
 

  void add_modes_tree_get_reverse(int level, int first_elem, 
          const std::string & file_out,
          ModePropagatorGenerator &mpg, const TimeGrid &tgrid,
          const TruncatedSVD & trunc, int intermediate_sweep_n,
          double dict_zero, int verbosity);

  void add_modes_tree_reverse(const std::string & file_out, 
          ModePropagatorGenerator &mpg, const TimeGrid &tgrid, 
          const TruncatedSVD & trunc, int intermediate_sweep_n,
          double dict_zero, int verbosity);

  std::string dims(const std::string & file);

 }
}
#endif
