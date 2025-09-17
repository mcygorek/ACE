#ifndef ACE_PROCESS_TENSOR_BUFFER_DEFINED_H
#define ACE_PROCESS_TENSOR_BUFFER_DEFINED_H

#include "ProcessTensorForward.hpp"
#include "ProcessTensorElement.hpp"
#include "DummyException.hpp"
#include "TempFileName.hpp"
#include "ProcessTensorBufferSpec.hpp"
#include "DiagBB.hpp"
#include "TimeGrid.hpp"
#include "LiouvilleTools.hpp"
#include "ModePropagatorGenerator.hpp"
#include "ReaderBasics.hpp"
#include "TruncationLayout.hpp"
#include "PreloadHint.hpp"
#include <vector>
#include <string>
#include <thread>
#include <future>

namespace ACE{
/**
Buffers blocks of PT elements.

For now, we consider blocks with equal number Nel of elements.
Each block is stored in a single file, which is small enough to be 
stored a memory but large enough to justify the overhead (~ 100MB - 1GB)
There will be a header file organizing the blocks.
    
Within this concept, we realize several different options:

a) explicit header file + block files (see above)
b) is_temporary=true: like a) but files are deleted on destruction of object
c) single block file & fname_header="": fall back to single, monolithic PT file, which is read upon initialization and remains in memory until closed.
d) use_file=false: don't use any file, store everything in RAM

*/

class ProcessTensorBuffer: public ProcessTensorBufferSpec, public ProcessTensorForward{
public:
  struct ShiftExtend{   //For combination of PT with different lengths:
    int shift_second;   //Start of second PT is shifted wrt. first PT
    int truncate_at;    //Overall length (first PT will be expanded by 1s). 
                        //If < 0, length corresponds to shifted second PT 
    int sweep_more;     //for partial sweep combining two PTs: how many time steps extra are swept. <-1: all
    ShiftExtend() : shift_second(0), truncate_at(-1), sweep_more(-1){
    }
  };

  //from ProcessTensorForward: 
  //imported: int n; n_tot;
  virtual const ProcessTensorElement * current();


  std::vector<ProcessTensorElement> buffer;
  int current_block; //which block is currently loaded. -1 is none.
 
  //pre-loaded: note: this won't be written, so use it as if it was read-only
  std::vector<ProcessTensorElement> preload;
  int preload_block; 
  std::future<bool> preload_future;
  bool preload_lock;

  //similar construction for asynchronous writes
  std::vector<ProcessTensorElement> write_buffer;
  int write_buffer_block;
  std::future<bool> write_buffer_future;
  bool write_buffer_lock;

//basic functions:
  
  ProcessTensorElement & get(int n, PreloadHint hint=NoPreload);
  inline ProcessTensorElement & operator[](int n){ return get(n);}

  const ProcessTensorElement & peek(int n); //like get but does not read/write buffer from/to files. Uses 'preload' if necessary.

  inline int get_blocksize()const{ return blocksize; }
  inline int get_nr_blocks()const{ 
    if(!use_file || blocksize<0){ return 1;}
    else{ return (n_tot+(blocksize-1))/blocksize; }
  }

  inline void clear_buffer(){ 
    std::vector<ProcessTensorElement>().swap(buffer);
  }
  inline void clear_preload(){ 
    std::vector<ProcessTensorElement>().swap(preload);
  }
  inline void clear_write_buffer(){ 
    std::vector<ProcessTensorElement>().swap(write_buffer);
  }

  void check_buffer_bounds(int i)const;
  void check_preload_bounds(int i)const;
  void check_write_buffer_bounds(int i)const;
  void print_dims(std::ostream &ofs=std::cout, bool print_both=false);
  void check_consistency();
  void clear();
  void push_back(const ProcessTensorElement &templ);
  void create(int n, const ProcessTensorElement &templ=ProcessTensorElement());
  void append(int n, const ProcessTensorElement &templ=ProcessTensorElement());
  void resize(int n, const ProcessTensorElement &templ=ProcessTensorElement());
  
  virtual void dict_expand(const ReadPT_struct &readPT);

  void copy_content(ProcessTensorBuffer & other);
  void copy_content(const std::string &filename);
  void copy_content(const ReadPT_struct &readPT);
  void copy_read_only(ProcessTensorBuffer & other);
  void delete_files();

  void read_block(int bl, PreloadHint hint=NoPreload );
  void read_block_preload(int bl);

  static bool request_async_preload_fct(ProcessTensorBuffer *PTB, int bl);
  void request_async_preload(int bl, PreloadHint hint);
  void wait_preload();
  void wait_write();

  static bool write_release_block_async_fct(
                 std::vector<ProcessTensorElement> &wbuf, 
                 std::string filename, bool debug_);

  void write_release_block(int bl);
  static bool can_read(const std::string &filename);
  void read(const std::string &filename, bool ro=false); //construction from header file
  void read(const ReadPT_struct &readPT, bool ro=false);
  void read_header();
  void write_header();
  void write_all();

  void print_info(std::ostream &os=std::cout)const;

// operations
  void set_trivial(int n_max, int sysdim);
  void calculate_closures();
  void distribute_weights();

  void sweep_forward(const TruncatedSVD &trunc, int verbosity,
                     int range_start=0, int range_end=-1);

  void sweep_backward(const TruncatedSVD &trunc, int verbosity,
                     int range_start=0, int range_end=-1);

  void sweep_pair_forward(const TruncatedSVD &trunc, int verbosity);
  void sweep_pair_backward(const TruncatedSVD &trunc, int verbosity);

  void join_and_sweep_forward(ProcessTensorBuffer & PTB2, 
                     const TruncatedSVD &trunc, int verbosity,
                     ShiftExtend shift_extend=ShiftExtend(), 
                     bool alternate=false);

  //implements symmetric Trotter; PTB2 must have time steps dt/2
  void join_symmetric_and_sweep_forward(ProcessTensorBuffer & PTB2, 
                     const TruncatedSVD &trunc, int verbosity,
                     ShiftExtend shift_extend=ShiftExtend()); 

  void join_and_sweep_backward(ProcessTensorBuffer & PTB2, 
                     const TruncatedSVD &trunc, int verbosity,
                     ShiftExtend shift_extend=ShiftExtend());

//  void join_select_and_sweep_forward(ProcessTensorBuffer & PTB2, 
//                     const TruncatedSVD &trunc, int verbosity,
//                     ShiftExtend shift_extend=ShiftExtend());

  void join_select_and_sweep_backward(ProcessTensorBuffer & PTB2, 
                     const TruncatedSVD &trunc_select, 
                     const TruncatedSVD &trunc_bw, 
                     int verbosity,
                     ShiftExtend shift_extend=ShiftExtend());

 
  void expand_DiagBB(DiagBB &diagBB, double dict_zero=0.);
  void apply_HilbertSpaceRotation(const HilbertSpaceRotation &hs_rot, double dict_zero=0.);

  bool split_inner(ProcessTensorBuffer &second, int center);
  bool split_inner_and_sweep_fbf(ProcessTensorBuffer &second, const TruncatedSVD &trunc, int center, int verbosity);
  void additive_join(ProcessTensorBuffer & PTB2);
//  void additive_join_and_sweep_backward(ProcessTensorBuffer & PTB2, 
//                                const TruncatedSVD &trunc, int verbosity);

  void sweep_intermediate_or_final_start_forward(const TruncationLayout &trunc, int cur_line, int max_lines, int verbosity, int range_start=0, int range_end=-1);
  void sweep_intermediate_or_final_start_backward(const TruncationLayout &trunc, int cur_line, int max_lines, int verbosity, int range_start=0, int range_end=-1);

  void set_from_DiagBB_single_line(DiagBB &diagBB, double dt, int n);
  void set_from_DiagBB_single_line_auto(DiagBB &diagBB, double dt, int n_tot, int n_mem, const TruncatedSVD & trunc, bool print_info=false, bool reverse=false);

  void set_from_DiagBB_fw(DiagBB &diagBB, const TimeGrid &tgrid,
                           TruncationLayout trunc, 
                           double dict_zero, int verbosity,
                           int stop_at_row=-1);
  void set_from_DiagBB(DiagBB &diagBB, const TimeGrid &tgrid,
                           TruncationLayout trunc, 
                           double dict_zero, int verbosity,
                           int stop_at_row=-1);
  void set_from_DiagBB_select(DiagBB &diagBB, const TimeGrid &tgrid,
                           TruncationLayout trunc, 
                           double dict_zero, int verbosity,
                           int stop_at_row=-1);
  void set_from_DiagBB_log(DiagBB &diagBB, const TimeGrid &tgrid,
                           TruncationLayout trunc, 
                           double dict_zero, int verbosity, 
                           int stop_at_row=-1, int start_at_row=-1,
                           bool copy_compress_second=false);


  void set_from_ModePropagator(ModePropagator &mprop, const TimeGrid &tgrid, 
                                                           double dict_zero);

  void add_modes(ModePropagatorGenerator &mpg, const TimeGrid &tgrid, 
                 TruncationLayout trunc, 
                 double dict_zero, int verbosity);

  void add_modes_firstorder(ModePropagatorGenerator &mpg, const TimeGrid &tgrid, 
                 TruncationLayout trunc, 
                 double dict_zero, int verbosity);

  void add_modes_select(ModePropagatorGenerator &mpg, const TimeGrid &tgrid, 
                 TruncationLayout trunc, 
                 double dict_zero, int verbosity);

  void add_modes_tree_get(int level, int max_level, int first_elem, 
          ModePropagatorGenerator &mpg, const TimeGrid &tgrid,
          const TruncationLayout & trunc, 
          double dict_zero, int verbosity);

  void add_modes_tree(ModePropagatorGenerator &mpg, const TimeGrid &tgrid, 
          TruncationLayout trunc, 
          double dict_zero, int verbosity);

  void set_from_modes_tree(ModePropagatorGenerator &mpg, const TimeGrid &tgrid, 
          TruncationLayout trunc, 
          double dict_zero, int verbosity);

  void set_from_coarse_grain(ProcessTensorBuffer &PTB, int coarse_grain);
// initializers

  void clean_up();
  void initialize(const ProcessTensorBufferSpec & spec=ProcessTensorBufferSpec());

  void set_new_single_file(const std::string &fname);
  
  void set_new_file(const std::string &fname, int blocksize_);
  
  void set_new_temporary(int blocksize_);
  
  void set_new_temporary(const ProcessTensorBufferSpec &other);
 
  void set_preload_none();
  void set_write_buffer_none();

  void expand_outer(int Nfront, int Nback);

  ProcessTensorBuffer(const std::string &fname, bool ro=false){
    preload_lock=write_buffer_lock=false;
    read(fname, ro);
  }
  ProcessTensorBuffer(const ReadPT_struct &readPT, bool ro=false){
    preload_lock=write_buffer_lock=false;
    read(readPT, ro);
  }
  ProcessTensorBuffer(const ProcessTensorBufferSpec & spec=ProcessTensorBufferSpec()){
    preload_lock=write_buffer_lock=false;
//    preload_future=std::async( []()->bool{return true;});
    
    set_preload_none();
    set_write_buffer_none();
    initialize(spec);
  }
  ~ProcessTensorBuffer(){
    clean_up();
  }
};

}//namespace
#endif
