#ifndef ACE_PROCESS_TENSOR_BUFFER_SPEC_DEFINED_H
#define ACE_PROCESS_TENSOR_BUFFER_SPEC_DEFINED_H

#include "ProcessTensorForward.hpp"
#include "ProcessTensorElement.hpp"
#include "DummyException.hpp"
#include "TempFileName.hpp"
#include <vector>
#include <string>

namespace ACE{
/**
Specifications for Buffers blocks of PT elements.

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
//specifies how a new buffer shall be created 
struct ProcessTensorBufferSpec{
  bool use_file; // write to files or keep everything in memory;
  bool is_temporary;
  bool use_single_file; //if only one block is used: Use one monolithic file instead of separate header and content files
  bool read_only; //suppress writing to files
  bool use_async_write; //switches on asynchronous writing (not much advantage foud so far)
  int blocksize;  // blocksize <0 encodes the use of a single file, 
                  //           =0 uninitialized
  std::string fname_header;        //name of header file

  std::string get_fname(int n);

  inline bool use_multiple_files()const{
    return use_file && !is_temporary && !use_single_file;
  }

  inline void copy(const ProcessTensorBufferSpec &other){
    use_file=other.use_file;
    is_temporary=other.is_temporary;
    use_single_file=other.use_single_file;
    read_only=other.read_only;
    use_async_write=other.use_async_write;
    blocksize=other.blocksize;
    fname_header=other.fname_header;
  }

  inline ProcessTensorBufferSpec clone_with_name(const std::string &fname)const{
    ProcessTensorBufferSpec new_spec(*this);
    new_spec.fname_header=fname;
    return new_spec;
  }
  inline ProcessTensorBufferSpec clone_read_only()const{
    ProcessTensorBufferSpec new_spec(*this);
    new_spec.read_only=true;
    return new_spec;
  }

  ProcessTensorBufferSpec & operator=(const ProcessTensorBufferSpec &other){
    copy(other);
    return *this;
  }
  ProcessTensorBufferSpec(const ProcessTensorBufferSpec &other){
    copy(other);
  }
  ProcessTensorBufferSpec() : use_file(false), is_temporary(false),
                              use_single_file(true), read_only(false), 
                              use_async_write(false),
                              blocksize(-1), fname_header(){
  }
  ProcessTensorBufferSpec(const std::string & fheader, int blocksize_=-1, bool temp=false)
   : use_file(true), is_temporary(temp), use_single_file(blocksize_<=0), 
     read_only(false), use_async_write(false), 
     blocksize(blocksize_), fname_header(fheader) {
  }
};


}//namespace
#endif
