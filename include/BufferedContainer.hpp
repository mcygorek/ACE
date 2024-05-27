#ifndef ACE_BUFFERED_CONTAINER_DEFINED_H
#define ACE_BUFFERED_CONTAINER_DEFINED_H
#include <string>
#include <vector>
#include "BufferedElement.hpp"
#include "PreloadHint.hpp"
#include <future>
#include <iostream>

namespace ACE{

/** 
A std::vector-like container that supports storage of blocks of equal size 
in files. Supported is also buffered pre-loading from files.

class T should implement BufferedElement
*/

template <typename T> class BufferedContainer{
public:
  int blocksize;  // blocksize <0 encodes the use of a single block, which is always kept in memory
  enum class ON_EXIT {WriteOnDestruction, DeleteOnDestruction} on_exit;
  bool read_only; //suppress writing to files
  bool was_modified;  // ...between write events. Don't need to write block before loading new one when block was unchanged.
  std::string fname_header;  //file name
  std::string magicString;

  int n_tot;
  std::vector<T> buffer;
  int current_block; //which block is currently loaded. -1 is none.

  //pre-loaded: note: this won't be written, so use it as if it was read-only
  std::vector<T> preload;
  int preload_block;
  std::future<bool> preload_future;
  bool preload_lock;

//basic functions
  inline size_t size()const{return n_tot;}
  inline int get_blocksize()const{ return blocksize; }
  inline int get_nr_blocks()const{
    if(blocksize<=0){ return 1;}
    else{ return (n_tot+(blocksize-1))/blocksize; }
  }
  inline void clear_buffer(){std::vector<T>().swap(buffer);}
  inline void clear_preload(){std::vector<T>().swap(preload);}
  void set_preload_none();

  void check_buffer_bounds(int i)const;
  void check_preload_bounds(int i)const;

  std::string get_fname(int n);

//accessors
  T & get(int n, PreloadHint hint=NoPreload);
  const T & get_ro(int n, PreloadHint hint=NoPreload);
  inline T & operator[](int n){ return get(n);}
  const T & peek(int n); //like get but does not read/write buffer from/to files. Uses 'preload' if necessary.

  void clear();
  void push_back(const T &templ);
  void set_read_only(bool ro);

  void read_block_from_file(int bl, std::vector<T> & buf);
  void read_block(int bl, PreloadHint hint=NoPreload );
  void read_block_preload(int bl);

  static bool request_async_preload_fct(BufferedContainer<T> *PTB, int bl);
  void request_async_preload(int bl, PreloadHint hint);
  void wait_preload();

  void write_release_block(int bl=-1);
  void read(const std::string &filename, const std::string &magic="", bool ro=false); //construction from header file
  void read_header();
  void write_header();
  void write_all();

  void copy_content(BufferedContainer<T> &other);

  void print_info(std::ostream &os=std::cout)const;
  void resize(int sz);
  void initialize(const std::string &fname_h="", int bs=-1);
  void delete_files();
  void clean_up();

/* NOTE: should reserve BufferedContainer(string) for reading from file
  BufferedContainer(const std::string &fname_h="", int bs=-1, int ntot=0) : preload_lock(false){
    initialize(fname_h, bs);
    if(ntot>0)resize(ntot);
  }
*/
 BufferedContainer() : preload_lock(false){
   initialize();
 }
 BufferedContainer(const std::string &fname_h, const std::string & magic="") 
  : preload_lock(false){
   read(fname_h, magic);
 }
  ~BufferedContainer(){
    clean_up();
  }
};

}//namespace
#endif
