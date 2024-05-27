#include <string>
#include <vector>
#include <iostream>
#include <cstdio>
#include "DummyException.hpp"
#include "BufferedElement.hpp"
#include "ProcessTensorElement.hpp"
#include "BufferedContainer.hpp"
#include "BinaryReader.hpp"


#ifdef DEBUG_BUFFERED_CONTAINER_ALL
#define DEBUG_BUFFERED_CONTAINER
#endif

namespace ACE{

template <typename T> void BufferedContainer<T>::set_preload_none(){
  wait_preload();
  clear_preload();
  preload_block=-1;
  preload_lock=false;
}

template <typename T> void BufferedContainer<T>::check_buffer_bounds(int n)const{
  if(n<0||n>=(int)buffer.size()){
    std::cerr<<"BufferedContainer::check_buffer_bounds: Out of bounds: "<<n<<"/"<<buffer.size()<<std::endl;
    std::cerr<<"info: "; print_info(std::cerr); std::cerr<<std::endl;
    throw DummyException();
  }

}
template <typename T> void BufferedContainer<T>::check_preload_bounds(int n)const{
  if(n<0||n>=(int)preload.size()){
    std::cerr<<"BufferedContainer::check_preload_bounds: Out of bounds: "<<n<<"/"<<preload.size()<<std::endl;
    std::cerr<<"info: "; print_info(std::cerr); std::cerr<<std::endl;
    throw DummyException();
  }
}

template <typename T> std::string BufferedContainer<T>::get_fname(int n){

  if(fname_header==""){
    std::cerr<<"BufferedContainer: get_fname(n) called with fname_header=''!"<<std::endl;
    throw DummyException();
  }
  if(n<0){
    std::cerr<<"BufferedContainer: get_fname(n) called with n="<<n<<"<0!"<<std::endl;
    throw DummyException();
  }
  return fname_header+"_"+std::to_string(n);
}

//accessors
template <typename T>  T &  BufferedContainer<T>::get(int n, PreloadHint hint){
#ifdef DEBUG_BUFFERED_CONTAINER_ALL
  std::cout<<"DEBUG_BC["<<fname_header<<"]: get("<<n<<")"<<std::endl;
#endif
  if(read_only){
    std::cerr<<"BufferedContainer::get("<<n<<") called for a read-only Container!"<<std::endl;
    throw DummyException();
  }

  was_modified=true;
  if(n<0||n>=n_tot){
    std::cerr<<"BufferedContainer::get: Out of bounds: "<<n<<"/"<<n_tot<<std::endl;
    std::cerr<<"info: "; print_info(std::cerr); std::cerr<<std::endl;
    throw DummyException();
  } 
  if(blocksize<1){ 
    check_buffer_bounds(n);
    return buffer[n];
  }
  int bl=n/blocksize;
  int in_bl=n%blocksize;
  
  read_block(bl, hint);
  check_buffer_bounds(in_bl);
  return buffer[in_bl];
}

template <typename T>  const T &  BufferedContainer<T>::peek(int n){

  if(n<0||n>=n_tot){
    std::cerr<<"BufferedContainer::peek: Out of bounds: "<<n<<"/"<<n_tot<<std::endl;
    throw DummyException();
  }
  if(blocksize<0){ 
    check_buffer_bounds(n);
    return buffer[n];
  }
  int bl=n/blocksize;
  int in_bl=n%blocksize;

  if(bl==current_block){
    check_buffer_bounds(in_bl);
    return buffer[in_bl];
  }

  if(bl!=preload_block){
    read_block_preload(bl);
  }
  check_preload_bounds(in_bl);
  return preload[in_bl];
}

template <typename T>  const T &  BufferedContainer<T>::get_ro(int n, PreloadHint hint){

  bool was_modified_bck=was_modified;
  bool read_only_bck=read_only;
  read_only=false;
  T & Tref=get(n, hint);
  read_only=read_only_bck;
  was_modified=was_modified_bck;
  return Tref;
}

template <typename T> void  BufferedContainer<T>::set_read_only(bool ro){
  read_only=ro;
}
template <typename T> void BufferedContainer<T>::clear(){
  if(blocksize<=0){
    clear_buffer();
  }else{
    set_preload_none();
    current_block=-1;
    clear_buffer();
    delete_files();
    was_modified=false;
  }
  n_tot=0;
}

template <typename T> void BufferedContainer<T>::push_back(const T &templ){
#ifdef DEBUG_BUFFERED_CONTAINER
  std::cout<<"DEBUG_BC["<<fname_header<<"]: push_back (applied to n_tot="<<n_tot<<")"<<std::endl;
#endif
  if(blocksize<=0){
    buffer.push_back(templ);
    n_tot++;
  }else{
    int nr_blocks=get_nr_blocks();
    int space_in_last_block=nr_blocks*blocksize-n_tot;
    if(space_in_last_block>0){
      read_block(nr_blocks-1);
      buffer.push_back(templ);
      n_tot++;
      was_modified=true;
    }else{ //start a new block
      set_preload_none();
      if(current_block>=0 && was_modified)write_release_block(current_block);
      clear_buffer();
      buffer.reserve(blocksize);
      buffer.push_back(templ);
      n_tot++;
      current_block=nr_blocks;
      was_modified=true;
    }
    if(fname_header!="")write_header();
  }
}

template <typename T> void BufferedContainer<T>::
            read_block_from_file(int bl, std::vector<T> & buf){
#ifdef DEBUG_BUFFERED_CONTAINER_ALL
  std::cout<<"DEBUG_BC["<<fname_header<<"]: read_block_from_file("<<bl<<")"<<std::endl;
#endif

  int local_blocksize=blocksize;
  if(blocksize>0 && bl==get_nr_blocks()-1){
    local_blocksize=n_tot-bl*blocksize;
  }


  std::vector<T>().swap(buf);
  buf.resize(local_blocksize);
  std::string fname=get_fname(bl);
 
  std::ifstream ifs(fname.c_str());
  if(!ifs.good()){std::cerr<<"BufferedContainer<T>::read_block_from_file: file '"<<fname<<"' incomplete!"<<std::endl;}

  for(int i=0; i<local_blocksize; i++){
    buf[i].read_binary(ifs);
    if(!ifs.good()){std::cerr<<"BufferedContainer<T>::read_block_from_file: file '"<<fname<<"' incomplete!"<<std::endl;}
  }

/*
  //check total size:
  if(blocksize>=0){
    if(buf.size()>blocksize){
      std::cerr<<"BufferedContainer::read_block("<<bl<<"): file '"<<fname<<"': buf.size()>blocksize() ("<<buf.size()<<" vs. "<<blocksize<<")!"<<std::endl;
      throw DummyException();
    }
    if(bl*blocksize+buf.size()>n_tot){
      std::cerr<<"BufferedContainer::read_block("<<bl<<"): file '"<<fname<<"': bl*blocksize+buffer.size()>n_tot ("<<bl*blocksize+buf.size()<<" vs. "<<n_tot<<")!"<<std::endl;
      throw DummyException();
    }
  }
*/
}

template <typename T> void BufferedContainer<T>::read_block(int bl, PreloadHint hint){
#ifdef DEBUG_BUFFERED_CONTAINER_ALL
  std::cout<<"DEBUG_BC["<<fname_header<<"]: read_block("<<bl<<")"<<std::endl;
#endif
  if(current_block==bl)return;  //nothing to do
  if(current_block>=0 && was_modified){
    write_release_block(current_block);
  }

  //check if preloaded
  wait_preload();
  if(preload_block>=0 && preload_block==bl){
    clear_buffer();
    buffer.swap(preload);
    current_block=bl;
    preload_block=-1;
    was_modified=false;

    request_async_preload(bl, hint);
    return;
  }

  clear_buffer();
  read_block_from_file(bl, buffer);
  current_block=bl;
  was_modified=false;
  
  request_async_preload(bl, hint);
}


template <typename T> void BufferedContainer<T>::read_block_preload(int bl){
#ifdef DEBUG_BUFFERED_CONTAINER_ALL
  std::cout<<"DEBUG_BC["<<fname_header<<"]: read_block_preload("<<bl<<")"<<std::endl;
#endif

  if(blocksize<0)return;
  if(preload_block==bl)return;  //nothing to do

  read_block_from_file(bl, preload);
  preload_block=bl;
}

template <typename T> bool BufferedContainer<T>::
  request_async_preload_fct(BufferedContainer<T> *PTB, int bl){

  PTB->read_block_preload(bl);
  return true;
}

template <typename T> void BufferedContainer<T>::request_async_preload(int bl, PreloadHint hint){

  wait_preload();
  switch(hint){
    case ForwardPreload:
      if(bl+1<get_nr_blocks()){ 
        preload_lock=true;
        preload_future=std::async(request_async_preload_fct, this, bl+1);
      }
      break;
    case BackwardPreload:
      if(bl-1>=0){
        preload_lock=true;
        preload_future=std::async(request_async_preload_fct, this, bl-1);
      }
      break;
  }
}

template <typename T> void BufferedContainer<T>::wait_preload(){
  if(preload_lock){
    preload_future.get();
    preload_lock=false;
  }
}

template <typename T> void BufferedContainer<T>::write_release_block(int bl){
  if(bl<0)bl=current_block;
#ifdef DEBUG_BUFFERED_CONTAINER
  std::cout<<"DEBUG_BC["<<fname_header<<"]: write_release_block("<<bl<<")"<<std::endl;
#endif
  if(read_only){
    std::cerr<<"BufferedContainer::write_release_block: won't write when 'read_only' is set!"<<std::endl;
    throw DummyException();
  }
  if(fname_header==""){
    std::cerr<<"BufferedContainer::write_release_block: file name is set to ''!"<<std::endl;
    throw DummyException();
  }

  std::string fname=get_fname(bl);
  std::ofstream os(fname.c_str());
  for(size_t i=0; i<buffer.size(); i++){
    buffer[i].write_binary(os);
  }
  current_block=-1;
}

template <typename T> void BufferedContainer<T>::read(const std::string &filename, const std::string &magic, bool ro){
  clean_up();
  if(magic!=""){
    magicString=magic;
  }
  initialize();
  set_read_only(ro);
  fname_header=filename;
  read_header();
}

template <typename T> void BufferedContainer<T>::read_header(){
#ifdef DEBUG_BUFFERED_CONTAINER
  std::cout<<"DEBUG_BC["<<fname_header<<"]: read_header"<<std::endl;
#endif

  if(fname_header==""){
    std::cerr<<"BufferedContainer::read_header: fname_header==''!"<<std::endl;
    throw DummyException();
  }
  set_preload_none();

  std::ifstream ifs(fname_header.c_str());
  //check if file exists and can be read:
  if(!ifs.good()){
    std::cerr<<"BufferedContainer::read_header: cannot open file '"<<fname_header<<"'!"<<std::endl;
    throw DummyException();
  }
  std::string magic=binary_read_fixedSizeString(ifs, magicString.size());
  if(magic!=magicString){
    std::cerr<<"BufferedContainer::read_header: '"<<fname_header<<"': magic='"<<magic<<"' (should be '"<<magicString<<"')!"<<std::endl;
    throw DummyException();
  }
  ifs>>n_tot;
  ifs>>blocksize;
//  n_tot=binary_read_int(ifs, fname_header+" n_tot");
//  blocksize=binary_read_int(ifs, fname_header+" blocksize");
 
  if(!ifs.good()){
    std::cerr<<"BufferedContainer::read_header: file '"<<fname_header<<"' incomplete!"<<std::endl;
    throw DummyException();
  }
}

template <typename T> void BufferedContainer<T>::write_header(){
#ifdef DEBUG_BUFFERED_CONTAINER
  std::cout<<"DEBUG_BC["<<fname_header<<"]: write_header: '"<<fname_header<<"'"<<std::endl;
#endif

  if(fname_header==""){
    std::cerr<<"BufferedContainer::write_header: fname_header=''!"<<std::endl;
    throw DummyException();
  }
  std::ofstream ofs(fname_header.c_str());
  binary_write_fixedSizeString(ofs, magicString.size(), magicString);

  ofs<<std::endl;
  ofs<<n_tot<<std::endl;
  ofs<<blocksize<<std::endl;
//  binary_write_int(ofs, n_tot);
//  binary_write_int(ofs, blocksize);
}

template <typename T> void BufferedContainer<T>::write_all(){
#ifdef DEBUG_BUFFERED_CONTAINER
  std::cout<<"DEBUG_BC["<<fname_header<<"]: write_all"<<std::endl;
#endif
//std::cout<<"write_all: start"<<std::endl;
//print_info(); std::cout<<std::endl;

  if(fname_header=="" || read_only)return;
  if(blocksize<=0){
    blocksize=n_tot;
    current_block=0;
  }

  if(current_block>=0)write_release_block();
  write_header();
  set_preload_none();
//std::cout<<"write_all: end"<<std::endl;
}


template <typename T> void BufferedContainer<T>::copy_content(BufferedContainer<T> &other){
#ifdef DEBUG_BUFFERED_CONTAINER
  std::cout<<"DEBUG_BC["<<fname_header<<"]: copy_content"<<std::endl;
#endif

  resize(other.n_tot);
#ifdef DEBUG_BUFFERED_CONTAINER
  std::cout<<"DEBUG_BC["<<fname_header<<"]: copy_content: resize successful"<<std::endl;
#endif
  for(int n=0; n<n_tot; n++){
    get(n)=other.get_ro(n);
  }
}

template <typename T> void BufferedContainer<T>::print_info(std::ostream &os)const{

  os<<"fname_header='"<<fname_header<<"'";
  os<<" n_tot="<<n_tot;
  os<<" buffer.size()="<<(int)buffer.size();
  os<<" blocksize="<<blocksize;
  os<<" on_exit=";
  os<<std::string((on_exit==ON_EXIT::WriteOnDestruction)?"WriteOnDestruction":"");
  os<<std::string((on_exit==ON_EXIT::DeleteOnDestruction)?"DeleteOnDestruction":"");
  os<<" read_only="<<std::string((read_only)?"true":"false");
  os<<" was_modified="<<std::string((was_modified)?"true":"false");
  os<<" preload_block="<<preload_block;
  os<<" current_block="<<current_block;
}

template <typename T> void BufferedContainer<T>::resize(int sz){
#ifdef DEBUG_BUFFERED_CONTAINER
  std::cout<<"DEBUG_BC["<<fname_header<<"]: resize("<<sz<<")"<<std::endl;
#endif
  if(sz==n_tot){ //nothing to do
    return;
  }
  if(read_only){
    std::cerr<<"BufferedContainer::resize: can't resize when 'read_only'!"<<std::endl;
    throw DummyException();
  }
  if(blocksize<=0 || fname_header==""){
    n_tot=sz;
    buffer.resize(sz);
    was_modified=true;
    return;
  }
 
  set_preload_none();

  //remove or add blocks?
  int old_n_tot=n_tot;
  int old_nr_blocks=get_nr_blocks();
  n_tot=sz;
  int nr_blocks=get_nr_blocks();

//std::cout<<"TEST1"<<std::endl;
  for(int i=nr_blocks; i<old_nr_blocks; i++){ //remove excess blocks
    std::string fname=get_fname(i);
    std::remove(fname.c_str());
  }

  if(old_nr_blocks<nr_blocks){
//std::cout<<"TEST2"<<std::endl;
    if(old_nr_blocks>0){ //fill up last block of old container
      read_block(old_nr_blocks-1);
      buffer.resize(blocksize);
      was_modified=true;
      write_release_block();
    }
    for(int i=old_nr_blocks; i<nr_blocks-1; i++){ //add new empty blocks
//std::cout<<"TEST2: i="<<i<<std::endl;
//      read_block(i);
      buffer.resize(blocksize);
      current_block=i;
      was_modified=true;
      write_release_block();
    }
  }

  //fix the last block:
  if(nr_blocks<1){ //resize was called with 0
//std::cout<<"TEST3"<<std::endl;
    buffer.clear();
    was_modified=true;
  }else{
    if(old_nr_blocks>=nr_blocks){ //last block existed
//std::cout<<"TEST4"<<std::endl;
      read_block(nr_blocks-1);
    }else{
//std::cout<<"TEST5"<<std::endl;
      current_block=nr_blocks-1;
    }
//std::cout<<"TEST6"<<std::endl;
    buffer.resize(sz-(nr_blocks-1)*blocksize);
    was_modified=true;
  }

//std::cout<<"TEST7"<<std::endl;
  write_header();
//std::cout<<"TEST8"<<std::endl;
}

template <typename T> void BufferedContainer<T>::initialize(const std::string &fname_h, int bs){
  blocksize=bs; 
  on_exit=ON_EXIT::WriteOnDestruction;
  read_only=false;
  was_modified=false;
  fname_header=fname_h;
  n_tot=0;
  clear_buffer();
  current_block=-1;
  set_preload_none();
}

template <typename T> void BufferedContainer<T>::delete_files(){
  if(blocksize<=0 || fname_header=="" || read_only)return;
  int nr_blocks=get_nr_blocks();
  for(int i=0; i<nr_blocks; i++){
    std::string fname=get_fname(i);
    std::remove(fname.c_str());
  }
  std::remove(fname_header.c_str());
}

template <typename T> void BufferedContainer<T>::clean_up(){
bool debug=false;
if(debug){std::cout<<"BufferedContainer::clean_up called."<<std::endl;}
  if(read_only){ 
if(debug){std::cout<<"read_only was set"<<std::endl;}
    return;
  }
  if(on_exit==ON_EXIT::WriteOnDestruction){
if(debug){std::cout<<"WriteOnDestruction was set"<<std::endl;}
    write_all();
  }else if(on_exit==ON_EXIT::DeleteOnDestruction){
if(debug){std::cout<<"DeleteOnDestruction was set"<<std::endl;}
    delete_files();
  }
}


//to prepare object files during compile time:
template class BufferedContainer<BufferedInt>;
template class BufferedContainer<ProcessTensorElement>;

}//namespace
