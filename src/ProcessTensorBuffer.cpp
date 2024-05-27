#include "ProcessTensorBuffer.hpp"
#include "ProcessTensorForward.hpp"
#include "ProcessTensorStream_ro.hpp"
#include "ProcessTensorStream_wo.hpp"
#include "ProcessTensor.hpp"
#include "BinaryReader.hpp"
#include "DummyException.hpp"
#include "AddPT.hpp"
#include "TruncationLayout.hpp"
#include <fstream>
#include <string>

namespace ACE{

void ProcessTensorBuffer::check_buffer_bounds(int n)const{
  if(n<0||n>=(int)buffer.size()){
    std::cerr<<"ProcessTensorBuffer::check_buffer_bounds: Out of bounds: "<<n<<"/"<<buffer.size()<<std::endl;
    std::cerr<<"info: "; print_info(std::cerr); std::cerr<<std::endl;
    throw DummyException();
  }
}
void ProcessTensorBuffer::check_preload_bounds(int n)const{
  if(n<0||n>=(int)preload.size()){
    std::cerr<<"ProcessTensorBuffer::check_preload_bounds: Out of bounds: "<<n<<"/"<<preload.size()<<std::endl;
    std::cerr<<"info: "; print_info(std::cerr); std::cerr<<std::endl;
    throw DummyException();
  }
}
void ProcessTensorBuffer::check_write_buffer_bounds(int n)const{
  if(n<0||n>=(int)write_buffer.size()){
    std::cerr<<"ProcessTensorBuffer::check_write_buffer_bounds: Out of bounds: "<<n<<"/"<<write_buffer.size()<<std::endl;
    std::cerr<<"info: "; print_info(std::cerr); std::cerr<<std::endl;
    throw DummyException();
  }
}
void ProcessTensorBuffer::print_dims(std::ostream &ofs, bool print_both){
  for(int n=0; n<n_tot; n++){
    ProcessTensorElement &e=get(n, ForwardPreload);
    if(n==0)ofs<<e.M.dim_d1;
    ofs<<" "<<e.M.dim_d2;
    if(print_both && n<n_tot-1){
      const ProcessTensorElement &e2=peek(n+1);
      ofs<<"/"<<e2.M.dim_d1;
    }
  }
}
void ProcessTensorBuffer::check_consistency(){
  int dim_d2_last=1;
  for(int n=0; n<n_tot; n++){
    ProcessTensorElement &e=get(n, ForwardPreload);
    if(e.M.dim_d1!=dim_d2_last){
      std::cerr<<"ProcessTensorBuffer::check_consistency: e.M.dim_d1!=dim_d2_last ("<<e.M.dim_d1<<" vs. "<<dim_d2_last<<") at n="<<n<<std::endl;
      std::cerr<<"dimensions: "; print_dims(std::cerr); std::cerr<<std::endl;
      throw DummyException();
    }
    dim_d2_last=e.M.dim_d2;
  }
}

ProcessTensorElement & ProcessTensorBuffer::get(int n, PreloadHint hint){
//std::cout<<"get("<<n<<") (fname_header='"<<fname_header<<"')."<<std::endl;
  if(n<0||n>=n_tot){
    std::cerr<<"ProcessTensorBuffer::get: Out of bounds: "<<n<<"/"<<n_tot<<std::endl;
    std::cerr<<"info: "; print_info(std::cerr); std::cerr<<std::endl;
    throw DummyException();
  }
  if(!use_file){
    check_buffer_bounds(n);
    return buffer[n];
  }
  if(blocksize<0){ //max. one file; Then, buffer is expected to be loaded
    check_buffer_bounds(n);
    return buffer[n];
  }
  int bl=n/blocksize;
  int in_bl=n%blocksize;

  read_block(bl, hint);
  check_buffer_bounds(in_bl);
  return buffer[in_bl];
}

const ProcessTensorElement & ProcessTensorBuffer::peek(int n){
  if(n<0||n>=n_tot){
    std::cerr<<"ProcessTensorBuffer::peek: Out of bounds: "<<n<<"/"<<n_tot<<std::endl;
    throw DummyException();
  }
  if(!use_file){
    check_buffer_bounds(n);
    return buffer[n];
  }
  if(blocksize<0){ //max. one file; Then, buffer is expected to be loaded
    check_buffer_bounds(n);
    return buffer[n];
  }
  int bl=n/blocksize;
  int in_bl=n%blocksize;

  if(bl==current_block){
    check_buffer_bounds(in_bl);
    return buffer[in_bl];
  }

if(false){std::cout<<"called peek("<<n<<"): bl="<<bl<<" not current_block="<<current_block<<" "; print_info(); std::cout<<std::endl;}

  wait_preload();
  if(bl!=preload_block){
if(false){std::cout<<"bl!=preload_block ("<<bl<<" vs. "<<preload_block<<")."<<std::endl;}
    read_block_preload(bl);
  }
  check_preload_bounds(in_bl);
  return preload[in_bl];
}


void ProcessTensorBuffer::clear(){
  if(blocksize<=0){  
    clear_buffer();
  }else{
    set_preload_none();
    wait_write();
    clear_write_buffer();
    current_block=-1;
    clear_buffer();
    delete_files();
  }
  n_tot=0;
}
void ProcessTensorBuffer::push_back(const ProcessTensorElement &templ){
  bool debug=false;
if(debug)std::cout<<"PUSH_BACK START, blocksize="<<blocksize<<std::endl;
  if(blocksize<=0){ 
    buffer.push_back(templ);
    n_tot++;
  }else{
    if(!use_file || use_single_file || fname_header==""){
      std::cerr<<"ProcessTensorBuffer::push: blocksize>0 but !use_file || fname_header==\"\"!"<<std::endl;
      throw DummyException();
    }
    
    int nr_blocks=get_nr_blocks();
if(debug){std::cout<<"INFO: "; print_info(); std::cout<<" nr_blocks="<<nr_blocks<<std::endl;}
    int space_in_last_block=nr_blocks*blocksize-n_tot;
    if(space_in_last_block>0){
if(debug)std::cout<<"space_in_last_block>0"<<std::endl;
      read_block(nr_blocks-1);
      buffer.push_back(templ);
      n_tot++;
    }else{ //start a new block
if(debug)std::cout<<"space_in_last_block<=0"<<std::endl;
      if(current_block>=0)write_release_block(current_block);
      buffer.push_back(templ);
      n_tot++;
      current_block=nr_blocks;
    }
  }
  if(use_file){
    write_header();
  }
if(debug)std::cout<<"PUSH_BACK END"<<std::endl;
}

void ProcessTensorBuffer::create(int n_new, const ProcessTensorElement &templ){
  if(blocksize<=0){  
    buffer.resize(n_new, templ);
    n_tot=n_new;
  }else{
    if(!use_file || use_single_file || fname_header==""){
      std::cerr<<"ProcessTensorBuffer::create: blocksize>0 but !use_file || fname_header==\"\"!"<<std::endl;
      throw DummyException();
    }
    current_block=-1;
    n_tot=n_new;
    int nr_blocks=(n_tot+blocksize-1)/blocksize;
    for(int bl=0; bl<nr_blocks; bl++){
      int bs = (bl==nr_blocks-1) ? n_tot-bl*blocksize : blocksize;
      buffer.resize(bs, templ);
      write_release_block(bl);
    }
  }
  if(use_file){
    write_header();
  }
}
void ProcessTensorBuffer::append(int n, const ProcessTensorElement &templ){
  if(n==0)return;
  if(n<0){
    std::cerr<<"ProcessTensorBuffer::append: n<0!"<<std::endl;
    throw DummyException();
  }
  if(n_tot==0)create(n, templ);
  for(int i=0; i<n; i++){
    push_back(templ);
  }
}

void ProcessTensorBuffer::resize(int n_new, const ProcessTensorElement &templ){
  if(n_new<0){
    std::cerr<<"ProcessTensorBuffer::resize: n_new<0!"<<std::endl;
    throw DummyException();
  } 
  clear();
  create(n_new, templ);
}

void ProcessTensorBuffer::dict_expand(const ReadPT_struct & readPT){
  for(int n=0; n<n_tot; n++){
    get(n, ForwardPreload).accessor.dict_expand(readPT);
  }
}

void ProcessTensorBuffer::copy_content(ProcessTensorBuffer & other){
  clear();
  for(int n=0; n<other.n_tot; n++){
    push_back(other.get(n));
  }
}
void ProcessTensorBuffer::copy_content(const std::string &filename){
  ProcessTensorBuffer other(filename, true);
  copy_content(other);
}
void ProcessTensorBuffer::copy_content(const ReadPT_struct &readPT){
  ProcessTensorBuffer other(readPT, true);
  copy_content(other);
}

void ProcessTensorBuffer::copy_read_only(ProcessTensorBuffer & other){
  ProcessTensorBufferSpec::copy(other);
  read_only=true;
  n_tot=other.n_tot;
  buffer=other.buffer;
  current_block=other.current_block;
  set_preload_none();
  set_write_buffer_none();
}

void ProcessTensorBuffer::copy_InfluenceFunctional_OD(const InfluenceFunctional_OD &IF){
  clear();
  int n_tot=IF.size();
  for(int n=0; n<n_tot; n++){
/*
    ProcessTensorElement e;
    e.clearNF();
    e.M=IF.a[n];
    e.closure=IF.c[n];
    e.env_ops.ops=IF.env_ops[n];
    e.accessor.set_from_dict(IF.dict);
    push_back(e);
*/
    push_back(ProcessTensorElement(IF,n));
  }
}

void ProcessTensorBuffer::delete_files(){
  bool debug=false;
  if(debug){std::cout<<"delete files: info: "; print_info(); std::cout<<std::endl;}

  if(read_only || !use_file)return;
  int nr_blocks=get_nr_blocks();
  for(int i=0; i<nr_blocks; i++){
    std::string fname=get_fname(i);
    if(debug)print_file_exists(fname);
    std::remove(fname.c_str());
  }
  std::remove(fname_header.c_str());
}

void ProcessTensorBuffer::read_block(int bl, PreloadHint hint){
  if(current_block==bl)return;  //nothing to do
  if(current_block>=0){
    write_release_block(current_block);
  }

  //check if preloaded
  wait_preload();
  if(preload_block>=0 && preload_block==bl){
if(false)std::cout<<"found preloaded block: "<<bl<<std::endl;
    clear_buffer();
    buffer.swap(preload);
    current_block=bl;
    preload_block=-1;

    request_async_preload(bl, hint);
    return;
  }

  std::string fname=get_fname(bl);

  //Don't reinvent the wheel: reuse implementation in 'ProcessTensor':
  clear_buffer();
  {ProcessTensor PT_tmp(fname); buffer.swap(PT_tmp.elements);}
  current_block=bl;
  request_async_preload(bl,hint); 

  //check total size:
  if(blocksize<0){
    n_tot=buffer.size();
  }else{
    if(buffer.size()>blocksize){
      std::cerr<<"ProcessTensorBuffer::read_block("<<bl<<"): file '"<<fname<<"': buffer.size()>blocksize() ("<<buffer.size()<<" vs. "<<blocksize<<")!"<<std::endl;
      throw DummyException();
    }
    if(bl*blocksize+buffer.size()>n_tot && !read_only){ //check for read_only: inequality broken if ro copy is read after rw copy is modified!
      std::cerr<<"ProcessTensorBuffer::read_block("<<bl<<"): file '"<<fname<<"': bl*blocksize+buffer.size()>n_tot ("<<bl*blocksize+buffer.size()<<" vs. "<<n_tot<<")!"<<std::endl;
      throw DummyException();
    }
  }
}
void ProcessTensorBuffer::read_block_preload(int bl){
  if(blocksize<0)return;
  if(preload_block==bl)return;  //nothing to do
  std::string fname=get_fname(bl);

if(false)std::cout<<"preload block: "<<bl<<std::endl;
  //Don't reinvent the wheel: reuse implementation in 'ProcessTensor':
  std::vector<ProcessTensorElement>().swap(preload);
  {ProcessTensor PT_tmp(fname); preload.swap(PT_tmp.elements);}
  preload_block=bl;

  //check total size:
  if(blocksize<0){
    return;
  }else{
    if(preload.size()>blocksize){
      std::cerr<<"ProcessTensorBuffer::read_block_preload("<<bl<<"): file '"<<fname<<"': preload.size()>blocksize() ("<<preload.size()<<" vs. "<<blocksize<<")!"<<std::endl;
      throw DummyException();
    }
    if(bl*blocksize+preload.size()>n_tot){
      std::cerr<<"ProcessTensorBuffer::read_block_preload("<<bl<<"): file '"<<fname<<"': bl*blocksize+preload.size()>n_tot ("<<bl*blocksize+preload.size()<<" vs. "<<n_tot<<")!"<<std::endl;
      throw DummyException();
    }
  }
}
  

bool ProcessTensorBuffer::request_async_preload_fct(
                                 ProcessTensorBuffer *PTB, int bl){
  PTB->read_block_preload(bl);
  return true;
}
void ProcessTensorBuffer::request_async_preload(int bl, PreloadHint hint){
  wait_preload();

  switch(hint){
    case ForwardPreload:
      if(bl+1<get_nr_blocks()){ 
        if(write_buffer_block==bl+1){
          wait_write();
        }
        preload_lock=true;
        preload_future=std::async(request_async_preload_fct, this, bl+1);
      }
      break;
    case BackwardPreload:
      if(bl-1>=0){
        if(write_buffer_block==bl-1){
          wait_write();
        }
        preload_lock=true;
        preload_future=std::async(request_async_preload_fct, this, bl-1);
      }
      break;
    default: 
      break;
  }
}
void ProcessTensorBuffer::wait_preload(){
  if(preload_lock){
    preload_future.get();
    preload_lock=false;
  }
} 
void ProcessTensorBuffer::wait_write(){
  if(use_async_write && write_buffer_lock){
    write_buffer_future.get();
    write_buffer_lock=false;
  }
} 

bool ProcessTensorBuffer::write_release_block_async_fct(
                                   std::vector<ProcessTensorElement> &wbuf,
                                   std::string filename, bool debug_){
  ProcessTensor PT_tmp;
  PT_tmp.elements.swap(wbuf);
  PT_tmp.write_binary(filename);
  if(debug_)std::cout<<"writing block to file '"<<filename<<" finished"<<std::endl;
  return true;
}
void ProcessTensorBuffer::write_release_block(int bl){
//std::cout<<"write_release_block("<<bl<<") called."<<std::endl;
  if(read_only||!use_file)return;

  bool debug=false;
//  bool debug=true;
  if(!use_single_file && (bl<0 || bl>=get_nr_blocks())){
    std::cerr<<"ProcessTensorBuffer::write_release_block: fname_header '"<<fname_header<<"': block out of bounds ("<<bl<<"/"<<get_nr_blocks()<<")!"<<std::endl;
    throw DummyException();
  }
  if(use_single_file)bl=0;

  if(debug)std::cout<<"writing block "<<bl<<" to file '"<<get_fname(bl)<<" started"<<std::endl;

  if(use_async_write){
    wait_write();
    clear_write_buffer();
    write_buffer_block=bl;
    write_buffer_lock=true;

    write_buffer.swap(buffer);
    current_block=-1;
    write_header();
    write_buffer_future=std::async(std::launch::async, 
          write_release_block_async_fct, std::ref(write_buffer), get_fname(bl), debug);
  }else{
    ProcessTensor PT_tmp;
    PT_tmp.elements.swap(buffer);
    PT_tmp.write_binary(get_fname(bl));
if(debug)std::cout<<"wrote block "<<bl<<" to file '"<<get_fname(bl)<<std::endl;
    current_block=-1;
    write_header();
  }
}

void ProcessTensorBuffer::read(const std::string &filename, bool ro){
  fname_header=filename;
  read_only=ro;
  use_file=true;
  current_block=-1;
  set_preload_none();
  set_write_buffer_none();
  read_header();
//  was_modified=false;
}
void ProcessTensorBuffer::read(const ReadPT_struct &readPT, bool ro){
  read(readPT.fname, ro);
  dict_expand(readPT);
}

bool ProcessTensorBuffer::can_read(const std::string & fname){

  std::string magic=read_first_bytes(fname,4);
  if(magic=="PT__"||magic=="PTr_" ||magic=="PTBH")return true;
  return false;
}

void ProcessTensorBuffer::read_header(){
  bool debug=false;
  if(!use_file || fname_header==""){
    std::cerr<<"ProcessTensorBuffer::read_header: 'fname_header' not set!"<<std::endl;
    throw DummyException();
  }

  set_preload_none();

  std::ifstream ifs(fname_header.c_str());
  //check if file exists and can be read:
  if(!ifs.good()){
    std::cerr<<"Cannot open PT header file '"<<fname_header<<"'!"<<std::endl;
    throw DummyException();
  }

  //check is single monolithic file without explicit header:
  std::string magic=binary_read_fixedSizeString(ifs, 4);
  if(magic=="PT__"||magic=="PTr_"){
    use_single_file=true;
    blocksize=-1;
    ifs.close();
    read_block(0);
    n_tot=buffer.size();
    return;    
  }

  //remaining case: explicit header file
  if(magic!="PTBH"){
    std::cerr<<"File '"<<fname_header<<"' does not match magic number for a valid process tensor file!"<<std::endl; 
    throw DummyException();
  }

  is_temporary=false;
  use_single_file=false;
  current_block=-1;
  
  n_tot=binary_read_int(ifs, fname_header+" n_tot");
if(debug)std::cout<<"header file '"<<fname_header<<"': n_tot="<<n_tot<<std::endl;

  blocksize=binary_read_int(ifs, fname_header+" blocksize");
if(debug)std::cout<<"header file '"<<fname_header<<"': blocksize="<<blocksize<<std::endl;

  int nr_blocks=get_nr_blocks();
  //check if all files exist:
  for(int i=0; i<nr_blocks; i++){
    std::string fname2=get_fname(i);
    std::ifstream ifs2(fname2.c_str());
    if(!ifs.good()){
      std::cerr<<"Cannot open PT content file '"<<fname2<<"'!"<<std::endl;
      throw DummyException();
    }
  } 
}

void ProcessTensorBuffer::write_header(){
  if(!use_file||use_single_file||read_only)return;
  std::ofstream ofs(fname_header.c_str());
  if(!ofs.good()){
    std::cerr<<"ProcessTensorBuffer::write_header: fname_header='"<<fname_header<<"': !ofs.good()!"<<std::endl;
    throw DummyException();
  }

  binary_write_fixedSizeString(ofs, 4, "PTBH");
  binary_write_int(ofs, n_tot);
  binary_write_int(ofs, blocksize);
}
void ProcessTensorBuffer::write_all(){
  if(!use_file)return;
  if(use_single_file){
//    write_release_block(0); <- don't release !
//    if(false)//<-force write all
    if(!read_only){//<-force write all
      ProcessTensor PT_tmp;
      PT_tmp.elements.swap(buffer);
      PT_tmp.write_binary(get_fname(0));
      PT_tmp.elements.swap(buffer);
    }
    return;
  }
  if(current_block>=0){
    write_release_block(current_block);
    write_header();
    return;
  }
}


void ProcessTensorBuffer::print_info(std::ostream &os)const{
  os<<"use_file="; if(use_file)os<<"true";else os<<"false"; 
  os<<" is_temporary="; if(is_temporary)os<<"true";else os<<"false"; 
  os<<" use_single_file="; if(use_single_file)os<<"true";else os<<"false"; 
  os<<" read_only="; if(read_only)os<<"true";else os<<"false"; 
  os<<" use_async_write="; if(use_async_write)os<<"true";else os<<"false"; 
  os<<" blocksize="<<blocksize;
  os<<" fname_header='"<<fname_header<<"'";
  os<<" n_tot="<<n_tot;
  os<<" preload_block="<<preload_block;
  os<<" current_block="<<current_block;
}

//Operations:
void ProcessTensorBuffer::set_trivial(int n_max, int sysdim){
  if(n_max<1){
    std::cerr<<"ProcessTensorBuffer::set_trivial: n_max<1!"<<std::endl;
    throw DummyException();
  }
  if(sysdim<2){
    std::cerr<<"ProcessTensorBuffer::set_trivial: sysdim<2!"<<std::endl;
    throw DummyException();
  }
  ProcessTensorElement element;
  element.set_trivial(sysdim);
  resize(n_max, element);

  clear_preload();
  preload_block=-1;
}
void ProcessTensorBuffer::calculate_closures(){
  if(n_tot<1)return;

  ProcessTensorElement last=get(n_tot-1);
  last.calculate_closure(NULL);
  get(n_tot-1)=last;
  for(int n=n_tot-2; n>=0; n--){
    get(n).calculate_closure(&last);
    last=get(n);
  }
}

void ProcessTensorBuffer::sweep_forward(const TruncatedSVD &trunc, 
                                        int verbosity,
                                        int range_start, int range_end){

  if(range_start<0)range_start=0;
  if(range_end<0||range_end>n_tot)range_end=n_tot;
  if(range_start==range_end){
     return; 
  }
  if(range_start>range_end){
    std::cerr<<"ProcessTensorBuffer::sweep_forward: range_start>range_end!"<<std::endl; 
    throw DummyException();
  }
  if(verbosity>0){
    std::cout<<"sweep_forward: range: ["<<range_start<<".."<<range_end<<"["<<std::endl;
  }
  
  int maxdim_in=0, maxdim_out=0, maxdim_at=0;
  PassOn pass_on;
  for(int n=range_start; n<range_end; n++){
    ProcessTensorElement & element = get(n, ForwardPreload);
    if(n==range_start){
      pass_on=PassOn(element.M.dim_d1);
      maxdim_in=maxdim_out=element.M.dim_d1;
    }
    if(element.M.dim_d2>maxdim_in)maxdim_in=element.M.dim_d2;

    element.sweep_forward(trunc, pass_on, (n==range_end-1));

    if(element.M.dim_d2>maxdim_out){
      maxdim_out=element.M.dim_d2;
      maxdim_at=n;
    }
  }
  if(verbosity>0)std::cout<<"Maxdim at n="<<maxdim_at<<": "<<maxdim_in<<" -> "<<maxdim_out<<std::endl;
}
void ProcessTensorBuffer::sweep_backward(const TruncatedSVD &trunc,
                                        int verbosity,
                                        int range_start, int range_end){

  if(range_start<0)range_start=0;
  if(range_end<0||range_end>n_tot)range_end=n_tot;
  if(range_start==range_end){
     return; 
  }
  if(range_start>range_end){
    std::cerr<<"ProcessTensorBuffer::sweep_backward: range_start>range_end!"<<std::endl; 
    throw DummyException();
  }
  if(verbosity>0){
    std::cout<<"sweep_backward: range: ["<<range_start<<".."<<range_end<<"["<<std::endl;
  }
  
  PassOn pass_on;
  int maxdim_in=0, maxdim_out=0, maxdim_at=0;
  for(int n=range_end-1; n>=range_start; n--){
    ProcessTensorElement & element = get(n, BackwardPreload);
    if(n==range_end-1){
      pass_on=PassOn(element.M.dim_d2);
      maxdim_in=maxdim_out=element.M.dim_d2;
      maxdim_at=n;
    }
    if(element.M.dim_d1>maxdim_in)maxdim_in=element.M.dim_d1;
    element.sweep_backward(trunc, pass_on, (n==range_start));
    if(element.M.dim_d1>maxdim_out){
      maxdim_out=element.M.dim_d1;
      maxdim_at=n;
    }
  }

  if(verbosity>0)std::cout<<"Maxdim at n="<<maxdim_at<<": "<<maxdim_in<<" -> "<<maxdim_out<<std::endl;
}

void ProcessTensorBuffer::sweep_pair_forward(const TruncatedSVD &trunc, int verbosity){
  if(verbosity>0)std::cout<<"sweep_pair_forward"<<std::endl;
  for(int n=0; n<get_n_tot()-1; n++){
    if(verbosity>1)std::cout<<"n="<<n<<std::endl;
    ProcessTensorElement &e=get(n, ForwardPreload);
    ProcessTensorElement e2=peek(n+1);
    e.sweep_pair_forward(e2, trunc);
    get(n+1, ForwardPreload)=e2;
  }
}
void ProcessTensorBuffer::sweep_pair_backward(const TruncatedSVD &trunc, int verbosity){
  if(verbosity>0)std::cout<<"sweep_pair_backward"<<std::endl;
  for(int n=(int)get_n_tot()-1; n>0; n--){
    if(verbosity>1)std::cout<<"n="<<n<<std::endl;
    ProcessTensorElement &e=get(n, BackwardPreload);
    ProcessTensorElement e2=peek(n-1);
    e.sweep_pair_backward(e2, trunc);
    get(n-1, BackwardPreload)=e2;
  }
}

void ProcessTensorBuffer::join_and_sweep_forward(
                                ProcessTensorBuffer & PTB2,
                                const TruncatedSVD &trunc, int verbosity,
                                ShiftExtend shift_extend, bool alternate){

  PassOn pass_on;

//  bool debug=false;

  int n_tot_old=n_tot;

  if(n_tot_old<=shift_extend.shift_second){ //no overlap in this case
    std::cerr<<"ProcessTensorBuffer::join_and_sweep_forward: ";
    std::cerr<<"PT_ro1.size()<=shift_extend.shift_second!"<<std::endl;
    throw DummyException();
  }

  //determine total length of result
  int n_tot_new=PTB2.n_tot+shift_extend.shift_second;
  if(shift_extend.truncate_at>=0){
    if(shift_extend.truncate_at>n_tot_new){ //no truncation required
       shift_extend.truncate_at=-1; 
    }else{
      n_tot_new=shift_extend.truncate_at;
    }
  } 
  //does PT1 have to be extendend?
  int extend_first=n_tot_new-n_tot_old;
  //does PT1 have to be truncated?
  int n_trunc1=0;
  if(extend_first<0){
    n_trunc1=-extend_first;
    extend_first=0;
  }

  if(verbosity>0){
    std::cout<<"join_and_sweep_forward: range: ["<<shift_extend.shift_second<<".."<<n_tot_new<<"["<<std::endl;
    std::cout<<""<<fname_header<<"','"<<PTB2.fname_header<<std::endl;
  }
  int maxdim_in1=0, maxdim_in2=0, maxdim_out=0, maxdim_at=0;

  //first n=shift_extend elements are already set. Combine loop:
  for(int n=shift_extend.shift_second; n<n_tot_old-n_trunc1; n++){
    ProcessTensorElement & e = get(n, ForwardPreload);
    ProcessTensorElement & e2 = PTB2.get(n-shift_extend.shift_second, ForwardPreload);

    if(e.M.dim_d1>maxdim_in1)maxdim_in1=e.M.dim_d2;
    if(e2.M.dim_d1>maxdim_in2)maxdim_in2=e2.M.dim_d2;

    if(n==n_tot_old-n_trunc1-1){
      ProcessTensorElement e3=e2;
      e3.close_off(); 
//      e.join_thisfirst(e3);
      e.join(e3, (n%2==1) && (alternate) );
    }else{
//      e.join_thisfirst(e2);
      e.join(e2, (n%2==1) && (alternate));
    }

    if(n==shift_extend.shift_second){
      pass_on=PassOn(e.M.dim_d1);
    }
    e.sweep_forward(trunc, pass_on, (n==n_tot_new-1));
    if(e.M.dim_d2>maxdim_out){
      maxdim_out=e.M.dim_d2;
      maxdim_at=n;
    }
  }

  //write overhanging elements of shiftet PT_ro2 (length: extend_first) )
  if(extend_first>0){
    for(int n=n_tot_old; n<n_tot_new; n++){
      push_back(PTB2.get(n-shift_extend.shift_second, ForwardPreload));
      ProcessTensorElement & e = get(n);
      e.sweep_forward(trunc, pass_on, (n==n_tot_new-1));
    }
  }
  get(n_tot_new-1).close_off();
  
  if(verbosity>0)std::cout<<"Maxdim at n="<<maxdim_at<<": "<<maxdim_in1<<","<<maxdim_in2<<" -> "<<maxdim_out<<std::endl;
}

void ProcessTensorBuffer::join_symmetric_and_sweep_forward(
                                ProcessTensorBuffer & PTB2,
                                const TruncatedSVD &trunc, int verbosity,
                                ShiftExtend shift_extend){

  PassOn pass_on;

//  bool debug=false;

  int n_tot_old=n_tot;

  if(n_tot_old<=shift_extend.shift_second){ //no overlap in this case
    std::cerr<<"ProcessTensorBuffer::join_symmetric_and_sweep_forward: ";
    std::cerr<<"PT_ro1.size()<=shift_extend.shift_second!"<<std::endl;
    throw DummyException();
  }
  if(PTB2.n_tot%2!=0){
    std::cerr<<"ProcessTensorBuffer::join_symmetric_and_sweep_forward: ";
    std::cerr<<"PTB2 must have a number of elements divisible by 2!"<<std::endl;
    throw DummyException();
  }

  //determine total length of result
  int n_tot_new=PTB2.n_tot/2+shift_extend.shift_second;
  if(shift_extend.truncate_at>=0){
    if(shift_extend.truncate_at>n_tot_new){ //no truncation required
       shift_extend.truncate_at=-1; 
    }else{
      n_tot_new=shift_extend.truncate_at;
    }
  } 
  //does PT1 have to be extendend?
  int extend_first=n_tot_new-n_tot_old;
  //does PT1 have to be truncated?
  int n_trunc1=0;
  if(extend_first<0){
    n_trunc1=-extend_first;
    extend_first=0;
  }

  if(verbosity>0){
    std::cout<<"join_symmetric_and_sweep_forward: range: ["<<shift_extend.shift_second<<".."<<n_tot_new<<"["<<std::endl;
    std::cout<<""<<fname_header<<"','"<<PTB2.fname_header<<std::endl;
  }
  int maxdim_in1=0, maxdim_in2=0, maxdim_out=0, maxdim_at=0;

  //first n=shift_extend elements are already set. Combine loop:
  for(int n=shift_extend.shift_second; n<n_tot_old-n_trunc1; n++){
    int n2=n-shift_extend.shift_second;

    ProcessTensorElement & e = get(n, ForwardPreload);
    ProcessTensorElement & e2 = PTB2.get(2*n2, ForwardPreload);
    ProcessTensorElement & e3 = PTB2.get(2*n2+1, ForwardPreload);

    if(e.M.dim_d2>maxdim_in1)maxdim_in1=e.M.dim_d2;
    if(e3.M.dim_d2>maxdim_in2)maxdim_in2=e3.M.dim_d2;

    if(n==n_tot_old-n_trunc1-1){
      ProcessTensorElement e4=e3;
      e4.close_off(); 
      e.join_symmetric(e2,e4);
    }else{
      e.join_symmetric(e2,e3);
    }

    if(n==shift_extend.shift_second){
      pass_on=PassOn(e.M.dim_d1);
    }
    e.sweep_forward(trunc, pass_on, (n==n_tot_new-1));
    if(e.M.dim_d2>maxdim_out){
      maxdim_out=e.M.dim_d2;
      maxdim_at=n;
    }
  }

  //write overhanging elements of shiftet PT_ro2 (length: extend_first) )
  if(extend_first>0){
    std::cerr<<"ProcessTensorBuffer::join_symmetric_and_sweep_forward: ";
    std::cerr<<"extend_first>0 NOT SUPPORTED YET!"<<std::endl;
    throw DummyException();
/*
    for(int n=n_tot_old; n<n_tot_new; n++){
      int n2=n-shift_extend.shift_second;
      push_back(PTB2.get(n-shift_extend.shift_second, ForwardPreload));
      ProcessTensorElement & e = get(n);
      e.sweep_forward(trunc, pass_on, (n==n_tot_new-1));
    }
*/
  }

  get(n_tot_new-1).close_off();
  
  if(verbosity>0)std::cout<<"Maxdim at n="<<maxdim_at<<": "<<maxdim_in1<<","<<maxdim_in2<<" -> "<<maxdim_out<<std::endl;
}


void ProcessTensorBuffer::join_and_sweep_backward(
                                ProcessTensorBuffer & PTB2,
                                const TruncatedSVD &trunc, int verbosity,
                                ShiftExtend shift_extend){

  bool sweep_overhang = true;
  PassOn pass_on;

//  bool debug=false;

  int n_tot_old=n_tot;

  if(n_tot_old<=shift_extend.shift_second){ //no overlap in this case
    std::cerr<<"ProcessTensorBuffer::join_and_sweep_backward: ";
    std::cerr<<"PT_ro1.size()<=shift_extend.shift_second!"<<std::endl;
    throw DummyException();
  }

  //determine total length of result
  int n_tot_new=PTB2.n_tot+shift_extend.shift_second;
  if(shift_extend.truncate_at>=0){
    if(shift_extend.truncate_at>n_tot_new){
       shift_extend.truncate_at=-1;
    }else{
      n_tot_new=shift_extend.truncate_at;
    }
  } 
  //does PT1 have to be extendend?
  int extend_first=n_tot_new-n_tot_old;
  int n_trunc1=0;
  if(extend_first<0){
    n_trunc1=-extend_first;
    extend_first=0;
  }

  if(verbosity>0){
    std::cout<<"join_and_sweep_backward: range: ["<<shift_extend.shift_second<<".."<<n_tot_new<<"["<<std::endl;
    std::cout<<fname_header<<"','"<<PTB2.fname_header<<std::endl;
  }
  int maxdim_in1=0, maxdim_in2=0, maxdim_out=0;

  //write overhanging elements of shiftet PT_ro2 (length: extend_first) )
  if(extend_first>0){
    append(extend_first);
    for(int n=n_tot_new-1; n>=n_tot_old; n--){
      ProcessTensorElement & e=get(n, BackwardPreload);
      e=PTB2.get(n-shift_extend.shift_second, BackwardPreload);
      if(n==n_tot_new-1){
        e.close_off();
        if(sweep_overhang){
          pass_on.set(e.M.dim_d2);
        }
      }
      if(sweep_overhang){
        e.sweep_backward(trunc, pass_on, (n==0));
      }
    }
  }
  //overlapping elements:
  for(int n=n_tot_old-1-n_trunc1; n>=shift_extend.shift_second; n--){
    ProcessTensorElement & e = get(n, BackwardPreload);
    ProcessTensorElement & e2 = PTB2.get(n-shift_extend.shift_second, BackwardPreload);
    if(e.M.dim_d1>maxdim_in1)maxdim_in1=e.M.dim_d1;
    if(e2.M.dim_d1>maxdim_in2)maxdim_in2=e2.M.dim_d1;

    if(n==n_tot_old-1-n_trunc1 && extend_first<=0){
      ProcessTensorElement e2_tmp=e2; 
      e2_tmp.close_off();
      e.join_thisfirst(e2_tmp);
    }else{
      e.join_thisfirst(e2);
    }

    if(n==n_tot_old-1-n_trunc1 && !(extend_first>0 && sweep_overhang)){
      pass_on=PassOn(e.M.dim_d2);
    }
    e.sweep_backward(trunc, pass_on, (n==shift_extend.shift_second));
    if(e.M.dim_d1>maxdim_out)maxdim_out=e.M.dim_d1;
  }
  if(verbosity>0)std::cout<<"Maxdim: "<<maxdim_in1<<","<<maxdim_in2<<" -> "<<maxdim_out<<std::endl;
}


void ProcessTensorBuffer::join_select_and_sweep_backward(
                                ProcessTensorBuffer & PTB2,
                                const TruncatedSVD &trunc_select, 
                                const TruncatedSVD &trunc_bw, 
                                int verbosity,
                                ShiftExtend shift_extend){

  int sweep_more_low=shift_extend.sweep_more;
  bool subsequent_sweep = true;
  PassOn pass_on;

  bool debug=false;//true;
if(debug)std::cout<<"jnsbw: started"<<std::endl;

  int n_tot_old=n_tot;

  if(n_tot_old<=shift_extend.shift_second){ //no overlap in this case
    std::cerr<<"ProcessTensorBuffer::join_select_and_sweep_backward: ";
    std::cerr<<"PT_ro1.size()<=shift_extend.shift_second!"<<std::endl;
    throw DummyException();
  }

  //determine total length of result
  int n_tot_new=PTB2.n_tot+shift_extend.shift_second;
  if(shift_extend.truncate_at>=0){
    if(shift_extend.truncate_at>n_tot_new){
       shift_extend.truncate_at=-1;
    }else{
      n_tot_new=shift_extend.truncate_at;
    }
  } 
  //does PT1 have to be extendend?
  int extend_first=n_tot_new-n_tot_old;
  int n_trunc1=0;
  if(extend_first<0){
    n_trunc1=-extend_first;
    extend_first=0;
  }
std::cout<<"extend_first="<<extend_first<<std::endl;

  //if sweep_more!=0: lowest index of element touched:
  int sweep_more_limit=shift_extend.shift_second;
  if(subsequent_sweep && (sweep_more_low!=0)){
    sweep_more_limit=0;
    if(sweep_more_low>0 && shift_extend.shift_second-sweep_more_low>0){
      sweep_more_limit=shift_extend.shift_second-sweep_more_low;
    }
  }

  if(verbosity>0){
    std::cout<<"join_select_and_sweep_backward: range: ["<<sweep_more_limit<<".."<<n_tot_new<<"["<<std::endl;
    std::cout<<fname_header<<"','"<<PTB2.fname_header<<std::endl;
  }
  int maxdim_in1=0, maxdim_in2=0, maxdim_out=0, maxdim_at=-1;

  if(debug){std::cout<<"n_tot_old="<<n_tot_old<<" n_tot_new="<<n_tot_new<<" extend_first="<<extend_first<<" shift_second="<<shift_extend.shift_second<<std::endl;}

  if(debug)std::cout<<"jnsbw: extend"<<std::endl;
  //write overhanging elements of shiftet PT_ro2 (length: extend_first) )
  if(extend_first>0){
    append(extend_first);
    if(debug){std::cout<<"after append: ";print_info();std::cout<<std::endl;}
    for(int n=n_tot_new-1; n>=n_tot_old; n--){
try{
      if(debug)std::cout<<"jnsbw: extend: n="<<n<<std::endl;
      ProcessTensorElement & element=get(n, BackwardPreload);
      if(n==n_tot_new-1){
        element = PTB2.get(n-shift_extend.shift_second, BackwardPreload);
        element.close_off();
//        pass_on=PassOn(element.M.dim_d2);
      }else{
        element = PTB2.get(n-shift_extend.shift_second, BackwardPreload);
      }

      if(subsequent_sweep){ //prepare pass_on and compress
        if(n==n_tot_new-1){
          pass_on=PassOn(element.M.dim_d2);
        }
        if(shift_extend.sweep_more==0){ //don't compress, update pass_on
          pass_on=PassOn(element.M.dim_d1);
        }else if(shift_extend.sweep_more<0){ //compress all
          element.sweep_backward(trunc_bw, pass_on, (n==0));
        }else if(shift_extend.sweep_more>0){//compress some
          if(n==n_tot_old-1+shift_extend.sweep_more){
            pass_on=PassOn(element.M.dim_d2);
          }
          if(n<=n_tot_old-1+shift_extend.sweep_more){ 
            element.sweep_backward(trunc_bw, pass_on, (n==0));
          }
        }
      }
}catch(std::exception &e){
  std::cerr<<"called from join_select_and_sweep_backward (ext); n="<<n<<std::endl;
  throw e;
}
    }
  }

  if(debug)std::cout<<"jnsbw: prepare overlapping"<<std::endl;
  ProcessTensorElement element1=get(n_tot_old-1-n_trunc1, BackwardPreload);
  ProcessTensorElement element2=PTB2.get(n_tot_old-1-n_trunc1-shift_extend.shift_second, BackwardPreload);
  if(extend_first==0){  //previous loop was not triggered: close PTs
    element1.close_off();
    element2.close_off();
  }

  ProcessTensorElement next_element1, next_element2;  
  SelectIndices k_list_left;
  SelectIndices k_list_right = 
             element1.get_forwardNF_selected_indices(element2, trunc_select);

  if(subsequent_sweep && (!extend_first)){
    pass_on=PassOn(k_list_right.size());
  }

  int maxdim_pre_select=k_list_right.size();
  maxdim_in1=element1.M.dim_d2;
  maxdim_in2=element2.M.dim_d2;
  maxdim_out=k_list_right.size();


  if(debug)std::cout<<"jnsbw: start overlapping"<<std::endl;
  int last_overlapping=n_tot_old-1-n_trunc1;
  for(int n=last_overlapping; n>=shift_extend.shift_second; n--){
try{
    if(verbosity>=2){
      std::cout<<"join_select_and_sweep_backward: n="<<n<<"/"<<n_tot_new<<std::endl;
    }
    if(n==shift_extend.shift_second){
      k_list_left.set_full(element1.M.dim_d1, element2.M.dim_d1);
    }else{
      next_element1=peek(n-1);
      next_element2=PTB2.peek(n-1-shift_extend.shift_second);
      k_list_left=next_element1.get_forwardNF_selected_indices(next_element2, trunc_select);

    }

// before actual selection:
    if(maxdim_pre_select<k_list_left.size())maxdim_pre_select=k_list_left.size();
    if(element1.M.dim_d1>maxdim_in1)maxdim_in1=element1.M.dim_d1;
    if(element2.M.dim_d1>maxdim_in2)maxdim_in2=element2.M.dim_d1;

//std::cout<<"before join_selected n="<<n<<std::endl;
#ifdef JNSBW_AVERAGE
    element1.join_average_selected(element2, k_list_left, k_list_right);
#elif defined(JNSBW_NOSYM)
    element1.join_selected(0, element2, k_list_left, k_list_right);
#else
    element1.join_selected(n, element2, k_list_left, k_list_right);
#endif
//std::cout<<"after join_selected n="<<n<<std::endl;

    if(n==last_overlapping){ //modify PassOn at interface.
/*
std::cout<<"element1.M.dim_d1="<<element1.M.dim_d1<<std::endl;
std::cout<<"element1.M.dim_d2="<<element1.M.dim_d2<<std::endl;
std::cout<<"element2.M.dim_d1="<<element2.M.dim_d1<<std::endl;
std::cout<<"element2.M.dim_d2="<<element2.M.dim_d2<<std::endl;
std::cout<<"k_list_right.size()="<<k_list_right.size()<<std::endl;
std::cout<<"pass_on.P.rows()="<<pass_on.P.rows()<<" pass_on.P.cols()="<<pass_on.P.cols()<<std::endl;
*/
      PassOn p2=pass_on;
      pass_on.P=Eigen::MatrixXcd::Zero(k_list_right.size(), p2.P.cols());
      pass_on.Pinv=Eigen::MatrixXcd::Zero(p2.Pinv.rows(),k_list_right.size());
      for(int i=0; i<k_list_right.size(); i++){
        for(int c=0; c<p2.P.cols(); c++){
          pass_on.P(i, c)=p2.P(k_list_right[i].second, c);
          pass_on.Pinv(c, i)=p2.Pinv(c, k_list_right[i].second);
        }
      }
//std::cout<<"pass_on.P.rows()="<<pass_on.P.rows()<<" pass_on.P.cols()="<<pass_on.P.cols()<<std::endl;
    }


//backward sweep of selected block
    if(shift_extend.sweep_more!=0){
      int fin=0;
      if(shift_extend.sweep_more>0 && shift_extend.shift_second>shift_extend.sweep_more){
        fin=shift_extend.shift_second-shift_extend.sweep_more;
      }
      element1.sweep_backward(trunc_bw, pass_on, (n==fin));
    }else if(subsequent_sweep){
      element1.sweep_backward(trunc_bw, pass_on, (n==shift_extend.shift_second));
    }
//std::cout<<"after sweep backward n="<<n<<std::endl;

    if(element1.M.dim_d1>maxdim_out){
      maxdim_out=element1.M.dim_d1;
      maxdim_at=n;
    }

//prepare for next element
    get(n)=element1;
    k_list_right=k_list_left;
    element1=next_element1;
    element2=next_element2;
   
}catch(DummyException &e){
  std::cerr<<"called from join_select_and_sweep_backward (overlap); n="<<n<<std::endl;
  throw e;
} 
  }
  
//handle overhanging lower end:
  if(debug)std::cout<<"jnsbw: compress low end"<<std::endl;
  for(int n=shift_extend.shift_second-1; n>=sweep_more_limit; n--){
    get(n).sweep_backward(trunc_bw, pass_on, (n==sweep_more_limit));
  }
 
  
  if(verbosity>0){
    std::cout<<"Maxdim at preselection: "<<maxdim_pre_select<<std::endl;
    std::cout<<"Maxdim at n="<<maxdim_at<<": "<<maxdim_in1<<","<<maxdim_in2<<" -> "<<maxdim_out<<std::endl;
  }
  if(debug)std::cout<<"jnsbw: done"<<std::endl;
}


void ProcessTensorBuffer::expand_DiagBB(DiagBB &diagBB, double dict_zero){
  if(diagBB.sys_dim()==diagBB.get_dim() && !diagBB.hs_rot.used()){ 
     return; //nothing to do!
  }
  
  if(diagBB.sys_dim()!=diagBB.get_dim()){
    std::cout<<"Expanding DiagBB from diagBB.get_dim()="<<diagBB.get_dim()<<" to diagBB.sys_dim()="<<diagBB.sys_dim()<<std::endl;
  }

  for(int n=0; n<n_tot; n++){
    if(diagBB.sys_dim()!=diagBB.get_dim()){

      get(n, ForwardPreload).expand_DiagBB(diagBB);
    }
    if(diagBB.hs_rot.used()){
      get(n, ForwardPreload).apply_HilbertSpaceRotation(diagBB.hs_rot, dict_zero);
    }
  }
}

void ProcessTensorBuffer::apply_HilbertSpaceRotation(const HilbertSpaceRotation &hs_rot, double dict_zero){
  if(!hs_rot.used())return;
  for(int n=0; n<n_tot; n++){
    get(n, ForwardPreload).apply_HilbertSpaceRotation(hs_rot, dict_zero);
  }
}


bool ProcessTensorBuffer::split_inner(ProcessTensorBuffer &second, int center){

//returns false if center dimension is less than 2

/* Split MPO into two: 
- copy current MPO to "second"
- take "center" element
- split inner indices in two (first try: alternating), first part (even) stays with 'this', second part goes to 'second'
- both MPOs: sweep forward from center, sweep backward, sweep forward again

=> We hope to end up with two smaller MPOs that can be separately multiplied
to another MPO. The results can then be added (with adds the dimensions, not
multiplies them).
*/
   
  if(center<0||center>=get_n_tot()-1){
    std::cerr<<"ProcessTensorBuffer::split_inner: center<0||center>=get_n_tot()-1 ("<<center<<" vs. "<<get_n_tot()-1<<std::endl;
    throw DummyException();
  }

  ProcessTensorElement E11=get(center);
  if(E11.M.dim_d2<2){  //Need to generate a zero element w.r.t addition
    second.clear();
    second.resize(get_n_tot(), ProcessTensorElement(E11.get_N()));
    second.get(0).M(0,0,0)=0;
    return false;
  }
  ProcessTensorElement &E12=get(center+1);

  second.copy_content(*this);
  ProcessTensorElement E21=second.get(center);
  ProcessTensorElement &E22=second.get(center+1);

  bool ret=AddPT::split_alternating(E11, E12, E21, E22);
  get(center)=E11;
  second.get(center)=E21;
  return ret;
}

bool ProcessTensorBuffer::split_inner_and_sweep_fbf(ProcessTensorBuffer &second, const TruncatedSVD &trunc, int center, int verbosity){


  bool can_split=split_inner(second, center);
  if(!can_split){
    if(verbosity>0){
      std::cout<<"Can't split PTB at center="<<center<<std::endl;
    }
    return false;
  }

  if(verbosity>0){
    std::cout<<"Split at center: "<<center<<std::endl;
  }
  sweep_forward(trunc, verbosity, center+1);
  second.sweep_forward(trunc, verbosity, center+1);

  sweep_backward(trunc, verbosity);
  second.sweep_forward(trunc, verbosity);

  sweep_forward(trunc, verbosity);
  second.sweep_forward(trunc, verbosity);
  return true;
}

void ProcessTensorBuffer::additive_join(ProcessTensorBuffer & PTB2){
  if(get_n_tot()!=PTB2.get_n_tot()){
    std::cerr<<"ProcessTensorBuffer::additive_join: get_n_tot()!=PTB2.get_n_tot()!"<<std::endl;
    throw DummyException();
  }
  if(get_n_tot()<2){
    std::cerr<<"ProcessTensorBuffer::additive_join: get_n_tot()<2!"<<std::endl;
    throw DummyException();
  }

  AddPT::add_head(get(0,ForwardPreload), PTB2.get(0,ForwardPreload));
  for(int n=1; n<get_n_tot()-1; n++){
    AddPT::add(get(n,ForwardPreload), PTB2.get(n,ForwardPreload));  
  }
  AddPT::add_tail( get(get_n_tot()-1,ForwardPreload), 
                   PTB2.get(get_n_tot()-1,ForwardPreload) );

}

//void ProcessTensorBuffer::additive_join_and_sweep_backward(
//        ProcessTensorBuffer & PTB2, const TruncatedSVD &trunc, int verbosity){
//}

void ProcessTensorBuffer::sweep_intermediate_or_final_start_forward(const TruncationLayout &trunc, int cur_line, int max_lines, int verbosity, int range_start, int range_end){

    bool is_final= ( cur_line==max_lines-1 );

    int iter=is_final ? trunc.get_final_sweep_n()
                      : trunc.get_intermediate_sweep_n();
    //use same TruncatedSVD as in sweep before, except for last final sweep
    TruncatedSVD trunc_tmp=trunc.get_backward(cur_line, max_lines);
    TruncatedSVD trunc_final = (!is_final || trunc.final_sweep_half) 
                                ? trunc_tmp : trunc.get_final_sweep();

//std::cout<<"TEST: PTB::siorfsf: iter="<<iter<<" trunc.get_intermediate_sweep_n()="<<trunc.get_intermediate_sweep_n()<<" trunc.get_final_sweep_n()="<<trunc.get_final_sweep_n()<<std::endl;
//if(is_final){std::cout<<"is_final=true"<<std::endl;}else{std::cout<<"is_final=false"<<std::endl;}

    if(is_final){range_start=0; range_end=-1;}
    for(int loop=0; loop<iter; loop++){
      if(verbosity>0){
        if(is_final){
          std::cout<<"final loop="<<loop<<"/"<<iter<<std::endl;
        }else{
          std::cout<<"intermediate loop="<<loop<<"/"<<iter<<std::endl;
        }
      }
      if(verbosity>0){trunc_tmp.print_info();std::cout<<std::endl;}
      sweep_forward(trunc_tmp, verbosity, range_start, range_end);

      TruncatedSVD trunc_tmp2=trunc_tmp;
      if(loop==iter-1){ trunc_tmp2=trunc_final; }
      if(verbosity>0){trunc_tmp2.print_info();std::cout<<std::endl;}
      sweep_backward(trunc_tmp2, verbosity, range_start, range_end);
    }

    if(is_final && trunc.final_sweep_half){
      if(verbosity>0){std::cout<<"final forward sweep"<<std::endl;}
      trunc_final=trunc.get_final_sweep();
      if(verbosity>0){trunc_final.print_info();std::cout<<std::endl;}
      sweep_forward(trunc_final, verbosity, range_start, range_end);
    }
}

void ProcessTensorBuffer::sweep_intermediate_or_final_start_backward(const TruncationLayout &trunc, int cur_line, int max_lines, int verbosity, int range_start, int range_end){

    bool is_final= ( cur_line==max_lines-1 );

    int iter=is_final ? trunc.get_intermediate_sweep_n()
                      : trunc.get_final_sweep_n();
    //use same TruncatedSVD as in sweep before, except for last final sweep
    TruncatedSVD trunc_tmp=trunc.get_forward(cur_line, max_lines);
    TruncatedSVD trunc_final = (is_final && !trunc.final_sweep_half) 
                                ? trunc_tmp : trunc.get_final_sweep();


    if(is_final){range_start=0; range_end=-1;}
    for(int loop=0; loop<iter; loop++){
      if(verbosity>0){
        if(loop==iter-1){
          std::cout<<"final loop="<<loop<<"/"<<iter<<std::endl;
        }else{
          std::cout<<"intermediate loop="<<loop<<"/"<<iter<<std::endl;
        }
      }
      if(verbosity>0){trunc_tmp.print_info();std::cout<<std::endl;}
      sweep_backward(trunc_tmp, verbosity, range_start, range_end);

      TruncatedSVD trunc_tmp2=trunc_tmp;
      if(loop==iter-1){ trunc_tmp2=trunc_final; }
      if(verbosity>0){trunc_tmp2.print_info();std::cout<<std::endl;}
      sweep_forward(trunc_tmp2, verbosity, range_start, range_end);
    }

    if(is_final && trunc.final_sweep_half){
      if(verbosity>0){std::cout<<"final backward sweep"<<std::endl;}
      trunc_final=trunc.get_final_sweep();
      if(verbosity>0){trunc_final.print_info();std::cout<<std::endl;}
      sweep_backward(trunc_final, verbosity, range_start, range_end);
    }
}

void ProcessTensorBuffer::set_from_DiagBB_single_line(DiagBB &diagBB, double dt, int n){

  if(!diagBB.is_set_up()){
    std::cerr<<"ProcessTensorBuffer::set_from_DiagBB_single_line: DiagBB not set up!"<<std::endl;
    throw DummyException();
  }
  int N=diagBB.get_dim();  //sys_dim();
  int NL=N*N;
  if(n<1){
    std::cerr<<"ProcessTensorBuffer::set_from_DiagBB_single_line: n<1!"<<std::endl;
    throw DummyException();
  }
  resize(n);

  ProcessTensorElement e;
  if(n==1){  //special case: matrix dimensions NL,1,1
    e.clear();
    e.accessor.dict.set_default_diag(N);
    e.closure=Eigen::VectorXcd::Ones(1);
    e.env_ops.ops=std::vector<Eigen::VectorXcd>(1, e.closure);
    e.M.resize(NL,1,1);
    Eigen::MatrixXcd expS=diagBB.calculate_expS(0,dt);
    for(int i=0; i<NL; i++){
      e.M(i,0,0)=expS(i,i);
    }
    get(0)=e;

  }else{
    e.clear();
    e.accessor.dict.set_default_diag(N);
    e.closure=Eigen::VectorXcd::Ones(NL);
    e.env_ops.ops=std::vector<Eigen::VectorXcd>(1, e.closure);
    e.M.resize(NL,1,NL);
    e.M.set_zero();
    Eigen::MatrixXcd expS=diagBB.calculate_expS(0,dt);
    for(int i=0; i<NL; i++){
        e.M(i, 0, i)=expS(i,i);
    }
    get(0)=e;

    for(int k=1; k<n-1; k++){
      e.clear();
      e.accessor.dict.set_default_diag(N);
      e.closure=Eigen::VectorXcd::Ones(NL);
      e.env_ops.ops=std::vector<Eigen::VectorXcd>(1, e.closure);
      e.M.resize(NL,NL,NL);
      e.M.set_zero();
      expS=diagBB.calculate_expS(k,dt);
      for(int i=0; i<NL; i++){
        for(int j=0; j<NL; j++){
          e.M(i, j, j)=expS(i,j);
        }
      }
      get(k)=e;
    }

    e.clear();
    e.accessor.dict.set_default_diag(N);
    e.closure=Eigen::VectorXcd::Ones(1);
    e.env_ops.ops=std::vector<Eigen::VectorXcd>(1, e.closure);
    e.M.resize(NL,NL,1);
    e.M.set_zero();
    expS=diagBB.calculate_expS(n-1,dt);
    for(int i=0; i<NL; i++){
      for(int j=0; j<NL; j++){
        e.M(i, j, 0)=expS(i,j);
      }
    }
    get(n-1)=e;
  }
  std::cout<<"info after single line: "; print_info(); std::cout<<std::endl;
}

void ProcessTensorBuffer::set_from_DiagBB_single_line_auto(DiagBB &diagBB, double dt, int n_tot, int n_mem, const TruncatedSVD & trunc, bool print_info, bool reverse){

//std::cout<<"TEST: dt="<<dt<<" n_tot="<<n_tot<<" n_mem="<<n_mem<<std::endl;
  if(n_mem==-2){ //encodes automatic detection
    set_from_DiagBB_single_line(diagBB, dt, n_tot);
    if(get_n_tot()<1){
      std::cerr<<"ProcessTensorBuffer::set_from_DiagBB_single_line_auto: PTB_line.get_n_tot()<1!"<<std::endl;
      throw DummyException();
    }

    if(reverse){
      sweep_backward(trunc, 0);
      sweep_forward(trunc, 0);
    }else{
      sweep_forward(trunc, 0);
      sweep_backward(trunc, 0);
    }

    int n_cut=-1;
    if(get(n_tot-1).M.dim_d1!=1){
      n_cut=n_tot;
    }else{
      for(int n=n_tot-1; n>=0; n--){
        if(get(n).M.dim_d1!=1){
          n_cut=n+1;
          break;
        } 
      }
      calculate_closures();
      if(n_cut>0)
      get(n_cut-1).close_off();
      {
        ProcessTensorBuffer tmp_buf; tmp_buf.copy_content(*this);
        clear();
        for(int n=0; n<n_cut; n++){
          push_back(tmp_buf.get(n));
        }
      }
    }
    if(print_info){
      std::cout<<"Automatic memory truncation: "<<n_tot<<" -> "<<get_n_tot()<<" time steps; memory time: "<<get_n_tot()*dt<<std::endl;
    }
  }else{
    set_from_DiagBB_single_line(diagBB, dt, n_mem);
    if(print_info){
      if(n_mem>=n_tot){
        std::cout<<"No memory truncation"<<std::endl;
      }else{
        std::cout<<"Memory truncation: "<<n_mem<<" of "<<n_tot<<std::endl;
      }
    }
#ifdef DIAGBB_COMPRESS_SINGLE_LINE
  {  //Compress single line before combination
    if(print_info){
      std::cout<<"Compressing single line before contraction"<<std::endl;
    }
    sweep_forward(trunc, 0);
    sweep_backward(trunc, 0);
  }
#endif
  }
}

void ProcessTensorBuffer::set_from_DiagBB_fw(
                    DiagBB &diagBB, const TimeGrid &tgrid, 
                    TruncationLayout trunc, 
                    double dict_zero, int verbosity, int stop_at_row){

  trunc.keep = diagBB.get_dim();

//  bool do_check=false;
  bool debug=false;
//  bool debug=true;

  double dt = tgrid.dt;
  int n_tot_new = tgrid.n_tot;
  int n_mem = tgrid.n_mem;
  if(stop_at_row<1){
    stop_at_row=n_tot_new;
  }
  if(verbosity>1)trunc.print_info();

  if(verbosity>2){
    std::cout<<"n_tot="<<n_tot_new<<" n_mem="<<n_mem<<std::endl;
    std::cout<<"Calculating single line..."<<std::endl;
  }

  ProcessTensorBuffer PTB_line;
  PTB_line.set_from_DiagBB_single_line(diagBB, dt, n_mem);
#ifdef DIAGBB_COMPRESS_SINGLE_LINE
  PTB_line.sweep_forward(trunc.get_forward(0,stop_at_row), verbosity);
  PTB_line.sweep_backward(trunc.get_backward(0,stop_at_row), verbosity);
#endif

  copy_content(PTB_line);
  for(int line=1; line<stop_at_row; line++){
    ShiftExtend shift_extend;
    shift_extend.shift_second=line;
    shift_extend.truncate_at = n_tot_new;

    if(verbosity>0){
      std::cout<<"line: "<<line<<"/"<<n_tot_new<<std::endl;   
    }
    if(verbosity>0 || debug){
      std::cout<<"PTB1: "; print_info(); std::cout<<std::endl;
      if(debug){std::cout<<"line: "; PTB_line.print_info(); std::cout<<std::endl;}
    }


    {
      TruncatedSVD trunc_tmp=trunc.get_forward(line-1, stop_at_row-1);
      if(trunc.use_QR)trunc_tmp.use_QR=true;
      if(verbosity>0 ||debug){trunc_tmp.print_info();std::cout<<std::endl;}

      join_and_sweep_forward(PTB_line, trunc_tmp, verbosity, shift_extend);
    }
    {
      TruncatedSVD trunc_tmp=trunc.get_backward(line-1, stop_at_row-1);
      if(verbosity>0 ||debug){trunc_tmp.print_info();std::cout<<std::endl;}

      sweep_backward(trunc_tmp, verbosity, shift_extend.shift_second+1, -1);

      if(debug){std::cout<<"after PTB1 info: "; print_info(); std::cout<<std::endl;}
    }

    sweep_intermediate_or_final_start_forward(trunc, line-1, stop_at_row-1, verbosity);

if(debug)std::cout<<"line: "<<line<<" done."<<std::endl;   
  }
  apply_HilbertSpaceRotation(diagBB.hs_rot, dict_zero);
if(debug)std::cout<<"from_DiagBB_fw: done."<<std::endl;
}

void ProcessTensorBuffer::set_from_DiagBB(
                    DiagBB &diagBB, const TimeGrid &tgrid, 
                    TruncationLayout trunc, 
                    double dict_zero, int verbosity, int stop_at_row){

  trunc.keep = diagBB.get_dim();

//  bool do_check=false;
  bool debug=false;
//  bool debug=true;

  double dt = tgrid.dt;
  int n_tot_new = tgrid.n_tot;
  int n_mem = tgrid.n_mem;
  if(stop_at_row<1){
    stop_at_row=n_tot_new;
  }
  if(verbosity>1)trunc.print_info();

  if(verbosity>2){
    std::cout<<"n_tot="<<n_tot_new<<" n_mem="<<n_mem<<std::endl;
    std::cout<<"Calculating single line..."<<std::endl;
  }

  ProcessTensorBuffer PTB_line;
  PTB_line.set_from_DiagBB_single_line(diagBB, dt, n_mem);
  n_mem=PTB_line.get_n_tot();
#ifdef DIAGBB_COMPRESS_SINGLE_LINE
  PTB_line.sweep_backward(TruncatedSVD(), verbosity);
  PTB_line.sweep_forward(trunc.get_first(), verbosity);
#endif

  copy_content(PTB_line);
  for(int line=1; line<stop_at_row; line++){
    ShiftExtend shift_extend;
    shift_extend.shift_second=line;
    shift_extend.truncate_at = n_tot_new;

    if(verbosity>0){
      std::cout<<"line: "<<line<<"/"<<n_tot_new<<std::endl;   
    }
    if(verbosity>0 || debug){
      std::cout<<"PTB1: "; print_info(); std::cout<<std::endl;
      if(debug){std::cout<<"line: "; PTB_line.print_info(); std::cout<<std::endl;}
    }

    {
      TruncatedSVD trunc_tmp=trunc.get_backward(line-1, stop_at_row-1);
      if(trunc.use_QR)trunc_tmp.use_QR=true;
      if(verbosity>0 ||debug){trunc_tmp.print_info();std::cout<<std::endl;}

      join_and_sweep_backward(PTB_line, trunc_tmp, verbosity, shift_extend);

      if(debug){std::cout<<"after PTB1 info: "; print_info(); std::cout<<std::endl;}
    }

    {
      TruncatedSVD trunc_tmp=trunc.get_forward(line-1, stop_at_row-1);
      if(trunc.use_QR)trunc_tmp.use_QR=true;
      if(verbosity>0 ||debug){trunc_tmp.print_info();std::cout<<std::endl;}

      sweep_forward(trunc_tmp, verbosity, shift_extend.shift_second, -1);
    }

    sweep_intermediate_or_final_start_backward(trunc, line-1, stop_at_row-1, verbosity);
if(debug)std::cout<<"line: "<<line<<" done."<<std::endl;   
  }
  if(verbosity>0){std::cout<<"from_diagBB_finalizing"<<std::endl;}
  expand_DiagBB(diagBB, dict_zero);
}

void ProcessTensorBuffer::set_from_DiagBB_select(
                    DiagBB &diagBB, const TimeGrid &tgrid, 
                    TruncationLayout trunc, 
                    double dict_zero, int verbosity, int stop_at_row){

  trunc.keep = diagBB.get_dim();

//  bool do_check=false;
  bool debug=false;
//  bool debug=true;

  double dt = tgrid.dt;
  int n_tot_new = tgrid.n_tot;
  int n_mem = tgrid.n_mem;
  if(stop_at_row<1){
    stop_at_row=n_tot_new;
  }
  if(verbosity>1)trunc.print_info();

  if(verbosity>2){
    std::cout<<"n_tot="<<n_tot_new<<" n_mem="<<n_mem<<std::endl;
    std::cout<<"Calculating single line..."<<std::endl;
  }

  ProcessTensorBuffer PTB_line;
  PTB_line.set_from_DiagBB_single_line(diagBB, dt, n_mem);
  PTB_line.sweep_forward(TruncatedSVD(), verbosity);

  copy_content(PTB_line);

  for(int line=1; line<stop_at_row; line++){
    ShiftExtend shift_extend;
    shift_extend.shift_second=line;
    shift_extend.truncate_at = n_tot_new;

    if(verbosity>0){
      std::cout<<"line: "<<line<<"/"<<n_tot_new<<std::endl;   
    }

    {
      TruncatedSVD trunc_tmp=trunc.get_forward(line-1, stop_at_row-1);
      if(verbosity>0 ||debug){trunc_tmp.print_info();std::cout<<std::endl;}
      if(shift_extend.shift_second-2<0){
        sweep_forward(trunc_tmp, verbosity, 0);
      }else{
        sweep_forward(trunc_tmp, verbosity, shift_extend.shift_second-2);
      }
    }

    if(verbosity>0 || debug){
      std::cout<<"PTB1: "; print_info(); std::cout<<std::endl;
      if(debug){std::cout<<"line: "; PTB_line.print_info(); std::cout<<std::endl;}
    }

    TruncatedSVD trunc_select=trunc.get_backward(line-1, stop_at_row-1, true);
    TruncatedSVD trunc_bw=trunc.get_backward(line-1, stop_at_row-1);
    if(trunc.use_QR)trunc_bw.use_QR=true;
    if(verbosity>0){
      trunc_bw.print_info();std::cout<<std::endl;
    }  
    join_select_and_sweep_backward(PTB_line, trunc_select, trunc_bw, verbosity, shift_extend);
    if(verbosity>0){
    }


    if(shift_extend.shift_second>0){
      sweep_backward(trunc_bw, verbosity, shift_extend.shift_second-1, shift_extend.shift_second+1);
    }


if(debug){std::cout<<"after PTB1 info: "; print_info(); std::cout<<std::endl;}
    sweep_intermediate_or_final_start_forward(trunc, line-1, stop_at_row-1, verbosity, (line>n_mem)?line-n_mem:0, -1);
    
if(debug)std::cout<<"line: "<<line<<" done."<<std::endl;   
  }
  expand_DiagBB(diagBB, dict_zero);
if(debug)std::cout<<"from_DiagBB_select: done."<<std::endl;
}


void ProcessTensorBuffer::set_from_DiagBB_log(
                    DiagBB &diagBB, const TimeGrid &tgrid, 
                    TruncationLayout trunc, double dict_zero, 
                    int verbosity, int stop_at_row, int start_at_row,
                    bool copy_compress_second){
  trunc.keep = diagBB.get_dim();
std::cout<<"set_from_DiagBB_log: "; trunc.print_info(); std::cout<<std::endl;

//  bool do_check=false;
  bool debug=false;

  double dt = tgrid.dt;
  int n_tot_new = tgrid.n_tot;
  if(n_tot_new<1){
    std::cerr<<"set_from_DiagBB_log: n_tot_new<1!"<<std::endl;
    throw DummyException();
  }
  int n_mem = tgrid.n_mem;

  if(stop_at_row<1){
    stop_at_row=n_tot_new;
  }
  if(stop_at_row < 1 || ( (stop_at_row & (stop_at_row-1)) != 0 )){
    std::cerr<<"set_from_DiagBB_log: number of rows "<<stop_at_row<<" not a power of 2!"<<std::endl;
    throw DummyException();
  }
   

  if(verbosity>2){
    std::cout<<"n_tot="<<n_tot_new<<" n_mem="<<n_mem<<std::endl;
    std::cout<<"Calculating single line..."<<std::endl;
  }

//  int LOGN=0; while(pow(2,LOGN)<stop_at_row)LOGN++;
  int LOGN=0; while(pow(2,LOGN)<n_tot_new)LOGN++;

  if( start_at_row <= 0 ){
//  set_from_DiagBB_single_line(diagBB, dt, n_mem);
    TruncatedSVD this_trunc; this_trunc.keep=trunc.keep;
    set_from_DiagBB_single_line_auto(diagBB, dt, n_tot_new, n_mem, this_trunc, true);
    n_mem=get_n_tot();
    start_at_row=1;
  }else{
    if( start_at_row != round(pow(2., round(log(start_at_row)/log(2.)))) ){
      std::cerr<<"ProcessTensorBuffer::set_from_DiagBB_log: start_at_row="<<start_at_row<<" is not a power of 2!"<<std::endl;
      throw DummyException();
    }
    n_mem=get_n_tot();
  }

  for(int line=start_at_row; line<stop_at_row; line*=2){
    int logline=round(log(line)/log(2.));
    ShiftExtend shift_extend;
    shift_extend.shift_second=line;
    shift_extend.truncate_at = n_tot_new;

#if defined(DIAGBB_SWEEP_MORE)
    shift_extend.sweep_more = DIAGBB_SWEEP_MORE;
    std::cout<<"shift_extend.sweep_more="<<shift_extend.sweep_more<<std::endl;
#else
    if(n_mem<n_tot_new){
      shift_extend.sweep_more = n_mem;
      std::cout<<"shift_extend.sweep_more="<<shift_extend.sweep_more<<std::endl;
    }
#endif

    if(verbosity>0){
      std::cout<<"line: "<<line<<"/"<<n_tot_new<<std::endl;   
    }

/* if using memory cutoff: 
   sweep back touches elements [line-sweep_more; line+mem+sweep_more[ 
   => previous run only changed [line/2-sweep_more; line/2+mem+sweep_more[
*/
    if(shift_extend.sweep_more>0){
      int range_start=line/2-shift_extend.sweep_more;
      if(range_start<0)range_start=0;
      int range_end=line/2+n_mem+shift_extend.sweep_more;
      if(range_end<0||range_end>get_n_tot())range_end=-1;

      TruncatedSVD trunc_tmp=trunc.get_forward(logline, LOGN);
      if(verbosity>0 ||debug){trunc_tmp.print_info();std::cout<<std::endl;}

      sweep_forward(trunc_tmp, verbosity, range_start, range_end);
//      sweep_forward(trunc_tmp, verbosity);

    }else{
      TruncatedSVD trunc_tmp=trunc.get_forward(logline, LOGN);
      if(verbosity>0 ||debug){trunc_tmp.print_info();std::cout<<std::endl;}

      sweep_forward(trunc_tmp, verbosity);
    }

    ProcessTensorBuffer PTB2;
    if(copy_compress_second){
      PTB2.set_new_temporary(*this);

      PTB2.clear();
      int n_tot_new2=n_tot_new-line;
      if(this->n_tot<n_tot_new2){n_tot_new2=this->n_tot;}
      for(int n=0; n<n_tot_new2; n++){
         PTB2.push_back(this->get(n));
      }
      PTB2.get(n_tot_new2-1).close_off();

      TruncatedSVD trunc_tmp=trunc.get_backward(logline, LOGN);
      if(verbosity>0 ||debug){std::cout<<"PTB2.n_tot="<<PTB2.n_tot<<std::endl;
                              trunc_tmp.print_info();std::cout<<std::endl;}

      sweep_backward(trunc_tmp, verbosity);
      sweep_forward(trunc_tmp, verbosity);
    }else{
      PTB2.copy_read_only(*this);
//      std::cout<<"PTB2: "; PTB2.print_info(); std::cout<<std::endl;
    }

    TruncatedSVD trunc_select=trunc.get_backward(logline,LOGN, true);
    TruncatedSVD trunc_bw=trunc.get_backward(logline,LOGN);
    if(trunc.use_QR)trunc_bw.use_QR=true;
    if(verbosity>0){
      std::cout<<"PTB1: "; print_info(); std::cout<<std::endl;
      std::cout<<"trunc_select: ";
      trunc_select.print_info(); std::cout<<std::endl;
      std::cout<<"trunc_bw: ";
      trunc_bw.print_info(); std::cout<<std::endl;
    }
    { 
      join_select_and_sweep_backward(PTB2, trunc_select, trunc_bw, verbosity, shift_extend);
    }


  sweep_intermediate_or_final_start_forward(trunc, logline, LOGN, verbosity);  

  }

  if(verbosity>0){std::cout<<"finalizing"<<std::endl;}
  expand_DiagBB(diagBB, dict_zero);

if(debug)std::cout<<"from_DiagBB_log: done."<<std::endl;
//print_dims();std::cout<<std::endl;

}



void ProcessTensorBuffer::set_from_ModePropagator(
            ModePropagator &mprop, const TimeGrid &tgrid, double dict_zero){

  clear(); 
  int n_max=tgrid.n_tot;

//  int N=mprop.get_N_system();
  int N_mode=mprop.get_N_mode();
 
  ProcessTensorElement element;
  for(int n=0; n<n_max; n++){
    element.set_from_ModePropagator(mprop, tgrid.get_t(n), tgrid.get_dt(n), dict_zero);
    if(n==0){
      Eigen::VectorXcd bath_init=H_Matrix_to_L_Vector(mprop.get_bath_init());
      element.M.inner_multiply_left(bath_init.transpose());
    }
    if(n==n_max-1){
      Eigen::VectorXcd Tr=H_Matrix_to_L_Vector(Eigen::MatrixXcd::Identity(N_mode,N_mode));
      element.closure=Eigen::VectorXcd::Ones(1);
      element.env_ops.set_ill_defined();
      element.M.inner_multiply_right(Tr);
    }
    push_back(element);
  } 
}

void ProcessTensorBuffer::add_modes(
          ModePropagatorGenerator &mpg, const TimeGrid &tgrid,
          TruncationLayout trunc, 
          double dict_zero, int verbosity){
  
  int n_max=tgrid.n_tot;
  if(n_max!=n_tot){
    std::cerr<<"ProcessTensorBuffer::add_modes: n_max!=n_tot!"<<std::endl;
    throw DummyException();
  }
  if(n_max<1){
    return;
  }
  int N=mpg.get_N();
  if(N!=get(0).get_N()){
    std::cerr<<"ProcessTensorBuffer::add_modes: N!=get(0).get_N()!"<<std::endl;
    throw DummyException();
  }
  if(verbosity>0){
    std::cout<<"Calculating PT for Generator '"<<mpg.name()<<"'"<<std::endl;
  } 

  TimeGrid tgrid2=tgrid.construct_half_dt();
  
//std::cout<<"-----"<<std::endl<<"TEST "<<trunc.base_threshold<<std::endl<<"-----"<<std::endl;
  trunc.print_info();std::cout<<std::endl;
 
  for(int k=mpg.first(); k<mpg.get_N_modes(); k=mpg.next(k)){
    if(verbosity>0){
      std::cout<<"Mode "<<k<<"/"<<mpg.get_N_modes()<<std::endl;
    }
    ModePropagatorPtr mpp=mpg.getModePropagator(k);

    ProcessTensorBuffer PTB;
    PTB.set_new_temporary(*this);
    PTB.set_from_ModePropagator(*mpp.get(), tgrid2, dict_zero);

    TruncatedSVD trunc_fw=trunc.get_forward(k, mpg.get_N_modes());
    if(trunc.use_QR)trunc_fw.use_QR=true;
    TruncatedSVD trunc_bw=trunc.get_backward(k, mpg.get_N_modes());

    if(verbosity>0){trunc_fw.print_info();std::cout<<std::endl;}
    join_symmetric_and_sweep_forward(PTB, trunc_fw, verbosity);

    if(verbosity>0){trunc_bw.print_info();std::cout<<std::endl;}
    sweep_backward(trunc_bw, verbosity);
   
    sweep_intermediate_or_final_start_forward(trunc, k, mpg.get_N_modes(), verbosity);
  }
}

void ProcessTensorBuffer::add_modes_firstorder(
          ModePropagatorGenerator &mpg, const TimeGrid &tgrid,
          TruncationLayout trunc, 
          double dict_zero, int verbosity){
  
  int n_max=tgrid.n_tot;
  if(n_max!=n_tot){
    std::cerr<<"ProcessTensorBuffer::add_modes: n_max!=n_tot!"<<std::endl;
    throw DummyException();
  }
  if(n_max<1){
    return;
  }
  int N=mpg.get_N();
  if(N!=get(0).get_N()){
    std::cerr<<"ProcessTensorBuffer::add_modes: N!=get(0).get_N()!"<<std::endl;
    throw DummyException();
  }
  if(verbosity>0){
    std::cout<<"Calculating PT for Generator '"<<mpg.name()<<"'"<<std::endl;
  }
  
  trunc.print_info();
 
  for(int k=mpg.first(); k<mpg.get_N_modes(); k=mpg.next(k)){
    if(verbosity>0){
      std::cout<<"Mode "<<k<<"/"<<mpg.get_N_modes()<<std::endl;
    }
    ModePropagatorPtr mpp=mpg.getModePropagator(k);

    ProcessTensorBuffer PTB;
    PTB.set_new_temporary(*this);
    PTB.set_from_ModePropagator(*mpp.get(), tgrid, dict_zero);

    TruncatedSVD trunc_fw=trunc.get_forward(k, mpg.get_N_modes());
    if(trunc.use_QR)trunc_fw.use_QR=true;
    TruncatedSVD trunc_bw=trunc.get_backward(k, mpg.get_N_modes());

    if(verbosity>0){trunc_fw.print_info();std::cout<<std::endl;}
    join_and_sweep_forward(PTB, trunc_fw, verbosity);

    if(verbosity>0){trunc_bw.print_info();std::cout<<std::endl;}
    sweep_backward(trunc_bw, verbosity);
   
    sweep_intermediate_or_final_start_forward(trunc, k, mpg.get_N_modes(), verbosity);
  }
}

void ProcessTensorBuffer::add_modes_select(
          ModePropagatorGenerator &mpg, const TimeGrid &tgrid,
          TruncationLayout trunc, 
          double dict_zero, int verbosity){
  
  int n_max=tgrid.n_tot;
  if(n_max!=n_tot){
    std::cerr<<"ProcessTensorBuffer::add_modes: n_max!=n_tot!"<<std::endl;
    throw DummyException();
  }
  if(n_max<1){
    return;
  }
  int N=mpg.get_N();
  if(N!=get(0).get_N()){
    std::cerr<<"ProcessTensorBuffer::add_modes: N!=get(0).get_N()!"<<std::endl;
    throw DummyException();
  }
  if(verbosity>0){
    std::cout<<"Calculating PT for Generator '"<<mpg.name()<<"'"<<std::endl;
  }

  for(int k=mpg.first(); k<mpg.get_N_modes(); k=mpg.next(k)){
    if(verbosity>0){
      std::cout<<"Mode "<<k<<"/"<<mpg.get_N_modes()<<std::endl;
    }
    ModePropagatorPtr mpp=mpg.getModePropagator(k);

    ProcessTensorBuffer PTB;
    PTB.set_new_temporary(*this);
    PTB.set_from_ModePropagator(*mpp.get(), tgrid, dict_zero);

    TruncatedSVD trunc_fw=trunc.get_forward(k, mpg.get_N_modes());
    TruncatedSVD trunc_bw=trunc.get_backward(k, mpg.get_N_modes());
    if(trunc.use_QR)trunc_bw.use_QR=true;
    TruncatedSVD trunc_select=trunc.get_backward(k, mpg.get_N_modes(),true);

    if(verbosity>0){trunc_fw.print_info();std::cout<<std::endl;}
    sweep_forward(trunc_fw, verbosity);
    PTB.sweep_forward(trunc_fw, verbosity);

    if(verbosity>0){trunc_bw.print_info();std::cout<<std::endl;}
    join_select_and_sweep_backward(PTB, trunc_select, trunc_bw, verbosity);
   
    sweep_intermediate_or_final_start_forward(trunc, k, mpg.get_N_modes(), verbosity);
  }
}

void ProcessTensorBuffer::add_modes_tree_get(int level, int max_level, 
          int first_elem, ModePropagatorGenerator &mpg, const TimeGrid &tgrid,
          const TruncationLayout & trunc, 
          double dict_zero, int verbosity){

  if(level==0 && first_elem==0){
    trunc.print_info(); std::cout<<std::endl;
  }

  if(level==0){
    if(first_elem>=mpg.get_N_modes()){
      std::cerr<<"add_modes_tree_get: level="<<level<<"/"<<max_level<<" first_elem="<<first_elem<<": first_element>=mpg.get_N_modes()="<<mpg.get_N_modes()<<"!"<<std::endl;
      throw DummyException();
    }
    
    if(verbosity>0){
      std::cout<<"--------------------------------------------"<<std::endl;
      std::cout<<"level: "<<level<<"/"<<max_level<<" first_elem: "<<first_elem<<"/"<<pow(2,max_level-1-level)<<std::endl;
      std::cout<<"--------------------------------------------"<<std::endl;
    }
    if(mpg.skip_list[first_elem]){
      set_trivial(tgrid.n_tot, mpg.get_N());
    }else{
      ModePropagatorPtr mpp=mpg.getModePropagator(first_elem);
      set_from_ModePropagator(*mpp.get(), tgrid, dict_zero);

      TruncatedSVD trunc_bw;
//      trunc_bw.use_QR=true;
      if(verbosity>0){
        std::cout<<"Sweep backward"<<std::endl;
        trunc_bw.print_info(); std::cout<<std::endl;
      }
      sweep_backward(trunc_bw, verbosity);
    }
  }else{  
    clear();
    add_modes_tree_get(level-1, max_level, 2*first_elem, mpg, tgrid, trunc, dict_zero, verbosity);

    ProcessTensorBuffer PTB2;
    PTB2.set_new_temporary(*this);
    PTB2.add_modes_tree_get(level-1, max_level, 2*first_elem+1, mpg, tgrid, trunc, dict_zero, verbosity);

    if(verbosity>0){
      std::cout<<"--------------------------------------------"<<std::endl;
      std::cout<<"level: "<<level<<"/"<<max_level<<" first_elem: "<<first_elem<<"/"<<pow(2,max_level-1-level)<<std::endl;
      std::cout<<"--------------------------------------------"<<std::endl;
    }

    TruncatedSVD trunc_fw=trunc.get_forward(level, max_level);
    if(verbosity>0){
      std::cout<<"sweep first forward"<<std::endl;
      trunc_fw.print_info(); std::cout<<std::endl;
    }
    sweep_forward(trunc_fw, verbosity);

    if(verbosity>0){
      std::cout<<"sweep second forward"<<std::endl;
      trunc_fw.print_info(); std::cout<<std::endl;
    }
    PTB2.sweep_forward(trunc_fw, verbosity);

    if(verbosity>0)std::cout<<"join and sweep backward"<<std::endl;
    TruncatedSVD trunc_select=trunc.get_backward(level, max_level, true);
    TruncatedSVD trunc_bw=trunc.get_backward(level, max_level);
    if(trunc.use_QR)trunc_bw.use_QR=true;
    trunc_select.print_info(); std::cout<<std::endl;
    trunc_bw.print_info(); std::cout<<std::endl;
    join_select_and_sweep_backward(PTB2, trunc_select, trunc_bw, verbosity);

#ifdef TEST_SWEEP_PAIR
std::cout<<"TEST_SWEEP_PAIR is set!"<<std::endl;
//for(int n=0; n<get_n_tot(); n++){std::cout<<"n="<<n<<": "; get(n).printNF();}
    TruncatedSVD trunc_pair=trunc_fw;
    sweep_pair_forward(trunc_pair, verbosity);
//for(int n=0; n<get_n_tot(); n++){std::cout<<"n="<<n<<": "; get(n).printNF();}
    sweep_pair_backward(trunc_pair, verbosity);
//    calculate_closures();
//    sweep_backward(trunc_fw, verbosity);
#endif
    sweep_intermediate_or_final_start_forward(trunc, level, max_level, verbosity);
  }
}

void ProcessTensorBuffer::set_from_modes_tree(
          ModePropagatorGenerator &mpg, const TimeGrid &tgrid,
          TruncationLayout trunc,
          double dict_zero, int verbosity){

  int N_modes=mpg.get_N_modes();
  if(N_modes<1)return;

  int N_hierarchy=round(log((double)N_modes)/log(2.));
  if(N_modes!=pow(2, N_hierarchy-1)){
    if(N_modes<pow(2, N_hierarchy-1)){
      mpg.zero_pad(pow(2, N_hierarchy-1));
    }else{
      mpg.zero_pad(pow(2, N_hierarchy));
    }
    N_modes=mpg.get_N_modes();
  }
/*
  int N_hierarchy=1;
  {int N_shift=N_modes;
    while(N_shift>1){
      N_hierarchy++;
      N_shift=N_shift>>1;
    }
  }
  if(N_modes!=pow(2, N_hierarchy-1)){
    std::cerr<<"ProcessTensorBuffer::add_modes_tree: N_modes="<<N_modes<<" != 2^"<<N_hierarchy-1<<"="<<pow(2, N_hierarchy-1)<<" => not a power of 2!"<<std::endl;
    throw DummyException();
  } 
*/
  if(verbosity>0)std::cout<<"N_modes="<<N_modes<<"=2^"<<N_hierarchy-1<<std::endl;

  add_modes_tree_get(N_hierarchy-1, N_hierarchy, 0, mpg, tgrid, trunc, 
                          dict_zero, verbosity);
 
}
void ProcessTensorBuffer::add_modes_tree(
          ModePropagatorGenerator &mpg, const TimeGrid &tgrid,
          TruncationLayout trunc, 
          double dict_zero, int verbosity){

  if(n_tot != tgrid.n_tot){
    std::cerr<<"ProcessTensorBuffer::add_modes_tree: n_tot != tgrid.n_tot ("<<n_tot<<" vs. "<<tgrid.n_tot<<")!"<<std::endl;
    throw DummyException();
  }

  TruncatedSVD trunc_tmp=trunc.get_base();
  if(verbosity>0){
    trunc_tmp.print_info(); std::cout<<std::endl;
  }
  sweep_forward(trunc_tmp, verbosity);
  ProcessTensorBuffer PTB2;
  PTB2.set_new_temporary(*this);
  PTB2.set_from_modes_tree(mpg, tgrid, trunc, dict_zero, verbosity);

  if(verbosity>0){
    std::cout<<"------------------------"<<std::endl;
    std::cout<<"combine_tree: last step:"<<std::endl;
    std::cout<<"------------------------"<<std::endl;
  }


  PTB2.sweep_forward(trunc_tmp, verbosity);

  TruncatedSVD trunc_select=trunc.get_backward(0,0,true);
  TruncatedSVD trunc_bw=trunc.get_backward(0,0);
  if(trunc.use_QR)trunc_bw.use_QR=true;
  if(verbosity>0){
    trunc_select.print_info(); std::cout<<std::endl;
  }
  join_select_and_sweep_backward(PTB2, trunc_select, trunc_bw, verbosity);
 
  sweep_intermediate_or_final_start_forward(trunc, 0, 1, verbosity);
}

void ProcessTensorBuffer::set_from_coarse_grain(ProcessTensorBuffer &PTB, int coarse_grain){

  if(PTB.get_n_tot()<1){
    std::cerr<<"PTB.get_n_tot()<1!"<<std::endl;
    throw DummyException();
  }
  if(PTB.get_n_tot()%coarse_grain!=0 &&coarse_grain<1){
    std::cerr<<"PTB.get_n_tot()\%coarse_grain!=0"<<std::endl;
    throw DummyException();
  }


  int n_tot_old=PTB.get_n_tot();
  int n_tot_new=PTB.get_n_tot()/coarse_grain;
  resize(n_tot_new);

  for(int n2=0; n2<n_tot_new; n2++){
    ProcessTensorElement &e2=get(n2,ForwardPreload);
    e2=PTB.get(n2*coarse_grain,ForwardPreload);

    for(int c=1; c<coarse_grain; c++){
      int n=n2*coarse_grain+c;

      ProcessTensorElement &e=PTB.get(n,ForwardPreload);
      IF_OD_Dictionary dict_first=e2.accessor.dict;
      IF_OD_Dictionary dict_second=e.accessor.dict;
      IF_OD_Dictionary dict_new=dict_first; dict_new.join(dict_second);
      int NL=dict_new.get_NL();
      std::vector<std::vector<int> > rev_new=dict_new.get_reverse_beta();
 
      MPS_Matrix M(dict_new.get_reduced_dim(), e2.M.dim_d1, e.M.dim_d2);
      M.set_zero();
      for(int k=0; k<NL; k++){
        for(int j=0; j<NL; j++){
          int i_ind=dict_first.beta[j*NL+k]; if(i_ind<0)continue;
          for(int i=0; i<NL; i++){
            int i_ind2=dict_second.beta[i*NL+j]; if(i_ind2<0)continue;
            int i_ind3=dict_new.beta[i*NL+k]; if(i_ind3<0)continue;
            if(rev_new[i_ind3][0]!=i*NL+k)continue; //only modify first occurance
            for(int d1=0; d1<e2.M.dim_d1; d1++){
              for(int d2=0; d2<e2.M.dim_d2; d2++){
                for(int d3=0; d3<e.M.dim_d2; d3++){
                  M(i_ind3, d1, d3)+= e2.M(i_ind, d1, d2)*e.M(i_ind2, d2, d3);
                }
              }
            }
          }
        }
      }
      e2.M.swap(M);
      e2.accessor.dict=dict_new;
      e2.closure=e.closure;
      e2.env_ops=e.env_ops;
      e2.forwardNF=e.forwardNF;
    }
  } 
}

void ProcessTensorBuffer::clean_up(){
//  std::cout<<"clean up: info: "; print_info(); std::cout<<std::endl;
//bool debug=true;
bool debug=false;
if(debug){std::cout<<"CLEAN UP: ";if(read_only)std::cout<<"read_only=true";else std::cout<<"read_only=false"; std::cout<<std::endl;}

  if(use_file){
    if(is_temporary){ 
if(debug)std::cout<<"temporary"<<std::endl;
      delete_files();
    }else if(read_only){
if(debug)std::cout<<"read_only"<<std::endl;
    }else if(use_single_file){
if(debug)std::cout<<"single file"<<std::endl;
      write_release_block(0);
    }else if(current_block>=0){
if(debug)std::cout<<"current_block="<<current_block<<std::endl;
      write_release_block(current_block);
      write_header();
      wait_write();
    }
  }
}


//Initializers
void ProcessTensorBuffer::initialize(const ProcessTensorBufferSpec & spec){
  ProcessTensorBufferSpec::copy(spec);
  n_tot=0; n=0;
  current_block=-1;
  clear_preload();
  preload_block=-1;
  clear_buffer();
  if(is_temporary && fname_header==""){
    TempFileName tmpname; fname_header=tmpname; tmpname.fname="";
  }
}
void ProcessTensorBuffer::set_new_single_file(const std::string &fname){
  if(fname==""){
    std::cerr<<"ProcessTensorBuffer::set_new_single_file: fname is empty!"<<std::endl;
    throw DummyException();
  }
  use_file=true;
  use_single_file=true;
  is_temporary=false;
  blocksize=-1;
  fname_header=fname;
  clear_preload();
  clear_write_buffer();
  current_block=-1;
}
void ProcessTensorBuffer::set_new_file(const std::string &fname, int blocksize_){
  if(fname==""){
    std::cerr<<"ProcessTensorBuffer::set_new_file: fname is empty!"<<std::endl;
    throw DummyException();
  }
  if(blocksize_<1){
    set_new_single_file(fname);
    return;
  }
  initialize();
  blocksize=blocksize_;
  use_file=true;
  is_temporary=false;
  use_single_file=false;
  blocksize=blocksize_;
  fname_header=fname;
  clear_preload();
  clear_write_buffer();
  write_all();
}

void ProcessTensorBuffer::set_new_temporary(int blocksize_){
  std::string tmpfile;
  {TempFileName tmpname; tmpfile=tmpname; tmpname.fname="";}

  set_new_file(tmpfile, blocksize_);
  is_temporary=true;
  clear_preload();
  clear_write_buffer();
  current_block=-1;
  write_all();
}

void ProcessTensorBuffer::set_new_temporary(const ProcessTensorBufferSpec &other){
  ProcessTensorBufferSpec::copy(other);
  if(use_file){
    TempFileName tmpname; fname_header=tmpname; tmpname.fname="";
  }
  is_temporary=true;
  clear_preload();
  clear_write_buffer();
  current_block=-1;
  write_all();
}

void ProcessTensorBuffer::set_preload_none(){
  wait_preload();
  clear_preload();
  preload_block=-1;
  preload_lock=false;
}
void ProcessTensorBuffer::set_write_buffer_none(){
  if(use_async_write){
    wait_write();
  }
  clear_write_buffer();
  write_buffer_block=-1;
  write_buffer_lock=false;
}
void ProcessTensorBuffer::expand_outer(int N_front, int N_back){
  if(N_front<2 && N_back<2){
    std::cerr<<"ProcessTensorBuffer::expand_outer: N_front<2 && N_back<2: nothing to extend"<<std::endl;
    throw DummyException();
  }
  for(int n=0; n<get_n_tot(); n++){
    ProcessTensorElement &e=get(n, ForwardPreload);
    if(N_front>1)e.expand_space_front(N_front);
    if(N_back>1)e.expand_space_back(N_back);
  }
  
}

//Inherited from ProcessTensorForward:

const ProcessTensorElement * ProcessTensorBuffer::current(){
  return &get(ProcessTensorForward::n, ForwardPreload);
}

}
