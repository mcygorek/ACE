#include "ProcessTensorRepeat.hpp"
#include "DummyException.hpp"
#include "BinaryReader.hpp"

namespace ACE{

ProcessTensorElement & ProcessTensorRepeat::get(int n_, PreloadHint hint){
  if(n_<(int)initial.size()){
    return initial.get(n_, hint);
  }
  if(repeated.size()<1){
    std::cerr<<"ProcessTensorRepeat::get("<<n_<<"): Repeated block is empty!"<<std::endl;
    throw DummyException();
  }
  return repeated.get((n_-(int)initial.size()) % ((int)repeated.size()), hint);
}

const ProcessTensorElement & ProcessTensorRepeat::get_ro(int n_, PreloadHint hint){
  if(n_<(int)initial.size()){
    return initial.get_ro(n_, hint);
  }
  if(repeated.size()<1){
    std::cerr<<"ProcessTensorRepeat::get("<<n_<<"): Repeated block is empty!"<<std::endl;
    throw DummyException();
  }
  return repeated.get_ro((n_-(int)initial.size()) % ((int)repeated.size()), hint);
}

int ProcessTensorRepeat::get_N_system(){
  int this_n=ProcessTensorForward::n; 
  if(this_n>=n_tot){ this_n=0; }
  if(n_tot<1){
    std::cerr<<"ProcessTensorRepeat::get_N_system(): n_tot<1!"<<std::endl;
    throw DummyException();
  }
  const ProcessTensorElement &e=get_ro(this_n);
  return e.get_N();
}


int ProcessTensorRepeat::get_n_mem(const TimeGrid &tgrid){
  int n_mem_tmp=tgrid.n_mem;
  if(n_mem_tmp<1)n_mem_tmp=tgrid.n_tot;
  int n_mem=pow(2, ceil(log((double)n_mem_tmp)/log(2.)) );
  return n_mem;
}

bool ProcessTensorRepeat::fall_back_to_DiagBB_log(const TimeGrid &tgrid, int verbosity){
  if(2*get_n_mem(tgrid) > tgrid.n_tot){
    if(verbosity>0){
      std::cout<<"2*n_mem="<<2*get_n_mem(tgrid)<<" > n_tot="<<tgrid.n_tot<<" => falling back to use_Gaussian_log"<<std::endl;
    }
    return true;
  }else{
    if(verbosity>0){
      std::cout<<"2*n_mem="<<2*get_n_mem(tgrid)<<" <= n_tot="<<tgrid.n_tot<<" => use_Gaussian_repeat"<<std::endl;
    }
    return false;
  }
}

bool ProcessTensorRepeat::is_read_only(bool val)const{
  return (initial.read_only && repeated.read_only);
}
void ProcessTensorRepeat::set_read_only(bool val){
  initial.set_read_only(val);
  repeated.set_read_only(val);
}
void ProcessTensorRepeat::set_specs(const std::string &write_PT, int blocksize){
  initial.clear(); repeated.clear();
  set_read_only(false);
  initial.was_modified=repeated.was_modified=false;
  initial.blocksize=repeated.blocksize=blocksize;
  if(write_PT==""){
    if(blocksize<1){
      initial.fname_header=repeated.fname_header="";
    }else{
      initial.fname_header=TempFileName().set_noremove(true).get();
      repeated.fname_header=TempFileName().set_noremove(true).get();
      initial.on_exit=repeated.on_exit=BufferedContainer<ProcessTensorElement>::ON_EXIT::DeleteOnDestruction;
    }
  }else{
    initial.fname_header=write_PT+std::string("_initial");
    repeated.fname_header=write_PT+std::string("_repeated");
    initial.on_exit=BufferedContainer<ProcessTensorElement>::ON_EXIT::WriteOnDestruction;
    repeated.on_exit=BufferedContainer<ProcessTensorElement>::ON_EXIT::WriteOnDestruction;
  }
}

void ProcessTensorRepeat::calculate(DiagBB &diagBB, const TimeGrid &tgrid,
          TruncationLayout trunc, double dict_zero, int verbosity, bool use_log){

  trunc.keep = diagBB.get_dim();
  int n_mem=get_n_mem(tgrid);

  TimeGrid tgrid2(2*n_mem, tgrid.dt, 0);//tgrid.ta);
  tgrid2.n_mem=n_mem;
  
  if(verbosity>0){
    std::cout<<"Calculating ProcessTensorRepeat with n_mem="<<n_mem<<std::endl;
  }

  ProcessTensorBuffer PTB; 
  if(initial.blocksize>0){
    PTB.set_new_temporary(initial.blocksize);
  }

/*
  //modify TruncationLayout to use the consistant threshold range
  TruncationLayout trunc2=trunc;
  {  
    int LOGN=0; while(pow(2,LOGN)<n_mem)LOGN++;
    //Derivation: set 1/r^(1-i/N) = 1/(r')^(1-i/(N-1))*f
    double f=1./pow(trunc.threshold_range_factor, 1./(LOGN));
    trunc2.base_threshold*=f;
    trunc2.final_sweep_threshold*=f;
    trunc2.threshold_range_factor=pow(trunc.threshold_range_factor, (LOGN-1.)/((double)LOGN));
  }
*/

  TruncatedSVD trunc_last=trunc.get_base();
  if(use_log){
//    PTB.set_from_DiagBB_log(diagBB, tgrid2, trunc2, dict_zero, verbosity, n_mem);
    PTB.set_from_DiagBB_log(diagBB, tgrid2, trunc, dict_zero, verbosity, n_mem);
  }else{
//    PTB.set_from_DiagBB(diagBB, tgrid2, trunc2, dict_zero, verbosity, n_mem);
    PTB.set_from_DiagBB(diagBB, tgrid2, trunc, dict_zero, verbosity, n_mem);
    PTB.sweep_backward(trunc_last, verbosity);
  }

  if(verbosity>0){
    std::cout<<"Composing periodic part"<<std::endl;
    trunc_last.print_info(); std::cout<<std::endl;
  }
  PTB.sweep_forward(trunc_last, verbosity);
  ProcessTensorBuffer PTB2;
  PTB2.copy_read_only(PTB);
  ProcessTensorBuffer::ShiftExtend shift_extend;
  shift_extend.shift_second=n_mem;
  shift_extend.truncate_at =3*n_mem;
  shift_extend.sweep_more=0;

  TruncatedSVD trunc_select=trunc.get_backward(-1,-1,true);
  TruncatedSVD trunc_bw=trunc.get_backward(-1,-1);
  if(trunc.use_QR)trunc_bw.use_QR=true;
  if(verbosity>0){
    trunc_select.print_info(); std::cout<<std::endl;
    trunc_bw.print_info(); std::cout<<std::endl;
  }
  PTB.join_select_and_sweep_backward(PTB2, trunc_select, trunc_bw, verbosity, shift_extend);
  PTB.sweep_forward(trunc_last, verbosity, n_mem, 2*n_mem);
  PTB.calculate_closures();
 
  if(verbosity>0){
    std::cout<<"n_mem="<<n_mem<<" PTB.get_n_tot()="<<PTB.get_n_tot()<<std::endl;
  }


  initial.resize(n_mem);
  for(int i=0; i<n_mem; i++){ 
    initial.get(i,ForwardPreload).swap(PTB.get(i,ForwardPreload));
  }

  repeated.resize(n_mem);
  for(int i=0; i<n_mem; i++){ 
    repeated.get(i,ForwardPreload).swap(PTB.get(i+n_mem,ForwardPreload));
  }
  n=0;

  sweep_final_start_backward(trunc, verbosity);
}

void ProcessTensorRepeat::sweep_backward(TruncatedSVD &trunc, bool flip_last, bool sweep_initial, int verbosity){

bool debug=false;
  if(repeated.size()<1 || initial.size()<1){
    std::cerr<<"ProcessTensorRepeat::sweep_backward: PTR is empty!"<<std::endl;
    throw DummyException();
  }

  int last_dim_d2=repeated.get((int)repeated.size()-1, BackwardPreload).M.dim_d2;
  PassOn pass_on=PassOn(last_dim_d2);
  int maxdim_in=last_dim_d2;
  int maxdim_out=last_dim_d2;
  int maxdim_at=repeated.size();
  for(int n=(int)repeated.size()-1; n>0; n--){
   try{
    ProcessTensorElement & element = repeated.get(n, BackwardPreload);

//std::cout<<"before: n="<<n<<" element.M.dim_d1="<<element.M.dim_d1<<" element.M.dim_d2="<<element.M.dim_d2<<std::endl;
    if(element.M.dim_d1>maxdim_in)maxdim_in=element.M.dim_d1;
if(debug){std::cout<<"repeated["<<n<<"].sweep_backward"<<std::endl;}
    element.sweep_backward(trunc, pass_on, (n==0));
    if(element.M.dim_d1>maxdim_out){
      maxdim_out=element.M.dim_d1;
      maxdim_at=n+initial.size();
    }
//std::cout<<"after: n="<<n<<" element.M.dim_d1="<<element.M.dim_d1<<" element.M.dim_d2="<<element.M.dim_d2<<std::endl;
 
   }catch(DummyException &e){
     std::cerr<<"called by ProcessTensorRepeat::sweep_backward: repeated["<<n<<"]"<<std::endl;
     throw e;
   }
  }

 try{
  ProcessTensorElement & element = repeated.get(0);

  if(element.M.dim_d1>maxdim_in)maxdim_in=element.M.dim_d1;
  element.sweep_backward(trunc, pass_on, !flip_last);
  if(element.M.dim_d1>maxdim_out){
    maxdim_out=element.M.dim_d1;
    maxdim_at=n+initial.size();
  }
  if(!flip_last){
//    if(verbosity>0)std::cout<<"Maxdim at n="<<maxdim_at<<": "<<maxdim_in<<" -> "<<maxdim_out<<std::endl;
    pass_on=PassOn(element.M.dim_d1);
  }else{ 
  //sweep first in "repeated" pass_on to both: last in "repeated" initial"
    PassOn pass_on2=pass_on; //keep original for sweeping "initial"
if(debug){std::cout<<"flip_last: sweep_backward"<<std::endl;}
    repeated.get((int)repeated.size()-1).sweep_backward(trunc, pass_on2, true);
    initial.get((int)initial.size()-1).sweep_backward(trunc, pass_on, true);
  }
 }catch(DummyException &e){
   std::cerr<<"called by ProcessTensorRepeat::sweep_backward: repeated["<<0<<"]"<<std::endl;
   throw e;
 }

  if(sweep_initial){
    pass_on=PassOn(initial.get((int)initial.size()-1).M.dim_d2);
    for(int n=(int)initial.size()-1; n>=0; n--){
     try{
      ProcessTensorElement & element = initial.get(n, BackwardPreload);

      if(element.M.dim_d1>maxdim_in)maxdim_in=element.M.dim_d1;
if(debug){std::cout<<"initial["<<n<<"].sweep_backward"<<std::endl;}
      element.sweep_backward(trunc, pass_on, (n==0));
      if(element.M.dim_d1>maxdim_out){
        maxdim_out=element.M.dim_d1;
        maxdim_at=n;
      }
     }catch(DummyException &e){
       std::cerr<<"called by ProcessTensorRepeat::sweep_backward: initial["<<n<<"]"<<std::endl;
       throw e;
     }
    }
  }
  if(verbosity>0)std::cout<<"Maxdim at n="<<maxdim_at<<": "<<maxdim_in<<" -> "<<maxdim_out<<std::endl;
}

void ProcessTensorRepeat::sweep_forward(TruncatedSVD &trunc, bool flip_last, bool sweep_initial, int verbosity){
  if(repeated.size()<1 || initial.size()<1){
    std::cerr<<"ProcessTensorRepeat::sweep_forward: PTR is empty!"<<std::endl;
    throw DummyException();
  }

/*
  if(flip_last){ //<-Would this even make sense?
    std::cerr<<"ProcessTensorRepeat::sweep_forward: flip_last not implemented!"<<std::endl;
    throw DummyException();
  }
*/
  int first_dim_d1 = sweep_initial ? initial.get(0, ForwardPreload).M.dim_d1
                                   : repeated.get(0, ForwardPreload).M.dim_d1;
  PassOn pass_on=PassOn(first_dim_d1);
  int maxdim_in=first_dim_d1;
  int maxdim_out=first_dim_d1;
  int maxdim_at=sweep_initial ? 0 : initial.size();
  
  if(sweep_initial){
    for(int n=0; n<initial.size(); n++){
     try{
      ProcessTensorElement & element = initial.get(n, ForwardPreload);

      if(element.M.dim_d1>maxdim_in)maxdim_in=element.M.dim_d1;
      element.sweep_forward(trunc, pass_on, (n==(int)initial.size()-1) && !flip_last);
      if(element.M.dim_d1>maxdim_out){
        maxdim_out=element.M.dim_d1;
        maxdim_at=n;
      }
     }catch(DummyException &e){
       std::cerr<<"called by ProcessTensorRepeat::sweep_forward: initial["<<n<<"]"<<std::endl;
       throw e;
     } 
    }
  }

  if(flip_last){
    try{
      //keep pass_on for the first step in iteration over "repeated"
      //have to act on the backward side of the last element of "repeated"
      PassOn pass_on2=pass_on.get_reversed(); 
      repeated.get((int)repeated.size()-1).sweep_backward(trunc, pass_on2, true);
    }catch(DummyException &e){
      std::cerr<<"called by ProcessTensorRepeat::sweep_forward: flip_last"<<std::endl;
      throw e;
    } 
  }
 
  for(int n=0; n<repeated.size(); n++){
    try{
      ProcessTensorElement & element = repeated.get(n, ForwardPreload);

      if(element.M.dim_d1>maxdim_in)maxdim_in=element.M.dim_d1;
      element.sweep_forward(trunc, pass_on, (!flip_last) && (n==(int)repeated.size()-1));
      if(element.M.dim_d1>maxdim_out){
        maxdim_out=element.M.dim_d1;
        maxdim_at=n+initial.size();
      }
    }catch(DummyException &e){
      std::cerr<<"called by ProcessTensorRepeat::sweep_forward: repeated["<<n<<"]"<<std::endl;
      throw e;
    } 
  }

  if(verbosity>0)std::cout<<"Maxdim at n="<<maxdim_at<<": "<<maxdim_in<<" -> "<<maxdim_out<<std::endl;
}


void ProcessTensorRepeat::sweep_final_start_backward(const TruncationLayout &trunc, int verbosity){
  TruncatedSVD trunc_tmp=trunc.get_base();
  for(int loop=0; loop<trunc.get_final_sweep_n(); loop++){
    if(verbosity>0){std::cout<<"final loop="<<loop<<"/"<<trunc.get_final_sweep_n()<<" backward sweep"<<std::endl;}
    if(verbosity>0){trunc_tmp.print_info();std::cout<<std::endl;}
    sweep_backward(trunc_tmp, false, false, verbosity);

    if(loop==trunc.get_final_sweep_n()-1 && !trunc.final_sweep_half){
      trunc_tmp=trunc.get_final_sweep();
    }
    if(verbosity>0){std::cout<<"final loop="<<loop<<"/"<<trunc.get_final_sweep_n()<<" forward sweep"<<std::endl;}
    if(verbosity>0){trunc_tmp.print_info();std::cout<<std::endl;}
    sweep_forward(trunc_tmp, false, false, verbosity);
  }
  if(trunc.final_sweep_half){
    if(verbosity>0){std::cout<<"final backward sweep"<<std::endl;}
    trunc_tmp=trunc.get_final_sweep();
    if(verbosity>0){trunc_tmp.print_info();std::cout<<std::endl;}
    sweep_backward(trunc_tmp, true, true, verbosity);
//    sweep_backward(trunc_tmp, false, true, verbosity);
//    sweep_backward(trunc_tmp, false, false, verbosity);
  }
}

void ProcessTensorRepeat::dict_expand(const ReadPT_struct &readPT){
  for(int n=0; n<initial.n_tot; n++){
    initial.get(n, ForwardPreload).accessor.dict_expand(readPT);
  }
  for(int n=0; n<repeated.n_tot; n++){
    repeated.get(n, ForwardPreload).accessor.dict_expand(readPT);
  }
}

bool ProcessTensorRepeat::can_read(const std::string &fname_h){
//std::cout<<"TEST: '"<<fname_h<<"'"<<std::endl;
  std::string magic;
  magic=read_first_bytes(fname_h+std::string("_initial"), 4);
//std::cout<<"magic: '"<<magic<<"'"<<std::endl;
  if(magic!="PTRI")return false;
  magic=read_first_bytes(fname_h+std::string("_repeated"), 4);
//std::cout<<"magic: '"<<magic<<"'"<<std::endl;
  if(magic!="PTRR")return false;

  return true;
}

void ProcessTensorRepeat::read(const std::string &fname_h, bool ro){
 try{
  initial.read(fname_h+std::string("_initial"), "PTRI", ro);
 }catch(DummyException &e){
  std::cerr<<"while reading ProcessTensorRepeat::initial"<<std::endl;
  throw e;
 }
 try{
  repeated.read(fname_h+std::string("_repeated"), "PTRR", ro);
 }catch(DummyException &e){
  std::cerr<<"while reading ProcessTensorRepeat::repeated"<<std::endl;
  throw e;
 }
 set_read_only(ro);
}

void ProcessTensorRepeat::set_trivial(int Nsys){
  set_specs("",-1);
  initial.clear();
  initial.push_back(ProcessTensorElement(Nsys));
  repeated.clear();
  repeated.push_back(ProcessTensorElement(Nsys));
  n=0;
}

void ProcessTensorRepeat::print_info(std::ostream &os)const{
  os<<"initial.size()="<<initial.size()<<" repeated.size()="<<repeated.size()<<std::endl;
}


}
