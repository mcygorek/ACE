#include "ProcessTensorForwardList.hpp"
#include "ProcessTensorStream_ro.hpp"
#include "ProcessTensorStream.hpp"
#include "ProcessTensorBuffer.hpp"
#include "ProcessTensorRepeat.hpp"
#include "PT_infinite.hpp"
#include <iostream>
#include "otimes.hpp"
#include "LiouvilleTools.hpp"
#include "MPG_Selector.hpp"
#include "TempFileName.hpp"
#include "Reader.hpp"
#include "InitialState.hpp"
#include "Timings.hpp"

namespace ACE{

void ProcessTensorForwardList::complain_if_null()const{
  for(size_t i=0; i<size(); i++){
    if(!list[i]){
      std::cerr<<"ProcessTensorForwardList: list["<<i<<"]: PT == NULL!"<<std::endl; 
      throw DummyException();
    }
  }
}

void ProcessTensorForwardList::print_info()const{
  std::cout<<"PT list: size: "<<size()<<std::endl;
}

std::vector<const ProcessTensorElement *> ProcessTensorForwardList::current_list(){
  complain_if_null();
  std::vector<const ProcessTensorElement *> ret;
  for(size_t i=0; i<size(); i++){
    ret.push_back(list[i]->current());
  }
  return ret;
}

void ProcessTensorForwardList::reset(){
  for(size_t i=0; i<list.size(); i++){
    list[i]->reset();
  }
}
bool ProcessTensorForwardList::done()const{
  for(size_t i=0; i<list.size(); i++){
    if(list[i]->done())return true;
  }
  return false;
}
void ProcessTensorForwardList::load_next(){
  for(size_t i=0; i<list.size(); i++){
    list[i]->load_next();
  }
}

std::shared_ptr<ProcessTensorForward> ProcessTensorForwardList::PTptr_from_file(const std::string &fname, bool read_only){

  std::shared_ptr<ProcessTensorForward> ret;

  if(ProcessTensorRepeat::can_read(fname)){
    std::cout<<"read as ProcessTensorRepeat"<<std::endl;
    ret=std::shared_ptr<ProcessTensorForward>(new ProcessTensorRepeat(fname));
    if(read_only){
      dynamic_cast<ProcessTensorRepeat*>(ret.get())->set_read_only(read_only);
    }
  }else{
    std::cout<<"read as ProcessTensorBuffer"<<std::endl;
    ret=std::shared_ptr<ProcessTensorForward>(new ProcessTensorBuffer(fname));
    if(read_only){
      dynamic_cast<ProcessTensorBuffer*>(ret.get())->read_only=read_only;
    }
  }
  return ret;
}

void ProcessTensorForwardList::add_PT(const ReadPT_struct &expand){
    std::cout<<"add_PT: '"<<expand.fname<<"'";
    if(expand.expand_front>1){std::cout<<" expand_front="<<expand.expand_front;}
    if(expand.expand_back>1){std::cout<<" expand_back="<<expand.expand_back;}
    std::cout<<std::endl;

    list.push_back(PTptr_from_file(expand.fname,true));
    for(size_t j=temp_expand.size(); j<list.size(); j++){
      temp_expand.push_back(ReadPT_struct());
    }
    temp_expand[list.size()-1]=expand;
}
void ProcessTensorForwardList::add_PT(Parameters &param){
  std::vector<std::vector<std::string> > add_PT_svv = param.get("add_PT");

  for(size_t i=0; i<add_PT_svv.size(); i++){
    bool complain=false;
    if(add_PT_svv[i].size()<1){
      std::cerr<<"USAGE: add_PT FILENAME [expand_dim_front] [expand_dim_back]!"<<std::endl;
      throw DummyException();
    }
    ReadPT_struct expand(add_PT_svv[i][0]);
    if(add_PT_svv[i].size()>1){
      expand.expand_front=readDouble(add_PT_svv[i][1],"add_PT FILENAME [expand_dim_front]");
    }
    if(add_PT_svv[i].size()>2){
      expand.expand_back=readDouble(add_PT_svv[i][2],"add_PT FILENAME [expand_dim_front] [expand_dim_back]");
    }

    add_PT(expand);
/*
    std::cout<<"add_PT: '"<<expand.fname<<"'";
    if(expand.expand_front>1){std::cout<<" expand_front="<<expand.expand_front;}
    if(expand.expand_back>1){std::cout<<" expand_back="<<expand.expand_back;}
    std::cout<<std::endl;

    list.push_back(PTptr_from_file(add_PT[i][0],true));
    for(size_t j=temp_expand.size(); j<list.size(); j++){
      temp_expand.push_back(ReadPT_struct());
    }
    temp_expand[list.size()-1]=expand;
*/
  }
}

void ProcessTensorForwardList::setup2(Parameters &param, std::vector<std::shared_ptr<ModePropagatorGenerator> > & initial_mpgs, int setdim, bool print_timings){

 time_point time1=now();
 try{
  list.clear();
  temp_expand.clear();
  if(param.is_specified("multi_PT")){
    std::cerr<<"Parameter 'multi_PT' is now replaced by 'add_PT'!"<<std::endl;
    throw DummyException();
  }
  if(param.is_specified("read_PT")){
    std::cerr<<"Parameter 'read_PT' is now replaced by 'initial_PT' (rw) or 'add_PT' (ro)!"<<std::endl;
    throw DummyException();
  }

  std::vector<std::shared_ptr<ModePropagatorGenerator> > mpgs=MPG_Selector(param);
  for(size_t i=0; i<initial_mpgs.size(); i++){
    mpgs.push_back(initial_mpgs[i]);
  }


  TimeGrid tgrid(param);
  TruncationLayout trunc(param);
  std::string initial_PT = param.get_as_string("initial_PT");  
  std::string write_PT = param.get_as_string("write_PT");
  double dict_zero = param.get_as_double("dict_zero",-1.);
  int verbosity = 1;
  int buffer_blocksize = param.get_as_int("buffer_blocksize",-1);
  bool use_Gaussian_repeat_JP = param.get_as_bool("use_Gaussian_repeat_JP",
     param.get_as_bool("use_Gaussian_periodic_JP",false));
  bool use_Gaussian_repeat = param.get_as_bool("use_Gaussian_repeat",
     param.get_as_bool("use_Gaussian_periodic",false)) | use_Gaussian_repeat_JP;
  int use_Gaussian_log_hybrid = param.get_as_size_t("use_Gaussian_log_hybrid",0);
  bool use_Gaussian_log = param.get_as_bool("use_Gaussian_log",
           param.get_as_bool("use_Gaussian_DnC",
            param.get_as_bool("use_Gaussian_divide_and_conquer",false))
       ) | use_Gaussian_repeat | (use_Gaussian_log_hybrid>0);

  bool use_Gaussian_infinite = param.get_as_bool("use_Gaussian_infinite");

  bool use_Gaussian_fw = param.get_as_bool("use_Gaussian_fw");
  bool use_Gaussian_select = param.get_as_bool("use_Gaussian_select",false);
  bool use_Gaussian = param.get_as_bool("use_Gaussian",use_Gaussian_fw|use_Gaussian_log|use_Gaussian_select|use_Gaussian_repeat|use_Gaussian_infinite);

  int intermediate_sweep_n = param.get_as_size_t("intermediate_sweep_n", 0);

  int final_sweep_n = param.get_as_size_t("final_sweep_n", 0);

  if(initial_PT!="" && use_Gaussian){
    std::cerr<<"Error: Parameter 'read_PT' incompatible with 'use_Gaussian true'. Use 'add_PT' instead!"<<std::endl;
    throw DummyException();
  }
 
  //setup expansion structure for initial_PT:
  ReadPT_struct initial_PT_struct(initial_PT);
  if(initial_PT!=""){
    std::vector<std::string> row=param.get_row("initial_PT",0);
    if(row.size()>1){
      initial_PT_struct.expand_front=readDouble(row[1],"initial_PT: expand_front");
    }
    if(row.size()>2){
      initial_PT_struct.expand_back=readDouble(row[2],"initial_PT: expand_back");
    }
  }


  //Difficult to modify a ProcessTensorRepeat after it's calculated => Only read, don't modify (except for expanding system dimension)
 try{
  if(initial_PT!="" && ProcessTensorRepeat::can_read(initial_PT)){
    list.push_back(std::shared_ptr<ProcessTensorForward>(new ProcessTensorRepeat(initial_PT)));
//    dynamic_cast<ProcessTensorRepeat*>(list.back().get())->read_only=true;
    //Modify permanently:
    temp_expand.push_back(ReadPT_struct(initial_PT));
    if(initial_PT_struct.have_to_expand()){
      list.back()->dict_expand(initial_PT_struct);
    }
    
    add_PT(param);

    time_point time2=now();
    std::cout<<"runtime for setting up PT: "<<time_diff(time2-time1)<<"ms"<<std::endl;
    return;
  }


  bool need_do_work=(mpgs.size()>0 || use_Gaussian);
  if(!need_do_work){
    if(initial_PT!=""){
      list.push_back(PTptr_from_file(initial_PT,false));

      //for initial_PT, expansion is not temporary but permanent
      temp_expand.push_back(ReadPT_struct(initial_PT));
      if(initial_PT_struct.have_to_expand()){
        list.back()->dict_expand(initial_PT_struct);
      }

    }else{ //neither calculating any PT nor reading initial_PT
      int Nsys = InitialState(param).rho.rows();
      list.push_back(std::shared_ptr<ProcessTensorForward>(new ProcessTensorRepeat(Nsys)));
      temp_expand.push_back(ReadPT_struct());
    }

    add_PT(param);

    time_point time2=now();
    std::cout<<"runtime for setting up PT: "<<time_diff(time2-time1)<<"ms"<<std::endl;
    return; 
  }
 }catch(DummyException &e){
   std::cerr<<"while reading PTs; ";
   throw e;
 }

  DiagBB diagBB;  //setup diagBB:
  if(use_Gaussian){
   try{
    std::string Gaussian_prefix=param.get_as_string("Gaussian_prefix", "Boson");
    diagBB.setup(param, Gaussian_prefix);
   }catch(DummyException &e){
    std::cerr<<"while setting up diagBB; ";
    throw e;
   }
    if((!param.is_specified("t_mem"))&&(!param.is_specified("n_mem")) 
       && param.is_specified("mem_threshold")){
      double mem_threshold=param.get_as_double("mem_threshold");
      int n_mem_est=diagBB.estimate_memory_length(tgrid.n_mem, tgrid.dt, mem_threshold, true);
      tgrid.n_mem=n_mem_est;
    }
  }

  if(use_Gaussian_infinite){
    list.push_back(PT_infinite(param, diagBB));

  }else if(use_Gaussian_repeat){ 
    std::cout<<"use_Gaussian_repeat=true"<<std::endl;
    list.push_back(std::shared_ptr<ProcessTensorForward>(new ProcessTensorRepeat()));
    ProcessTensorRepeat *PTR = dynamic_cast<ProcessTensorRepeat*>(list.back().get()); 
    PTR->set_specs(write_PT, buffer_blocksize);
    PTR->calculate(diagBB, tgrid, trunc, dict_zero, verbosity, !use_Gaussian_repeat_JP);
  }else{  // use blockwise processing of PT elements, i.e. ProcessTensorBuffer:

    //need temporary files when blocksize finite:
    bool did_use_temporary=false;
    if(write_PT=="" && !(buffer_blocksize<0)){ 
      temp_file_list.push_back(TempFileName().set_noremove(true).get());
      write_PT=temp_file_list.back();
      did_use_temporary=true;
    }

    list.push_back(std::shared_ptr<ProcessTensorForward>(new ProcessTensorBuffer()));
    ProcessTensorBuffer *PTB = dynamic_cast<ProcessTensorBuffer*>(list.back().get()); 
    if(write_PT!=""){
      PTB->set_new_file(write_PT, buffer_blocksize);
      PTB->use_async_write=param.get_as_bool("use_async_write",false);
    }

    if(use_Gaussian){
      if(use_Gaussian_log){ 
        bool copy_compress_second=param.get_as_bool("copy_compress_second",false);
        if(use_Gaussian_log_hybrid>0){
          PTB->set_from_DiagBB(diagBB, tgrid, trunc, dict_zero, verbosity, use_Gaussian_log_hybrid);
          PTB->set_from_DiagBB_log(diagBB, tgrid, trunc, dict_zero, verbosity, -1, use_Gaussian_log_hybrid,copy_compress_second);
        }else{
          PTB->set_from_DiagBB_log(diagBB, tgrid, trunc, dict_zero, verbosity);
        }
      }else if(use_Gaussian_select){
        PTB->set_from_DiagBB_select(diagBB, tgrid, trunc, dict_zero, verbosity);
      }else if(use_Gaussian_fw){
        PTB->set_from_DiagBB_fw(diagBB, tgrid, trunc, dict_zero, verbosity);
      }else{
        PTB->set_from_DiagBB(diagBB, tgrid, trunc, dict_zero, verbosity);
      }
    }else{
      //read/set initial PT:
      bool use_combine_tree = param.get_as_bool("use_combine_tree");
      bool use_combine_tree_randomized = param.get_as_bool("use_combine_tree_randomized");
      bool use_combine_select = param.get_as_bool("use_select");
      bool use_combine_alternate = param.get_as_bool("use_combine_alternate");
      int sysdim=mpgs[0]->get_N();


      std::shared_ptr<CompressionTree> TTree;
      std::shared_ptr<CompressionTree> TTree_inv;
      int TTree_at;
      std::string TTree_filename;
      if(param.is_specified("calculate_CompressionTree")){
        TTree_filename = param.get_as_string_check("calculate_CompressionTree",0,0);
        TTree_at = param.get_as_size_t("calculate_CompressionTree",tgrid.n_tot/2,0,1);
        if(TTree_at <1 || TTree_at > tgrid.n_tot-3 ){
          std::cerr<<"TTree_at="<<TTree_at<<" out of bounds!"<<std::endl;
          throw DummyException();
        }
        TTree = std::make_shared<CompressionTree>(1);
        TTree_inv = std::make_shared<CompressionTree>(1);
        std::cout<<"Prepare CompressionTree at site "<<TTree_at<<" and select file "<<TTree_filename<<std::endl;
      }


      if(initial_PT!=""){
        if(TTree){
          std::cerr<<"Combining initial_PT with CompressionTree not implemented yet!"<<std::endl;
          throw DummyException();
        }
        PTB->copy_content(initial_PT_struct);
      }else{
        PTB->set_trivial(tgrid.n_tot, sysdim);
        PTB->TTree_at=TTree_at;
        PTB->TTree=TTree;
        PTB->TTree_inv=TTree_inv;
      }

      for(int i=0; i<(int)mpgs.size(); i++){
        if(use_combine_tree_randomized){
          PTB->set_from_modes_tree_randomized(*mpgs[i].get(), tgrid, trunc, \
                   param.get_as_double("randomized_chi_a", 0), \
                   param.get_as_double("randomized_chi_b", 16), \
                   dict_zero, verbosity);
        }else if(use_combine_tree){
          if(initial_PT=="" && i==0){
std::cout<<"using: set_from_modes_tree"<<std::endl;
            PTB->set_from_modes_tree(*mpgs[i].get(), tgrid, trunc, dict_zero, verbosity);
          }else{
std::cout<<"using: add_modes_tree"<<std::endl;
            PTB->add_modes_tree(*mpgs[i].get(), tgrid, trunc, dict_zero, verbosity);
          }
        }else if(use_combine_select){
          PTB->add_modes_select(*mpgs[i].get(), tgrid, trunc, dict_zero, verbosity);
        }else if(use_combine_alternate){
          PTB->add_modes_firstorder(*mpgs[i].get(), tgrid, trunc, dict_zero, verbosity, true);
        }else{
          PTB->add_modes(*mpgs[i].get(), tgrid, trunc, dict_zero, verbosity);
        }
      }


      if(TTree){  
        TTree->write(TTree_filename);
        TTree->print_info_file(TTree_filename+"_info");
      }
      if(TTree_inv){  
        TTree_inv->write(TTree_filename+"_inv");
        TTree_inv->print_info_file(TTree_filename+"_inv_info");
      }
    }
    
    PTB->write_all(); if(PTB->use_async_write)PTB->wait_write();
    PTB->read_only=true;
    if(did_use_temporary && PTB->use_multiple_files()){ 
      for(int bl=0; bl<PTB->get_nr_blocks(); bl++){
        temp_file_list.push_back(PTB->get_fname(bl));
      }
    } 
  }

  add_PT(param);

 }catch(DummyException &e){
   std::cerr<<"called by ProcessTensorForwardList::setup"<<std::endl;
   throw e;
 }
  time_point time2=now();
  std::cout<<"runtime for setting up PT: "<<time_diff(time2-time1)<<"ms"<<std::endl;
}

void ProcessTensorForwardList::read(const std::string &fname){

  list.push_back(std::shared_ptr<ProcessTensorForward>(new ProcessTensorBuffer(fname)));
//  list.push_back(std::shared_ptr<ProcessTensorForward>(new ProcessTensorStream_ro(fname, true)));
}

/*
void ProcessTensorForwardList::propagate_select(Eigen::MatrixXcd & state, const TruncatedSVD &trunc){
  std::cerr<<"ProcessTensorForwardList::propagate_select: NOT IMPLEMENTED YET!"<<std::endl; 
  throw DummyException();
}
*/

Eigen::VectorXcd ProcessTensorForwardList::get_rho_reduced(const Eigen::MatrixXcd & state){
  if(done()){
    std::cerr<<"ProcessTensorForwardList::get_rho_reduced: 'done' was set!"<<std::endl;
    throw DummyException();
  }
  std::vector<const ProcessTensorElement *> elements = current_list();

  Eigen::VectorXcd closure=Eigen::VectorXcd::Ones(1);
  for(int i=0; i<(int)elements.size(); i++){
    closure = Vector_otimes(closure, elements[i]->closure);   
  }
  return state * closure;
}

std::vector<std::complex<double> > ProcessTensorForwardList::get_env_reduced(
    const Eigen::MatrixXcd & state, const Which_Env_Ops_List & which_env_ops){

  if(done()){
    std::cerr<<"ProcessTensorForwardList::get_env_reduced: 'done' was set!"<<std::endl;
    throw DummyException();
  }
  std::vector<const ProcessTensorElement *> elements = current_list();

  std::vector<std::complex<double> > env_reduced;

  for(const Which_Env_Ops & w : which_env_ops){
    Eigen::VectorXcd sysop = H_Matrix_to_L_Vector(w.A.transpose());
    if(sysop.rows()!=state.rows()){
      std::cerr<<"Env_Ops: sysop.rows()!=state.rows() ("<<sysop.rows()<<" vs. "<<state.rows()<<std::endl;
      throw DummyException();
    }
  
    if(w.i<0||w.i>=(int)elements.size()){
      std::cerr<<"Env_Ops: PT out of bounds ("<<w.i<<" vs. "<<elements.size()<<")!"<<std::endl;
      throw DummyException();
    }
    if(w.o<0||w.o>=(int)elements[w.i]->env_ops.size()){
      std::cerr<<"Env_Ops: operator out of bounds ("<<w.o<<" vs. "<<elements[w.i]->env_ops.size()<<")!"<<std::endl;
      throw DummyException();
    }

    Eigen::VectorXcd closure=Eigen::VectorXcd::Ones(1);
    for(int j=0; j<w.i; j++){
      closure = Vector_otimes(closure, elements[j]->closure);
    }
    closure = Vector_otimes(closure, elements[w.i]->env_ops[w.o]);
    for(int j=w.i+1; j<(int)elements.size(); j++){
      closure = Vector_otimes(closure, elements[j]->closure);
    }

    if(state.cols() != closure.rows()){
      std::cerr<<"ProcessTensorForwardList::get_env_reduced: state.cols() != closure.rows()! ("<<state.cols()<<" vs. "<<closure.rows()<<")!"<<std::endl;
      throw DummyException();
    }
    
    env_reduced.push_back((sysop.transpose() * state * closure)(0));
  }
  return env_reduced;
}


void ProcessTensorForwardList::propagate(Eigen::MatrixXcd & state, bool reverse_order){
  std::vector<const ProcessTensorElement *> elements = current_list();

  int d1tot=1; 
  for(size_t i=0; i<elements.size(); i++)d1tot*=elements[i]->M.dim_d1;
  
  int d2tot=1; 
  for(size_t i=0; i<elements.size(); i++)d2tot*=elements[i]->M.dim_d2;

  if(state.cols()!=d1tot){
    std::cerr<<"ProcessTensorForwardList::propagate: state.cols()!=d1tot ("<<state.cols()<<" vs. "<<d1tot<<")!"<<std::endl;
    throw DummyException();
  }
  
  if(!reverse_order){
    for(size_t i=0; i<elements.size(); i++){
      int dim1_front=1;
      for(int j=0; j<i; j++)dim1_front*=elements[j]->M.dim_d2;
      
      if(i<temp_expand.size()){
        elements[i]->propagate(state, dim1_front, temp_expand[i]);
      }else{
        elements[i]->propagate(state, dim1_front);
      }
    }
  }else{
    for(int i=elements.size()-1; i>=0; i--){
      int dim1_front=1;
      for(int j=0; j<i; j++)dim1_front*=elements[j]->M.dim_d1;

      if(i<temp_expand.size()){
        elements[i]->propagate(state, dim1_front, temp_expand[i]);
      }else{
        elements[i]->propagate(state, dim1_front);
      }     
    }
  }
}

}//namespace
