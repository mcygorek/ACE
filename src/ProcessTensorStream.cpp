#include "ProcessTensorStream_ro.hpp"
#include "ProcessTensorStream_wo.hpp"
#include "ProcessTensorStream.hpp"
#include "ProcessTensor.hpp"
#include "LiouvilleTools.hpp"
#include "ModePropagatorGenerator.hpp"
#include "DiagBB.hpp"
#include "TimeGrid.hpp"
#include "TempFileName.hpp"
#include "DummyException.hpp"


namespace ACE{

int ProcessTensorStream::get_length(const std::string & file){ 
  std::ifstream ifs(file.c_str());
  return ProcessTensor::read_header(ifs, std::string("get_length('")+file+std::string("')")).first;
}

void ProcessTensorStream::check_consistency(const std::string & file){
  ProcessTensorStream_ro PT_ro(file, true);
 
  ProcessTensorElement e=PT_ro.get_first();
  e.check_consistency(file+" n=0");
  if(e.M.dim_d1!=1){
    std::cerr<<"ProcessTensorStream::check_consistency: file '"<<file<<"': ";
    std::cerr<<"first dim_d1!=1 ("<<e.M.dim_d1<<")!"<<std::endl;
    throw DummyException();
  }
  int max_abs=e.M.max_element_abs();
  int dim_d2=e.M.dim_d2;
  for(int n=1; n<(int)PT_ro.size(); n++){
    e=PT_ro.get(n);
    e.check_consistency(file+" n="+std::to_string(n));
    if(e.M.max_element_abs()>max_abs || !std::isfinite(e.M.max_element_abs())){
      max_abs=e.M.max_element_abs();
    }
    if(e.M.dim_d1!=dim_d2){
      std::cerr<<"ProcessTensorStream::check_consistency: file '"<<file<<"': ";
      std::cerr<<"n="<<n<<": previous.dim_d2!=dim_d1 ("<<dim_d2<<","<<e.M.dim_d1<<")!"<<std::endl;
      throw DummyException();
    }   
    dim_d2=e.M.dim_d2;
  }
  if(e.M.dim_d2!=1){
    std::cerr<<"ProcessTensorStream::check_consistency: file '"<<file<<"': ";
    std::cerr<<"last dim_d2!=1 ("<<e.M.dim_d2<<")!"<<std::endl;
    throw DummyException();
  }
  std::cout<<"Check consistency: max_abs="<<max_abs<<std::endl;
}

void ProcessTensorStream::set_trivial(const std::string & file_out, 
                                             int n_max, int sysdim){

  if(n_max<1){
    std::cerr<<"ProcessTensorStream::set_trivial: n_max<1!"<<std::endl;
    throw DummyException();
  }
  if(sysdim<2){
    std::cerr<<"ProcessTensorStream::set_trivial: sysdim<2!"<<std::endl;
    throw DummyException();
  }
 
  ProcessTensorStream_wo PT_wo(file_out, n_max, false);
  ProcessTensorElement element;
  element.set_trivial(sysdim);
  for(int n=0; n<n_max; n++){
    PT_wo.put(element);
  }
}

void ProcessTensorStream::copy(const std::string & file_out,
                               const std::string & file_in){

  ProcessTensorStream_ro PT_ro(file_in, false);
  ProcessTensorStream_wo PT_wo(file_out, PT_ro.size(), false);

  for(int n=0; n<(int)PT_ro.size(); n++){
    PT_wo.put(PT_ro.get(n));
  }
}

void ProcessTensorStream::calculate_closures(const std::string & file_out,
                                             const std::string & file_in){

  ProcessTensorStream_ro PT_ro(file_in, true);
  ProcessTensorStream_wo PT_wo(file_out, PT_ro.size(), true);

  ProcessTensorElement last=PT_ro.get_last();
  last.calculate_closure(NULL);
  PT_wo.put(last);
  for(int n=(int)PT_ro.size()-2; n>=0; n--){
    ProcessTensorElement e=PT_ro.get(n);
    e.calculate_closure(&last);
    PT_wo.put(e);
    last=e;
  }
}

void ProcessTensorStream::sweep_forward(const std::string & file_out, 
                                        const std::string & file_in,
                                        const TruncatedSVD &trunc, 
                                        int verbosity,
                                        int range_start, int range_end){

  ProcessTensorStream_ro PT_ro(file_in, true);  
  ProcessTensorStream_wo PT_wo(file_out, PT_ro.size(), false);  
  if(range_start<0)range_start=0;
  if(range_end<0)range_end=PT_ro.size();
  if(range_start==range_end){
     copy(file_out, file_in);
     return; 
  }
  if(range_start>range_end){
    std::cerr<<"ProcessTensorStream::sweep_forward: range_start>range_end!"<<std::endl; 
    throw DummyException();
  }
  if(verbosity>0){
    std::cout<<"sweep_forward: '"<<file_in<<"' -> '"<<file_out<<"'"<<std::endl;
    std::cout<<"range: ["<<range_start<<".."<<range_end<<"["<<std::endl;
  }
  
  for(int n=0; n<range_start; n++){
    PT_wo.put(PT_ro.get(n));
  }

  int maxdim_in=0, maxdim_out=0;
  PassOn pass_on;
  for(int n=range_start; n<range_end; n++){
    ProcessTensorElement element = PT_ro.get(n);
    if(n==range_start){
      pass_on=PassOn(element.M.dim_d1);
      maxdim_in=maxdim_out=element.M.dim_d1;
    }
    if(element.M.dim_d2>maxdim_in)maxdim_in=element.M.dim_d2;

    element.sweep_forward(trunc, pass_on, (n==range_end-1));

    if(element.M.dim_d2>maxdim_out)maxdim_out=element.M.dim_d2;
    PT_wo.put(element);
  }
  for(int n=range_end; n<(int)PT_ro.size(); n++){
    PT_wo.put(PT_ro.get(n));
  }
  if(verbosity>0)std::cout<<"Maxdim: "<<maxdim_in<<" -> "<<maxdim_out<<std::endl;
}
void ProcessTensorStream::sweep_backward(const std::string & file_out, 
                                        const std::string & file_in,
                                        const TruncatedSVD &trunc,
                                        int verbosity,
                                        int range_start, int range_end){

  ProcessTensorStream_ro PT_ro(file_in, true);  
  ProcessTensorStream_wo PT_wo(file_out, PT_ro.size(), true);  
  if(range_start<0)range_start=0;
  if(range_end<0)range_end=PT_ro.size();
  if(range_start==range_end){
     copy(file_out, file_in);
     return; 
  }
  if(range_start>range_end){
    std::cerr<<"ProcessTensorStream::sweep_backward: range_start>range_end!"<<std::endl; 
    throw DummyException();
  }
  if(verbosity>0){
    std::cout<<"sweep_backward: '"<<file_in<<"' -> '"<<file_out<<"'"<<std::endl;
    std::cout<<"range: ["<<range_start<<".."<<range_end<<"["<<std::endl;
  }
  
  for(int n=(int)PT_ro.size()-1; n>=range_end; n--){
    PT_wo.put(PT_ro.get(n));
  }  

  PassOn pass_on;
  int maxdim_in=0, maxdim_out=0;
  for(int n=range_end-1; n>=range_start; n--){
    ProcessTensorElement element = PT_ro.get(n);
    if(n==range_end-1){
      pass_on=PassOn(element.M.dim_d2);
      maxdim_in=maxdim_out=element.M.dim_d2;
    }
    if(element.M.dim_d1>maxdim_in)maxdim_in=element.M.dim_d1;
    element.sweep_backward(trunc, pass_on, (n==range_start));
    if(element.M.dim_d1>maxdim_out)maxdim_out=element.M.dim_d1;
    PT_wo.put(element);
  }

  for(int n=range_start-1; n>=0; n--){
    PT_wo.put(PT_ro.get(n));
  }
  if(verbosity>0)std::cout<<"Maxdim: "<<maxdim_in<<" -> "<<maxdim_out<<std::endl;
}

void ProcessTensorStream::join_select_and_sweep_backward(
            const std::string & file_out, 
            const std::string & file_in1, const std::string & file_in2,
            const TruncatedSVD &trunc, int verbosity,
            ShiftExtend shift_extend){

  bool subsequent_sweep_full = false;
  bool subsequent_sweep = subsequent_sweep_full | true;
  PassOn pass_on;

  ProcessTensorStream_ro PT_ro1(file_in1, true);  
  ProcessTensorStream_ro PT_ro2(file_in2, true);  

  //determine total length of result
  int n_tot=PT_ro2.size()+shift_extend.shift_second;
  if(shift_extend.truncate_at>=0){
    if(shift_extend.truncate_at>n_tot){
       shift_extend.truncate_at=-1;
    }else{
      n_tot=shift_extend.truncate_at;
    }
  } 
  if(n_tot<=0){
    std::cerr<<"join_select_and_sweep_backward called with n_tot="<<n_tot<<"<=0!"<<std::endl;
    throw DummyException();
  }
  if(PT_ro1.size()<=shift_extend.shift_second){ //no overlap in this case
    std::cerr<<"ProcessTensorStream::join_select_and_sweep_backward: ";
    std::cerr<<"PT_ro1.size()<=shift_extend.shift_second!"<<std::endl;
    throw DummyException();
  }

  int extend_first=n_tot-PT_ro1.size(); 
  int n_trunc1=0;
  if(extend_first<0){
    n_trunc1=-extend_first;
    extend_first=0;
  }


  ProcessTensorStream_wo PT_wo(file_out, n_tot, true); 

  if(verbosity>0)std::cout<<"join_select_and_sweep_backward: '"<<file_in1<<"','"<<file_in2<<"' -> '"<<file_out<<"'"<<std::endl;
  int maxdim_in1=0, maxdim_in2=0, maxdim_out=0;

if(verbosity>0){
std::cout<<"n_tot="<<n_tot;
std::cout<<" PT_ro1.size()="<<PT_ro1.size();
std::cout<<" PT_ro2.size()="<<PT_ro2.size();
std::cout<<" shift_second="<<shift_extend.shift_second;
std::cout<<" truncate_at="<<shift_extend.truncate_at;
std::cout<<" extend_first="<<extend_first;
std::cout<<" n_trunc1="<<n_trunc1<<std::endl;
}

  //write overhanging elements of shiftet PT_ro2 (length: extend_first) )

  if(extend_first>0){
    ProcessTensorElement element2;
    for(int n=n_tot-1; n>=(int)PT_ro1.size(); n--){
      if(n==n_tot-1){
        element2 = PT_ro2.get_as_last(n-shift_extend.shift_second);
        pass_on=PassOn(element2.M.dim_d2);
      }else{
        element2 = PT_ro2.get(n-shift_extend.shift_second);
      }

      if(subsequent_sweep_full){
        element2.sweep_backward(trunc, pass_on, (n==0));
      }
      PT_wo.put(element2);
    }
  }

  ProcessTensorElement element1, element2;
  if(extend_first>0){  //previous loop was triggered: PT2 open
    element1=PT_ro1.get_last();
    element2=PT_ro2.get(PT_ro1.size()-1-shift_extend.shift_second);
  }else{ //may need to close off both PT2
    element1=PT_ro1.get_as_last(n_tot-1);
    element2=PT_ro2.get_as_last(n_tot-1-shift_extend.shift_second);
  }
  ProcessTensorElement next_element1, next_element2;  
  
  SelectIndices k_list_left;
  SelectIndices k_list_right = 
             element1.get_forwardNF_selected_indices(element2, trunc);

  if(subsequent_sweep && (!(subsequent_sweep_full && extend_first))){
    pass_on=PassOn(k_list_right.size());
  }

  maxdim_in1=element1.M.dim_d2;
  maxdim_in2=element2.M.dim_d2;
  maxdim_out=k_list_right.size();


  for(int n=(int)PT_ro1.size()-1-n_trunc1; n>=shift_extend.shift_second; n--){
    if(n==shift_extend.shift_second){
      k_list_left.set_full(element1.M.dim_d1, element2.M.dim_d1);
    }else{
      next_element1=PT_ro1.get(n-1);
      next_element2=PT_ro2.get(n-1-shift_extend.shift_second);
      k_list_left=next_element1.get_forwardNF_selected_indices(next_element2, trunc);
    }
    if(element1.M.dim_d1>maxdim_in1)maxdim_in1=element1.M.dim_d1;
    if(element2.M.dim_d1>maxdim_in2)maxdim_in2=element2.M.dim_d1;

//std::cout<<"before join_selected n="<<n<<std::endl;
    element1.join_selected(n, element2, k_list_left, k_list_right);
//std::cout<<"after join_selected n="<<n<<std::endl;

    if(subsequent_sweep_full){
      element1.sweep_backward(trunc, pass_on, (n==0));
    }else if(subsequent_sweep){
      element1.sweep_backward(trunc, pass_on, (n==shift_extend.shift_second));
    }
//std::cout<<"after sweep backward n="<<n<<std::endl;

    if(element1.M.dim_d1>maxdim_out)maxdim_out=element1.M.dim_d1;

    PT_wo.put(element1);
    k_list_right=k_list_left;
    element1=next_element1;
    element2=next_element2;
  }
  for(int n=shift_extend.shift_second-1; n>=0; n--){
    element1=PT_ro1.get(n);
    if(subsequent_sweep_full){
      element1.sweep_backward(trunc, pass_on, (n==0));
    }
    PT_wo.put(element1);
  }
  
  if(verbosity>0)std::cout<<"Maxdim: "<<maxdim_in1<<","<<maxdim_in2<<" -> "<<maxdim_out<<std::endl;
}

void ProcessTensorStream::join_select_and_sweep_forward(
            const std::string & file_out, 
            const std::string & file_in1, const std::string & file_in2,
            const TruncatedSVD &trunc, int verbosity,
            ShiftExtend shift_extend){

  ProcessTensorStream_ro PT_ro1(file_in1, true);  
  ProcessTensorStream_ro PT_ro2(file_in2, true);  

  //determine total length of result
  int n_tot=PT_ro2.size()+shift_extend.shift_second;
  if(shift_extend.truncate_at>=0){
    if(shift_extend.truncate_at>n_tot){
       shift_extend.truncate_at=-1;
    }else{
      n_tot=shift_extend.truncate_at;
    }
  } 
  if(n_tot<=0){
    std::cerr<<"join_select_and_sweep_forward called with n_tot="<<n_tot<<"<=0!"<<std::endl;
    throw DummyException();
  }
  if(PT_ro1.size()<=shift_extend.shift_second){ //no overlap in this case
    std::cerr<<"ProcessTensorStream::join_select_and_sweep_forward: ";
    std::cerr<<"PT_ro1.size()<=shift_extend.shift_second!"<<std::endl;
    throw DummyException();
  }

  int extend_first=n_tot-PT_ro1.size(); 
  int n_trunc1=0;
  if(extend_first<0){
    n_trunc1=-extend_first;
    extend_first=0;
  }

 {
  ProcessTensorStream_wo PT_wo(file_out, n_tot, false); 

  if(verbosity>0)std::cout<<"join_select_and_sweep_forward: '"<<file_in1<<"','"<<file_in2<<"' -> '"<<file_out<<"'"<<std::endl;
  int maxdim_in1=0, maxdim_in2=0, maxdim_out=0;

if(verbosity>0){
std::cout<<"n_tot="<<n_tot;
std::cout<<" PT_ro1.size()="<<PT_ro1.size();
std::cout<<" PT_ro2.size()="<<PT_ro2.size();
std::cout<<" shift_second="<<shift_extend.shift_second;
std::cout<<" truncate_at="<<shift_extend.truncate_at;
std::cout<<" extend_first="<<extend_first;
std::cout<<" n_trunc1="<<n_trunc1<<std::endl;
}
 
  for(int n=0; n<shift_extend.shift_second; n++){
    PT_wo.put(PT_ro1.get(n));
  }

  ProcessTensorElement element1, element2;
  element1=PT_ro1.get(shift_extend.shift_second);
  element2=PT_ro2.get_first();
  if(shift_extend.shift_second == n_tot-1){
    element1.close_off();
    element2.close_off();
  }

  ProcessTensorElement next_element1, next_element2;  
  
  SelectIndices k_list_left = 
             element1.get_backwardNF_selected_indices(element2, trunc);
  SelectIndices k_list_right;  
  PassOn pass_on(k_list_left.size());

  maxdim_in1=element1.M.dim_d1;
  maxdim_in2=element2.M.dim_d1;
  maxdim_out=k_list_left.size();

  for(int n=shift_extend.shift_second; n<PT_ro1.size()-n_trunc1; n++){
    if(n==PT_ro1.size()-1-n_trunc1){
      k_list_right.set_full(element1.M.dim_d2, element2.M.dim_d2);
    }else{
      next_element1=PT_ro1.get(n+1);
      next_element2=PT_ro2.get(n+1-shift_extend.shift_second);
      if(n+1==n_tot-1){
        next_element1.close_off();
        next_element2.close_off();
      }
      k_list_right=next_element1.get_backwardNF_selected_indices(next_element2, trunc);
    }
    if(element1.M.dim_d2>maxdim_in1)maxdim_in1=element1.M.dim_d2;
    if(element2.M.dim_d2>maxdim_in2)maxdim_in2=element2.M.dim_d2;

    element1.join_selected(n, element2, k_list_left, k_list_right);
    element1.sweep_forward(trunc, pass_on, (n==(int)PT_ro1.size()-1));

    if(element1.M.dim_d2>maxdim_out)maxdim_out=element1.M.dim_d2;

    PT_wo.put(element1);
    k_list_left=k_list_right;
    element1=next_element1;
    element2=next_element2;
  }

  //overhang of PT_ro2:
  if(extend_first>0){
    for(int n=PT_ro1.size(); n<n_tot; n++){
      if(n==n_tot-1){
        PT_wo.put(PT_ro2.get_as_last(n-shift_extend.shift_second));
      }else{
        PT_wo.put(PT_ro2.get(n-shift_extend.shift_second));
      }
    }
  }
  
  if(verbosity>0)std::cout<<"Maxdim: "<<maxdim_in1<<","<<maxdim_in2<<" -> "<<maxdim_out<<std::endl;

 }
 // std::cout<<ProcessTensorStream::dims(file_out)<<std::endl;
}

void ProcessTensorStream::set_from_DiagBB_single_line(
                    const std::string & file_out,
                    DiagBB &diagBB, double dt, int n){
  int N=diagBB.sys_dim();
  int NL=N*N;
//std::cout<<"test single line: NL="<<NL<<std::endl;

  if(n<1){
    std::cerr<<"ProcessTensorStream::set_from_DiagBB_single_line: n<1!"<<std::endl; 
    throw DummyException();
  }

//  TempFileName tmpname(file_out);
  std::string tmpname=file_out;
 {
  ProcessTensorStream_wo PT_wo(tmpname, n, false);
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
    PT_wo.put(e);

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
    PT_wo.put(e);
//std::cout<<"i="<<0<<":"<<std::endl; e.print_debug();

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
      PT_wo.put(e);
//std::cout<<"i="<<k<<":"<<std::endl; e.print_debug();
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
    PT_wo.put(e);
//std::cout<<"i="<<n-1<<":"<<std::endl; e.print_debug();
  }

/*
  e.clear();
  e.accessor.dict.set_default_diag(N);
  e.closure=Eigen::VectorXcd::Ones(1);
  e.env_ops.ops=std::vector<Eigen::VectorXcd>(1, e.closure);
  e.M.resize(NL,1,1);
  for(int i=0; i<NL; i++){
    e.M(i,0,0)=1;
  }
  for(int i=n; i<n; i++){
    PT_wo.put(e);
  }
*/

 }//close PT_wo for file tmpname
 
//  calculate_closures(file_out, tmpname);
}

void ProcessTensorStream::set_from_DiagBB(const std::string & file_out,
                    DiagBB &diagBB, const TimeGrid &tgrid, 
                    const TruncatedSVD & trunc, int intermediate_sweep_n,
                    double dict_zero, int verbosity){

  double dt=tgrid.dt;
  int n_tot=tgrid.n_tot;
  int n_mem=tgrid.n_mem;

  TempFileName tmpline(file_out);
  if(verbosity>0){
    std::cout<<"n_tot="<<n_tot<<" n_mem="<<n_mem<<std::endl;
    std::cout<<"Calculating single line..."<<std::endl;
  }
  set_from_DiagBB_single_line(file_out, diagBB, dt, n_mem);
  sweep_forward(tmpline, file_out, trunc, verbosity);
  copy(file_out, tmpline);

  if(verbosity>2){
    std::cout<<"line: dims: "<<ProcessTensorStream::dims(file_out)<<std::endl;
  }


  TempFileName cur(file_out);

  for(int line=1; line<n_tot; line++){
    ShiftExtend shift_extend;
    shift_extend.shift_second=line;
    shift_extend.truncate_at = n_tot;

    if(verbosity>0){
      std::cout<<"line: "<<line<<std::endl;   
    }
    sweep_forward(cur, file_out, trunc, verbosity);

    join_select_and_sweep_backward(file_out, cur, tmpline, trunc, verbosity, shift_extend);

    if(verbosity>2){
      std::cout<<"dims: "<<ProcessTensorStream::dims(file_out)<<std::endl;
    }
    for(int loop=0; loop<intermediate_sweep_n; loop++){
      std::cout<<"intermediate loop="<<loop<<std::endl;
      sweep_forward(cur, file_out, trunc, verbosity);
      sweep_backward(file_out, cur, trunc, verbosity);
      if(verbosity>2){
        std::cout<<"dims: "<<ProcessTensorStream::dims(file_out)<<std::endl;
      }
    }
  }
}

void ProcessTensorStream::set_from_DiagBB_reverse(const std::string & file_out,
                    DiagBB &diagBB, const TimeGrid &tgrid, 
                    const TruncatedSVD & trunc, int intermediate_sweep_n,
                    double dict_zero, int verbosity){

  double dt=tgrid.dt;
  int n_tot=tgrid.n_tot;
  int n_mem=tgrid.n_mem;

  TempFileName tmpline(file_out);
  if(verbosity>0){
    std::cout<<"n_tot="<<n_tot<<" n_mem="<<n_mem<<std::endl;
    std::cout<<"Calculating single line..."<<std::endl;
  }
  set_from_DiagBB_single_line(file_out, diagBB, dt, n_mem);
  sweep_backward(tmpline, file_out, trunc, verbosity);
  copy(file_out, tmpline);

  if(verbosity>2){
    std::cout<<"line: dims: "<<ProcessTensorStream::dims(file_out)<<std::endl;
  }


  TempFileName cur(file_out);

  for(int line=1; line<n_tot; line++){
    ShiftExtend shift_extend;
    shift_extend.shift_second=line;
    shift_extend.truncate_at = n_tot;

    if(verbosity>0){
      std::cout<<"line: "<<line<<std::endl;   
    }
    sweep_backward(cur, file_out, trunc, verbosity);

    join_select_and_sweep_forward(file_out, cur, tmpline, trunc, verbosity, shift_extend);

    if(verbosity>2){
      std::cout<<"dims: "<<ProcessTensorStream::dims(file_out)<<std::endl;
    }
    for(int loop=0; loop<intermediate_sweep_n; loop++){
      std::cout<<"intermediate loop="<<loop<<std::endl;
      sweep_backward(cur, file_out, trunc, verbosity);
      sweep_forward(file_out, cur, trunc, verbosity);
      if(verbosity>2){
        std::cout<<"dims: "<<ProcessTensorStream::dims(file_out)<<std::endl;
      }
    }
  }
}

void ProcessTensorStream::set_from_DiagBB_log(
                    const std::string & file_out,
                    DiagBB &diagBB, const TimeGrid &tgrid, 
                    TruncatedSVD trunc, int intermediate_sweep_n,
                    double dict_zero, int verbosity){

  trunc.keep = diagBB.get_dim();

  bool do_check=false;

  double dt = tgrid.dt;
  int n_tot = tgrid.n_tot;
  int n_mem = tgrid.n_mem;

  if(verbosity>2){
    std::cout<<"n_tot="<<n_tot<<" n_mem="<<n_mem<<std::endl;
    std::cout<<"Calculating single line..."<<std::endl;
  }
  set_from_DiagBB_single_line(file_out, diagBB, dt, n_mem);

  if(verbosity>2){
    std::cout<<"line: dims: "<<ProcessTensorStream::dims(file_out)<<std::endl;
  }
  if(do_check){
    std::cout<<"consistency check after single line"<<std::endl;
    ProcessTensorStream::check_consistency(file_out);
  }


  TempFileName cur(file_out);
  for(int line=1; line<n_tot; line*=2){
    ShiftExtend shift_extend;
    shift_extend.shift_second=line;
    shift_extend.truncate_at = n_tot;

    if(verbosity>0){
      std::cout<<"line: "<<line<<std::endl;   
    }

    sweep_forward(cur, file_out, trunc, verbosity);

    if(do_check){
      std::cout<<"line="<<line<<" consistency check after sweep_forward"<<std::endl;
      ProcessTensorStream::check_consistency(cur);
    }

    join_select_and_sweep_backward(file_out, cur, cur, trunc, verbosity, shift_extend);
    if(do_check){
      std::cout<<"line="<<line<<" consistency check after join_select_and_sweep_backward"<<std::endl;
      ProcessTensorStream::check_consistency(file_out);
    }

    if(verbosity>2){
      std::cout<<"dims: "<<ProcessTensorStream::dims(file_out)<<std::endl;
    }
    for(int loop=0; loop<intermediate_sweep_n; loop++){
      std::cout<<"intermediate loop="<<loop<<std::endl;
      int range_start=0, range_end=-1;
      if(loop==0){
        if(line>n_mem)range_start=line-n_mem; 
      }
      sweep_forward(cur, file_out, trunc, verbosity, range_start, range_end);
      if(do_check){
        std::cout<<"line="<<line<<" loop="<<loop<<" consistency check after sweep_forward"<<std::endl;
        ProcessTensorStream::check_consistency(cur);
      }

      sweep_backward(file_out, cur, trunc, verbosity, range_start, range_end);
      if(do_check){
        std::cout<<"line="<<line<<" loop="<<loop<<" consistency check after sweep_backward"<<std::endl;
        ProcessTensorStream::check_consistency(file_out);
      }
      if(verbosity>2){
        std::cout<<"dims: "<<ProcessTensorStream::dims(file_out)<<std::endl;
      }
    }
  }
//  copy(cur, file_out);
//  calculate_closures(file_out, cur);
}

void ProcessTensorStream::set_from_DiagBB_log_reverse(
                    const std::string & file_out,
                    DiagBB &diagBB, const TimeGrid &tgrid, 
                    const TruncatedSVD & trunc, int intermediate_sweep_n,
                    double dict_zero, int verbosity){

  double dt=tgrid.dt;
  int n_tot=tgrid.n_tot;
  int n_mem=tgrid.n_mem;

  if(verbosity>2){
    std::cout<<"n_tot="<<n_tot<<" n_mem="<<n_mem<<std::endl;
    std::cout<<"Calculating single line..."<<std::endl;
  }
  set_from_DiagBB_single_line(file_out, diagBB, dt, n_mem);

  if(verbosity>2){
    std::cout<<"line: dims: "<<ProcessTensorStream::dims(file_out)<<std::endl;
  }


  TempFileName cur(file_out);
  for(int line=1; line<n_tot; line*=2){
    ShiftExtend shift_extend;
    shift_extend.shift_second=line;
    shift_extend.truncate_at = n_tot;

    if(verbosity>0){
      std::cout<<"line: "<<line<<std::endl;   
    }
    sweep_backward(cur, file_out, trunc, verbosity);

    join_select_and_sweep_forward(file_out, cur, cur, trunc, verbosity, shift_extend);

    if(verbosity>2){
      std::cout<<"dims: "<<ProcessTensorStream::dims(file_out)<<std::endl;
    }
    for(int loop=0; loop<intermediate_sweep_n; loop++){
      std::cout<<"intermediate loop="<<loop<<std::endl;
      sweep_backward(cur, file_out, trunc, verbosity);
      sweep_forward(file_out, cur, trunc, verbosity);
      if(verbosity>2){
        std::cout<<"dims: "<<ProcessTensorStream::dims(file_out)<<std::endl;
      }
    }
  }
}


void ProcessTensorStream::set_from_ModePropagator(const std::string & file_out,
               ModePropagator &mprop, const TimeGrid &tgrid, double dict_zero){

  int n_max=tgrid.n_tot;
  ProcessTensorStream_wo PT_wo(file_out, n_max, false);

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
    PT_wo.put(element);
  } 
}

void ProcessTensorStream::add_modes(
          const std::string & file_out, 
          ModePropagatorGenerator &mpg, const TimeGrid &tgrid,
          const TruncatedSVD & trunc, int intermediate_sweep_n,
          double dict_zero, int verbosity){
  
  int n_max=tgrid.n_tot;
  if(n_max<1){
    std::cerr<<"ProcessTensorStream::set_from_ModePropagatorGenerator: n_max<1!"<<std::endl;
    throw DummyException();
  }
  int N=mpg.get_N();
  if(N<2){
    std::cerr<<"ProcessTensorStream::set_from_ModePropagatorGenerator: N<2!"<<std::endl;
    throw DummyException();
  }
  
  if(verbosity>0){
    std::cout<<"Calculating PT for Generator '"<<mpg.name()<<"'"<<std::endl;
  }

//  set_trivial(file_out, n_max, N);
//  copy(file_out, file_in);
  TempFileName tmpname(file_out);
  TempFileName tmpname2((std::string)tmpname);

  for(int k=mpg.first(); k<mpg.get_N_modes(); k=mpg.next(k)){
    if(verbosity>0){
      std::cout<<"Mode "<<k<<"/"<<mpg.get_N_modes()<<std::endl;
    }
    ModePropagatorPtr mpp=mpg.get_ModePropagator(k);

    set_from_ModePropagator(tmpname, *mpp.get(), tgrid, dict_zero);

    sweep_forward(tmpname2, tmpname, trunc, verbosity);
    sweep_forward(tmpname, file_out, trunc, verbosity);

    join_select_and_sweep_backward(file_out, tmpname, tmpname2, trunc, verbosity);
   
    for(int loop=0; loop<intermediate_sweep_n; loop++){
      sweep_forward(tmpname, file_out, trunc, verbosity);
      sweep_backward(file_out, tmpname, trunc, verbosity);
    }
  }
}

void ProcessTensorStream::add_modes_reverse(
          const std::string & file_out, 
          ModePropagatorGenerator &mpg, const TimeGrid &tgrid,
          const TruncatedSVD & trunc, int intermediate_sweep_n,
          double dict_zero, int verbosity){
  
  int n_max=tgrid.n_tot;
  if(n_max<1){
    std::cerr<<"ProcessTensorStream::set_from_ModePropagatorGenerator: n_max<1!"<<std::endl;
    throw DummyException();
  }
  int N=mpg.get_N();
  if(N<2){
    std::cerr<<"ProcessTensorStream::set_from_ModePropagatorGenerator: N<2!"<<std::endl;
    throw DummyException();
  }
  
  if(verbosity>0){
    std::cout<<"Calculating PT for Generator '"<<mpg.name()<<"'"<<std::endl;
  }

//  set_trivial(file_out, n_max, N);
//  copy(file_out, file_in);
  TempFileName tmpname(file_out);
  TempFileName tmpname2((std::string)tmpname);

  for(int k=mpg.first(); k<mpg.get_N_modes(); k=mpg.next(k)){
    if(verbosity>0){
      std::cout<<"Mode "<<k<<"/"<<mpg.get_N_modes()<<std::endl;
    }
    ModePropagatorPtr mpp=mpg.get_ModePropagator(k);

    set_from_ModePropagator(tmpname, *mpp.get(), tgrid, dict_zero);

    sweep_backward(tmpname2, tmpname, trunc, verbosity);
    sweep_backward(tmpname, file_out, trunc, verbosity);

    join_select_and_sweep_forward(file_out, tmpname, tmpname2, trunc, verbosity);
   
    for(int loop=0; loop<intermediate_sweep_n; loop++){
      sweep_backward(tmpname, file_out, trunc, verbosity);
      sweep_forward(file_out, tmpname, trunc, verbosity);
    }
  }
}


void ProcessTensorStream::add_modes_tree_get(
          int level, int first_elem, 
          const std::string & file_out, 
          ModePropagatorGenerator &mpg, const TimeGrid &tgrid,
          const TruncatedSVD & trunc, int intermediate_sweep_n,
          double dict_zero, int verbosity){
                                             
  if(level==0){
    if(first_elem>=mpg.get_N_modes()){
      std::cerr<<"add_modes_tree_get: level="<<level<<" first_elem="<<first_elem<<": first_element>=mpg.get_N_modes()="<<mpg.get_N_modes()<<"!"<<std::endl;
      throw DummyException();
    }
    
    if(verbosity>0){
      std::cout<<"--------------------------------------------"<<std::endl;
      std::cout<<"level: "<<level<<" first_elem: "<<first_elem<<std::endl;
      std::cout<<"--------------------------------------------"<<std::endl;
    }
    if(mpg.skip_list[first_elem]){
      set_trivial(file_out, tgrid.n_tot, mpg.get_N());
    }else{
      ModePropagatorPtr mpp=mpg.get_ModePropagator(first_elem);
      set_from_ModePropagator(file_out, *mpp.get(), tgrid, dict_zero);
    }
  }else{  
    TempFileName tmpname0;
    add_modes_tree_get(level-1, 2*first_elem, tmpname0, mpg, tgrid, trunc, intermediate_sweep_n, dict_zero, verbosity);

    TempFileName tmpname1;
    add_modes_tree_get(level-1, 2*first_elem+1, tmpname1, mpg, tgrid, trunc, intermediate_sweep_n, dict_zero, verbosity);

    if(verbosity>0){
      std::cout<<"--------------------------------------------"<<std::endl;
      std::cout<<"level: "<<level<<" first_elem: "<<first_elem<<std::endl;
      std::cout<<"--------------------------------------------"<<std::endl;
    }
    if(verbosity>0)std::cout<<"sweep first forward"<<std::endl;
    copy(file_out, tmpname0);
    sweep_forward(tmpname0, file_out, trunc, verbosity);
    if(verbosity>0)std::cout<<"sweep second forward"<<std::endl;
    copy(file_out, tmpname1);
    sweep_forward(tmpname1, file_out, trunc, verbosity);
    if(verbosity>0)std::cout<<"join and sweep backward"<<std::endl;
    join_select_and_sweep_backward(file_out, tmpname0, tmpname1, trunc, verbosity);
  
    for(int loop=0; loop<intermediate_sweep_n; loop++){
      if(verbosity>0)std::cout<<"sweep intermediate "<<loop<<"/"<<intermediate_sweep_n<<" forward"<<std::endl;
      sweep_forward(tmpname0, file_out, trunc, verbosity);
      if(verbosity>0)std::cout<<"sweep intermediate "<<loop<<"/"<<intermediate_sweep_n<<" backward"<<std::endl;
      sweep_backward(file_out, tmpname0, trunc, verbosity);
    }
  }
}

void ProcessTensorStream::add_modes_tree(
          const std::string & file_out, 
          ModePropagatorGenerator &mpg, const TimeGrid &tgrid,
          const TruncatedSVD & trunc, int intermediate_sweep_n,
          double dict_zero, int verbosity){

  int N_modes=mpg.get_N_modes();
  if(N_modes<1)return;

  int N_hierarchy=1;
  {int N_shift=N_modes;
    while(N_shift>1){
      N_hierarchy++;
      N_shift=N_shift>>1;
    }
  }
  if(N_modes!=pow(2, N_hierarchy-1)){
    std::cerr<<"ProcessTensorStream::add_modes_tree: N_modes="<<N_modes<<" != 2^"<<N_hierarchy-1<<"="<<pow(2, N_hierarchy-1)<<" => not a power of 2!"<<std::endl;
    throw DummyException();
  } 
  if(verbosity>0)std::cout<<"N_modes="<<N_modes<<"=2^"<<N_hierarchy-1<<std::endl;

  if(ProcessTensorStream_ro(file_out).size() != tgrid.n_tot){
    std::cerr<<"ProcessTensorStream::add_modes_tree: ProcessTensorStream_ro(file_out).size() != tgrid.n_tot ("<<ProcessTensorStream_ro(file_out).size()<<" vs. "<<tgrid.n_tot<<")!"<<std::endl;
    throw DummyException();
  }

  TempFileName tmpname0;
  sweep_forward(tmpname0, file_out, trunc, verbosity);

  TempFileName tmpname1;
  add_modes_tree_get(N_hierarchy-1, 0, file_out, mpg, tgrid, trunc, 
                                 intermediate_sweep_n, dict_zero, verbosity);
 
  if(verbosity>0){
    std::cout<<"------------------------"<<std::endl;
    std::cout<<"combine_tree: last step:"<<std::endl;
    std::cout<<"------------------------"<<std::endl;
  }
  sweep_forward(tmpname1, file_out, trunc, verbosity);

  join_select_and_sweep_backward(file_out, tmpname0, tmpname1, trunc, verbosity);
 
  for(int loop=0; loop<intermediate_sweep_n; loop++){
    if(verbosity>0)std::cout<<"sweep intermediate "<<loop<<"/"<<intermediate_sweep_n<<" forward"<<std::endl;
    sweep_forward(tmpname0, file_out, trunc, verbosity);
    if(verbosity>0)std::cout<<"sweep intermediate "<<loop<<"/"<<intermediate_sweep_n<<" backward"<<std::endl;
    sweep_backward(file_out, tmpname0, trunc, verbosity);
  }
}



void ProcessTensorStream::add_modes_tree_get_reverse(
          int level, int first_elem, 
          const std::string & file_out, 
          ModePropagatorGenerator &mpg, const TimeGrid &tgrid,
          const TruncatedSVD & trunc, int intermediate_sweep_n,
          double dict_zero, int verbosity){
                                             
  if(level==0){
    if(first_elem>=mpg.get_N_modes()){
      std::cerr<<"add_modes_tree_get: level="<<level<<" first_elem="<<first_elem<<": first_element>=mpg.get_N_modes()="<<mpg.get_N_modes()<<"!"<<std::endl;
      throw DummyException();
    }
    
    if(verbosity>0){
      std::cout<<"--------------------------------------------"<<std::endl;
      std::cout<<"level: "<<level<<" first_elem: "<<first_elem<<std::endl;
      std::cout<<"--------------------------------------------"<<std::endl;
    }
    if(mpg.skip_list[first_elem]){
      set_trivial(file_out, tgrid.n_tot, mpg.get_N());
    }else{
      ModePropagatorPtr mpp=mpg.get_ModePropagator(first_elem);
      set_from_ModePropagator(file_out, *mpp.get(), tgrid, dict_zero);
    }
  }else{  
    TempFileName tmpname0;
    add_modes_tree_get_reverse(level-1, 2*first_elem, tmpname0, mpg, tgrid, trunc, intermediate_sweep_n, dict_zero, verbosity);

    TempFileName tmpname1;
    add_modes_tree_get_reverse(level-1, 2*first_elem+1, tmpname1, mpg, tgrid, trunc, intermediate_sweep_n, dict_zero, verbosity);

    if(verbosity>0){
      std::cout<<"--------------------------------------------"<<std::endl;
      std::cout<<"level: "<<level<<" first_elem: "<<first_elem<<std::endl;
      std::cout<<"--------------------------------------------"<<std::endl;
    }
    if(verbosity>0)std::cout<<"sweep first backward"<<std::endl;
    copy(file_out, tmpname0);
    sweep_backward(tmpname0, file_out, trunc, verbosity);
    if(verbosity>0)std::cout<<"sweep second backward"<<std::endl;
    copy(file_out, tmpname1);
    sweep_backward(tmpname1, file_out, trunc, verbosity);
    if(verbosity>0)std::cout<<"join and sweep forward"<<std::endl;
    join_select_and_sweep_forward(file_out, tmpname0, tmpname1, trunc, verbosity);
  
    for(int loop=0; loop<intermediate_sweep_n; loop++){
      if(verbosity>0)std::cout<<"sweep intermediate "<<loop<<"/"<<intermediate_sweep_n<<" backward"<<std::endl;
      sweep_backward(tmpname0, file_out, trunc, verbosity);
      if(verbosity>0)std::cout<<"sweep intermediate "<<loop<<"/"<<intermediate_sweep_n<<" forward"<<std::endl;
      sweep_forward(file_out, tmpname0, trunc, verbosity);
    }
  }
}

void ProcessTensorStream::add_modes_tree_reverse(
          const std::string & file_out, 
          ModePropagatorGenerator &mpg, const TimeGrid &tgrid,
          const TruncatedSVD & trunc, int intermediate_sweep_n,
          double dict_zero, int verbosity){

  int N_modes=mpg.get_N_modes();
  if(N_modes<1)return;

  int N_hierarchy=1;
  {int N_shift=N_modes;
    while(N_shift>1){
      N_hierarchy++;
      N_shift=N_shift>>1;
    }
  }
  if(N_modes!=pow(2, N_hierarchy-1)){
    std::cerr<<"ProcessTensorStream::add_modes_tree: N_modes="<<N_modes<<" != 2^"<<N_hierarchy-1<<"="<<pow(2, N_hierarchy-1)<<" => not a power of 2!"<<std::endl;
    throw DummyException();
  } 
  if(verbosity>0)std::cout<<"N_modes="<<N_modes<<"=2^"<<N_hierarchy-1<<std::endl;

  if(ProcessTensorStream_ro(file_out).size() != tgrid.n_tot){
    std::cerr<<"ProcessTensorStream::add_modes_tree: ProcessTensorStream_ro(file_out).size() != tgrid.n_tot ("<<ProcessTensorStream_ro(file_out).size()<<" vs. "<<tgrid.n_tot<<")!"<<std::endl;
    throw DummyException();
  }

  TempFileName tmpname0;
  sweep_backward(tmpname0, file_out, trunc, verbosity);

  TempFileName tmpname1;
  add_modes_tree_get_reverse(N_hierarchy-1, 0, file_out, mpg, tgrid, trunc, 
                                 intermediate_sweep_n, dict_zero, verbosity);
 
  if(verbosity>0){
    std::cout<<"------------------------"<<std::endl;
    std::cout<<"combine_tree: last step:"<<std::endl;
    std::cout<<"------------------------"<<std::endl;
  }
  sweep_backward(tmpname1, file_out, trunc, verbosity);

  join_select_and_sweep_forward(file_out, tmpname0, tmpname1, trunc, verbosity);
 
  for(int loop=0; loop<intermediate_sweep_n; loop++){
    if(verbosity>0)std::cout<<"sweep intermediate "<<loop<<"/"<<intermediate_sweep_n<<" backward"<<std::endl;
    sweep_backward(tmpname0, file_out, trunc, verbosity);
    if(verbosity>0)std::cout<<"sweep intermediate "<<loop<<"/"<<intermediate_sweep_n<<" forward"<<std::endl;
    sweep_forward(file_out, tmpname0, trunc, verbosity);
  }
}



std::string ProcessTensorStream::dims(const std::string & file){
  std::stringstream ss;
  ProcessTensorStream_ro PT_ro(file);
  for(size_t i=0; i<PT_ro.size(); i++){
    if(i==0)ss<<PT_ro.get(0).M.dim_d1;
    ss<<" "<<PT_ro.get(i).M.dim_d2;
//    if(i<PT_ro.size()-1)ss<<"/"<<PT_ro.get(i+1).M.dim_d1;
  }
  return ss.str();
}

}//namespace
