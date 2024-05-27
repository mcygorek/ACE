#include "Parameters.hpp"
#include "Reader.hpp"
#include "ProcessTensorBuffer.hpp"
#include "DummyException.hpp"

using namespace ACE;

int main(int args, char ** argv){
 try{
//  Parameters param(args, argv, true, false);
  Parameters param(args, argv);

  std::string outfile=param.get_as_string("outfile","PT_analyse.out");

  std::string read_PT=param.get_as_string_check("read_PT");


  std::string print_SVD=param.get_as_string("print_SVD","");
  int nr_print_SVD=param.get_as_int("print_SVD",-1,0,1);
/*
  DECISION: don't change PT at all. Use binaries like PTB_sweep_forward instead!
  //force users to acknowledge if the PT is to be modified:
  bool sweep_forward=param.get_as_bool("sweep_forward",false);
  bool sweep_backward=param.get_as_bool("sweep_backward",false);
  if((sweep_forward || sweep_backward) && !print_SVD){
    std::cerr<<"no reason to sweep if 'print_SVD' is not set!"<<std::endl;
    exit(1);
  }
*/

  //DON'T MODIFY PT
  ProcessTensorBuffer PTB(read_PT);
  PTB.print_info(); std::cout<<std::endl;
  PTB.read_only=true;
  if(PTB.get_n_tot()<1){
    std::cout<<"Process Tensor empty!"<<std::endl;
    return 0;
  }

  {
    const ProcessTensorElement & e=PTB.get(0, ForwardPreload);
    std::cout<<"System dimension: "<<e.get_N()<<std::endl;
    std::cout<<"First dictionary: ";e.accessor.dict.print_beta();std::cout<<std::endl;
    std::cout<<"First element is "; if(!e.is_forwardNF())std::cout<<"not ";
    std::cout<<"in forward normal form"<<std::endl;
    std::cout<<"First element is "; if(!e.is_backwardNF())std::cout<<"not ";
    std::cout<<"in backward normal form"<<std::endl;
    std::cout<<"First element env_ops size: "<<e.env_ops.size()<<std::endl;
  }


  std::vector<std::vector<int> > dims; 
  int maxdim=0, maxdim_at=0;
  for(int n=0; n<PTB.get_n_tot(); n++){
    const ProcessTensorElement & e=PTB.get(n, ForwardPreload);
    dims.push_back(std::vector<int>{e.M.dim_i, e.M.dim_d1, e.M.dim_d2});
    if(n==0)std::cout<<dims.back()[1];
    std::cout<<" "<<dims.back()[2]<<std::flush;
    if(maxdim<dims.back()[2]){
      maxdim=dims.back()[2];
      maxdim_at=n;
    }
  }
  std::cout<<std::endl;
  std::cout<<"Maxdim "<<maxdim<<" at "<<maxdim_at<<std::endl;


  if(print_SVD!=""){
    std::ofstream ofs(print_SVD);
    if(nr_print_SVD<0){
      nr_print_SVD=PTB.get_n_tot()/2-1;
      if(nr_print_SVD<0)nr_print_SVD=0;//could happen if n_tot==1
    }
    if(param.get_as_bool("SVD_at_maxdim",false)){
      nr_print_SVD=maxdim_at;
    }

    Eigen::VectorXd SVs;
    if(PTB.get(nr_print_SVD).is_forwardNF()){
      std::cout<<"PT is in forward NF"<<std::endl;
      SVs=PTB.get(nr_print_SVD).forwardNF;

    }else if(PTB.get(nr_print_SVD).is_backwardNF()){
      std::cout<<"PT is in backward NF"<<std::endl;
      SVs=PTB.get(nr_print_SVD).backwardNF;
   
    }else{
      std::cout<<"PT is neither in forward NF nor in backward NF."<<std::endl;
      return 1;
    }
    for(int i=0; i<SVs.rows(); i++){
      ofs<<SVs(i)<<std::endl;
    }
    std::cout<<"Singular values of PT element "<<nr_print_SVD<<" written to '"<<print_SVD<<"'."<<std::endl;
  }
 }catch (DummyException &e){
  return 1;
 }
  return 0;
}
