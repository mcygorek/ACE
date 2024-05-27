#include "Parameters.hpp"
#include "Reader.hpp"
#include "ProcessTensorRepeat.hpp"
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

  ProcessTensorRepeat PTR(read_PT,true);
  PTR.print_info(); std::cout<<std::endl;
  if(PTR.initial.size()<1 || PTR.repeated.size()<1){
    std::cout<<"Process Tensor empty!"<<std::endl;
    return 0;
  }

  {
    const ProcessTensorElement & e=PTR.get_ro(0, ForwardPreload);
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
  for(int n=0; n<PTR.initial.size(); n++){
    const ProcessTensorElement & e=PTR.initial.get_ro(n, ForwardPreload);
    dims.push_back(std::vector<int>{e.M.dim_i, e.M.dim_d1, e.M.dim_d2});
    if(n==0)std::cout<<dims.back()[1];
    std::cout<<" "<<dims.back()[2]<<std::flush;
    if(maxdim<dims.back()[2]){
      maxdim=dims.back()[2];
      maxdim_at=n;
    }
  }
  std::cout<<std::endl;
  for(int n=0; n<PTR.repeated.size(); n++){
    const ProcessTensorElement & e=PTR.repeated.get_ro(n, ForwardPreload);
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
      nr_print_SVD=PTR.initial.size()+PTR.repeated.size()/2-1;
      if(nr_print_SVD<0)nr_print_SVD=0;//could happen if n_tot==1
    }
    if(param.get_as_bool("SVD_at_maxdim",false)){
      nr_print_SVD=maxdim_at;
    }

    Eigen::VectorXd SVs;
    if(PTR.get_ro(nr_print_SVD).is_forwardNF()){
      std::cout<<"PT is in forward NF"<<std::endl;
      SVs=PTR.get_ro(nr_print_SVD).forwardNF;

    }else if(PTR.get_ro(nr_print_SVD).is_backwardNF()){
      std::cout<<"PT is in backward NF"<<std::endl;
      SVs=PTR.get_ro(nr_print_SVD).backwardNF;
   
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
