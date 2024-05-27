#include "ProcessTensorRepeat.hpp"
#include "ProcessTensorForwardList.hpp"
#include "Parameters.hpp"
#include "DummyException.hpp"
#include "TimeGrid.hpp"
#include "Simulation_PT.hpp"
#include "OutputPrinter.hpp"
#include <sstream>
#include <complex>
#include <Eigen/Eigenvalues>


using namespace ACE;
int main(int args, char** argv){

  Parameters param(args, argv, true);
  param.add_to("use_Gaussian_repeat", "true");

  TimeGrid tgrid(param);
  FreePropagator fprop(param);
  int N=fprop.get_dim(); int NL=N*N;

  ProcessTensorForwardList PT(param, N);
  if(PT.list.size()<1){
    std::cerr<<"No ProcessTensor specified!"<<std::endl;
    exit(1);
  }
  if(PT.list.size()>1){
    std::cerr<<"PT.list.size()="<<PT.list.size()<<">1: not yet implemented!"<<std::endl;
    exit(1);
  }
  ProcessTensorRepeat *PTR = dynamic_cast<ProcessTensorRepeat*>(PT.list[0].get()); 
  if(PTR==NULL){
    std::cerr<<"ProcessTensor: not a 'ProcessTensorRepeat'!"<<std::endl;
    exit(1);
  }

  int n_init=PTR->initial.size();
  int n_rep=PTR->repeated.size();
  if(n_rep<1){
    std::cerr<<"PTR->repeated.size()<1!"<<std::endl;
  }
/*
  for(int n=0; n<n_rep; n++){
    std::cout<<"PTR->repeated.get_ro("<<n<<").M.dim_d1="<<PTR->repeated.get_ro(n).M.dim_d1<<", M.dim_d2="<<PTR->repeated.get_ro(n).M.dim_d2<<std::endl;
  }
*/

  int D=PTR->repeated.get_ro(0).M.dim_d1;
  {
    int d2=PTR->repeated.get_ro(n_rep-1).M.dim_d2;
    if(D!=d2){
      std::cerr<<"repeated: first M.dim_d1="<<D<<"!=last M.dim_d2="<<d2<<"!"<<std::endl,
      exit(1);
    }
  }
  std::cout<<"Inner dimension at interface between repeated blocks: "<<D<<std::endl;
  Simulation_PT sim(param);

  TimeGrid tgrid2; tgrid2.set_default(n_rep, tgrid.dt, 0); //n_init*tgrid.dt);
 
  Eigen::MatrixXcd BigProp=Eigen::MatrixXcd::Zero(NL*D,NL*D);
  //construct effective Liouvillian using a basis of system density matrices:
  for(int b=0; b<NL; b++){
    for(int d=0; d<D; d++){
      Eigen::MatrixXcd state=Eigen::MatrixXcd::Zero(NL, D);
      state(b,d)=1.;
  
      PT.reset(); PT.list[0]->n=n_init; // set current n
      for(int n=0; n<n_rep; n++){
        sim.propagate_state(state, n, tgrid2, fprop, PT);

//        printer.print(n+1, tgrid.get_t(n+1), PT.get_rho_reduced(state), 
//                       PT.get_env_reduced(state, printer.which_env_ops));
        PT.load_next();
      }
      for(int b2=0; b2<NL; b2++){
        for(int d2=0; d2<D; d2++){ 
          BigProp(b2*D+d2, b*D+d)=state(b2,d2);
        }
      }
    }
  }
  std::cout<<"Propagator calculated."<<std::endl;
  
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solv(BigProp);
  Eigen::VectorXcd EV=solv.eigenvalues();

  std::string dump_EV=param.get_as_string("dump_EV");
  if(dump_EV!=""){ 
    std::ofstream ofs(dump_EV.c_str());
    for(int i=0; i<EV.rows(); i++){
      ofs<<EV(i).real()<<" "<<EV(i).imag()<<std::endl;
    }
  }
  
  std::cout<<std::endl<<"Largest Eigenvalues:"<<std::endl;
  for(int i=0; i<5; i++){
    int j=EV.rows()-1-i;
    if(j>0){
      std::cout<<EV(j).real()<<" "<<EV(j).imag()<<std::endl;
    }
  }
 
  //get eq_state
  Eigen::MatrixXcd eq_state(NL, D);
  int which_EV=EV.rows()-1-param.get_as_size_t("which_EV",0);
  if(which_EV<0){
    std::cerr<<"which_EV="<<param.get_as_size_t("which_EV")<<">=EV.rows()="<<EV.rows()<<"!"<<std::endl; 
    exit(1);
  }
  for(int b=0; b<NL; b++){
    for(int d=0; d<D; d++){
      eq_state(b,d)=solv.eigenvectors()(b*D+d, which_EV);
    }
  }

  Output_Ops ops(param);
  PT.list[0]->n=n_init+n_rep-1;
  Eigen::MatrixXcd rho_reduced=L_Vector_to_H_Matrix(PT.get_rho_reduced(eq_state));
  //Fix free phase of Liouville space eigenvector using trace:
  { std::complex<double> c=rho_reduced.trace(); 
    rho_reduced*=std::conj(c)/abs(c)/abs(c);
  }
  std::cout<<"Density matrix:"<<std::endl<<rho_reduced<<std::endl;
//  std::cout<<"rho_reduced.rows()="<<rho_reduced.rows()<<" rho_reduced_cols()="<<rho_reduced.cols()<<std::endl;

  std::cout<<"Observables:"<<std::endl;
  for(const Eigen::MatrixXcd &op : ops.ops){
//    std::cout<<"op.rows()="<<op.rows()<<" op.cols()="<<op.cols()<<std::endl;
    std::complex<double> c=(op*rho_reduced).trace();
    std::cout<<c.real()<<" "<<c.imag()<<std::endl;
  }


  return 0;
}
