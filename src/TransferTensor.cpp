#include "TransferTensor.hpp"
//#include "OutputExtractor.hpp"
#include "Simulation_PT.hpp"
#include "LiouvilleTools.hpp"
#include "DummyException.hpp"
#include <memory>


namespace ACE {

std::vector<Eigen::MatrixXcd> TransferTensor::calculate_E(Propagator &prop, 
 ProcessTensorForwardList &PT, Simulation_PT &sim, const TimeGrid &tgrid,
 int verbosity){

  int N=prop.get_dim(); int NL=N*N;
  if(N<2){
    std::cerr<<"TransferTensor::calculate: N<2!"<<std::endl;
    throw DummyException();
  }
  if(tgrid.n_tot<1){
    std::cerr<<"TransferTensor::calculate: tgrid.n_tot<2!"<<std::endl;
    throw DummyException();
  }

  //Initialize E
  std::vector<Eigen::MatrixXcd> E(tgrid.n_tot, Eigen::MatrixXcd::Zero(NL,NL));
 
  //Calculate E
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      if(verbosity>0){
        std::cout<<"("<<i<<","<<j<<")/("<<N<<","<<N<<")"<<std::endl;
      }
      Eigen::MatrixXcd initial_rho=Eigen::MatrixXcd::Zero(N,N);
      initial_rho(i,j)=1;

      OutputPrinter extractor; 
      extractor.do_extract=true; extractor.start_extract=0;
      sim.run_std(prop, PT, initial_rho, tgrid, extractor);
//std::cout<<"test: extractor.rho_t.size()="<<extractor.rho_t.size()<<std::endl;
//std::cout<<"test: tgrid.n_tot="<<tgrid.n_tot<<std::endl;

      for(int l=0; l<tgrid.n_tot; l++){
        E[l].col(i*N+j)=extractor.rho_t[l+1];  //Note: rho_t[0]~Identity
      }
    }
  }
  return E;
}

void TransferTensor::set_T_from_E(const std::vector<Eigen::MatrixXcd> &E, int n_max){
  if(n_max<=0){
    n_max=E.size();
  }
  if(E.size()<1){
    std::cerr<<"TransferTensor::set_T_from_E: E.size()<1!"<<std::endl;
    throw DummyException();
  }
  int NL=E[0].rows();
  if(n_max>E.size()){
    std::cerr<<"TransferTensor::set_T_from_E: n_max>E.size()!"<<std::endl;
    throw DummyException();
  }
  
  //calculate T:
  T=std::vector<Eigen::MatrixXcd>(n_max, Eigen::MatrixXcd::Zero(NL,NL));
  for(int l=0; l<n_max; l++){
    T[l]=E[l];

    Eigen::MatrixXcd subtr=Eigen::MatrixXcd::Zero(NL,NL);
    for(int k=l-1; k>=0; k--){
      subtr+=T[k]*E[(l-1)-k];
    }
    T[l]-=subtr;
  }
}


void TransferTensor::set_LT_from_E(const std::vector<Eigen::MatrixXcd> &E, int verbosity){

   if(E.size()<2){
     std::cerr<<"TransferTensor::set_LT_from_E: E.size()<2!"<<std::endl;
     throw DummyException();
   }
   if(verbosity>0){ 
      std::cout<<"Long-time propagator: "<<std::endl;
   }

   Eigen::MatrixXcd E_last=E.back();
   Eigen::MatrixXcd E_penult=E[E.size()-2];
   int NL=E_last.rows();

   Eigen::JacobiSVD<Eigen::MatrixXcd> svd_last(E_last, Eigen::ComputeFullU | Eigen::ComputeFullV);
   Eigen::JacobiSVD<Eigen::MatrixXcd> svd_penult(E_penult, Eigen::ComputeFullU | Eigen::ComputeFullV);

   if(verbosity>0){ 
     std::cout<<"Eigenvalues:"<<std::endl;
     for(int i=0; i<NL; i++){
       std::cout<<svd_last.singularValues()(i)<<"/"<<svd_penult.singularValues()(i)<<"="<<svd_last.singularValues()(i)/svd_penult.singularValues()(i)<<std::endl;
     }
   }
   if(verbosity>1){
     std::cout<<"svd_last.matrixU().adjoint()*svd_penult.matrixU():"<<std::endl<<svd_last.matrixU().adjoint()*svd_penult.matrixU()<<std::endl;
     std::cout<<"svd_last.matrixV().adjoint()*svd_penult.matrixV():"<<std::endl<<svd_last.matrixV().adjoint()*svd_penult.matrixV()<<std::endl;
   }

   Eigen::VectorXcd inv_sigma(NL);
   for(int i=0; i<NL; i++)inv_sigma(i)=1./svd_penult.singularValues()(i);
   
   LT=Eigen::MatrixXcd::Zero(NL,NL);
   LT=E_last*svd_last.matrixV()*(inv_sigma.asDiagonal())*svd_last.matrixU().adjoint();
}

void TransferTensor::calculate(Propagator &prop, ProcessTensorForwardList &PT, 
                               Simulation_PT &sim, const TimeGrid &tgrid){

  std::vector<Eigen::MatrixXcd> E=calculate_E(prop, PT, sim, tgrid);
  set_T_from_E(E);
  set_LT_from_E(E);
}


void TransferTensor::propagate_initial(std::vector<Eigen::VectorXcd> &rho_t, const TimeGrid &tgrid, OutputPrinter &printer)const{

  if(rho_t.size()<1){
    std::cerr<<"TransferTensor::propagate_initial: rho_t.size()<1!"<<std::endl;
    throw DummyException();
  }
  if(T.size()<1){
    std::cerr<<"TransferTensor::propagate_initial: T.size()<1!"<<std::endl;
    throw DummyException();
  }
  if(rho_t.size()>T.size()){
    std::cerr<<"TransferTensor::propagate_initial: rho_t.size()>T.size()!"<<std::endl;
    throw DummyException();
  }

  int NL=rho_t[0].rows(); 
  if(tgrid.n_tot<1){
    std::cerr<<"TransferTensor::propagate: tgrid.n_tot<1!"<<std::endl;
    throw DummyException();
  }
  if(T[0].rows()!=NL){
    std::cerr<<"TransferTensor::propagate: T[0].rows()!=NL!"<<std::endl;
    throw DummyException();
  }

  //rho_t will contain rho at last time steps: rho_t[0]=rho(t_{m-1}), rho_t[1]=rho(t_{m-2}), ...
  //reduce to a single initial density matrix
  rho_t.resize(1);
  //zero-pad 
  rho_t.resize(T.size(), Eigen::VectorXcd::Zero(NL));
//  Eigen::VectorXcd rho=H_Matrix_to_L_Vector(initial_rho);

  printer.print(0, tgrid.get_t(0), rho_t[0]);

  int lim=T.size()-1; if(tgrid.n_tot-1<lim){lim=tgrid.n_tot-1;}
  for(int n=0; n<lim; n++){
    //calculate rho at step t_{n-1}:
    Eigen::VectorXcd rho=Eigen::VectorXcd::Zero(NL);
    for(int k=0; k<n+1; k++){
      rho+=T[k]*rho_t[k];
    }

    //print:
    double t_next=tgrid.get_t(n+1);
    printer.print(n+1, t_next, rho);

    //prepare rho_t for next step:
    for(int l=(int)n+1; l>0; l--){ 
      rho_t[l].swap(rho_t[l-1]);
    } 
    rho_t[0]=rho;
  }
}

void TransferTensor::propagate_LT(std::vector<Eigen::VectorXcd> &rho_t, const TimeGrid &tgrid, OutputPrinter &printer)const{

  if(rho_t.size()<1){
    std::cerr<<"TransferTensor::propagate_LT: rho_t.size()<1!"<<std::endl;
    throw DummyException();
  }

  int NL=rho_t[0].rows(); 
  if(rho_t[0].rows()!=NL ){
    std::cerr<<"TransferTensor::propagate_LT: rho_t[0].rows()!=NL!"<<std::endl;
    throw DummyException();
  }
  if(tgrid.n_tot<1){
    std::cerr<<"TransferTensor::propagate_LT: tgrid.n_tot<1!"<<std::endl;
    throw DummyException();
  }
  if(LT.rows()!=NL || LT.cols()!=NL){
    std::cerr<<"TransferTensor::propagate_LT: T[0].rows()!=NL!"<<std::endl;
    throw DummyException();
  }

  for(int n=0; n<tgrid.n_tot; n++){
    Eigen::VectorXcd rho=LT*rho_t[0];

    //print:
    double t_next=tgrid.get_t(n+1);
    printer.print(n+1, t_next, rho);

    //prepare rho_t for next step:
    for(int l=(int)rho_t.size()-1; l>0; l--){ 
      rho_t[l].swap(rho_t[l-1]);
    } 
    rho_t[0]=rho;
  }
}
void TransferTensor::propagate_bulk(std::vector<Eigen::VectorXcd> &rho_t, const TimeGrid &tgrid, OutputPrinter &printer)const{

  if(tgrid.n_tot<1){
    std::cerr<<"TransferTensor::propagate_bulk: tgrid.n_tot<1!"<<std::endl;
    throw DummyException();
  }
  if(rho_t.size()<1){
    std::cerr<<"TransferTensor::propagate_bulk: rho_t.size()<1!"<<std::endl;
    throw DummyException();
  }
  if(T.size()<1){
    std::cerr<<"TransferTensor::propagate_bulk: T.size()<1!"<<std::endl;
    throw DummyException();
  }
  if(rho_t.size()>T.size()){
    std::cerr<<"TransferTensor::propagate_bulk: rho_t.size()>T.size()!"<<std::endl;
    throw DummyException();
  }

  int NL=rho_t[0].rows(); 
  for(int l=0; l<rho_t.size(); l++){
    if(rho_t[l].rows()!=NL ){
      std::cerr<<"TransferTensor::propagate_bulk: rho_t["<<l<<"].rows()!=NL!"<<std::endl;
      throw DummyException();
    } 
    if(T[l].rows()!=NL){
      std::cerr<<"TransferTensor::propagate: T["<<l<<"].rows()!=NL!"<<std::endl;
      throw DummyException();
    }
  }

  for(int n=0; n<tgrid.n_tot; n++){
    //calculate rho at step t_{n-1}:
    Eigen::VectorXcd rho=Eigen::VectorXcd::Zero(NL);
    for(size_t k=0; k<rho_t.size(); k++){
      rho+=T[k]*rho_t[k];
    }

    //print:
    double t_next=tgrid.get_t(n+1);
    printer.print(n+1, t_next, rho);

    //prepare rho_t for next step:
    for(int l=(int)rho_t.size()-1; l>0; l--){ 
      rho_t[l].swap(rho_t[l-1]);
    } 
    rho_t[0]=rho;
  }
}

void TransferTensor::propagate(const Eigen::MatrixXcd &initial_rho, const TimeGrid &tgrid, OutputPrinter &printer)const{
  std::vector<Eigen::VectorXcd> rho_t(1, H_Matrix_to_L_Vector(initial_rho));

  propagate_initial(rho_t, tgrid, printer);
  if(tgrid.n_tot<=(int)T.size()-1){return;}

  TimeGrid tgrid_bulk=tgrid; 
  tgrid_bulk.ta+=tgrid.get_t(T.size()-1); tgrid_bulk.n_tot-=T.size()-1;
  if(use_LT){
    propagate_LT(rho_t, tgrid_bulk, printer);
  }else{
    propagate_bulk(rho_t, tgrid_bulk, printer);
  }
}

void TransferTensor::print_norms(std::ostream &os, double dt)const{
  for(size_t l=0; l<T.size(); l++){
    os<<l*dt<<" "<<(T[l].adjoint()*T[l]).trace().real()<<std::endl;
  }
}
void TransferTensor::print_norms(const std::string &fname, double dt)const{
  std::ofstream ofs(fname.c_str());
  print_norms(ofs, dt);
}

void TransferTensor::read(const std::string &fname){
  std::cerr<<"TransferTensor::read: NOT IMPLEMENTED YET!!!"<<std::endl;
}
 
void TransferTensor::write(const std::string &fname)const{
  std::cerr<<"TransferTensor::write: NOT IMPLEMENTED YET!!!"<<std::endl;
}

bool TransferTensor::get_use_TT(Parameters &param){
  return param.get_as_bool("use_TT",param.get_as_bool("use_LT",false));
}
int TransferTensor::get_TT_n_from(Parameters &param){
  TimeGrid tgrid(param);
  double TT_t_from=param.get_as_double("TT_t_from");   
  int TT_n_from=param.get_as_size_t("TT_n_from", tgrid.get_closest_n(TT_t_from));
  return TT_n_from;
}
int TransferTensor::get_TT_n_mem(Parameters &param, bool do_complain){
  TimeGrid tgrid(param);
  double TT_t_mem=param.get_as_double("TT_t_mem");
  int TT_n_mem=param.get_as_size_t("TT_n_mem", TT_t_mem/tgrid.dt);
  if(do_complain && TT_n_mem<1){
    std::cerr<<"TransferTensor::get_TT_n_mem: TT_n_mem="<<TT_n_mem<<"<1!"<<std::endl;
    throw DummyException();
  }
  return TT_n_mem;
}

Parameters TransferTensor::get_paramPT(Parameters &param){
  Parameters param2(param);

  if(!get_use_TT(param)){return param2;}

  TimeGrid tgrid(param);
  int TT_n_mem=get_TT_n_mem(param, true);
  int TT_n_from=get_TT_n_from(param);

  int lim=tgrid.n_tot;
  if(TT_n_from<lim){
    lim=TT_n_from;
  }
  if(TT_n_mem>lim){
    lim=TT_n_mem;
  }
  if(tgrid.n_tot < TT_n_mem){
    std::cerr<<"TransferTensor::get_paramPT: tgrid.n_tot="<<tgrid.n_tot<<" < TT_n_mem="<<TT_n_mem<<"!"<<std::endl;
    throw DummyException();
  }
  param2.override_param("te", tgrid.get_t(lim));
  return param2;
}

}//namespace
