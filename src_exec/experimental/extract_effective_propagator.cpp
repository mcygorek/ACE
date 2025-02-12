#include "PCH.hpp"
#include "GenericSimulation.hpp"
#include "ProcessTensorForwardList.hpp"
#include "TransferTensor.hpp"
#include "BinaryReader.hpp"
#include "DummyException.hpp"
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>
#include "LindbladMasterEquation.hpp"

using namespace ACE;

int main(int args, char **argv){
  Parameters param(args, argv, true);

  std::string read_from_file=param.get_as_string("read");
  TimeGrid tgrid(param);
  double print_threshold=param.get_as_double("print_threshold",0.);
  
  double t_start=param.get_as_double("t_start", tgrid.ta);
  int n_init=round( (t_start-tgrid.ta)/tgrid.dt );
  if(n_init<0){ 
    std::cout<<"t_start < tgrid.ta -> t_start ignored"<<std::endl;
    n_init=0;
  }
  int n_tot=tgrid.n_tot;
  int n_diff=n_tot-n_init;
  std::cout<<"n_init="<<n_init<<" n_diff="<<n_diff<<" n_tot="<<n_tot<<std::endl;
  if(n_diff<1){
    std::cerr<<"te-t_start<dt!"<<std::endl;
    throw DummyException();
  }

  Eigen::MatrixXcd E;
  if(read_from_file!=""){
    if(n_init>0){
      std::cerr<<"Please don't use 't_start' together with 'read'!"<<std::endl;
      throw DummyException(); 
    }
    std::cout<<"Reading effective propagator from file '"<<read_from_file<<"'."<<std::endl; 
 
    E=read_MatrixXcd(read_from_file);

  }else{
    std::string outfile=param.get_as_string("outfile","Propagator.prop");
    param.override_param("outfile","/dev/null");

    GenericSimulation sim(param);
    std::cout<<"Setting up Process Tensor..."<<std::endl;
    ProcessTensorForwardList PT(param, sim.sysdim);
    PT.print_info();

    std::cout<<"Calculating effective propagator..."<<std::endl;
    std::vector<Eigen::MatrixXcd> vecE=TransferTensor::calculate_E(sim.fprop, PT, sim.sim, tgrid);
   
    if(vecE.size()<n_tot){
      std::cerr<<"vecE.size()<n_tot!"<<std::endl;
      throw DummyException();
    }
    if(n_init>0){
      Eigen::MatrixXcd prevE=vecE[n_init-1];
//      E=prevE.fullPivLu().solve(vecE[ntot-1]);

      Eigen::JacobiSVD<Eigen::MatrixXcd> svd(prevE,Eigen::ComputeFullU|Eigen::ComputeFullV);
      int odim=svd.singularValues().rows();
      int dim=odim;
      std::cout<<"Inverting E(t_start, ta): SVD: "<<svd.singularValues().transpose()<<std::endl;
/*      for(int i=1; i<odim; i++){
        if(svd.singularValues()(i)<1e-9){
          dim=i; break;
        }
      }*/
      Eigen::VectorXd sigma_inv(dim);
      for(int i=0; i<dim; i++){
        sigma_inv(i)=1./svd.singularValues()(i);
      }
      std::cout<<"1/sigma_min="<<sigma_inv(dim-1)<<std::endl;
            
      Eigen::MatrixXcd Einv=svd.matrixV().block(0,0,odim,dim)*sigma_inv.asDiagonal()*svd.matrixU().adjoint().block(0,0,dim, odim);

      E=vecE[n_tot-1]*Einv;

    }else{
      E=vecE[n_tot-1];
    }
    write_MatrixXcd(E, outfile);
    std::cout<<"Effective propagator written to file '"<<outfile<<"'."<<std::endl; 
  }
//  std::cout<<"Propagator:"<<std::endl;
//  std::cout<<E<<std::endl;

  bool analyze_LME=param.get_as_bool("analyze_LME",false);
  double threshold=param.get_as_double("threshold",0.);
  if(analyze_LME){

    double t_log=tgrid.get_t(n_tot)-tgrid.get_t(n_init);
    Eigen::MatrixXcd L=E.log()/t_log;
    std::cout<<"t_log="<<t_log<<std::endl;

    //Fitting Liouvillian 
    LindbladMasterEquation lme(L, threshold, 1);
    lme.print(std::cout, print_threshold);
    std::string print_param=param.get_as_string("print_param");
    bool print_Markov=param.get_as_bool("print_Markov");
    if(print_param!=""){
      lme.print_param(print_param, tgrid, print_Markov, print_threshold);
    }

    //difference:
    Eigen::MatrixXcd Lrec=lme.construct_Liouvillian();

    std::cout<<"|L-Lrec|="<<(L-Lrec).norm()<<std::endl;
    std::cout<<"|L-Lrec|/|L|="<<(L-Lrec).norm()/L.norm()<<std::endl;
    std::cout<<"Tr(L-Lrec)="<<(L-Lrec).trace()<<std::endl;
  }

  return 0;
}
