#include "ACE.hpp"
#include "CompressedPropagator.hpp"
#include "TimeGrid.hpp"
#include "ProcessTensorElement.hpp"
#include "ProcessTensorRepeat.hpp"
#include "LiouvilleTools.hpp"

using namespace ACE;
using namespace ACE::CompressedPropagator;

PTE_rhoE0 calculate_center_Q_ACE_select_alternate(int n_alt,
                       std::shared_ptr<CompressionTree> TTree, 
                       std::shared_ptr<CompressionTree> TTree_inv,
                       std::shared_ptr<ModePropagatorGenerator> mpg, 
                       int mode, double ta, double dt, double dict_zero){

  std::cout<<"calculate_center_Q_ACE_select_alternate: mode="<<mode<<std::endl;

  if(mode<0){
    std::cerr<<"calculate_center_Q_ACE_select_alternate: mode<0!"<<std::endl;
    throw DummyException();
  }else if(mode==0){
    return calculate_single(TTree, TTree_inv, mpg->get_ModePropagator(mode), ta, dt, dict_zero);
  }

  PTE_rhoE0 e_rhoE0=calculate_center_Q_ACE_select_alternate(n_alt,TTree->first,\
                     TTree_inv->first, mpg, mode-1, ta, dt, dict_zero);

  PTE_rhoE0 e_rhoE0_2=calculate_single(TTree->second, TTree_inv->second, \
                            mpg->get_ModePropagator(mode), ta, dt, dict_zero);

  std::cout<<"first:d1="<<e_rhoE0.first.M.dim_d1;
  std::cout<<" first:d2="<<e_rhoE0.first.M.dim_d2;
  std::cout<<" second:d1="<<e_rhoE0_2.first.M.dim_d1;
  std::cout<<" second:d2="<<e_rhoE0_2.first.M.dim_d2;
  std::cout<<" TTree->S.size()="<<TTree->S.size();
  std::cout<<" TTree->S.min_dim_first()="<<TTree->S.min_dim_first();
  std::cout<<" TTree->S.min_dim_second()="<<TTree->S.min_dim_second();
  std::cout<<" TTree_inv->S.size()="<<TTree_inv->S.size();
  std::cout<<" TTree_inv->S.min_dim_first()="<<TTree_inv->S.min_dim_first();
  std::cout<<" TTree_inv->S.min_dim_second()="<<TTree_inv->S.min_dim_second();
  std::cout<<std::endl;

//std::cout<<"MARK0"<<std::endl;
  e_rhoE0.first.join_selected(n_alt,e_rhoE0_2.first, TTree_inv->S, TTree->S );
//std::cout<<"MARK1"<<std::endl;
  Eigen::VectorXcd rhoE0(TTree->S.size());
  for(int k=0; k<(int)TTree->S.size(); k++){
    rhoE0(k) = e_rhoE0.second(TTree->S[k].first)*e_rhoE0_2.second(TTree->S[k].second);
  }
//std::cout<<"MARK2"<<std::endl;
  e_rhoE0.second=rhoE0;
//std::cout<<"MARK3"<<std::endl;
 
  compress(e_rhoE0, TTree, TTree_inv);
//std::cout<<"MARK4"<<std::endl;


  std::cout<<"calculate_center_Q_ACE_select_alternate: done"<<std::endl;
  return e_rhoE0;
}
PTE_rhoE0 calculate_center_Q_ACE_select_thissecond(
                       std::shared_ptr<CompressionTree> TTree, 
                       std::shared_ptr<CompressionTree> TTree_inv,
                       std::shared_ptr<ModePropagatorGenerator> mpg, 
                       int mode, double ta, double dt, double dict_zero){
  return calculate_center_Q_ACE_select_alternate(0, TTree, TTree_inv, mpg, mode, ta, dt, dict_zero);
}
PTE_rhoE0 calculate_center_Q_ACE_select_thisfirst(
                       std::shared_ptr<CompressionTree> TTree, 
                       std::shared_ptr<CompressionTree> TTree_inv,
                       std::shared_ptr<ModePropagatorGenerator> mpg, 
                       int mode, double ta, double dt, double dict_zero){
  return calculate_center_Q_ACE_select_alternate(1, TTree, TTree_inv, mpg, mode, ta, dt, dict_zero);
}


int main(int args, char ** argv){
  Parameters param(args, argv, true);
 
  std::string infile=param.get_as_string_check("infile");

  std::shared_ptr<CompressionTree> TTree(new CompressionTree(infile));
  std::shared_ptr<CompressionTree> TTree_inv(new CompressionTree(infile+"_inv"));

  TTree->print_info();
  TTree_inv->print_info();

  std::string write_PT=param.get_as_string("write_PT","");
  if(write_PT!=""){
    std::cout<<"write_PT is set to '"<<write_PT<<"'"<<std::endl;
    TimeGrid tgrid(param);
    double dict_zero=param.get_as_double("dict_zero",-1.);
    std::vector<std::shared_ptr<ModePropagatorGenerator> > mpgs=MPG_Selector(param);
    if(mpgs.size()>0 && mpgs[0]){
      ProcessTensorRepeat PTR;
      PTR.set_specs(write_PT, param.get_as_double("buffer_blocksize",-1));
      PTR.initial.resize(1);
      PTR.repeated.resize(1);

      int N_modes=mpgs[0]->get_N_modes();
      std::cout<<"N_modes="<<N_modes<<std::endl;
      PTE_rhoE0 e_rhoE0;

      if(param.get_as_bool("use_combine_tree")){
        std::cerr<<"use_combine_tree NOT IMPLEMENTED YET!"<<std::endl;
        throw DummyException();

      }else if(param.get_as_bool("use_select", true)){

        bool use_symmetric_Trotter=param.get_as_bool("use_symmetric_Trotter", true);
        if(use_symmetric_Trotter){
std::cout<<"use_select and use_symmetric_Trotter"<<std::endl;
          e_rhoE0=calculate_center_Q_ACE_select_thissecond(TTree, TTree_inv,\
             mpgs[0], N_modes-1, tgrid.ta, tgrid.dt/2., dict_zero);
PTE_rhoE0 e_rhoE0_2=calculate_center_Q_ACE_select_thisfirst(TTree, TTree_inv,\
             mpgs[0], N_modes-1, tgrid.ta+tgrid.dt/2., tgrid.dt/2., dict_zero);
          e_rhoE0.first.join_thisfirst_sameinner(e_rhoE0_2.first);
        }else{
          if(param.get_as_bool("use_symmetric_thisfirst",false)){
            e_rhoE0=calculate_center_Q_ACE_select_thisfirst(TTree, TTree_inv,\
                    mpgs[0], N_modes-1, tgrid.ta, tgrid.dt, dict_zero);
          }else{
            e_rhoE0=calculate_center_Q_ACE_select_thissecond(TTree, TTree_inv,\
                    mpgs[0], N_modes-1, tgrid.ta, tgrid.dt, dict_zero);
          }         
        }
      }else{
        // must be sequential ACE
        e_rhoE0=CompressedPropagator::calculate_ACE_sequential(
                           TTree, TTree_inv,\
                           mpgs[0], N_modes-1, tgrid.ta, tgrid.dt, dict_zero);
      }
      std::cout<<"e_rhoE0.first.M.dim_d1="<<e_rhoE0.first.M.dim_d1;
      std::cout<<" e_rhoE0.first.M.dim_d2="<<e_rhoE0.first.M.dim_d2;
      std::cout<<" e_rhoE0.first.closure.rows()="<<e_rhoE0.first.closure.rows();
      std::cout<<" e_rhoE0.second.rows()="<<e_rhoE0.second.rows()<<std::endl;

      std::cout<<"e_rhoE0.second.transpose()*e_rhoE0.first.closure="<<e_rhoE0.second.transpose()*e_rhoE0.first.closure<<std::endl;

      PTR.initial.get(0)=e_rhoE0.first;
      PTR.initial.get(0).M.inner_multiply_left(e_rhoE0.second.transpose());
      std::cout<<"PTR.initial.get(0).M.dim_d1="<<PTR.initial.get(0).M.dim_d1;
      std::cout<<" PTR.initial.get(0).M.dim_d2="<<PTR.initial.get(0).M.dim_d2;
      std::cout<<std::endl;
      PTR.repeated.get(0)=e_rhoE0.first;


      std::string print_eigenvalues=param.get_as_string("print_eigenvalues","");
      if(print_eigenvalues!="" && print_eigenvalues!="/dev/null"){
        ProcessTensorElement &e=PTR.repeated.get(0);
        int NL=e.get_NL();
        int chi=e.M.dim_d2;
        Eigen::MatrixXcd Q=Eigen::MatrixXcd::Zero(NL*chi, NL*chi);
        for(int i2=0; i2<NL; i2++){
          for(int i1=0; i1<NL; i1++){
            int i_ind=e.accessor.dict.beta[i2*NL+i1]; if(i_ind<0)continue;
            for(int d2=0; d2<chi; d2++){
              for(int d1=0; d1<chi; d1++){
                Q(i2*chi+d2, i1*chi+d1)=e.M(i_ind, d1, d2);
              }
            }
          }
        }
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solv(Q);
        std::vector<std::complex<double> > evals(Q.rows());
        for(int i=0; i<Q.rows(); i++){evals[i]=solv.eigenvalues()(i);}
        std::ofstream ofs(print_eigenvalues); 
        for(int i=0; i<(int)evals.size(); i++){
          ofs<<std::real(evals[i])<<" "<<std::imag(evals[i])<<std::endl;
        }
      }
    }
  }

  return 0;
}
