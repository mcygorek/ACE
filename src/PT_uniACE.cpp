#include "ACE.hpp"
#include "PT_uniACE.hpp"
#include "ProcessTensorRepeat.hpp"
#include "Timings.hpp"

namespace ACE{

void PT_uniACE::set_from_ModePropagator(ModePropagator &mpg, const TimeGrid & tgrid, double dict_zero){

  e.set_from_ModePropagator(mpg, tgrid.get_t(0), tgrid.get_dt(0), dict_zero);
  bath_init = H_Matrix_to_L_Vector(mpg.get_bath_init());
  TTree=std::shared_ptr<CompressionTree>(new CompressionTree(bath_init.rows()));
}

void PT_uniACE::join_symmetric_and_sweep_forward(PT_uniACE &other, const TruncatedSVD & trunc, int n_steps){
 
  e.join_symmetric(other.e, other.e);
  bath_init = Vector_otimes(bath_init, other.bath_init);
  int d2_before=e.M.dim_d2;
  if(TTree){TTree->combine(other.TTree);}

  if(n_steps<2){
    std::cerr<<"PT_uniACE::jsasf: n_steps="<<n_steps<<" >= 2 required!"<<std::endl;
    throw DummyException();
  }
  std::cout<<"Join and sweep forward: "; trunc.print_info(std::cout); std::cout<<std::endl;
  VectorSpace VS; double VSthr=1e-8;
  ProcessTensorElement tmp=e;
  tmp.M.inner_multiply_left(bath_init.transpose());
  PassOn pass_on(1);
  tmp.sweep_forward(trunc, pass_on, false);
  VS.add_ortho_Mat_normalize(pass_on.P.adjoint(),VSthr);
  for(int n=1; n<n_steps-1; n++){
    tmp=e; 
    tmp.sweep_forward(trunc, pass_on, false);
    VS.add_ortho_Mat_normalize(pass_on.P.adjoint(),VSthr);
  }
  tmp=e; 
  tmp.sweep_forward(trunc, pass_on, true);
  Eigen::MatrixXcd A=tmp.M.get_Matrix_d1i_d2();
  TruncatedSVD_Return ret=trunc.compress(A);
  //Relevant subspace is now defined by ret.Vdagger
  int d2_after=ret.Vdagger.rows();
  std::cout<<"Maxdim at n="<<n_steps-1<<": "<<d2_before<<" -> "<<d2_after<<std::endl;
  VS.add_ortho_Mat(ret.Vdagger.adjoint(),VSthr);
  std::cout<<"VS.V.cols()="<<VS.V.cols()<<std::endl;
//ret.Vdagger=VS.V.adjoint();
  e.M.inner_multiply_right(ret.Vdagger.adjoint());
  e.M.inner_multiply_left(ret.Vdagger);
  bath_init = (bath_init.transpose()*ret.Vdagger.adjoint()).transpose();
  e.closure = ret.Vdagger * e.closure;
  for(Eigen::VectorXcd & op : e.env_ops.ops){
    op = ret.Vdagger * op;
  }
  //std::cout<<"sigma: "<<ret.sigma.transpose()<<std::endl;
  e.clearNF();
  e.forwardNF = ret.sigma;

  if(TTree){TTree->T = ret.Vdagger * TTree->T;}
  if(TTree){std::cout<<"TTree->T.dims: "<<TTree->T.rows()<<","<<TTree->T.cols()<<std::endl;}
}

void PT_uniACE::sweep_backward(const TruncatedSVD & trunc, int n_steps){
  if(n_steps<2){
    std::cerr<<"PT_uniACE::sweep_backward: n_steps="<<n_steps<<" >= 2 required!"<<std::endl;
    throw DummyException();
  }
  std::cout<<"Sweep backward: "; trunc.print_info(std::cout); std::cout<<std::endl;
  int d1_before=e.M.dim_d1;
  ProcessTensorElement tmp=e;
  tmp.close_off();
  PassOn pass_on(1);
  tmp.sweep_backward(trunc, pass_on, false);
  for(int n=1; n<n_steps-1; n++){
    tmp=e; 
    tmp.sweep_backward(trunc, pass_on, false);
  }
  tmp=e; //.join_thisfirst_sameinner(e);
  tmp.sweep_backward(trunc, pass_on, true);
  Eigen::MatrixXcd A=tmp.M.get_Matrix_d1_id2();
  TruncatedSVD_Return ret=trunc.compress(A);
  //Relevant subspace is now defined by ret.U
  int d1_after=ret.U.cols();
  std::cout<<"Maxdim at n="<<n_steps-1<<": "<<d1_before<<" -> "<<d1_after<<std::endl;
  e.M.inner_multiply_left(ret.U.adjoint());
  e.M.inner_multiply_right(ret.U);
  bath_init = (bath_init.transpose()*ret.U).transpose();
  e.closure = ret.U.adjoint() * e.closure;
  for(Eigen::VectorXcd & op : e.env_ops.ops){
    op = ret.U.adjoint() * op;
  }
  e.clearNF();
  e.backwardNF = ret.sigma;

  if(TTree){TTree->T = ret.U.adjoint() * TTree->T;}
  if(TTree){std::cout<<"TTree->T.dims: "<<TTree->T.rows()<<","<<TTree->T.cols()<<std::endl;}
  //std::cout<<"sigma: "<<ret.sigma.transpose()<<std::endl;
}

void PT_uniACE::sweep_backward_and_eliminate(const TruncatedSVD & trunc, int n_steps, bool stop_if_min_dim){
  if(n_steps<2){
    std::cerr<<"PT_uniACE::sweep_backward_and_eliminate: n_steps="<<n_steps<<" >= 2 required!"<<std::endl;
    throw DummyException();
  }
  if(!e.is_forwardNF()){
    std::cerr<<"PT_uniACE::sweep_backward_and_eliminate: e is not in forward normal form!"<<std::endl;
    throw DummyException();
  }
  std::cout<<"Sweep backward and eliminate: "; trunc.print_info(std::cout); std::cout<<std::endl;

  int d1_before=e.M.dim_d1;
  ProcessTensorElement tmp=e;
  tmp.close_off();
  PassOn pass_on(1);
  tmp.sweep_backward(trunc, pass_on, false);
  for(int n=1; n<n_steps; n++){
    tmp=e; 
    tmp.sweep_backward(trunc, pass_on, false);
    if(stop_if_min_dim && pass_on.P.rows()>d1_before){break;}
  }
  
  std::vector<std::pair<double, int> > sort_weights(e.forwardNF.size());
  for(size_t j=0; j<sort_weights.size(); j++){
    sort_weights[j].first = e.forwardNF[j]*pass_on.P.row(j).norm();
    sort_weights[j].second = j;
  } 
  sort(sort_weights.begin(), sort_weights.end(), [](const auto &a, const auto &b){return a.first>b.first;});

  Eigen::VectorXd weights(e.forwardNF.size());
  for(size_t j=0; j<sort_weights.size(); j++){
    weights[j] = sort_weights[j].first;
  }
  int d1_after=trunc.get_truncated_dim(weights);
  std::cout<<"Maxdim at n="<<n_steps-1<<": "<<d1_before<<" -> "<<d1_after<<std::endl;
  //std::cout<<weights.transpose()<<std::endl;

  
  Eigen::MatrixXcd T=Eigen::MatrixXcd::Zero(d1_after, d1_before);
  for(int j=0; j<d1_after; j++){
    T(j, sort_weights[j].second)=1;
  }

  e.M.inner_multiply_left(T);
  e.M.inner_multiply_right(T.transpose());
  bath_init = T * bath_init;
  e.closure = T * e.closure;
  for(Eigen::VectorXcd & op : e.env_ops.ops){
    op = T * op;
  }
  e.clearNF();
  e.backwardNF = weights.head(d1_after);

  if(TTree){TTree->T = T * TTree->T;}
  if(TTree){std::cout<<"TTree->T.dims: "<<TTree->T.rows()<<","<<TTree->T.cols()<<std::endl;}
  //std::cout<<"sigma: "<<ret.sigma.transpose()<<std::endl;
}

void PT_uniACE::add_modes(ModePropagatorGenerator &mpg,
                 const TimeGrid & tgrid, const TruncationLayout & trunc,
                 double dict_zero, int verbosity){
  std::cout<<"PT_uniACE::add_modes: mpg.get_N_modes()="<<mpg.get_N_modes()<<std::endl;

  TimeGrid tgrid2=tgrid.construct_half_dt();
  trunc.print_info();std::cout<<std::endl;

  if(mpg.get_N_modes()>0 && e.M.dim_i<2){
    e.set_trivial(mpg.get_N());
    bath_init = Eigen::VectorXcd::Ones(1);
  }
  for(int k=mpg.first(); k<mpg.get_N_modes(); k=mpg.next(k)){
    if(verbosity>0){
      std::cout<<"Mode "<<k<<"/"<<mpg.get_N_modes()<<std::endl;
    }
//    std::cout<<"bath_init: "<<bath_init.transpose()<<std::endl;
    ModePropagatorPtr mpp=mpg.get_ModePropagator(k);

    PT_uniACE other;
    other.set_from_ModePropagator(*mpp.get(), tgrid2, dict_zero);

    TruncatedSVD trunc_fw=trunc.get_forward(k, mpg.get_N_modes());
    if(trunc.use_QR)trunc_fw.use_QR=true;
    join_symmetric_and_sweep_forward(other, trunc_fw, tgrid.n_mem);

//    TruncatedSVD trunc_bw=trunc.get_backward(k, mpg.get_N_modes());
//    if(trunc.use_QR)trunc_bw.use_QR=true;
//    sweep_backward_and_eliminate(trunc_bw, tgrid.n_mem);

  }
}

std::shared_ptr<ProcessTensorRepeat> PT_uniACE::get_PTR(const std::string &write_PT) const{
  std::shared_ptr<ProcessTensorRepeat> PTR(new ProcessTensorRepeat(e.get_N()));
  PTR->set_specs(write_PT, -1);
  PTR->initial.resize(1);
  PTR->initial.get(0)=e;
  //e.M.print_dims(); std::cout<<std::endl;
  PTR->initial.get(0).M.inner_multiply_left(bath_init.transpose());
  PTR->repeated.resize(1);
  PTR->repeated.get(0)=e;

  return PTR;
}

} //namespace
