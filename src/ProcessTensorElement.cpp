#include "ProcessTensorElement.hpp"
#include "ModePropagator.hpp"
#include "otimes.hpp"
#include <iostream>
#include <fstream>
#include "BinaryReader.hpp"
#include "LiouvilleTools.hpp"
#include "DummyException.hpp"
#include "DiagBB.hpp"
#include "QRPinv_struct.hpp"


namespace ACE{

void ProcessTensorElement::check_consistency(const std::string & context)const{
  if(M.dim_d1<1 || M.dim_d2<1 || M.dim_i<1){
    std::cerr<<"ProcessTensorElement::check_consistency";
    if(context!="")std::cerr<<"(context: '"<<context<<"')";
    std::cerr<<": M.dim_d1<1 || M.dim_d2<1 || M.dim_i<1)!"<<std::endl;
    throw DummyException();
  }
  if(M.dim_d2!=closure.rows()){
    std::cerr<<"ProcessTensorElement::check_consistency";
    if(context!="")std::cerr<<"(context: '"<<context<<"')";
    std::cerr<<": M.dim_d2!=closure.rows() ("<<M.dim_d2<<" vs. "<<closure.rows()<<")!"<<std::endl;
    throw DummyException();
  }

//std::cout<<"TeSt0 context '"<<context<<"' "<<M.dim_i<<","<<M.dim_d1<<","<<M.dim_d2<<std::endl;
  if(false){
    double max=M.max_element_abs();
    if(!std::isfinite(max)){
      std::cerr<<"ProcessTensorElement::check_consistency";
      if(context!="")std::cerr<<"(context: '"<<context<<"')";
      std::cerr<<"M has non-finite element!"<<std::endl;
      throw DummyException();            
    }

    std::cout<<"check_consistency (context: '"<<context<<"'): max_abs="<<max<<std::endl;
    for(int r=0; r<closure.rows(); r++){
      if(std::isnan(abs(closure(r)))){
        std::cerr<<"ProcessTensorElement::check_consistency";
        if(context!="")std::cerr<<"(context: '"<<context<<"')";
        std::cerr<<"closure entry "<<r<<" is nan!"<<std::endl;
        throw DummyException();            
      }
    }
  }
//std::cout<<"TeSt1"<<std::endl;
}

void ProcessTensorElement::clearNF(){
  forwardNF=Eigen::VectorXd(0);
  backwardNF=Eigen::VectorXd(0);
}
bool ProcessTensorElement::is_forwardNF()const{
  return (forwardNF.size()>0 && forwardNF.size()==M.dim_d2);
}
bool ProcessTensorElement::is_backwardNF()const{
  return (backwardNF.size()>0 && backwardNF.size()==M.dim_d1);
}
void ProcessTensorElement::printNF(std::ostream &ofs)const{
  if(is_forwardNF() && is_backwardNF()){
    ofs<<"PT element is both forwardNF and backwardNF"<<std::endl;
  }else if(is_forwardNF()){
    ofs<<"PT element is forwardNF"<<std::endl;
  }else if(is_backwardNF()){
    ofs<<"PT element is backwardNF"<<std::endl;
  }else{
    ofs<<"PT element is neither forwardNF nor backwardNF"<<std::endl;
  }
}

void ProcessTensorElement::join_thisfirst(const ProcessTensorElement &other){
  accessor.join_thisfirst(other.accessor, M, other.M);
  closure=Vector_otimes(closure, other.closure);
  env_ops.join(other.env_ops);  

  clearNF();
}
void ProcessTensorElement::join_thissecond(const ProcessTensorElement &other){
  accessor.join_thissecond(other.accessor, M, other.M);
  closure=Vector_otimes(closure, other.closure);
  env_ops.join(other.env_ops);  

  clearNF();
}

void ProcessTensorElement::join_symmetric(
                         const ProcessTensorElement &other_left,
                         const ProcessTensorElement &other_right){

  accessor.join_symmetric(other_left.accessor, other_right.accessor,
                          M, other_left.M, other_right.M);

  closure=Vector_otimes(closure, other_right.closure);
  env_ops.join(other_right.env_ops);  

  clearNF();
}

void ProcessTensorElement::join_thisfirst_sameinner(const ProcessTensorElement &e2){
  if(M.dim_d2 != e2.M.dim_d1){
    std::cerr<<"ProcessTensorElement::join_thisfirst_sameinner: M.dim_d2 != e2.M.dim_d1!"<<std::endl;
    throw DummyException();
  }
  
  ProcessTensorElementAccessor::VVPI i_list=accessor.join_thisfirst_indices(e2.accessor);

  MPS_Matrix tmp(i_list.size(), M.dim_d1, e2.M.dim_d2);
  tmp.set_zero();
  for(int i=0; i<(int)i_list.size(); i++){
    for(int j=0; j<(int)i_list[i].size(); j++){
      for(int d1=0; d1<M.dim_d1; d1++){
        for(int d2=0; d2<M.dim_d2; d2++){
          for(int d3=0; d3<e2.M.dim_d2; d3++){
            tmp(i, d1, d3)+= M(i_list[i][j].first, d1, d2)*e2.M(i_list[i][j].second, d2, d3);
          }
        }
      }
    }
  }
  M.swap(tmp);
  closure=e2.closure;
  env_ops=e2.env_ops;
  forwardNF=e2.forwardNF;
}


SelectIndices ProcessTensorElement::get_forwardNF_selected_indices(
          const ProcessTensorElement &other, const TruncatedSVD &trunc)const{

  if(!is_forwardNF()){
    std::cerr<<"ProcessTensorElement::get_forwardNF_selected_indices: this MPS_Matrix is not in forward normal form!"<<std::endl;
    throw DummyException();
  }
  if(!other.is_forwardNF()){
    std::cerr<<"ProcessTensorElement::get_forwardNF_selected_indices: other MPS_Matrix not in forward normal form!"<<std::endl;
    throw DummyException();
  }  

  return trunc.get_select_indices(forwardNF, other.forwardNF);
}

SelectIndices ProcessTensorElement::get_backwardNF_selected_indices(
          const ProcessTensorElement &other, const TruncatedSVD &trunc)const{

  if(!is_backwardNF()){
    std::cerr<<"ProcessTensorElement::get_backwardNF_selected_indices: this MPS_Matrix is not in backward normal form!"<<std::endl;
    exit(1);
  }
  if(!other.is_backwardNF()){
    std::cerr<<"ProcessTensorElement::get_backwardNF_selected_indices: other MPS_Matrix not in backward normal form!"<<std::endl;
    exit(1);
  }  

  return trunc.get_select_indices(backwardNF, other.backwardNF);
}


void ProcessTensorElement::join_selected(
                              int n, const ProcessTensorElement &other, 
                              const SelectIndices & k_list_left, 
                              const SelectIndices & k_list_right){
//std::cout<<"tEsT0"<<std::endl;

  accessor.join_select_indices_alternate(n, M, other.M, other.accessor, k_list_left, k_list_right);
//std::cout<<"tEsT1"<<std::endl;

  //closure=k_list_right.Vector_otimes(closure, other.closure);
  Eigen::VectorXcd new_closure(k_list_right.size());
  for(int kr=0; kr<k_list_right.size(); kr++){
    new_closure(kr) = closure(k_list_right[kr].first)
                    * other.closure(k_list_right[kr].second);
  }
  closure=new_closure;

//std::cout<<"tEsT2"<<std::endl;
  env_ops.join_select_indices(other.env_ops, k_list_right);


//std::cout<<"tEsT3"<<std::endl;
  if(is_forwardNF() && other.is_forwardNF()){
    Eigen::VectorXd new_forwardNF(k_list_right.size());
    for(int k=0; k<k_list_right.size(); k++){
      new_forwardNF(k) = forwardNF(k_list_right[k].first)
                         * other.forwardNF(k_list_right[k].second);
    }
    clearNF();
    forwardNF=new_forwardNF;
  }else if(is_backwardNF() && other.is_backwardNF()){
    Eigen::VectorXd new_backwardNF(k_list_left.size());
    for(int k=0; k<k_list_left.size(); k++){
      new_backwardNF(k) = backwardNF(k_list_left[k].first)
                          * other.backwardNF(k_list_left[k].second);
    }
    clearNF();
    backwardNF=new_backwardNF;
/*  }else{
//    if(k_list_left.size()==M.dim_d1
    std::cerr<<"n="<<n<<" ";
    std::cerr<<"k_list_left.size()="<<k_list_left.size()<<" ";
    std::cerr<<"k_list_right.size()="<<k_list_right.size()<<" ";
    std::cerr<<"M.dim_d1="<<M.dim_d1<<" M.dim_d2="<<M.dim_d2<<" ";
    std::cerr<<"other.M.dim_d1="<<other.M.dim_d1<<" ";
    std::cerr<<"other.M.dim_d2="<<other.M.dim_d2<<std::endl;   
    std::cerr<<"ProcessTensorElement::join_selected: neither forwardNF nor backwardNF!"<<std::endl;
    throw DummyException();
*/
  }
//std::cout<<"tEsT4"<<std::endl;
}
void ProcessTensorElement::join_average_selected(
                              const ProcessTensorElement &other,
                              const SelectIndices & k_list_left,
                              const SelectIndices & k_list_right){

  ProcessTensorElement e_thissecond(*this);
  e_thissecond.join_selected(1, other, k_list_left, k_list_right);
  
  join_selected(0, other, k_list_left, k_list_right);

  MPS_Matrix tmp(M.dim_i, M.dim_d1, M.dim_d2);
  for(int i=0; i<M.dim_i; i++){
    for(int d1=0; d1<M.dim_d1; d1++){
      for(int d2=0; d2<M.dim_d2; d2++){
        tmp(i,d1,d2)=0.5*(M(i,d1,d2)+e_thissecond.M(i,d1,d2));
      }
    }
  }
  M.swap(tmp);
}


void ProcessTensorElement::sweep_forward(const TruncatedSVD &trunc, PassOn &pass_on, bool is_last){
  M.inner_multiply_left(pass_on.P);
  if(is_last){
    if(!is_forwardNF()){
      forwardNF=Eigen::VectorXd::Ones(M.dim_d2);
    }
    pass_on.set(M.dim_d2);
    return;
  }
  clearNF();

  Eigen::MatrixXcd A=M.get_Matrix_d1i_d2();
  pass_on=trunc.compress_forward(A, forwardNF);
  double keep=forwardNF(0);
  if(trunc.keep>0){
    keep=trunc.keep;
//  }else if(closure.size()==pass_on.P.cols()){
//    Eigen::VectorXcd new_closure=pass_on.P*closure;
//    std::cout<<"#"<<new_closure.norm()<<"#";
//    keep=new_closure.norm();
  }
  if(fabs(keep-1.)>1e-6){
    pass_on.P/=keep;
    pass_on.Pinv*=keep;
    A*=keep;
  }
  M.set_from_Matrix_d1i_d2(A, M.dim_i);

  if(closure.size()==pass_on.P.cols())closure=pass_on.P*closure;
  env_ops.process_forward(pass_on);
}

void ProcessTensorElement::sweep_backward(const TruncatedSVD &trunc, PassOn &pass_on, bool is_last){
  if(pass_on.P.rows()!=M.dim_d2){
    std::cerr<<"ProcessTensorElement::sweep_backward: pass_on.P.rows()="<<pass_on.P.rows()<<"!=M.dim_d2="<<M.dim_d2<<"!"<<std::endl;
    throw DummyException();
  }
  M.inner_multiply_right(pass_on.P);
  if(closure.size()==pass_on.Pinv.cols())closure=pass_on.Pinv*closure;
  env_ops.process_backward(pass_on);
  if(is_last){
    clearNF();
    if(!is_backwardNF()){
      backwardNF=Eigen::VectorXd::Ones(M.dim_d1);
    }
    pass_on.set(M.dim_d1);
    return;
  }
  clearNF();

  Eigen::MatrixXcd A=M.get_Matrix_d1_id2();
  pass_on=trunc.compress_backward(A, backwardNF);
  double keep=backwardNF(0);
  if(trunc.keep>0)keep=trunc.keep;
  if(fabs(keep-1.)>1e-6){
    pass_on.P/=keep;
    pass_on.Pinv*=keep;
    A*=keep;
  }
  M.set_from_Matrix_d1_id2(A, M.dim_i);
}

void ProcessTensorElement::sweep_forward_QR(const TruncatedSVD &trunc, PassOn &pass_on, bool is_last){
//NOTE: so far: no truncation implemented!

  M.inner_multiply_left(pass_on.P);
  if(is_last){
    if(!is_forwardNF()){
      forwardNF=Eigen::VectorXd::Ones(M.dim_d2);
    }
    pass_on.set(M.dim_d2);
    return;
  }
  clearNF();

  QRPinv_struct QRPinv(M.get_Matrix_d1i_d2());
  if(QRPinv.weights.size()<1){
    std::cerr<<"ProcessTensorElement::sweep_forward_QR: QRPinv.weights.size()<1!"<<std::endl;
    throw DummyException();
  }
  double keep=std::abs(QRPinv.weights(0));
  if(trunc.keep>0)keep=trunc.keep;
//std::cout<<"DEBUG: forward: keep="<<keep<<std::endl;

  M.set_from_Matrix_d1i_d2(QRPinv.Q*keep, M.dim_i);

  forwardNF=Eigen::VectorXd::Zero(QRPinv.weights.size());
  for(size_t i=0; i<forwardNF.size(); i++){
    forwardNF(i)=std::abs(QRPinv.weights(i));
  }

  pass_on.P=QRPinv.RPinv/keep;
  pass_on.Pinv=QRPinv.RPinv_inv*keep;
  if(closure.size()==pass_on.P.cols())closure=pass_on.P*closure;
  env_ops.process_forward(pass_on);
}

void ProcessTensorElement::sweep_backward_QR(const TruncatedSVD &trunc, PassOn &pass_on, bool is_last){
//NOTE: so far: no truncation implemented!
  if(pass_on.P.rows()!=M.dim_d2){
    std::cerr<<"ProcessTensorElement::sweep_backward: pass_on.P.rows()="<<pass_on.P.rows()<<"!=M.dim_d2="<<M.dim_d2<<"!"<<std::endl;
    throw DummyException();
  }
  M.inner_multiply_right(pass_on.P);
  if(closure.size()==pass_on.Pinv.cols())closure=pass_on.Pinv*closure;
  env_ops.process_backward(pass_on);
  if(is_last){
    clearNF();
    if(!is_backwardNF()){
      backwardNF=Eigen::VectorXd::Ones(M.dim_d1);
    }
    pass_on.set(M.dim_d1);
    return;
  }
  clearNF();

  QRPinv_struct QRPinv(M.get_Matrix_d1_id2().transpose());
  if(QRPinv.weights.size()<1){
    std::cerr<<"ProcessTensorElement::sweep_forward_QR: QRPinv.weights.size()<1!"<<std::endl;
    throw DummyException();
  }
  double keep=std::abs(QRPinv.weights(0));
  if(trunc.keep>0)keep=trunc.keep;
//std::cout<<"DEBUG: forward: keep="<<keep<<std::endl;
  
  M.set_from_Matrix_d1_id2(QRPinv.Q.transpose()*keep, M.dim_i);

  backwardNF=Eigen::VectorXd::Zero(QRPinv.weights.size());
  for(size_t i=0; i<forwardNF.size(); i++){
    backwardNF(i)=std::abs(QRPinv.weights(i));
  }

  pass_on.P=QRPinv.RPinv.transpose()/keep;
  pass_on.Pinv=QRPinv.RPinv_inv.transpose()*keep;

}

void ProcessTensorElement::sweep_pair_forward(ProcessTensorElement &e2, const TruncatedSVD &trunc){

  //define: M^{ij}_{dd'} = M^{i}_{dd''} M^{j}_{d'' d'} 
  //SVD:    M^{ij}_{dd'}-> U_{(d,i),k} sigma_{k} V^\dagger_{k, (j,d')}
  //update: M^{i}_{dd''} <- U_{(d,i),k} ; 
  //        M^{j}_{d'' d'} <- sigma_{k} V^\dagger_{k, (j,d')
  //Question: How to update env_ops?
  //Assume original matrix in backwardNF: (RS=Rescaling, e.g., largest SVD)
  //=>  M^{j}_{d'' d'} ~ RS*(\tilde{V}^\dagger)_{d'',(j,d')} 
  //->  T^{-1} removes old M^j and replaces it by sigma V^\dagger:
  //->  T^{-1}= sigma_{k} V^\dagger_{k, (j,d') (M^{j})^{-1}_{d' d''}
  //->  T^{-1}= sigma_{k} V^\dagger_{k, (j,d') \tilde{V}_{(j,d'),d''}/RS
  //->  T^{-1}= sigma_{k} V^\dagger_{k, (j,d') M^\dagger_{(j,d') d''}/RS^2
  //->  T= RS\tilde{V}^\dagger_{(j,d'),d''} V_{k, (j,d') sigma_{k}^{-1}
  //->  T=  M^{j}_{d'' d'} V_{(j,d'),k} sigma_{k}^{-1}
  //->  PassOn P=T^{-1}, Pinv=T

  if(!e2.is_backwardNF()){
    std::cerr<<"ProcessTensorElement::sweep_pair_forward requires e2 in backward normal form!"<<std::endl;
    throw DummyException();
  }
  if(e2.M.dim_d1!=M.dim_d2){
    std::cerr<<"ProcessTensorElement::sweep_pair_forward: e2.M.dim_d1!=M.dim_d2!"<<std::endl;
    throw DummyException();
  }


  //construct Mij
  Eigen::MatrixXcd A1=M.get_Matrix_d1i_d2();
  Eigen::MatrixXcd A2=e2.M.get_Matrix_d1_id2();
  Eigen::MatrixXcd tmp=A1*A2;

  //Extract rescaling of next e2 (used for inverse)  
  double RS=A2.row(0).norm();

  //SVD:
  TruncatedSVD_Return ret=trunc.compress(tmp); 
  clearNF();
  forwardNF=ret.sigma;
   
  double keep=forwardNF(0);
  if(trunc.keep>0)keep=trunc.keep;
  if(fabs(keep-1.)<=1e-6)keep=1.;

  //Adjust closures:
  PassOn pass_on;
  pass_on.P=ret.Vdagger*(A2.adjoint())/RS;
  pass_on.Pinv=pass_on.P.adjoint();
  pass_on.P=(ret.sigma/keep/RS).asDiagonal()*pass_on.P;
  Eigen::VectorXd sinv(ret.sigma.rows());

  for(int k=0; k<sinv.rows(); k++){sinv(k)=RS*keep/ret.sigma(k);}
  pass_on.Pinv=pass_on.Pinv*(sinv.asDiagonal());

  //Set matrices
  M.set_from_Matrix_d1i_d2(ret.U*keep, M.dim_i);
  e2.M.set_from_Matrix_d1_id2(ret.sigma.asDiagonal()*ret.Vdagger/keep, e2.M.dim_i);

  if(closure.size()==pass_on.P.cols()){
    closure=pass_on.P*closure;
  }
  env_ops.process_forward(pass_on);

/*
  if(M.dim_d2!=e2.M.dim_d1){
    std::cerr<<"M.dim_d2="<<M.dim_d2<<" e2.M.dim_d1="<<e2.M.dim_d1<<std::endl;
    throw DummyException();
  }
  if(closure.rows()!=M.dim_d2){
    std::cerr<<"closure.rows()="<<closure.rows()<<" != M.dim_d2="<<M.dim_d2<<std::endl;
    throw DummyException();
  }
  for(int i=0; i<env_ops.ops.size(); i++){
    if(env_ops.ops[i].rows()!=M.dim_d2){
      std::cerr<<"env_ops.ops["<<i<<"].rows()="<<env_ops.ops[i].rows()<<" != M.dim_d2="<<M.dim_d2<<std::endl;
      throw DummyException();
    }
  }
*/
}
void ProcessTensorElement::sweep_pair_backward(ProcessTensorElement &e2, const TruncatedSVD &trunc){

  //define: M^{ij}_{dd'} = M^{i}_{dd''} M^{j}_{d'' d'} 
  //SVD:    M^{ij}_{dd'}-> U_{(d,i),k} sigma_{k} V^\dagger_{k, (j,d')}
  //update: M^{i}_{dd''} <- U_{(d,i),k} sigma_{k} ; 
  //        M^{j}_{d'' d'} <- V^\dagger_{k, (j,d')
  //Question: How to update env_ops?
  //Assume original matrix in forwardNF: (RS=Rescaling, e.g., largest SVD)
  //=>  M^{i}_{d d''} ~ RS* \tilde{U}_{(d,i),d'')} 
  //->  T removes old M^i and replaces it by U*sigma:
  //->  T = (M^{i}_{d d''})^{-1} U_{(d'',i),k} sigma_{k} 
  //->  T = (\tilde{U}/RS) U_{(d'',i),k} \sigma_{k} 
  //->  T = (M^{i})^\dagger_{d d''}/RS^2  U_{(d'',i),k} sigma_{k} 
  //->  T^{-1} = sigma_{k}^{-1} U_{(d'',i),k}^\dagger M^{i}_{d d''}
  //->  PassOn P=T, Pinv=T^{-1}

  if(!e2.is_forwardNF()){
    std::cerr<<"ProcessTensorElement::sweep_pair_backward requires e2 in forward normal form!"<<std::endl;
    throw DummyException();
  }
  if(e2.M.dim_d2!=M.dim_d1){
    std::cerr<<"ProcessTensorElement::sweep_pair_backward: e2.M.dim_d2!=M.dim_d1!"<<std::endl;
    throw DummyException();
  }


  //construct Mij
  Eigen::MatrixXcd A1=e2.M.get_Matrix_d1i_d2();
  Eigen::MatrixXcd A2=M.get_Matrix_d1_id2();
  Eigen::MatrixXcd tmp=A1*A2;

  //Extract rescaling of next e2 (used for inverse)  
  double RS=A1.col(0).norm();

  //SVD:
  TruncatedSVD_Return ret=trunc.compress(tmp); 
  clearNF();
  backwardNF=ret.sigma;
   
  double keep=backwardNF(0);
  if(trunc.keep>0)keep=trunc.keep;
  if(fabs(keep-1.)<=1e-6)keep=1.;

  //Adjust closures:
  PassOn pass_on;
  pass_on.Pinv=ret.U.adjoint()*A1/RS;
  pass_on.P=pass_on.Pinv.adjoint();
//std::cout<<"TEST2.2"<<std::endl;
  pass_on.P=pass_on.P*(ret.sigma/keep/RS).asDiagonal();
//std::cout<<"TEST2.3"<<std::endl;
  Eigen::VectorXd sinv(ret.sigma.rows());

  for(int k=0; k<sinv.rows(); k++){sinv(k)=RS*keep/ret.sigma(k);}
  pass_on.Pinv=(sinv.asDiagonal())*pass_on.Pinv;
//std::cout<<"TEST3"<<std::endl;

  //Set matrices
  e2.M.set_from_Matrix_d1i_d2(ret.U*ret.sigma.asDiagonal()/keep, e2.M.dim_i);
  M.set_from_Matrix_d1_id2(ret.Vdagger*keep, M.dim_i);
//std::cout<<"TEST4"<<std::endl;

  if(e2.closure.size()==pass_on.Pinv.cols()){
    e2.closure=pass_on.Pinv*e2.closure;
  }
//std::cout<<"TEST5"<<std::endl;
  e2.env_ops.process_backward(pass_on);
//std::cout<<"TEST6"<<std::endl;

/*
  if(M.dim_d2!=e2.M.dim_d1){
    std::cerr<<"M.dim_d2="<<M.dim_d2<<" e2.M.dim_d1="<<e2.M.dim_d1<<std::endl;
    throw DummyException();
  }
  if(closure.rows()!=M.dim_d2){
    std::cerr<<"closure.rows()="<<closure.rows()<<" != M.dim_d2="<<M.dim_d2<<std::endl;
    throw DummyException();
  }
  for(int i=0; i<env_ops.ops.size(); i++){
    if(env_ops.ops[i].rows()!=M.dim_d2){
      std::cerr<<"env_ops.ops["<<i<<"].rows()="<<env_ops.ops[i].rows()<<" != M.dim_d2="<<M.dim_d2<<std::endl;
      throw DummyException();
    }
  }
*/
}


void ProcessTensorElement::clear(){
  accessor.set_trivial(0, M);
  closure=Eigen::VectorXcd::Ones(1);
  env_ops.ops.resize(1, closure);
  clearNF();
}

void ProcessTensorElement::set_trivial(int N_sys){  
  clear();
  accessor.set_trivial(N_sys, M); 
  closure=Eigen::VectorXcd::Ones(1);
  forwardNF=Eigen::VectorXd::Ones(1);
  backwardNF=Eigen::VectorXd::Ones(1);
}

void ProcessTensorElement::set_from_ModePropagator(ModePropagator &mprop, double ta, double dt, double dict_zero){
  clear();
  int N=mprop.get_N_system();
  int NL=N*N;
  int N_mode=mprop.get_N_mode();
//  int ML=N_mode*N_mode;

  mprop.update(ta, dt);
  M=mprop.get_MPS_Matrix();
  accessor.dict.set_default(N);
 
  IF_OD_Dictionary newdict=detect_dict(dict_zero);
  reduce_to_dict(newdict);

  closure=H_Matrix_to_L_Vector(Eigen::MatrixXcd::Identity(N_mode,N_mode));
  env_ops=mprop.env_ops;
}


IF_OD_Dictionary ProcessTensorElement::detect_dict(double dict_zero)const{
  IF_OD_Dictionary dict; 
  dict.detect(M, dict_zero);
  return dict;
}
void ProcessTensorElement::reduce_to_dict(const IF_OD_Dictionary &dict){
  dict.reduce_MPS_Matrix(M);
  accessor.dict=dict;
}
void ProcessTensorElement::expand_from_dict(){
  accessor.dict.expand_MPS_Matrix(M);
  accessor.dict.set_default(accessor.dict.get_N());
}
void ProcessTensorElement::expand_DiagBB(const DiagBB &diagBB){
  accessor.expand_DiagBB(diagBB);
}
void ProcessTensorElement::expand_space_front(int N_front){
  accessor.dict.expand_space_front(N_front);
}
void ProcessTensorElement::expand_space_back(int N_back){
  accessor.dict.expand_space_back(N_back);
}
void ProcessTensorElement::apply_HilbertSpaceRotation(const HilbertSpaceRotation &hs_rot, double dict_zero){
  accessor.apply_HilbertSpaceRotation(hs_rot, M, dict_zero);
}

void ProcessTensorElement::calculate_closure(const ProcessTensorElement *last){
  closure=Eigen::VectorXcd::Zero(M.dim_d2);

  if(last==NULL){
    if(closure.rows()!=1){
      std::cerr<<"ProcessTensorElement::calculate_closure: last but M.dim_d2="<<M.dim_d2<<"!=1!"<<std::endl;
      throw DummyException();
    }
    closure(0)=1;
  }else{
    if(last->M.dim_d2!=last->closure.rows()){
      std::cerr<<"ProcessTensorElement::calculate_closure: last->M.dim_d2!=last->closure.rows()!"<<std::endl;
      throw DummyException();
    }

    std::pair<double, std::vector<int> > ind_list = 
                                           last->accessor.closure_indices();
    const MPS_Matrix &Mlast = last->M;
    for(const int &i : ind_list.second){
      for(int d1=0; d1<Mlast.dim_d1; d1++){
        for(int d2=0; d2<Mlast.dim_d2; d2++){
          closure(d1) += ind_list.first * Mlast(i, d1, d2) * last->closure(d2);
        }
      }
    } 
  }
}
void ProcessTensorElement::close_off(){
  if(M.dim_d2==1)return;

  M.inner_multiply_right(closure);
  closure=Eigen::VectorXcd::Ones(1);
  env_ops.ops=std::vector<Eigen::VectorXcd>(env_ops.ops.size(), Eigen::VectorXcd::Ones(1)*1./0.);
//  env_ops.ops=std::vector<Eigen::VectorXcd>(env_ops.ops.size(), Eigen::VectorXcd::Zero(1));

  if(forwardNF.size()>0){
    forwardNF=Eigen::VectorXd::Ones(1);
  }
  if(backwardNF.size()>0){
    backwardNF=Eigen::VectorXd::Ones(1);
  }
}

void ProcessTensorElement::print_debug(std::ostream &ofs)const{
  ofs<<"dict: "; accessor.dict.print_beta(ofs); ofs<<std::endl;
  ofs<<"closure.rows(): "<<closure.rows()<<std::endl;
  ofs<<"env_ops: "; env_ops.print_debug(); ofs<<std::endl;
  ofs<<"forwardBF.rows(): "<<forwardNF.rows()<<std::endl;
  ofs<<"backwardBF.rows(): "<<backwardNF.rows()<<std::endl;
  ofs<<"M: "; M.print_dims(ofs);ofs<<std::endl;
}

void ProcessTensorElement::read_binary(std::istream &is){
  accessor.read_binary(is); 
  closure=binary_read_EigenMatrixXcd(is,"closure");
  env_ops.read_binary(is);
  forwardNF=binary_read_EigenMatrixXd(is,"forwardNF");
  backwardNF=binary_read_EigenMatrixXd(is,"backwardNF");
  M.read_binary(is);
}

void ProcessTensorElement::write_binary(std::ostream &os)const{
  accessor.write_binary(os);
  binary_write_EigenMatrixXcd(os, closure);
  env_ops.write_binary(os);
  binary_write_EigenMatrixXd(os,forwardNF);
  binary_write_EigenMatrixXd(os,backwardNF);
  M.write_binary(os);
}

/*
void ProcessTensorElement::swap(ProcessTensorElement &e){
  accessor.swap(e.accessor);
  closure.swap(e.closure);
  env_ops.swap(e.env_ops);
  forwardNF.swap(e.forwardNF);
  backwardNF.swap(e.backwardNF);
  M.swap(e.M);
}*/

}
