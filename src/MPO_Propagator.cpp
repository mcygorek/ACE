#include "MPO_Propagator.hpp"
#include "CCTEMPO.hpp"


namespace ACE{ 

void MPO_Propagator::propagate_comb(Comb &comb, const TimeGrid &tgrid, int step, const TruncatedSVD &trunc){
  for(size_t r=0; r<fprop.size(); r++){
    if(!fprop[r]){continue;}
    if(r>=comb.size()){
      std::cerr<<"MPO_Propagator::propagate_comb: r>=comb.size()!"<<std::endl;
      throw DummyException();
    }
    if(use_symmetric_Trotter){
      fprop[r]->update(tgrid.get_t(step-1), tgrid.get_dt(step-1)/2.);
    }else{
      fprop[r]->update(tgrid.get_t(step-1), tgrid.get_dt(step-1));
    }
    comb[r].f.multiply_d3(fprop[r]->M.transpose());
  }
  for(size_t i=0; i<sys_prop_mpo.size(); i++){
    comb.join_d3_d1_compress_d2_d4(sys_prop_mpo[i], trunc, 1 );
    comb.sweep_d4_d2(trunc);
  }
}
void MPO_Propagator::propagate_comb_second(Comb &comb, const TimeGrid &tgrid, int step, const TruncatedSVD &trunc){
  if(!use_symmetric_Trotter){ return; }
  for(int i=(int)sys_prop_mpo.size()-1; i>=0; i--){
    comb.join_d3_d1_compress_d2_d4(sys_prop_mpo[i], trunc, 1 );
    comb.sweep_d4_d2(trunc);
  }
  for(size_t r=0; r<fprop.size(); r++){
    if(!fprop[r]){continue;}
    if(r>=comb.size()){
      std::cerr<<"MPO_Propagator::propagate_comb: r>=comb.size()!"<<std::endl;
      throw DummyException();
    }
    //fprop[r]->update(tgrid.get_t(step-1), tgrid.get_dt(step-1)/2.);
    comb[r].f.multiply_d3(fprop[r]->M.transpose());
  }
}

void MPO_Propagator::setup(Parameters & param){
  use_symmetric_Trotter=param.get_as_bool("use_symmetric_Trotter", true);

  //CCTEMPO::add_HilbertSpace_Propagator_CCTEMPO(sys_prop_mpo, param, t);

  std::vector<int> dims=CCTEMPO::get_dimensions(param);
  int N_sites=dims.size();

  TimeGrid tgrid(param);
  TruncatedSVD trunc=CCTEMPO::get_SysPropMPO_trunc(param);

  //add local propagators
  fprop.clear(); fprop.resize(N_sites);
  for(int r=0; r<N_sites; r++){
    std::string key="S"+int_to_string(r);
    Parameters param2; param2.add_from_prefix(key, param);
    fprop[r]=std::make_shared<FreePropagator>(param2);
  }

  //add collective terms to Hamiltonian
  bool has_collective=false;
  MPO H_mpo=MPO_op::Zero(dims);
  //from single-excitation Manifold
  Eigen::MatrixXcd SXH = CCTEMPO::get_SX_Hamiltonian(param, tgrid.ta);
  if(SXH.norm()>1e-14){
    has_collective=true;
    for(int i=0; i<(int)dims.size(); i++){
      if(dims[i]<2)continue;
      MPO_op::add(H_mpo, MPO_op::OneBody(dims, i, SXH(i,i)*Operators(dims[i]).ketbra(1,1)));
      for(int j=0; j<(int)dims.size(); j++){
        if(dims[j]<2 || i==j)continue;
        MPO_op::add(H_mpo, MPO_op::TwoBody(dims, i, j, sqrt(SXH(i,j))*Operators(dims[i]).ketbra(1,0), sqrt(SXH(i,j))*Operators(dims[j]).ketbra(0,1)));
      }
    }
  }

  //explicit two-body terms
  for(int r=0; r<N_sites; r++){
    for(int r2=0; r2<N_sites; r2++){
      if(r==r2)continue;
      std::string key="S"+int_to_string(r)+"S"+int_to_string(r2);
      Parameters param2; param2.add_from_prefix(key, param);
      FreePropagator prop(param2);
      if(prop.const_H.norm()>1e-14){
        has_collective=true;
        if(prop.const_H.rows()!=dims[r]*dims[r2]){
          std::cerr<<"MPO_Propagator::setup: "<<key<<" wrong dimensions!"<<std::endl;
          throw DummyException();
        } 
        Eigen::MatrixXcd A(dims[r]*dims[r], dims[r2]*dims[r2]);
        for(int i1=0; i1<dims[r]; i1++){ for(int i2=0; i2<dims[r]; i2++){
            for(int j1=0; j1<dims[r2]; j1++){ for(int j2=0; j2<dims[r2]; j2++){
                A(i1*dims[r]+i2, j1*dims[r2]+j2) = 
                 prop.const_H(i1*dims[r2]+j1,i2*dims[r2]+j2);
        } } } }
        Eigen::JacobiSVD<Eigen::MatrixXcd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
        for(int i=0; i<svd.singularValues().rows(); i++){
          if(svd.singularValues()(i)<1e-15){ break; }
          Eigen::MatrixXcd A=L_Vector_to_H_Matrix(svd.matrixU().col(i));
          Eigen::MatrixXcd B=L_Vector_to_H_Matrix(svd.matrixV().col(i).conjugate());
std::cout<<"TwoBody: r="<<r<<" r2="<<r2<<" SV="<<svd.singularValues()(i)<<std::endl;
std::cout<<"A="<<A<<std::endl<<"B="<<B<<std::endl;
          MPO_op::add(H_mpo, MPO_op::TwoBody(dims, r, r2, 
            sqrt(svd.singularValues()(i))*A, sqrt(svd.singularValues()(i))*B));
        }
      }
    }
  }

  //if collective terms are given, 
  if(has_collective){
    //subtract local Hamiltonians and process with collective terms
    for(int r=0; r<N_sites; r++){
      if(fprop[r] && fprop[r]->const_H.norm()>1e-15){
        MPO_op::add(H_mpo, MPO_op::OneBody(dims, r, fprop[r]->const_H));
        fprop[r]->const_H = Eigen::MatrixXcd::Zero(dims[r],dims[r]);
      }
    }

    //exponentiate
    std::complex<double> dtfac = use_symmetric_Trotter
      ? std::complex<double>(0,-tgrid.dt/2./hbar_in_meV_ps)
      : std::complex<double>(0,-tgrid.dt/hbar_in_meV_ps);
    int N_Taylor=param.get_as_size_t("exp_Taylor",4);

    std::cout<<"H_mpo before: "; H_mpo.print_dims();
    trunc.print_info();std::cout<<std::endl;
    H_mpo.sweep_d2_d4(trunc);
    H_mpo.sweep_d4_d2(trunc);
    std::cout<<"H_mpo after: "; H_mpo.print_dims();

    MPO Hprop = MPO_op::exp_Taylor(H_mpo, dtfac, trunc, N_Taylor);
    std::cout<<"Hprop.print_dims(): "; Hprop.print_dims(); std::cout<<std::endl;
    sys_prop_mpo.push_back(MPO_op::H_to_L_forward(Hprop));
    sys_prop_mpo.push_back(MPO_op::H_to_L_backward(Hprop));
  }

  //remove local propagators if they don't do anything
  for(int r=0; r<N_sites; r++){
    if(fprop[r] && fprop[r]->does_nothing()){
      fprop[r].reset();
    }
  }
  for(int r=0; r<N_sites; r++){
    if(fprop[r]){ std::cout<<"fprop["<<r<<"]=true"<<std::endl; }
    else{        std::cout<<"fprop["<<r<<"]=false"<<std::endl; }
  }

  //add collective Lindbladians
  CCTEMPO::add_collective_decay(sys_prop_mpo, param);
  
  for(size_t i=0; i<sys_prop_mpo.size(); i++){
    std::cout<<"sys_prop_mpo["<<i<<"].print_dims(): "; 
    sys_prop_mpo[i].print_dims(); 
    std::cout<<std::endl;
  }

  //Read MPO_Apply_Operators
  apply_ops.setup(param);
}

}
