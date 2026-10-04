#include "FreePropagator.hpp"
#include "TimeGrid.hpp"
#include "FennaMatthewsOlson.hpp"
#include "MPO_op.hpp"
#include "MPO_state.hpp"
#include "Operators.hpp"

namespace ACE{
namespace CCTEMPO{


std::vector<int> get_dimensions(Parameters &param){
  int N_sites = param.get_as_size_t("N_sites");
  std::vector<int> list(N_sites, 2);
  if(param.is_specified("dimensions")){
    list=param.get_row_ints("dimensions");
    if(param.is_specified("N_sites") && list.size()!=N_sites){
      std::cerr<<"Mismatch: Size of 'dimensions' vs. 'N_sites'!"<<std::endl;
      throw DummyException();
    }
  }
//TODO: set my S?_initial 
  return list;
}
std::vector<int> get_dimensions_squared(Parameters &param){
  std::vector<int> list=get_dimensions(param);
  for(size_t i=0; i<list.size(); i++){
    list[i]=list[i]*list[i];
  }
  return list;
}

int get_N_sites(Parameters &param){
  if(param.is_specified("N_sites")){
    return (int)param.get_as_size_t("N_sites");
  }
  std::vector<int> dims = get_dimensions(param);
  if(dims.size()<1){
    std::cerr<<"get_N_sites: Need positive N_sites or dimensions!"<<std::endl;
    throw DummyException();
  }
  return dims.size();
}

TruncatedSVD get_SysProp_trunc(Parameters &param){
  TruncationLayout trunc_layout;
  Parameters param_bb;
  param_bb.add_from_prefix("SysProp", param);
  param_bb.add(param);
  trunc_layout.setup(param_bb); 
  return trunc_layout.get_forward(-1,-1);
}
TruncatedSVD get_SysPropMPO_trunc(Parameters &param){
  TruncationLayout trunc_layout;
  Parameters param_bb;
  param_bb.add_from_prefix("SysPropMPO", param);
  param_bb.add_from_prefix("SysProp", param);
  param_bb.add(param);
  trunc_layout.setup(param_bb); 
  return trunc_layout.get_forward(-1,-1);
}

FreePropagator get_SX_FreePropagator(Parameters &param){
  int N_sites=get_N_sites(param);

  FreePropagator fprop(param);
  if(fprop.get_dim()>1){
    if(fprop.get_dim()!=N_sites){
      std::cerr<<"single_excitation_manifold requires Hamiltonian to have dimension equals N_sites!"<<std::endl;
      throw DummyException();
    }else{
      fprop.set_dim(N_sites);
    }
  }
  if(param.get_as_bool("add_Hamiltonian_FMO_monomer")){
    fprop.add_Hamiltonian( H_FMO_monomer().block(0, 0, N_sites, N_sites) );
  }
  if(param.get_as_bool("add_Hamiltonian_FMO_trimer")){
    fprop.add_Hamiltonian( H_FMO_trimer().block(0, 0, N_sites, N_sites) );
  }
  return fprop;    
}
Eigen::MatrixXcd get_SX_Hamiltonian(Parameters &param, double t){
  return get_SX_FreePropagator(param).get_Htot(t);
}
Eigen::MatrixXcd get_SX_Hamiltonian_from_individual_terms(Parameters &param){
  std::vector<int> dims=CCTEMPO::get_dimensions(param);
  int N_sites=dims.size();
  Eigen::MatrixXcd H = Eigen::MatrixXcd::Zero(N_sites,N_sites);
  for(int r=0; r<N_sites; r++){
    std::string key="S"+int_to_string(r);
    Parameters param2; param2.add_from_prefix(key, param);
    FreePropagator fprop(param2);
    if(fprop.const_H.rows()>1){ H(r,r) = fprop.const_H(1,1); }
  }
  for(int r=0; r<N_sites; r++){
    for(int r2=0; r2<N_sites; r2++){
      if(r==r2)continue;
      std::string key="S"+int_to_string(r)+"S"+int_to_string(r2);
      Parameters param2; param2.add_from_prefix(key, param);
      FreePropagator prop(param2);
      if(prop.const_H.rows()==dims[r]*dims[r2]){
        H(r,r2)=prop.const_H(1*dims[r2]+0, 0*dims[r2]+1);
        H(r2,r)=prop.const_H(0*dims[r2]+1, 1*dims[r2]+0);
      }
    }
  }
  return H;
}


MPO get_initial_state(Parameters &param){
  std::cout<<"Setting up MPO initial state"<<std::endl;
  std::vector<int> dims=get_dimensions(param);
  std::vector<int> dims2=get_dimensions_squared(param);
  int N_sites = dims2.size();

  MPO state(N_sites);
  //Truncation for backbone: Let parameters "backbone_.." preceed the usual ones
  TruncationLayout trunc_layout;
  { Parameters param_bb;
    param_bb.add_from_prefix("backbone", param);
    param_bb.add(param);
    trunc_layout.setup(param_bb); }
  TruncatedSVD trunc=trunc_layout.get_forward(-1,-1);


  std::string print_initial_eigenvectors=param.get_as_string("print_initial_eigenvectors","");
  if(print_initial_eigenvectors!=""){ 
    TimeGrid tgrid(param);
    Eigen::MatrixXcd H=get_SX_Hamiltonian(param, tgrid.ta);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(H);;
    std::ofstream ofs(print_initial_eigenvectors);
    for(int i=0; i<H.rows(); i++){
      ofs<<solver.eigenvalues()(i);
      for(int j=0; j<H.rows(); j++){
        ofs<<" "<<solver.eigenvectors()(j,i).real();
        ofs<<" "<<solver.eigenvectors()(j,i).imag();
      }
      ofs<<std::endl;
    }
  }
  int initial_SX_eigenstate=param.get_as_int("initial_SX_eigenstate",-1);
  if(initial_SX_eigenstate>=0){
    if(initial_SX_eigenstate>=N_sites){
      std::cerr<<"CCTEMPO::get_initial_state: initial_SX_eigenstate>=N_sites"<<std::endl;
      throw DummyException();
    }
    TimeGrid tgrid(param);
    Eigen::MatrixXcd H=get_SX_Hamiltonian(param, tgrid.ta);
    { Eigen::MatrixXcd H_indiv=CCTEMPO::get_SX_Hamiltonian_from_individual_terms(param);
      if(H_indiv.norm()>1e-14){
        if(H.rows()==0){H=H_indiv;}
        else if(H.rows()!=H_indiv.rows()){
          std::cerr<<"initial_SX_eigenstate: H.rows()!=H_indiv.rows()!"<<std::endl; 
          throw DummyException();
        }
        H+=H_indiv;
      }
    }
    if(H.rows()==0){
      std::cerr<<"initial_SX_eigenstate: H.rows()==0!"<<std::endl; 
      throw DummyException();
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(H);
    std::cout<<"initial single-excitation manifold Hamiltonian eigenvalues: "<<solver.eigenvalues().transpose()<<std::endl;
 
    Eigen::VectorXcd phi=solver.eigenvectors().col(initial_SX_eigenstate);
    Eigen::MatrixXcd rho=phi*phi.adjoint();
    state=MPO_state::from_SingleExcitationManifold(rho);
    
    std::cout<<"initial state.print_dims(): "; state.print_dims(); 
    state.sweep_d2_d4(trunc);
    state.sweep_d4_d2(trunc);
    std::cout<<"initial state.print_dims(): "; state.print_dims(); std::cout<<std::endl;
    std::cout<<"dot: "<<MPO_state::dot(state, state)<<std::endl;

  }else{
    for(int r=0; r<N_sites; r++){
      std::string key="S"+int_to_string(r)+"_initial";
      Eigen::MatrixXcd default_init=Eigen::MatrixXcd::Zero(dims[r],dims[r]); 
      default_init(0,0)=1;
      Eigen::MatrixXcd rho_initial=param.get_as_operator(key,default_init);
      std::cout<<key<<": "<<std::endl<<rho_initial<<std::endl;
      int N=rho_initial.rows(); 
      if(N!=dims[r]){
        std::cerr<<"CCTEMPO::get_initial_state: N!=dims[r]"<<std::endl;
        throw DummyException();
      }
      state[r]=FourLeg(1,1,N*N,1);
      for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
          state[r](0,0,i*N+j,0)=rho_initial(i,j);
        }
      } 
    }  
  }

  int N_coeffs = param.get_nr_rows("add_initial_coeff");
  std::vector<int> ivec_zeros(state.size(),0);
  for(int i=0; i<N_coeffs; i++){
    std::vector<double> dvec = param.get_row_doubles("add_initial_coeff", i);
    if(dvec.size()!=N_sites+2){
      std::cerr<<"initial_state_MPO: add_initial_coeffs needs N_sites+2 parameters: Index list + real + imaginary value!"<<std::endl;
      throw DummyException();
    }
    std::vector<int> ivec(N_sites);
    for(int r=0; r<N_sites; r++){
      ivec[r]=dvec[r];
    }
    std::complex<double> value(dvec[N_sites], dvec[N_sites+1]);

    std::cout<<"dvec:"; for(int i=0; i<N_sites; i++){std::cout<<" "<<dvec[i];}
    std::cout<<std::endl; 
    std::cout<<"adding initial coefficients:";
    for(int i=0; i<N_sites; i++){std::cout<<" "<<ivec[i];} 
    std::cout<<": "<<value<<std::endl;

    MPO_op::add_coeff(state, ivec_zeros, ivec, value);
    state.sweep_d2_d4(trunc);
    state.sweep_d4_d2(trunc);
  }

  return state;
}

 

void add_HilbertSpace_Propagator_CCTEMPO(std::vector<MPO> & sys_prop, Parameters &param, double t){
  MPO Hprop;
  std::cout<<"Setting up HilbertSpace_Propagator"<<std::endl;

  std::vector<int> dims=CCTEMPO::get_dimensions(param);
  std::vector<int> dims2=CCTEMPO::get_dimensions_squared(param);
  int N_sites = dims2.size();

  bool use_symmetric_Trotter=param.get_as_bool("use_symmetric_Trotter", true);
  TimeGrid tgrid(param);
  bool do_something=false;
 
  TruncatedSVD trunc=CCTEMPO::get_SysPropMPO_trunc(param);
  
//  bool SX_Hamiltonian = param.get_as_bool("SX_Hamiltonian", false);
//  if(SX_Hamiltonian) //build N-level H, then make MPO
    Eigen::MatrixXcd SXH = get_SX_Hamiltonian(param, t);
    std::cout<<"Single-excitation manifold Hamiltonian: "<<std::endl<<SXH<<std::endl;
    if(SXH.norm()>1e-14){
       do_something=true;
       Eigen::MatrixXcd SXM = use_symmetric_Trotter 
            ? (SXH*std::complex<double>(0,-tgrid.dt/2./hbar_in_meV_ps)).exp()
            : (SXH*std::complex<double>(0,-tgrid.dt/hbar_in_meV_ps)).exp();

     if(false){
      Hprop = MPO_op::SingleExcitationManifold_to_MPO_H(SXM, trunc, dims);
     }else if(false){ //distribute exp
      Hprop = MPO_op::Identity(dims);
      for(int i=0; i<(int)dims.size(); i++){
//std::cout<<"i="<<i<<std::endl;
        if(dims[i]<2)continue;
        MPO_op::add(Hprop, MPO_op::OneBody(dims, i, (SXM(i,i)-1.)*Operators(dims[i]).ketbra(1,1)));

        for(int j=0; j<(int)dims.size(); j++){
//std::cout<<"j="<<j<<std::endl;
          if(dims[j]<2 || i==j)continue;
          MPO_op::add(Hprop, MPO_op::TwoBody(dims, i, j, sqrt(SXM(i,j))*Operators(dims[i]).ketbra(1,0), sqrt(SXM(i,j))*Operators(dims[j]).ketbra(0,1)));
        }
      }
/*      std::cout<<"test"<<std::endl;
      std::cout<<"0="<<0<<" <-> "<<MPO_op::evaluate(Hprop,{0,0,0,0,0,0,0},{0,0,0,0,0,0,0})<<std::endl;
      std::cout<<"SXM(0,0)="<<SXM(0,0)<<" <-> "<<MPO_op::evaluate(Hprop,{1,0,0,0,0,0,0},{1,0,0,0,0,0,0})<<std::endl;
      std::cout<<"SXM(0,1)="<<SXM(0,1)<<" <-> "<<MPO_op::evaluate(Hprop,{1,0,0,0,0,0,0},{0,1,0,0,0,0,0})<<std::endl;
      std::cout<<"SXM(1,0)="<<SXM(1,0)<<" <-> "<<MPO_op::evaluate(Hprop,{0,1,0,0,0,0,0},{1,0,0,0,0,0,0})<<std::endl;
      std::cout<<"SXM(0,6)="<<SXM(0,6)<<" <-> "<<MPO_op::evaluate(Hprop,{1,0,0,0,0,0,0},{0,0,0,0,0,0,1})<<std::endl;
      Hprop.print_dims(); */
      trunc.print_info();std::cout<<std::endl;
      Hprop.sweep_d2_d4(trunc);
      Hprop.sweep_d4_d2(trunc);
      Hprop.print_dims(); 
//      std::cout<<"SXM(0,6)="<<SXM(0,6)<<" <-> "<<MPO_op::evaluate(Hprop,{1,0,0,0,0,0,0},{0,0,0,0,0,0,1})<<std::endl;
     }else{
      MPO Hop = MPO_op::Zero(dims);
      for(int i=0; i<(int)dims.size(); i++){
//std::cout<<"i="<<i<<std::endl;
        if(dims[i]<2)continue;
        MPO_op::add(Hop, MPO_op::OneBody(dims, i, SXH(i,i)*Operators(dims[i]).ketbra(1,1)));

        for(int j=0; j<(int)dims.size(); j++){
//std::cout<<"j="<<j<<std::endl;
          if(dims[j]<2 || i==j)continue;
          MPO_op::add(Hop, MPO_op::TwoBody(dims, i, j, sqrt(SXH(i,j))*Operators(dims[i]).ketbra(1,0), sqrt(SXH(i,j))*Operators(dims[j]).ketbra(0,1)));
        }
      }
/*      std::cout<<"test"<<std::endl;
      std::cout<<"0="<<0<<" <-> "<<MPO_op::evaluate(Hop,{0,0,0,0,0,0,0},{0,0,0,0,0,0,0})<<std::endl;
      std::cout<<"SXH(0,0)="<<SXH(0,0)<<" <-> "<<MPO_op::evaluate(Hop,{1,0,0,0,0,0,0},{1,0,0,0,0,0,0})<<std::endl;
      std::cout<<"SXH(0,1)="<<SXH(0,1)<<" <-> "<<MPO_op::evaluate(Hop,{1,0,0,0,0,0,0},{0,1,0,0,0,0,0})<<std::endl;
      std::cout<<"SXH(1,0)="<<SXH(1,0)<<" <-> "<<MPO_op::evaluate(Hop,{0,1,0,0,0,0,0},{1,0,0,0,0,0,0})<<std::endl;
      std::cout<<"SXH(0,6)="<<SXH(0,6)<<" <-> "<<MPO_op::evaluate(Hop,{1,0,0,0,0,0,0},{0,0,0,0,0,0,1})<<std::endl;*/
      Hop.print_dims();
      trunc.print_info();std::cout<<std::endl;
      Hop.sweep_d2_d4(trunc);
      Hop.sweep_d4_d2(trunc);
      std::cout<<"Hop: ";
      Hop.print_dims(); 
      std::complex<double> dtfac = use_symmetric_Trotter 
            ? std::complex<double>(0,-tgrid.dt/2./hbar_in_meV_ps)
            : std::complex<double>(0,-tgrid.dt/hbar_in_meV_ps);
      int N_Taylor=4;
      Hprop = MPO_op::exp_Taylor(Hop, dtfac, trunc, N_Taylor);
      std::cout<<"Hprop: ";
      Hprop.print_dims(); 
     }

      std::cout<<"Single-excitation manifold Hamiltonian: Norm of difference: "<<MPO_op::compare_SingeExcitationManifold_MPO(Hprop, SXM)<<std::endl;
    }else{
      Hprop = MPO::Identity_d1_d3(dims);
    }

  //add single-body terms
  for(int r=0; r<N_sites; r++){
    std::string key="S"+int_to_string(r);
    Parameters param2; param2.add_from_prefix(key, param);
//    std::cout<<"Single-Body Parameters on site "<<r<<std::endl;
//    param2.print();
    FreePropagator fprop(param2);
    Eigen::MatrixXcd H = fprop.get_Htot(t);

    if(H.norm()>1e-14){
      do_something=true;
      std::cout<<"Single-body Hamiltonian at site "<<r<<":"<<std::endl;
      std::cout<<H<<std::endl;

      Eigen::MatrixXcd M = use_symmetric_Trotter 
        ? (H*std::complex<double>(0,-tgrid.dt/4./hbar_in_meV_ps)).exp()
        : (H*std::complex<double>(0,-tgrid.dt/2./hbar_in_meV_ps)).exp();
      Hprop[r].multiply_d1(M);
      Hprop[r].multiply_d3(M.transpose());
    }
  }

  if(do_something){
    std::cout<<"Hprop.print_dims(): "; Hprop.print_dims(); std::cout<<std::endl;
    sys_prop.push_back(MPO_op::H_to_L_forward(Hprop));
    sys_prop.push_back(MPO_op::H_to_L_backward(Hprop));
  }
}

void add_collective_decay(std::vector<MPO> & sys_prop, Parameters & param){
    double rate=param.get_as_double("add_collective_decay",0.);
    int N_Taylor=param.get_as_size_t("add_collective_decay",2,0,1);
    if(rate>0.){
      std::vector<int> dims2=CCTEMPO::get_dimensions_squared(param);
      int N_sites = dims2.size();
      TimeGrid tgrid(param);
      TruncatedSVD trunc = get_SysPropMPO_trunc(param);

      MPO L = MPO_op::Zero(dims2);
      Eigen::MatrixXcd sigma_m=Eigen::MatrixXcd::Zero(2,2);
      sigma_m(0,1)=1;
      Eigen::MatrixXcd sigma_p=Eigen::MatrixXcd::Zero(2,2);
      sigma_p(1,0)=1;
      
      for(int i=0; i<N_sites; i++){ //OneBody terms
        MPO tmp=MPO_op::OneBody_forward_backward(dims2, i, sigma_m, sigma_p);
        tmp.scale_first(rate);
        MPO_op::add(L, tmp);

        tmp=MPO_op::OneBody_forward(dims2, i, sigma_p*sigma_m);
        tmp.scale_first(-0.5*rate);
        MPO_op::add(L, tmp);

        tmp=MPO_op::OneBody_backward(dims2, i, sigma_p*sigma_m);
        tmp.scale_first(-0.5*rate);
        MPO_op::add(L, tmp);
      }
      for(int i=0; i<N_sites; i++){ //TwoBody terms
        for(int j=0; j<N_sites; j++){
          if(j==i)continue;
          MPO tmp=MPO_op::TwoBody_forward_backward(dims2, i, j, sigma_m, sigma_p);
          tmp.scale_first(rate);
          MPO_op::add(L, tmp);

          tmp=MPO_op::TwoBody_forward(dims2, i, j, sigma_p, sigma_m);
          tmp.scale_first(-0.5*rate);
          MPO_op::add(L, tmp);

          tmp=MPO_op::TwoBody_backward(dims2, i, j, sigma_p, sigma_m);
          tmp.scale_first(-0.5*rate);
          MPO_op::add(L, tmp);
        }
      }
      double dt=param.get_as_bool("use_symmetric_Trotter", true) ? tgrid.dt/2 : tgrid.dt;
std::cout<<"before sweep: "; L.print_dims();
      L.sweep_d2_d4(trunc);
      L.sweep_d4_d2(trunc);
std::cout<<"before exp_Taylor: "; L.print_dims();
      sys_prop.push_back(MPO_op::exp_Taylor(L, dt, trunc, N_Taylor));
std::cout<<"after exp_Taylor: "; sys_prop.back().print_dims();
    }
}

std::vector<MPO> get_sys_prop_mpo(Parameters &param, double t){
//TODO: propagation modes:
//      - Trotter: decompose propagator as a products of exp.
//      - propagate_Taylor: decompose Liouvillian as a sum of terms; then exponentiate using Taylor series. 
//      - Taylor is not necessarily linked to Hilbert/Liouville space. On-site Lindbladians are cheap. The point is that the individual terms shall have small bond dimensions.

  std::vector<MPO> sys_prop;
  CCTEMPO::add_HilbertSpace_Propagator_CCTEMPO(sys_prop, param, t);
  CCTEMPO::add_collective_decay(sys_prop, param);
  
  for(size_t i=0; i<sys_prop.size(); i++){
    std::cout<<"sys_prop_mpo["<<i<<"].print_dims(): "; 
    sys_prop[i].print_dims(); 
    std::cout<<std::endl;
  }

  return sys_prop;
}


}//namespace CCTEMPO
}//namespace ACE

