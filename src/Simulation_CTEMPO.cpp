#include "Simulation_CTEMPO.hpp"
#include "DummyException.hpp"
#include "Timings.hpp"
#include "Coupling_Groups.hpp"
#include "otimes.hpp"

using namespace ACE;


void Simulation_CTEMPO::initialize_PTB(ProcessTensorBuffer &PTB, int n_mem, const DiagBB & diagBB, const Eigen::VectorXcd & rho_reduced){

    int N = sqrt(rho_reduced.rows());
    int NL=N*N;
    int NdiagBB=diagBB.get_dim();
    int NLdiagBB=NdiagBB*NdiagBB;
    std::cout<<"diagBB.sys_dim()="<<diagBB.sys_dim()<<" ";
    std::cout<<"diagBB.get_dim()="<<diagBB.get_dim()<<std::endl;
    if(diagBB.sys_dim() != N){
      std::cerr<<"Simulation_CTEMPO::initialize_PTB: diagBB.sys_dim() != N ("<<diagBB.sys_dim()<<" vs. "<<N<<")!"<<std::endl;
      throw DummyException();
    }
 
    if(n_mem<1)n_mem=1;

//    ProcessTensorBuffer PTB; 
    //PTB.set_trivial(n_mem, N);
    {
      // work with dimension diagBB.get_dim(). 
      // Only at the end, expand using dict.expand_DiagBB(diagBB);
      // Expansion may only be needed just before merging with sys. prop.      
      ProcessTensorElement e;
      e.set_trivial(NdiagBB);
      e.accessor.dict.set_default_diag(NdiagBB); 
      e.accessor.dict.expand_DiagBB(diagBB); 
      e.M.set_from_Matrix_d1i_d2(Eigen::VectorXcd::Ones(NLdiagBB), NLdiagBB);
      e.env_ops.clear();
      PTB.resize(n_mem, e); 
    }
    { // set initial state
      ProcessTensorElement & e = PTB.get(0);
      e.accessor.dict.set_default_diag(N); // full set of outer bonds
      e.M.set_from_Matrix_d1i_d2(rho_reduced, NL); 
    }
//    return PTB;
}
std::vector<Eigen::MatrixXcd> Simulation_CTEMPO::initialize_b(int n_mem, DiagBB &diagBB, double dt){
  std::vector<Eigen::MatrixXcd> b(n_mem);
  for(int n=0; n<n_mem; n++){
    b[n]=diagBB.calculate_expS(n, dt);
  }
  return b;
}


void Simulation_CTEMPO::single_step(ProcessTensorBuffer & PTB, Propagator &prop,const std::vector<Eigen::MatrixXcd> &b, const TimeGrid &tgrid, int n, const TruncationLayout & trunc_layout)const{

      if(PTB.n_tot<2){
        std::cerr<<"Simulation_CTEMPO::single_step: PTB.n_tot<2!"<<std::endl;
        throw DummyException();
      }
      int n_mem=PTB.n_tot;
      if(b.size()<n_mem){
        std::cerr<<"Simulation_CTEMPO::single_step: b.size()<n_mem!"<<std::endl;
        throw DummyException();
      }
      
      int N; 
      IF_OD_Dictionary e1dict;
      {
        const ProcessTensorElement & e0 = PTB.get(0, ForwardPreload);
        N=e0.accessor.dict.get_N();
        const ProcessTensorElement & e1 = PTB.peek(1);
        e1dict=e1.accessor.dict;
      }
      int NL=N*N;
      int Ngrps2=e1dict.reduced_dim;

      int cur_length=n_mem;
      int new_length=n_mem;   //length of new MPO without PT
      if(n+n_mem > tgrid.n_tot){ 
        cur_length = tgrid.n_tot-n+1; 
        new_length = tgrid.n_tot-n; 
      }
      //std::cout<<"n="<<n<<" n_tot="<<tgrid.n_tot<<" cur_length="<<cur_length<<" new_length="<<new_length<<std::endl;

      //TruncatedSVD trunc = trunc_layout.get_base();
      TruncatedSVD trunc_fwd=trunc_layout.get_forward(n, tgrid.n_tot);
      trunc_fwd.keep = sqrt(Ngrps2);//diagBB.get_dim();
      std::cout<<"Forward sweep: ";trunc_fwd.print_info(); std::cout<<std::endl;


      using RowMatrixMap = Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >;
      using OuterStrideMap = Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 0, Eigen::OuterStride<Eigen::Dynamic> >;
      using InnerStrideMap = Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 0, Eigen::InnerStride<Eigen::Dynamic> >;
      using DynamicStrideMap = Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> >;

      // propagate first element:
      if(no_system_propagation){

        ProcessTensorElement & e0 = PTB.get(0, ForwardPreload);
        const ProcessTensorElement & e1 = PTB.peek(1);
        ProcessTensorElement Q0;
        Q0.set_trivial(N);
        Q0.accessor.dict.set_default_diag(N);
        Q0.closure = e1.closure;
        Q0.env_ops.clear(); // = e1.env_ops;
/*
        Q0.M.resize(NL, e0.M.dim_d1, e1.M.dim_d2*Ngrps2);
        Q0.M.set_zero();
        for(int a1=0; a1<NL; a1++){ 
          int g1=e1dict.beta[a1*NL+a1];//int g1=grp2[a1];
          if(g1<0)continue;
          for(int d1=0; d1<e1.M.dim_d2; d1++){
            for(int d_=0; d_<e0.M.dim_d1; d_++){
              for(int d0=0; d0<e1.M.dim_d1; d0++){
                Q0.M(a1, d_, d1*Ngrps2+g1) += b[0](g1,g1) * e1.M(g1, d0, d1) * e0.M(a1, d_, d0);
              }
            }
          }
        }      
        Q0.closure = Vector_otimes(e1.closure, Eigen::VectorXcd::Ones(Ngrps2));
        e0.swap(Q0);
*/
        Q0.M.resize(NL, e0.M.dim_d1, e1.M.dim_d2*Ngrps2);
        Q0.M.set_zero();
        for(int a1=0; a1<NL; a1++){ 
          int g1=e1dict.beta[a1*NL+a1]; if(g1<0)continue;

          DynamicStrideMap(Q0.M.mem+a1*Q0.M.dim_d2+g1, Q0.M.dim_d1, e1.M.dim_d2, Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>(Q0.M.dim_i*Q0.M.dim_d2, Ngrps2))= 
            b[0](g1,g1) *
            OuterStrideMap(e0.M.mem+a1*e0.M.dim_d2, e0.M.dim_d1, e0.M.dim_d2, Eigen::OuterStride<Eigen::Dynamic>(e0.M.dim_i*e0.M.dim_d2)) *
            OuterStrideMap(e1.M.mem+g1*e1.M.dim_d2, e1.M.dim_d1, e1.M.dim_d2, Eigen::OuterStride<Eigen::Dynamic>(e1.M.dim_i*e1.M.dim_d2));
        }
        e0.swap(Q0);
        e0.closure = Vector_otimes(e1.closure, Eigen::VectorXcd::Ones(Ngrps2));

/*
        Q0.M.resize(NL, e0.M.dim_d1, e1.M.dim_d2);
        Q0.M.set_zero();
        for(int a1=0; a1<NL; a1++){ 
          int g1=e1dict.beta[a1*NL+a1]; if(g1<0)continue;
          OuterStrideMap(Q0.M.mem+a1*Q0.M.dim_d2, Q0.M.dim_d1, Q0.M.dim_d2, Eigen::OuterStride<Eigen::Dynamic>(Q0.M.dim_i*Q0.M.dim_d2) ) = 
            OuterStrideMap(e0.M.mem+a1*e0.M.dim_d2, e0.M.dim_d1, e0.M.dim_d2, Eigen::OuterStride<Eigen::Dynamic>(e0.M.dim_i*e0.M.dim_d2)) *
            OuterStrideMap(e1.M.mem+g1*e1.M.dim_d2, e1.M.dim_d1, e1.M.dim_d2, Eigen::OuterStride<Eigen::Dynamic>(e1.M.dim_i*e1.M.dim_d2));
  
        }      
  
        e0.accessor.dict.set_default_diag(N);
        e0.env_ops.clear(); 
        e0.M.resize(Q0.M.dim_i, Q0.M.dim_d1, Q0.M.dim_d2*Ngrps2);
        e0.M.set_zero();
        for(int a1=0; a1<NL; a1++){ 
          int g1=e1dict.beta[a1*NL+a1]; if(g1<0)continue;

          DynamicStrideMap(e0.M.mem+a1*e0.M.dim_d2+g1, e0.M.dim_d1, e1.M.dim_d2, Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>(e0.M.dim_i*e0.M.dim_d2,Ngrps2)) = 
            b[0](g1,g1) * 
            OuterStrideMap(Q0.M.mem+a1*Q0.M.dim_d2, Q0.M.dim_d1, Q0.M.dim_d2, Eigen::OuterStride<Eigen::Dynamic>(Q0.M.dim_i*Q0.M.dim_d2));
        }
        e0.closure = Vector_otimes(e1.closure, Eigen::VectorXcd::Ones(Ngrps2));
*/
     }else{
        ProcessTensorElement & e0 = PTB.get(0, ForwardPreload);
        const ProcessTensorElement & e1 = PTB.peek(1);
        ProcessTensorElement Q0;
        Q0.set_trivial(N);
        Q0.accessor.dict.set_default_diag(N);
        Q0.M.resize(NL, e0.M.dim_d1, e0.M.dim_d2);
        Q0.closure = e0.closure;
        Q0.env_ops.clear();// = e0.env_ops;
        Q0.M.set_zero();
      
        prop.update(tgrid.get_t(n), tgrid.dt/2);
//        for(int d0=0; d0<e0.M.dim_d2; d0++){
//          for(int a1=0; a1<NL; a1++){ 
//            for(int a0=0; a0<NL; a0++){ 
//              for(int d1=0; d1<e0.M.dim_d1; d1++){
//                Q0.M(a1, d1, d0) += prop.M(a1, a0) * e0.M(a0, d1, d0);
//        } } } }  
        for(int d1=0; d1<e0.M.dim_d1; d1++){
          RowMatrixMap(Q0.M.mem+d1*NL*e0.M.dim_d2, NL, e0.M.dim_d2) = 
            prop.M * RowMatrixMap(e0.M.mem+d1*NL*e0.M.dim_d2, NL, e0.M.dim_d2);
  
        }
      
        ProcessTensorElement Qtmp;
        Qtmp.set_trivial(N);
        Qtmp.accessor.dict.set_default_diag(N);
        Qtmp.M.resize(NL, e0.M.dim_d1, e1.M.dim_d2);
        Qtmp.closure = e1.closure;
        Qtmp.env_ops.clear(); // = e1.env_ops;
        Qtmp.M.set_zero();
        for(int a1=0; a1<NL; a1++){ 
          int g1=e1dict.beta[a1*NL+a1];//int g1=grp2[a1];
          if(g1<0)continue;
//          for(int d1=0; d1<e1.M.dim_d2; d1++){
//            for(int d_=0; d_<e0.M.dim_d1; d_++){
//              for(int d0=0; d0<e1.M.dim_d1; d0++){
//                Qtmp.M(a1, d_, d1) += e1.M(g1, d0, d1) * Q0.M(a1, d_, d0);
//          } } }  
          OuterStrideMap(Qtmp.M.mem+a1*Qtmp.M.dim_d2, Qtmp.M.dim_d1, Qtmp.M.dim_d2, Eigen::OuterStride<Eigen::Dynamic>(Qtmp.M.dim_i*Qtmp.M.dim_d2) ) = 
            OuterStrideMap(Q0.M.mem+a1*Q0.M.dim_d2, Q0.M.dim_d1, Q0.M.dim_d2, Eigen::OuterStride<Eigen::Dynamic>(Q0.M.dim_i*Q0.M.dim_d2)) *
            OuterStrideMap(e1.M.mem+g1*e1.M.dim_d2, e1.M.dim_d1, e1.M.dim_d2, Eigen::OuterStride<Eigen::Dynamic>(e1.M.dim_i*e1.M.dim_d2));
  
        }      
      
        prop.update(tgrid.get_t(n)+tgrid.dt/2, tgrid.dt/2);

        Q0.M.resize(NL, Qtmp.M.dim_d1, Qtmp.M.dim_d2*Ngrps2);
        Q0.M.set_zero();
/*
        for(int d_=0; d_<e0.M.dim_d1; d_++){
          for(int d1=0; d1<e1.M.dim_d2; d1++){
            for(int a1=0; a1<NL; a1++){ 
              int g1=e1dict.beta[a1*NL+a1];//int g1=grp2[a1];
              if(g1<0)continue;
              for(int a_=0; a_<NL; a_++){ 
                Q0.M(a_, d_, d1*Ngrps2+g1) += prop.M(a_, a1) * b[0](g1,g1) * Qtmp.M(a1, d_, d1);
        } } } }
        Q0.closure = Vector_otimes(Qtmp.closure, Eigen::VectorXcd::Ones(Ngrps2));
        e0.swap(Q0);
*/
        for(int a1=0; a1<NL; a1++){ 
          int g1=e1dict.beta[a1*NL+a1]; if(g1<0)continue;

          DynamicStrideMap(Q0.M.mem+a1*Q0.M.dim_d2+g1, Qtmp.M.dim_d1, Qtmp.M.dim_d2, Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>(NL*Q0.M.dim_d2,Ngrps2)) = 
            b[0](g1,g1) * 
            OuterStrideMap(Qtmp.M.mem+a1*Qtmp.M.dim_d2, Qtmp.M.dim_d1, Qtmp.M.dim_d2, Eigen::OuterStride<Eigen::Dynamic>(Qtmp.M.dim_i*Qtmp.M.dim_d2));
        }
        e0.M.resize(NL, Q0.M.dim_d1, Q0.M.dim_d2);
        e0.M.set_zero();
        for(int d1=0; d1<Q0.M.dim_d1; d1++){
          RowMatrixMap(e0.M.mem+d1*NL*Q0.M.dim_d2, NL, Q0.M.dim_d2) = 
            prop.M * RowMatrixMap(Q0.M.mem+d1*NL*Q0.M.dim_d2, NL, Q0.M.dim_d2);
  
        }
        e0.closure = Vector_otimes(Qtmp.closure, Eigen::VectorXcd::Ones(Ngrps2));

      }

      PassOn pass_on(PTB.get(0, ForwardPreload).M.dim_d1);
      PTB.get(0, ForwardPreload).sweep_forward(trunc_fwd, pass_on, false);
     
      for(int l=1; l<cur_length-1; l++){
        ProcessTensorElement & e = PTB.get(l, ForwardPreload);
        const ProcessTensorElement & e1 = PTB.peek(l+1);
        e.M.resize(e1.M.dim_i, e1.M.dim_d1*Ngrps2, e1.M.dim_d2*Ngrps2);
        e.M.set_zero();
    
        Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> stride(Ngrps2*Ngrps2*Ngrps2*e1.M.dim_d2, Ngrps2);
        for(int i=0; i<Ngrps2; i++){
          for(int beta=0; beta<Ngrps2; beta++){
            DynamicStrideMap(e.M.mem+(beta*Ngrps2+i)*e1.M.dim_d2*Ngrps2+beta, e1.M.dim_d1, e1.M.dim_d2, stride) =  b[l](i, beta) * e1.M.get_Matrix_d1_d2(i);
//            for(int d1=0; d1<e1.M.dim_d1; d1++){
//              for(int d2=0; d2<e1.M.dim_d2; d2++){
//                e.M(i, d1*Ngrps2+beta, d2*Ngrps2+beta) = b[l](i,beta)*e1.M(i, d1, d2);
//              }
//            }
          } 
        }
        e.closure = Vector_otimes(e1.closure, Eigen::VectorXcd::Ones(Ngrps2));
        e.sweep_forward(trunc_fwd, pass_on, l==cur_length-2 && new_length!=cur_length);
      }
      if(new_length==cur_length){ //have to expand
        ProcessTensorElement & e = PTB.get(new_length-1);
        e.M.resize(e.M.dim_i, Ngrps2, 1);
        e.M.set_zero();
        for(int i=0; i<Ngrps2; i++){
          for(int beta=0; beta<Ngrps2; beta++){
            e.M(i, beta, 0) = b[new_length-1](i, beta);
          }
        }
//        e.closure=Eigen::VectorXcd::Ones(1);
        e.sweep_forward(trunc_fwd, pass_on, true);
      }else{
        PTB.get(new_length-1).close_off();
      }

      //std::cout<<"Dims before compression (i, d1, d2): "; for(int l=0; l<new_length; l++){PTB.get(l, ForwardPreload).M.print_dims();std::cout<<"; ";if(l>=9){std::cout<<"..."; break;}} std::cout<<std::endl;
      
//      PTB.sweep_forward(trunc_fwd, 1, 0, new_length);

      TruncatedSVD trunc_bwd=trunc_layout.get_backward(n, tgrid.n_tot);
      trunc_bwd.keep = sqrt(e1dict.reduced_dim);//diagBB.get_dim();
      std::cout<<"Backward sweep: ";trunc_bwd.print_info();std::cout<<std::endl;
      PTB.sweep_backward(trunc_bwd, 1, 0, new_length);

}

Eigen::MatrixXcd Simulation_CTEMPO::run(Propagator &prop, DiagBB &diagBB,
           const Eigen::MatrixXcd & initial_rho, const TimeGrid &tgrid,
           OutputPrinter &printer, const TruncationLayout &trunc_layout)const{

    int N=initial_rho.rows();
    int NL=N*N;
    if(N<2){
      std::cerr<<"Simulation_CTEMPO::run: N<2!"<<std::endl;
      throw DummyException();
    }
    if(prop.get_dim() != N){
      std::cerr<<"Simulation_CTEMPO: prop.get_dim() != N!"<<std::endl;
      throw DummyException();
    }

    int n_tot=tgrid.n_tot;
    if(n_tot<1){
      std::cerr<<"Simulation_CTEMPO::run: tgrid.n_tot<1!"<<std::endl;
      throw DummyException();
    }
    int n_mem=tgrid.n_mem;
    if(n_mem<2){
      std::cerr<<"Simulation_CTEMPO::run: tgrid.n_mem<1!"<<std::endl;
      throw DummyException();
    }

    //row of Influence Functional
    std::vector<Eigen::MatrixXcd> b=initialize_b(n_mem, diagBB, tgrid.dt);
//    Coupling_Groups_Liouville lgroups(diagBB.groups);
//    int Ngrps2=lgroups.get_Ngrps();
//    const std::vector<int> & grp2=lgroups.grp;


    Eigen::VectorXcd rho_reduced=H_Matrix_to_L_Vector(initial_rho);
    printer.print(0, tgrid.get_t(0), rho_reduced);

    ProcessTensorBuffer PTB;
    initialize_PTB(PTB, n_mem, diagBB, rho_reduced); 

    for(int n=0; n<tgrid.n_tot; n++){
      if(print_timesteps){
        std::cout<<"step: "<<n<<"/"<<tgrid.n_tot<<std::endl;
      }  

      single_step(PTB, prop, b, tgrid, n, trunc_layout);

      //std::cout<<"Dims after compression (i, d1, d2): "; for(int l=0; l<new_length; l++){PTB.get(l, ForwardPreload).M.print_dims();std::cout<<"; "; if(l>=9){std::cout<<"..."; break;}} std::cout<<std::endl;
      if(print_dims_file!="" && n==print_dims_step){
        std::ofstream ofs(print_dims_file);
        for(int l=0; l<PTB.n_tot-n; l++){
          ofs<<PTB.get(l, ForwardPreload).M.dim_d1<<" " \
             <<PTB.get(l, ForwardPreload).M.dim_d2<<std::endl;
        }
      }

      //extract rho_reduced:
      rho_reduced = PTB.get(0).M.get_Matrix_d1i_d2()*PTB.get(0).closure;

      //print output:
      double t_next=tgrid.get_t(n+1);
      printer.print(n+1, t_next, rho_reduced);

    }

    return L_Vector_to_H_Matrix(rho_reduced);
}

void Simulation_CTEMPO::setup(Parameters &param){
    print_timesteps=param.get_as_bool("print_timesteps",true);
    no_system_propagation=param.get_as_bool("no_system_propagation",false);
    print_dims_file=param.get_as_string("print_dims_file","");
    print_dims_step=param.get_as_int("print_dims_step",-1);
}

