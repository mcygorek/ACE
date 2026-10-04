#include "Simulation_TEMPO.hpp"
#include "DummyException.hpp"
#include "Timings.hpp"
#include "Coupling_Groups.hpp"
#include "otimes.hpp"

using namespace ACE;


void Simulation_TEMPO::initialize_PTB(ProcessTensorBuffer &PTB, int n_mem, const DiagBB & diagBB, const Eigen::VectorXcd & rho_reduced){

    int N = sqrt(rho_reduced.rows());
    int NL=N*N;
    int NdiagBB=diagBB.get_dim();
    int NLdiagBB=NdiagBB*NdiagBB;
    std::cout<<"diagBB.sys_dim()="<<diagBB.sys_dim()<<" ";
    std::cout<<"diagBB.get_dim()="<<diagBB.get_dim()<<std::endl;
    if(diagBB.sys_dim() != N){
      std::cerr<<"Simulation_TEMPO::initialize_PTB: diagBB.sys_dim() != N ("<<diagBB.sys_dim()<<" vs. "<<N<<")!"<<std::endl;
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
      //e.accessor.dict.set_default_diag(N); // full set of outer bonds
      e.accessor.dict.set_trivial(1); // rho_reduces encoded as vector wrt. d1.
      e.M.set_from_Matrix_d1i_d2(rho_reduced, 1); 
    }
//    return PTB;
}
std::vector<Eigen::MatrixXcd> Simulation_TEMPO::initialize_b(int n_mem, DiagBB &diagBB, double dt){
  std::vector<Eigen::MatrixXcd> b(n_mem);
  for(int n=0; n<n_mem; n++){
    b[n]=diagBB.calculate_expS(n, dt);
  }
  return b;
}


void Simulation_TEMPO::single_step(ProcessTensorBuffer & PTB, Propagator &prop,const std::vector<Eigen::MatrixXcd> &b, const TimeGrid &tgrid, int n, const TruncationLayout & trunc_layout)const{

      if(PTB.n_tot<2){
        std::cerr<<"Simulation_TEMPO::single_step: PTB.n_tot<2!"<<std::endl;
        throw DummyException();
      }
      int n_mem=PTB.n_tot;
      if(b.size()<n_mem){
        std::cerr<<"Simulation_TEMPO::single_step: b.size()<n_mem!"<<std::endl;
        throw DummyException();
      }
      
      int NL; 
      IF_OD_Dictionary e1dict;
      {
        const ProcessTensorElement & e0 = PTB.get(0, ForwardPreload);
        NL=e0.M.dim_d1;//accessor.dict.get_N();
        const ProcessTensorElement & e1 = PTB.peek(1);
        e1dict=e1.accessor.dict;
      }
      int Ngrps2=e1dict.reduced_dim;

      double keep=0.;//sqrt(Ngrps2);//diagBB.get_dim();

      int cur_length=n; 
      if(cur_length>n_mem){cur_length=n_mem;}
      int new_length=n+1;   
      if(new_length>n_mem){new_length=n_mem;}
//      std::cout<<"n="<<n<<" n_tot="<<tgrid.n_tot<<" cur_length="<<cur_length<<" new_length="<<new_length<<std::endl;

      //TruncatedSVD trunc = trunc_layout.get_base();
      TruncatedSVD trunc_fwd=trunc_layout.get_forward(n, tgrid.n_tot);
      trunc_fwd.keep = keep;
      std::cout<<"Forward sweep: ";trunc_fwd.print_info(); std::cout<<std::endl;


      using RowMatrixMap = Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >;
      using OuterStrideMap = Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 0, Eigen::OuterStride<Eigen::Dynamic> >;
      using InnerStrideMap = Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 0, Eigen::InnerStride<Eigen::Dynamic> >;
      using DynamicStrideMap = Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> >;
      using OuterStrideColMap = Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>, 0, Eigen::OuterStride<Eigen::Dynamic> >;

      ProcessTensorElement prev_e; //store previous e (expand from font)
      { // propagate first element(s):
        ProcessTensorElement & e0 = PTB.get(0, ForwardPreload);
        const ProcessTensorElement & e1 = PTB.peek(1);
        ProcessTensorElement Q0;
        Q0.set_trivial(Ngrps2);
        Q0.accessor.dict = e1.accessor.dict; 
        Q0.env_ops.clear();// = e0.env_ops;

 
        Eigen::MatrixXcd Mfirst, Msecond;
        if(!no_system_propagation){
          prop.update(tgrid.get_t(n), tgrid.dt/2);
          Mfirst=prop.M;
          prop.update(tgrid.get_t(n)+tgrid.dt/2, tgrid.dt/2);
          Msecond=prop.M;
        }
        if(n==0){
          Q0.M.resize(Ngrps2, NL, 1);
          Q0.M.set_zero();
          Q0.closure = Eigen::VectorXcd::Ones(1);
          Eigen::MatrixXcd rho0 = e0.M.get_Matrix_d1i_d2();
//std::cout<<"rho0.rows()="<<rho0.rows()<<" rho0.cols()="<<rho0.cols()<<std::endl;
//Q0.M.print_dims();std::cout<<std::endl;
          if(no_system_propagation){
            for(int i=0; i<NL; i++){ 
              int gi=e1dict.beta[i*NL+i];
              if(gi<0)continue;
                Q0.M(gi, i, 0) += b[0](gi,gi) * rho0(i);
            }
          }else{
            for(int d1=0; d1<NL; d1++){
              for(int i=0; i<NL; i++){ 
                int gi=e1dict.beta[i*NL+i];
                if(gi<0)continue;
                for(int d2=0; d2<NL; d2++){
                  Q0.M(gi, d1, 0) += Msecond(d1, i) * b[0](gi,gi) * Mfirst(i, d2) * rho0(d2);
            } } }
          }
          e0.swap(Q0);
          return;
        }

        //first site when n!=0:
        Q0.M.resize(Ngrps2, NL, NL*Ngrps2);
        Q0.M.set_zero();
        Q0.closure = Eigen::VectorXcd::Ones(Q0.M.dim_d2);

        if(no_system_propagation){
          for(int i=0; i<NL; i++){ 
            int gi=e1dict.beta[i*NL+i];
            if(gi<0)continue;
              Q0.M(gi, i, i*Ngrps2+gi) += b[0](gi,gi);
          } 
        }else{
          for(int d1=0; d1<NL; d1++){
            for(int i=0; i<NL; i++){ 
              int gi=e1dict.beta[i*NL+i];
              if(gi<0)continue;
              for(int d2=0; d2<NL; d2++){
                Q0.M(gi, d1, d2*Ngrps2+gi) += Msecond(d1, i) * b[0](gi,gi) * Mfirst(i, d2);
          } } }
        }
        prev_e.swap(e0);
        e0.swap(Q0);
      }

      //all others
      PassOn pass_on(PTB.get(0, ForwardPreload).M.dim_d1);
      PTB.get(0, ForwardPreload).sweep_forward(trunc_fwd, pass_on, false);

      for(int l=1; l<new_length; l++){
        ProcessTensorElement Q0;
        Q0.set_trivial(Ngrps2);
        Q0.accessor.dict = e1dict; 
        Q0.env_ops.clear();// = e0.env_ops;
        Q0.closure=Vector_otimes(prev_e.closure,Eigen::VectorXcd::Ones(Ngrps2));

#ifdef TEMPO_PASS_ON_MULTIPLY_SEPARATELY
        Q0.M.resize(Ngrps2, prev_e.M.dim_d1*Ngrps2, prev_e.M.dim_d2*Ngrps2);
        Q0.M.set_zero();
/*
        for(int d1=0; d1<prev_e.M.dim_d1; d1++){
          for(int d2=0; d2<prev_e.M.dim_d2; d2++){
            for(int gi=0; gi<Ngrps2; gi++){
              for(int gj=0; gj<Ngrps2; gj++){
                Q0.M(gj, d1*Ngrps2+gi, d2*Ngrps2+gi ) += b[l](gi,gj) * prev_e.M(gj, d1, d2);
        } } } }
*/
        for(int gj=0; gj<Ngrps2; gj++){
          for(int gi=0; gi<Ngrps2; gi++){
            DynamicStrideMap(Q0.M.mem+(gi*Ngrps2+gj)*Q0.M.dim_d2+gi, prev_e.M.dim_d1, prev_e.M.dim_d2, Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>(Ngrps2*Q0.M.dim_i*Q0.M.dim_d2, Ngrps2)) = 
            b[l](gi,gj) * prev_e.M.get_Matrix_d1_d2(gj);
            
        } }
#else
        Q0.M.resize(Ngrps2, pass_on.P.rows(), prev_e.M.dim_d2*Ngrps2);
        Q0.M.set_zero();
        for(int gi=0; gi<Ngrps2; gi++){
          for(int gj=0; gj<Ngrps2; gj++){
//            for(int r=0; r<pass_on.P.rows(); r++){
//              for(int d2=0; d2<prev_e.M.dim_d2; d2++){
//                for(int d1=0; d1<prev_e.M.dim_d1; d1++){
//                  Q0.M(gj, r, d2*Ngrps2+gi ) += pass_on.P(r,d1*Ngrps2+gi) *b[l](gi,gj) * prev_e.M(gj, d1, d2);
//        } } } } }

            DynamicStrideMap(Q0.M.mem+gj*Q0.M.dim_d2+gi, Q0.M.dim_d1, prev_e.M.dim_d2, Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>(Q0.M.dim_i*Q0.M.dim_d2, Ngrps2)) +=
            b[l](gi,gj) * OuterStrideColMap(pass_on.P.data()+gi*pass_on.P.rows(), pass_on.P.rows(), prev_e.M.dim_d1, Eigen::OuterStride<Eigen::Dynamic>(pass_on.P.rows()*Ngrps2)) * 
            OuterStrideMap(prev_e.M.mem+gj*prev_e.M.dim_d2, prev_e.M.dim_d1, prev_e.M.dim_d2, Eigen::OuterStride<Eigen::Dynamic>(prev_e.M.dim_i*prev_e.M.dim_d2));
        } }
        pass_on.P=Eigen::MatrixXcd(0,0);  //M already multiplied with P
#endif

        ProcessTensorElement & e = PTB.get(l, ForwardPreload);
        prev_e.swap(e);
        e.swap(Q0);
        e.sweep_forward(trunc_fwd, pass_on, l==new_length-1);
      }
      //close off last site
      if(cur_length==new_length){
        Eigen::VectorXcd cap=Eigen::VectorXcd::Zero(prev_e.M.dim_d1*Ngrps2);
        for(int d1=0; d1<prev_e.M.dim_d1; d1++){
          for(int gj=0; gj<Ngrps2; gj++){
            for(int gi=0; gi<Ngrps2; gi++){
              cap(d1*Ngrps2+gj)+=prev_e.M(gi, d1, 0);     
            }
          }
        } 
        ProcessTensorElement & e = PTB.get(new_length-1, ForwardPreload);
        e.M.inner_multiply_right(cap);
        e.closure=Eigen::VectorXcd::Ones(1);
      }else{
        ProcessTensorElement & e = PTB.get(new_length-1, ForwardPreload);
        MPS_Matrix M(e.M.dim_i, e.M.dim_d1, 1);
        M.set_zero();
        for(int d1=0; d1<e.M.dim_d1; d1++){
          for(int i=0; i<e.M.dim_i; i++){
            for(int d2=0; d2<e.M.dim_d2; d2++){
              M(i, d1, 0) += e.M(i, d1, d2);
        } } }
        e.M.swap(M);
        e.closure=Eigen::VectorXcd::Ones(1);
      }
 
//      std::cout<<"Dims before compression (i, d1, d2): "; for(int l=0; l<new_length; l++){PTB.get(l, ForwardPreload).M.print_dims();std::cout<<"; ";if(l>=9){std::cout<<"..."; break;}} std::cout<<std::endl;
//      PTB.sweep_forward(trunc_fwd, 1, 0, new_length);
      TruncatedSVD trunc_bwd=trunc_layout.get_backward(n, tgrid.n_tot);
      trunc_bwd.keep = keep;
      std::cout<<"Backward sweep: ";trunc_bwd.print_info();std::cout<<std::endl;
      PTB.sweep_backward(trunc_bwd, 1, 0, new_length);

}

Eigen::MatrixXcd Simulation_TEMPO::run(Propagator &prop,  DiagBB &diagBB,
           const Eigen::MatrixXcd & initial_rho, const TimeGrid &tgrid,
           OutputPrinter &printer, const TruncationLayout &trunc_layout)const{

    int N=initial_rho.rows();
    int NL=N*N;
    if(N<2){
      std::cerr<<"Simulation_TEMPO::run: N<2!"<<std::endl;
      throw DummyException();
    }
    if(prop.get_dim() != N){
      std::cerr<<"Simulation_TEMPO: prop.get_dim() != N!"<<std::endl;
      throw DummyException();
    }

    int n_tot=tgrid.n_tot;
    if(n_tot<1){
      std::cerr<<"Simulation_TEMPO::run: tgrid.n_tot<1!"<<std::endl;
      throw DummyException();
    }
    int n_mem=tgrid.n_mem;
    if(n_mem<2){
      std::cerr<<"Simulation_TEMPO::run: tgrid.n_mem<1!"<<std::endl;
      throw DummyException();
    }

    //row of Influence Functional
    //std::cout<<"n_mem="<<n_mem<<std::endl;
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
        for(int l=0; l<=n && l<PTB.n_tot; l++){
          ofs<<PTB.get(l, ForwardPreload).M.dim_d1<<" " \
             <<PTB.get(l, ForwardPreload).M.dim_d2<<std::endl;
        }
      }
      //extract rho_reduced:
      int length=n+1; if(length>n_mem){length=n_mem;}
      Eigen::VectorXcd cap=Eigen::VectorXcd::Ones(1);
      for(int l=length-1; l>=0; l--){
        const ProcessTensorElement & e = PTB.get(l,BackwardPreload);
        rho_reduced = Eigen::VectorXcd::Zero(e.M.dim_d1);
        for(int d1=0; d1<e.M.dim_d1; d1++){
          for(int i=0; i<e.M.dim_i; i++){
            for(int d2=0; d2<e.M.dim_d2; d2++){
              rho_reduced(d1)+=e.M(i,d1,d2)*cap(d2);
            }
          }
        }  
        cap=rho_reduced;
      }

      //print output:
      double t_next=tgrid.get_t(n+1);
      printer.print(n+1, t_next, rho_reduced);

    }

    return L_Vector_to_H_Matrix(rho_reduced);
}

void Simulation_TEMPO::setup(Parameters &param){
    print_timesteps=param.get_as_bool("print_timesteps",true);
    no_system_propagation=param.get_as_bool("no_system_propagation",false);
    print_dims_file=param.get_as_string("print_dims_file","");
    print_dims_step=param.get_as_int("print_dims_step",-1);
}

