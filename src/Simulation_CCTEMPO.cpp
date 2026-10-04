#include "Simulation_CCTEMPO.hpp"
#include "CCTEMPO.hpp"

namespace ACE{

Eigen::MatrixXcd Simulation_CCTEMPO::get_SX_Hamiltonian_for_eigenstates(double t){
  return SX_Hamiltonian_for_eigenstates;
}


void Simulation_CCTEMPO::run(MPO_Propagator &mpo_prop, CCTEMPO_IF_Nodes &nodes, const MPO & initial, const TimeGrid &tgrid, const TruncationLayout &trunc_layout, OutputPrinter_CCTEMPO &outp){
 
  Comb comb(initial, nodes, buffer_blocksize, buffer_filename);

  //TODO: check if all the sizes match!!!

  TruncatedSVD trunc=trunc_layout.get_base_line(0, tgrid.n_tot);

  mpo_prop.apply_ops.apply(0, comb, trunc);
  outp.print(0, tgrid.ta, comb);
  outp.print_eigenstate_occupations(tgrid.ta, get_SX_Hamiltonian_for_eigenstates(tgrid.ta), comb);

  for(int n=1; n<=tgrid.n_tot; n++){
    std::cout<<"step: "<<n<<std::endl;
    TruncatedSVD trunc=trunc_layout.get_base_line(n-1, tgrid.n_tot);
    //propagate system
    mpo_prop.propagate_comb(comb, tgrid, n, trunc);

    //propagate CTEMPO teeth
    for(int r=0; r<(int)comb.size(); r++){ 
      CombElement & e = comb.get(r, ForwardPreload);
      if(e.tooth && e.tooth->n_tot>1){
        std::cout<<"tooth: "<<r<<std::endl;
        e.tooth->get(0).M=FourLeg_to_MPS_Matrix(e.f);
 
        simC.single_step(*e.tooth, dummy_prop, nodes[r].b, tgrid, n-1, trunc_layout); 
        std::cout<<"tooth["<<r<<"] dims (i, d1, d2): "; for(int k=0; k<e.tooth->n_tot; k++){e.tooth->get(k, ForwardPreload).M.print_dims();std::cout<<"; "; if(k>=9){std::cout<<"..."; break;}} std::cout<<std::endl;

        e.f=MPS_Matrix_to_FourLeg(e.tooth->get(0).M,e.f.dim_d2);
        e.tooth->get(0).M.resize(0,0,0);
        e.closure_d1=e.tooth->get(0).closure;
        if(r<(int)comb.size()-1){
           Eigen::MatrixXcd P;
           std::cout<<"site "<<r<<" ";
           e.sweep_d2_d4(trunc, P, false, 1);
           comb.get(r+1, ForwardPreload).sweep_d2_d4(trunc, P, true, 1);
        }
      }else if(nodes[r].PT && nodes[r].PT->size()>0){
        e.closure_d1=nodes[r].PT->get_closure();
        FourLeg tmp(e.closure_d1.rows(), e.f.dim_d2, e.f.dim_d3, e.f.dim_d4);
//        std::cout<<"n="<<n<<" state["<<r<<"].print_dims(): "; e.f.print_dims(); std::cout<<std::endl;
//        std::cout<<"PT["<<r<<"].list[0].n="<<e.PT->list[0]->n<<std::endl;
        for(int d2=0; d2<e.f.dim_d2; d2++){
          for(int d4=0; d4<e.f.dim_d4; d4++){
            Eigen::MatrixXcd rho(e.f.dim_d3, e.f.dim_d1);
            for(int d3=0; d3<e.f.dim_d3; d3++){
              for(int d1=0; d1<e.f.dim_d1; d1++){
                rho(d3,d1)=e.f(d1,d2,d3,d4);
              } 
            }
            nodes[r].PT->propagate(rho,false);
            for(int d3=0; d3<tmp.dim_d3; d3++){
              for(int d1=0; d1<tmp.dim_d1; d1++){
                tmp(d1,d2,d3,d4)=rho(d3,d1);
              } 
            }
          } 
        }
        nodes[r].PT->load_next();
        e.f.swap(tmp);
//        std::cout<<"state["<<r<<"].print_dims(): "; state[r].print_dims(); std::cout<<std::endl;
      }
    }
    
    comb.sweep_d4_d2(trunc);

    mpo_prop.propagate_comb_second(comb, tgrid, n, trunc);

    //extract output
    //std::complex<double> n0=MPO_op::evaluate(state, closure, index0);
    //os<<tgrid.get_t(n)<<" "<<n0.real()<<" "<<n0.imag()<<std::endl;
    mpo_prop.apply_ops.apply(n, comb, trunc);
    outp.print(n, tgrid.get_t(n), comb);
    outp.print_eigenstate_occupations(tgrid.get_t(n), get_SX_Hamiltonian_for_eigenstates(tgrid.get_t(n)), comb);
    std::cout<<"comb.print_dims(): "; comb.print_dims(); std::cout<<std::endl;
  }
}


void Simulation_CCTEMPO::setup(Parameters &param){
  dummy_prop=FreePropagator();
  simC=Simulation_CTEMPO(param);
  simC.no_system_propagation=true;
  use_symmetric_Trotter=param.get_as_bool("use_symmetric_Trotter", true);

  buffer_blocksize=param.get_as_int("buffer_blocksize");
  buffer_filename=param.get_as_string("buffer_filename", "");

  FreePropagator SXprop=CCTEMPO::get_SX_FreePropagator(param);
  //TODO: include time-dependent individual terms?
  { Eigen::MatrixXcd H_indiv=CCTEMPO::get_SX_Hamiltonian_from_individual_terms(param);
    if(H_indiv.norm()>1e-14)SXprop.add_Hamiltonian(H_indiv); }
  SX_Hamiltonian_for_eigenstates = SXprop.const_H;

}
 
}//namespace
