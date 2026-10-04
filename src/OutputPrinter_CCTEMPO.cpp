#include "OutputPrinter_CCTEMPO.hpp"
#include "LiouvilleTools.hpp"
#include "MPO_state.hpp"
#include <iomanip>

namespace ACE{

void OutputPrinter_CCTEMPO::set_stream(const std::string &outfile, int precision){
    if(ofs){
      if(ofs->is_open())ofs->close();
      ofs.reset(nullptr);
    }
    if(outfile!="" && outfile!="/dev/null"){
      ofs.reset(new std::ofstream(outfile.c_str()));
      if(precision>0)*ofs<<std::setprecision(precision);
    } 
  }
void OutputPrinter_CCTEMPO::clear(){
    ofs.reset(nullptr);
    OneBody.clear();
    ofs_eigenstates.reset();
    eigenstate_components.clear();
    coeff_list.clear();
}

void OutputPrinter_CCTEMPO::print(int n, double t, const MPO & state, const std::vector<Eigen::VectorXcd> &closure){
    if(!ofs)return;
    if(!ofs->is_open()){
      std::cerr<<"OutputPrinter::print not set up!"<<std::endl;
      throw DummyException();
    }
  
    *ofs<<t;
    //OneBody terms:
    if(OneBody.size()>0){
      std::vector<Eigen::VectorXcd> cap_d2=MPO_op::get_cap_d2(state, closure);
      std::vector<Eigen::VectorXcd> cap_d4=MPO_op::get_cap_d4(state, closure);
      Eigen::VectorXcd rho;
      for(int i=0; i<(int)OneBody.size(); i++){
        int r=OneBody[i].first;
        if(r<0||r>=state.size()){
          std::cerr<<"OutputPrinter_CCTEMPO::print: r<0||r>=state.size()!"<<std::endl;
          throw DummyException();
        }
        if(i==0 || OneBody[i].first!=OneBody[i-1].first){
          rho = state[r].get_reduced_d3(closure[r], cap_d2[r], cap_d4[r]);
        }
        if( OneBody[i].second.rows() != rho.rows() ){
          std::cerr<<"OutputPrinter_CCTEMPO::print: OneBody[i].second.rows() != rho.rows()!"<<std::endl;
          throw DummyException();
        }
        std::complex<double> value = OneBody[i].second.transpose()*rho;
        *ofs<<" "<<value.real()<<" "<<value.imag();
      }
    }

    //coeff_list:
    for(int i=0; i<(int)coeff_list.size(); i++){
      std::complex<double> c=MPO_op::evaluate(state, closure, coeff_list[i]);
      *ofs<<" "<<c.real()<<" "<<c.imag();
    }
    *ofs<<std::endl;
}
void OutputPrinter_CCTEMPO::print(int n, double t, Comb &comb){
    if(!ofs)return;
    if(!ofs->is_open()){
      std::cerr<<"OutputPrinter::print not set up!"<<std::endl;
      throw DummyException();
    }
  
    *ofs<<t;
    //OneBody terms:
    if(OneBody.size()>0){
      std::vector<Eigen::VectorXcd> cap_d2=comb.get_cap_d2();
      std::vector<Eigen::VectorXcd> cap_d4=comb.get_cap_d4();
      Eigen::VectorXcd rho;
      for(int i=0; i<(int)OneBody.size(); i++){
        int r=OneBody[i].first;
        if(r<0||r>=comb.size()){
          std::cerr<<"OutputPrinter_CCTEMPO::print: r<0||r>=comb.size()!"<<std::endl;
          throw DummyException();
        }
        if(i==0 || OneBody[i].first!=OneBody[i-1].first){
          const CombElement & e= comb.get_ro(r, ForwardPreload);
          rho = e.f.get_reduced_d3(e.closure_d1, cap_d2[r], cap_d4[r]);
        }
        if( OneBody[i].second.rows() != rho.rows() ){
          std::cerr<<"OutputPrinter_CCTEMPO::print: OneBody[i].second.rows() != rho.rows()!"<<std::endl;
          throw DummyException();
        }
        std::complex<double> value = OneBody[i].second.transpose()*rho;
        *ofs<<" "<<value.real()<<" "<<value.imag();
      }
    }
    //collective operators:
    for(int i=0; i<(int)collective.size(); i++){
      if(collective[i].size()!=comb.size()){
          std::cerr<<"OutputPrinter_CCTEMPO::print: collective[i].size()!=comb.size()!"<<std::endl;
          throw DummyException();
      }

      Eigen::MatrixXcd cap=Eigen::MatrixXcd::Ones(1,1);
      for(size_t j=0; j<comb.size(); j++){
          const CombElement & e= comb.get_ro(j,ForwardPreload);

        if(collective[i][j].dim_d3!=1){
          std::cerr<<"OutputPrinter_CCTEMPO::print: collective[i][j].dim_d3!=1!"<<std::endl;
          throw DummyException();
        }
        if(collective[i][j].dim_d1!=e.f.dim_d3){
          std::cerr<<"OutputPrinter_CCTEMPO::print: collective[i][j].dim_d1!=comb[j].f.dim_d3!"<<std::endl;
          throw DummyException();
        }
        Eigen::MatrixXcd cap2 = Eigen::MatrixXcd::Zero(e.f.dim_d2,e.f.dim_d3*collective[i][j].dim_d4);
        for(int d2=0; d2<e.f.dim_d2; d2++){
          for(int d3=0; d3<e.f.dim_d3; d3++){
            for(int d4=0; d4<collective[i][j].dim_d4; d4++){
              for(int d=0; d<collective[i][j].dim_d2; d++){
                cap2(d2, d3*collective[i][j].dim_d4+d4) +=
                 cap(d2, d) * collective[i][j](d3, d, 0, d4);
              }
            }
          }
        }
        cap = Eigen::MatrixXcd::Zero(e.f.dim_d4, collective[i][j].dim_d4);
        for(int d4_=0; d4_<e.f.dim_d4; d4_++){
          for(int d4=0; d4<collective[i][j].dim_d4; d4++){
            for(int d2=0; d2<e.f.dim_d2; d2++){
              for(int d3=0; d3<e.f.dim_d3; d3++){
                for(int d1=0; d1<e.f.dim_d1; d1++){
                  cap(d4_, d4) += e.closure_d1(d1) * 
                                  cap2(d2, d3*collective[i][j].dim_d4+d4) * 
                                  e.f(d1, d2, d3, d4_);
        } } } } }
      }
      std::complex<double> c=cap(0,0);
      *ofs<<" "<<c.real()<<" "<<c.imag();
    } 

    //coeff_list:
    for(int i=0; i<(int)coeff_list.size(); i++){
      std::complex<double> c=comb.evaluate(coeff_list[i]);
      *ofs<<" "<<c.real()<<" "<<c.imag();
    }

    *ofs<<std::endl;
}


void OutputPrinter_CCTEMPO::print_eigenstate_occupations(double t, const Eigen::MatrixXcd & H, const MPO & state, const std::vector<Eigen::VectorXcd> & closure){
  if(ofs_eigenstates){
    int N=H.rows();
    if(state.size() != N){
      std::cerr<<"OutputPrinter::print_eigenstate_occupations: state.size() != N!"<<std::endl;
      throw DummyException();
    }
    if(closure.size()!=state.size()){
      std::cerr<<"OutputPrinter::print_eigenstate_occupations: closure.size()!=state.size()!"<<std::endl;
      throw DummyException();
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(H);

    *ofs_eigenstates<<t;
    //First N columns are eigenstate occupations
    for(int i=0; i<N; i++){ 
      Eigen::VectorXcd phi=solver.eigenvectors().col(i);
      Eigen::MatrixXcd rho=phi*phi.adjoint();
      MPO evec=MPO_state::from_SingleExcitationManifold(rho);
      MPO closed(state.size());
      for(int r=0; r<closed.size(); r++){
        closed[r]=state[r];
        closed[r].multiply_d1(closure[r]);
      }
      std::complex<double> c=MPO_state::dot(evec, closed);
      *ofs_eigenstates<<" "<<c.real();
    }
    //Next N columns are eigenvalues
    for(int i=0; i<N; i++){
      double e=solver.eigenvalues()(i);
      *ofs_eigenstates<<" "<<e;
    }
    //Finally, the selected eigenvectors as components
    for(int o=0; o<eigenstate_components.size(); ++o){
      if(eigenstate_components[o]>=N){
        std::cerr<<"OutputPrinter::print_eigenstate_occupations: eigenstate_components["<<o<<"]="<<eigenstate_components[o]<<" >= N="<<N<<std::endl;
        throw DummyException();
      }
      for(int i=0; i<N; i++){
        std::complex<double> c=solver.eigenvectors()(eigenstate_components[o],i);
        *ofs_eigenstates<<" "<<c.real()<<" "<<c.imag();
      }
    }
    *ofs_eigenstates<<std::endl;
  }
}
void OutputPrinter_CCTEMPO::print_eigenstate_occupations(double t, const Eigen::MatrixXcd & H, Comb & comb){
  if(ofs_eigenstates){
    int N=H.rows();
    if(comb.size() != N){
      std::cerr<<"OutputPrinter::print_eigenstate_occupations: comb.size() != N!"<<std::endl;
      throw DummyException();
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(H);

    //First N columns are eigenstate occupations
    MPO closed(comb.size());
    for(int r=0; r<comb.size(); r++){
      const CombElement & e = comb.get_ro(r, ForwardPreload);
      closed[r] = e.f;
      closed[r].multiply_d1( e.closure_d1 );
    }
    *ofs_eigenstates<<t;
    for(int i=0; i<N; i++){ 
      Eigen::VectorXcd phi=solver.eigenvectors().col(i);
      Eigen::MatrixXcd rho=phi*phi.adjoint();
      MPO evec=MPO_state::from_SingleExcitationManifold(rho);

      std::complex<double> c=MPO_state::dot(evec, closed);
      *ofs_eigenstates<<" "<<c.real();
    }
    //Next N columns are eigenvalues
    for(int i=0; i<N; i++){
      double e=solver.eigenvalues()(i);
      *ofs_eigenstates<<" "<<e;
    }
    //Finally, the selected eigenvectors as components
    for(int o=0; o<eigenstate_components.size(); ++o){
      if(eigenstate_components[o]>=N){
        std::cerr<<"OutputPrinter::print_eigenstate_occupations: eigenstate_components["<<o<<"]="<<eigenstate_components[o]<<" >= N="<<N<<std::endl;
        throw DummyException();
      }
      for(int i=0; i<N; i++){
        std::complex<double> c=solver.eigenvectors()(eigenstate_components[o],i);
        *ofs_eigenstates<<" "<<c.real()<<" "<<c.imag();
      }
    }
    *ofs_eigenstates<<std::endl;
  }
}


   
void OutputPrinter_CCTEMPO::setup(Parameters & param){
    clear();

    //setup stream
    std::string outfile = param.get_as_string("outfile", "CCTEMPO.out");
    int precision = param.get_as_int("set_precision",-1);
    set_stream(outfile, precision);

    //setup OneBody
    int N_sites = param.get_as_size_t_check("N_sites");
    for(int r=0; r<N_sites; r++){
      std::string key="S"+int_to_string(r)+"_add_Output";
      int rows=param.get_nr_rows(key);
      for(int i=0; i<rows; i++){
        Eigen::MatrixXcd op=param.get_as_operator(key, Eigen::MatrixXcd(), i, 0);
        OneBody.push_back(std::make_pair(r, H_Matrix_to_L_Vector(op.transpose())));
std::cout<<key<<": "<<std::endl<<L_Vector_to_H_Matrix(OneBody.back().second.transpose())<<std::endl;
      }
    } 
    //setup collective operators
    if(param.is_specified("add_Output_sum")){
      std::vector<std::vector<std::string> > entry=param.get("add_Output_sum");
      for(size_t i=0; i<entry.size(); i++){
        if(entry[i].size()<N_sites){
          std::cerr<<"add_Output {OP_1} {OP_2} ... {OP_N_sites}!"<<std::endl;
          throw DummyException();
        }
        std::vector<Eigen::MatrixXcd> mat_list(N_sites);
        for(size_t j=0; j<N_sites; j++){
          mat_list[j] = ReadExpression(entry[i][j]);
        }
        collective.push_back(MPO_op::OneBody_Sum(mat_list));
        MPO & mpo=collective.back();
        for(size_t j=0; j<N_sites; j++){
          FourLeg tmp(mpo[j].dim_d1*mpo[j].dim_d3, mpo[j].dim_d2, 1, mpo[j].dim_d4);
          for(int d1=0; d1<mpo[j].dim_d1; d1++){
            for(int d2=0; d2<mpo[j].dim_d2; d2++){
              for(int d3=0; d3<mpo[j].dim_d3; d3++){
                for(int d4=0; d4<mpo[j].dim_d4; d4++){
                  tmp(d1*mpo[j].dim_d3+d3, d2, 0, d4) = mpo[j](d1,d2,d3,d4);
          } } } }
          mpo[j]=tmp;
        }
      }
    }
    if(param.is_specified("add_Output_product")){
      std::vector<std::vector<std::string> > entry=param.get("add_Output_product");
      for(size_t i=0; i<entry.size(); i++){
        if(entry[i].size()<N_sites){
          std::cerr<<"add_Output {OP_1} {OP_2} ... {OP_N_sites}!"<<std::endl;
          throw DummyException();
        }
        std::vector<Eigen::MatrixXcd> mat_list(N_sites);
        for(size_t j=0; j<N_sites; j++){
          Eigen::MatrixXcd mat = ReadExpression(entry[i][j]);
          mat_list[j] = H_Matrix_to_L_Vector(mat.transpose()).transpose();
//std::cout<<"TEST: "<<mat_list[j]<<std::endl;
        }
        collective.push_back(MPO_op::Product(mat_list));
/*
        MPO & mpo=collective.back();
        for(size_t j=0; j<N_sites; j++){
          FourLeg tmp(mpo[j].dim_d1*mpo[j].dim_d3, mpo[j].dim_d2, 1, mpo[j].dim_d4);
          for(int d1=0; d1<mpo[j].dim_d1; d1++){
            for(int d2=0; d2<mpo[j].dim_d2; d2++){
              for(int d3=0; d3<mpo[j].dim_d3; d3++){
                for(int d4=0; d4<mpo[j].dim_d4; d4++){
                  tmp(d1*mpo[j].dim_d3+d3, d2, 0, d4) = mpo[j](d1,d2,d3,d4);
          } } } }
          mpo[j]=tmp;
        }
*/
      }
    }

    //setup eigenstates file
    std::string print_eigenstate_occupations = param.get_as_string("print_eigenstate_occupations", "");
    if(print_eigenstate_occupations!="" && print_eigenstate_occupations!="/dev/null"){
      ofs_eigenstates.reset(new std::ofstream(print_eigenstate_occupations.c_str()));
      int set_precision=param.get_as_int("set_precision",-1);
      if(set_precision>0)*ofs_eigenstates<<std::setprecision(set_precision);
    }else if(ofs_eigenstates){
      if(ofs_eigenstates->is_open())ofs_eigenstates->close();
      ofs_eigenstates.reset(nullptr);
    }
    eigenstate_components=param.get_all_size_t("print_eigenstate_component");


    //setup coeff_list
    int coeff_rows = param.get_nr_rows("Output_coeff");
    for(int i=0; i<coeff_rows; i++){
      std::vector<double> coeff = param.get_row_doubles("Output_coeff", i, N_sites);
      std::vector<int> coeff_i(N_sites);
      for(int j=0; j<N_sites; j++){
        coeff_i[j]=coeff[j];
      }
      coeff_list.push_back(coeff_i);
    }
}

}//namespace
