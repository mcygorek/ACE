#include "FreePropagator.hpp"
#include "Constants.hpp"
#include "Pulse_Selector.hpp"
#include "TimedepMatrix.hpp"
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>
#include <vector>
#include "Propagator.hpp"
#include "Parameters.hpp"
#include "Operators.hpp"
#include "HilbertSpaceRotation.hpp"
#include "MultitimeOp.hpp"
#include "TimeGrid.hpp"
#include "Reader.hpp"
#include "ComplexFunction_Interpolate.hpp"
#include "DummyException.hpp"

namespace ACE{

  //get Dimension of Hamiltonian
  int FreePropagator::get_dim()const{
    return const_H.rows();
  }
  //Either set dimension with a corresponding (zero) Hamiltonian matrix, or make sure adding a new term is consistent with existing dimension
  int FreePropagator::set_dim(int dim, const std::string &error_comment){
    if(!dim_fixed){
      const_H=Eigen::MatrixXcd::Zero(dim,dim);
      dim_fixed=true;
      return dim;
    }else{
      if(const_H.rows()!=dim || const_H.cols()!=dim){
        std::cerr<<"FreePropagator::set_dim: Mismatch with dimension of Hamiltonian!"<<std::endl;
        if(error_comment!=""){
          std::cerr<<"Comment: "<<error_comment<<std::endl;
        }
        throw DummyException();
      }
    
      return const_H.rows();
    }
  } 
  int FreePropagator::set_dim(const Eigen::MatrixXcd M, const std::string &error_comment){
    set_dim(M.cols(),error_comment);
    return set_dim(M.rows(), error_comment);
  }
  

  ///sanity checks
  void FreePropagator::check_dimensions()const{
    if(const_H.rows()!=const_H.cols()){
      std::cerr<<"Hamiltonian::check_dimensions: const_H.rows()!=const_H.cols()!"<<std::endl;
      throw DummyException();
    }
    for(size_t i=0; i<timedep_H.size(); i++){
      if(timedep_H[i]->get_dim()!=const_H.cols()){
        std::cerr<<"Hamiltonian::check_dimensions: timedep_H[i].get_dim()!=const_H.cols() for i="<<i<<"!"<<std::endl;
        throw DummyException();
      }
    } 

    for(size_t i=0; i<multitime_op.size(); i++){
      if(multitime_op[i].op.rows()!=const_H.rows()*const_H.rows()){
        std::cerr<<"Hamiltonian::check_dimensions: multitime_op[i].op.rows()!=const_H.rows()*const_H.rows() for i="<<i<<"!"<<std::endl;
        throw DummyException();
      }
    }
  }

  Eigen::MatrixXcd FreePropagator::get_Htot(double t)const{
    //constant part:
    Eigen::MatrixXcd matH=const_H;
    //add time-dependent parts:
    for(size_t i=0; i<timedep_H.size(); i++){
      matH += timedep_H[i]->f(t);
    }
    return matH;
  }
  Eigen::MatrixXcd FreePropagator::get_H_forward(double t)const{
    //constant part:
    Eigen::MatrixXcd matH=Eigen::MatrixXcd::Zero(get_dim(),get_dim());
    //add time-dependent parts:
    for(size_t i=0; i<timedep_H_forward.size(); i++){
      matH += timedep_H_forward[i]->f(t);
    }
    return matH;
  }
  Eigen::MatrixXcd FreePropagator::get_H_backward(double t)const{
    //constant part:
    Eigen::MatrixXcd matH=Eigen::MatrixXcd::Zero(get_dim(),get_dim());
    //add time-dependent parts:
    for(size_t i=0; i<timedep_H_backward.size(); i++){
      matH += timedep_H_backward[i]->f(t);
    }
    return matH;
  }

  /* returns a matrix in Liouville space whose ".exp()" gives the time evolution operator \mathcal{M} in Liouville space */
  Eigen::MatrixXcd FreePropagator::Hamiltonian_to_Liouville_Generator(const Eigen::MatrixXcd &H, double dt){
    return forward_Hamiltonian_to_Liouville_Generator(H,dt) 
             + backward_Hamiltonian_to_Liouville_Generator(H,dt);
  }

  Eigen::MatrixXcd FreePropagator::forward_Hamiltonian_to_Liouville_Generator(const Eigen::MatrixXcd &H, double dt){

    int dim=H.rows();
    std::complex<double> gamma(0.,-dt/hbar_in_meV_ps);

    Eigen::MatrixXcd gen=Eigen::MatrixXcd::Zero(dim*dim, dim*dim);

    for(int i0=0; i0<dim; i0++){
      for(int i1=0; i1<dim; i1++){
        for(int i2=0; i2<dim; i2++){
            gen(i0*dim+i1, i2*dim+i1) += gamma*H(i0,i2);
        }
      }
    }
    return gen;
  }
  Eigen::MatrixXcd FreePropagator::backward_Hamiltonian_to_Liouville_Generator(const Eigen::MatrixXcd &H, double dt){

    int dim=H.rows();
    std::complex<double> gamma(0.,-dt/hbar_in_meV_ps);

    Eigen::MatrixXcd gen=Eigen::MatrixXcd::Zero(dim*dim, dim*dim);

    for(int i0=0; i0<dim; i0++){
      for(int i1=0; i1<dim; i1++){
        for(int i3=0; i3<dim; i3++){
          gen(i0*dim+i1, i0*dim+i3) += -gamma*std::conj(H(i1,i3));
        }
      }
    } 
    return gen;
  }

  bool FreePropagator::has_to_recalculate(double dt)const{
    if((precalculated 
        && fabs(precalculated_dt-dt)<1e-16 
        && timedep_H.size()<1  
        && timedep_H_forward.size()<1
        && timedep_H_backward.size()<1
       ) || never_update){
      return false;
    }else{
      return true;
    }
  }

  Eigen::MatrixXcd FreePropagator::Total_Generator(double t, double dt){
    size_t dim=get_dim();
    Eigen::MatrixXcd gen=Eigen::MatrixXcd::Zero(dim*dim, dim*dim);
    if(nonH.rows()>1)gen=nonH*dt;
    gen+=Hamiltonian_to_Liouville_Generator(get_Htot(t),  dt);
    return gen; 
  }

  void FreePropagator::update_single(double t, double dt){
    Eigen::MatrixXcd gen=Total_Generator(t,dt);
    if(timedep_H_forward.size()>0){
      gen+=forward_Hamiltonian_to_Liouville_Generator(get_H_forward(t),dt);
    }
    if(timedep_H_backward.size()>0){
      gen+=backward_Hamiltonian_to_Liouville_Generator(get_H_backward(t),dt);
    }
    M=gen.exp();
    M=hs_rot.apply_Liouville(M);
  }

  void FreePropagator::update(double t, double dt){
#ifdef DEBUG
std::cout<<"DEBUG: FreePropagator: update called: t="<<t<<std::endl;
//if(precalculated)std::cout<<"DEBUG: precalculated=true";
//else std::cout<<"DEBUG: precalculated=false";
//std::cout<<" precalculated_dt-dt="<<precalculated_dt-dt<<" timedep_H.size()"<<std::endl;
#endif

    if(has_to_recalculate(dt)){
      if(Nintermediate<1){ 
        update_single(t, dt);
      }else{
        double dt2=dt/(Nintermediate+1);
        int NL=get_dim()*get_dim();
        Eigen::MatrixXcd N=Eigen::MatrixXcd::Identity(NL,NL);
        for(int i=0; i<Nintermediate+1; i++){
          update_single(t+i*dt2, dt2);
          N=M*N;
        }
        M=N;
      }
      precalculated=true;
      precalculated_dt=dt;
    }else{
#ifdef DEBUG
std::cout<<"DEBUG: skipping calculation (precalculated)"<<t<<std::endl;
#endif
       //skip calculation
    }
 
//apply multitime operators: 
//(handle case of two operators at equal times with apply_operator=false in second loop)
    for(int o=(int)multitime_op.size()-1; o>=0; o--){
      if((!multitime_op[o].apply_before) && multitime_op[o].is_now(t,dt)){
        multitime_op[o].apply(M);
        precalculated=false;
      }
    }
    for(size_t o=0; o<multitime_op.size(); o++){
      if(multitime_op[o].apply_before && multitime_op[o].is_now(t,dt)){
        multitime_op[o].apply(M);
        precalculated=false;
      }
    }
  }

  void FreePropagator::set_Hamiltonian(const Eigen::MatrixXcd & H){
    const_H=H;
    precalculated=false;
    dim_fixed=true;
    check_dimensions();
  }

  void FreePropagator::add_Hamiltonian(const Eigen::MatrixXcd & H){
    if(H.rows()!=H.cols()){
      std::cerr<<"FreePropagator::add_Hamiltonian: H not square!"<<std::endl;
      throw DummyException();
    }
    set_dim(H);
    const_H+=H;
    precalculated=false;
  }

  void FreePropagator::add_Pulse(const std::pair<std::vector<double>,std::vector<std::complex<double> > > & shape, const Eigen::MatrixXcd &A){

    ComplexFunctionPtr ptr=std::make_shared<ComplexFunction_Interpolate>();
    ComplexFunction_Interpolate* p = static_cast<ComplexFunction_Interpolate*>(ptr.get());
    int length=shape.first.size();
    p->val.resize(length);
    for(int i=0; i<length; i++){
      p->val[i].first=shape.first[i];
      p->val[i].second= (i<shape.second.size() ? shape.second[i] : 0.);
    }
    add_Pulse(ptr, A);
  }

  void FreePropagator::add_Pulse(ComplexFunctionPtr &f, const Eigen::MatrixXcd &A){
    if(A.rows()!=A.cols()){
        std::cerr<<"FreePropagator::add_Pulse: Operator not square!"<<std::endl;
        throw DummyException();
    }
    set_dim(A,"FreePropagator::add_Pulse: Mismatch in dimensions!");
    timedep_H.push_back(std::make_shared<TimedepMatrix_Const_Shape_Hermitian>(f, A));
  }
  void FreePropagator::add_forward_Pulse(ComplexFunctionPtr &f, const Eigen::MatrixXcd &A){
    if(A.rows()!=A.cols()){
        std::cerr<<"FreePropagator::add_forward_Pulse: Operator not square!"<<std::endl;
        throw DummyException();
    }
    set_dim(A,"FreePropagator::add_forward_Pulse: Mismatch in dimensions!");
    timedep_H_forward.push_back(std::make_shared<TimedepMatrix_Const_Shape_Hermitian>(f, A));
  }  
  void FreePropagator::add_backward_Pulse(ComplexFunctionPtr &f, const Eigen::MatrixXcd &A){
    if(A.rows()!=A.cols()){
        std::cerr<<"FreePropagator::add_backward_Pulse: Operator not square!"<<std::endl;
        throw DummyException();
    }
    set_dim(A,"FreePropagator::add_backward_Pulse: Mismatch in dimensions!");
    timedep_H_backward.push_back(std::make_shared<TimedepMatrix_Const_Shape_Hermitian>(f, A));
  }
  void FreePropagator::add_TimedepMatrix(TimedepMatrixPtr & mp){
    set_dim(mp->get_dim(), "FreePropagator::add_Timedep: mp");
    timedep_H.push_back(mp);
  }

  void FreePropagator::add_Lindblad(double gamma, const Eigen::MatrixXcd & L, const Eigen::MatrixXcd & Ldag_){
    if(nonH.rows()==0){
      nonH=Eigen::MatrixXcd::Zero(L.rows()*L.rows(), L.cols()*L.cols());
    }
    if(L.rows()!=L.cols()){
      std::cerr<<"FreePropagator::add_Lindblad: L.rows()!=L.cols()"<<std::endl;
      throw DummyException();
    }

    Eigen::MatrixXcd Ldag=Ldag_;
    if(Ldag.rows()<1){ 
      Ldag=L.adjoint();
    }else{
      if(Ldag.rows()!=Ldag.cols()){
        std::cerr<<"FreePropagator::add_Lindblad: Ldag.rows()!=Ldag.cols()"<<std::endl;
        throw DummyException();
      }
      if(L.rows()!=Ldag.cols()){
        std::cerr<<"FreePropagator::add_Lindblad: L.rows()!=Ldag.cols()"<<std::endl;
        throw DummyException();
      }
    }
   
    size_t dim=set_dim(L, "FreePropagator::add_Lindblad: L.rows()!=H.get_dim()");
    
    double g=gamma;
    Eigen::MatrixXcd L2=Ldag*L;

    for(size_t i0=0; i0<dim; i0++){
      for(size_t i1=0; i1<dim; i1++){
        for(size_t i2=0; i2<dim; i2++){
          for(size_t i3=0; i3<dim; i3++){
            nonH(i0*dim+i1, i2*dim+i3) += g*L(i0,i2)*Ldag(i3,i1);
          }
        }
      }
    }
    for(size_t i0=0; i0<dim; i0++){
      for(size_t i1=0; i1<dim; i1++){
        for(size_t i2=0; i2<dim; i2++){
          nonH(i0*dim+i1, i2*dim+i1) += -0.5*g*L2(i0,i2);
        }
      }
    }
    for(size_t i0=0; i0<dim; i0++){
      for(size_t i1=0; i1<dim; i1++){
        for(size_t i3=0; i3<dim; i3++){
          nonH(i0*dim+i1, i0*dim+i3) += -0.5*g*L2(i3,i1);
        }
      }
    }
    precalculated=false;
  }

  void FreePropagator::add_diagonal_loss(double gamma, int index){
    int dim=get_dim();
    if(index<0||index>=dim){
      std::cerr<<"FreePropagator_nonHamiltonian::add_diagonal_loss: index<0||index>=dim"<<std::endl;
      throw DummyException();
    }
    for(int k=0; k<dim; k++){
      nonH( index*dim+k , index*dim+k) -= gamma/2.;
    }
    for(int k=0; k<dim; k++){
      nonH( k*dim+index , k*dim+index) -= gamma/2.;
    }
    precalculated=false;
  }
  
  void FreePropagator::add_MultitimeOp(double t, Eigen::MatrixXcd m1, Eigen::MatrixXcd m2, bool apply_before){
    set_dim(m1); set_dim(m2); 
    multitime_op.push_back(MultitimeOp(t, m1, m2, apply_before));
  }

  void FreePropagator::print_pulses(Parameters &param){
    std::string file=param.get_as_string("print_pulses");
    if(file=="")return;

//std::cout<<"PULSES: "<<timedep_H.size()<<" "<<timedep_H_forward.size()<<" "<<timedep_H_backward.size()<<std::endl;
    if(timedep_H.size()>0){
      TimeGrid tgrid(param);
      std::ofstream ofs(file.c_str());
      for(int ti=0; ti<=tgrid.n_tot; ti++){
        double t=tgrid.get_t(ti);
        ofs<<t;
        Eigen::MatrixXcd mat=timedep_H[0]->f(t);
        for(size_t i=1; i<timedep_H.size(); i++){
          mat+=timedep_H[i]->f(t);
        }
        for(int i=0; i<mat.rows(); i++){
          for(int j=0; j<mat.cols(); j++){
            ofs<<" "<<mat(i,j).real()<<" "<<mat(i,j).imag();
          }
        } 
        ofs<<std::endl;
      }
    }
    if(timedep_H_forward.size()>0){
      TimeGrid tgrid(param);
      std::string file2=file+"_forward";
      std::ofstream ofs(file2.c_str());
      for(int ti=0; ti<=tgrid.n_tot; ti++){
        double t=tgrid.get_t(ti);
        ofs<<t;
        Eigen::MatrixXcd mat=timedep_H_forward[0]->f(t);
        for(size_t i=1; i<timedep_H_forward.size(); i++){
          mat+=timedep_H_forward[i]->f(t);
        }
        for(int i=0; i<mat.rows(); i++){
          for(int j=0; j<mat.cols(); j++){
            ofs<<" "<<mat(i,j).real()<<" "<<mat(i,j).imag();
          }
        } 
        ofs<<std::endl;
      }
    } 
    if(timedep_H_backward.size()>0){
      TimeGrid tgrid(param);
      std::string file2=file+"_backward";
      std::ofstream ofs(file2.c_str());
      for(int ti=0; ti<=tgrid.n_tot; ti++){
        double t=tgrid.get_t(ti);
        ofs<<t;
        Eigen::MatrixXcd mat=timedep_H_backward[0]->f(t);
        for(size_t i=1; i<timedep_H_backward.size(); i++){
          mat+=timedep_H_backward[i]->f(t);
        }
        for(int i=0; i<mat.rows(); i++){
          for(int j=0; j<mat.cols(); j++){
            ofs<<" "<<mat(i,j).real()<<" "<<mat(i,j).imag();
          }
        } 
        ofs<<std::endl;
      }
    }
  }

  void FreePropagator::setup(Parameters &param, int setdim){
    initialize();
    if(setdim>=0)set_dim(setdim);
    Nintermediate=param.get_as_size_t("Nintermediate", 0);

    //Special case: magnetic field of spin 1/2
    if(param.is_specified("TLS_B_eff")){
      set_dim(2);
      std::vector<double> dvec=param.get_row_doubles("TLS_B_eff", 0, 3);
      const_H+= dvec[0]*0.5*sigma_x()
              + dvec[1]*0.5*sigma_y()
              + dvec[2]*0.5*sigma_z();
    }
  
    //Pulses:
    std::vector<std::vector<std::string> >  svv;

    if(param.is_specified("add_Hamiltonian")){
      for(int i=0; i<param.get_nr_rows("add_Hamiltonian"); i++){
        std::string str=param.get_as_single_string("add_Hamiltonian", i);
std::cout<<"add_Hamiltonian: "<<str<<std::endl;
        add_Hamiltonian(ReadExpression(str));
      }

      std::cout<<"Hamiltonian: "<<std::endl<<const_H<<std::endl;
    }

    if(param.is_specified("add_Lindblad")){
      std::vector<std::vector<std::string> > ssv=param.get("add_Lindblad");
      for(size_t i=0; i<ssv.size(); i++){
        if(ssv[i].size()<2){
          std::cerr<<"Error: add_Lindblad needs 2 arguments: rate L [Ldagger]!"<<std::endl;
          throw DummyException();
        }
        double rate=readDouble(ssv[i][0],"add_Lindblad _rate_ L [Ldagger]");
        Eigen::MatrixXcd M=ReadExpression(ssv[i][1]);

        Eigen::MatrixXcd Mdag=M.adjoint();
        if(ssv[i].size()>2)Mdag=ReadExpression(ssv[i][2]);
        add_Lindblad(rate, M, Mdag);
      }
    }
    if(param.is_specified("apply_Operator")){
      std::cerr<<"Parameter 'apply_Operator' is deprecated. Please use 'apply_Operator_left_right  TIME  {OP_LEFT}  {OP_RIGHT}'!"<<std::endl;
      throw DummyException();
    }
    if(param.is_specified("apply_Operator_left_right")){
      for(int i=0; i<param.get_nr_rows("apply_Operator_left_right"); i++){
        param.complain_if_row_shorter("apply_Operator_left_right", 3, i, "USAGE: apply_Operator_left_right TIME {OP_LEFT}  {OP_RIGHT}");
        double time=param.get_as_double("apply_Operator_left_right", 0, i, 0);
        Eigen::MatrixXcd ML=param.get_as_operator_check("apply_Operator_left_right", i, 1);
        Eigen::MatrixXcd MR=param.get_as_operator_check("apply_Operator_left_right", i, 2);

        bool apply_before=param.get_as_bool("apply_Operator", false, i, 3);
        add_MultitimeOp(time, ML, MR, apply_before);
      }
    }
    if(param.is_specified("apply_Operator_left")){
      for(int i=0; i<param.get_nr_rows("apply_Operator_left"); i++){
        param.complain_if_row_shorter("apply_Operator_left", 2, i, "USAGE: apply_Operator_left TIME OPERATOR");
        double time=param.get_as_double("apply_Operator_left", 0, i, 0);
        Eigen::MatrixXcd M=param.get_as_operator_check("apply_Operator_left", i, 1);

        bool apply_before=param.get_as_bool("apply_Operator_left", false, i, 2);
        add_MultitimeOp(time, M, Eigen::MatrixXcd::Identity(M.rows(), M.cols()), apply_before);
      }
    }
    if(param.is_specified("apply_Operator_right")){
      for(int i=0; i<param.get_nr_rows("apply_Operator_right"); i++){
        param.complain_if_row_shorter("apply_Operator_right", 2, i, "USAGE: apply_Operator_left TIME OPERATOR");
        double time=param.get_as_double("apply_Operator_right", 0, i, 0);
        Eigen::MatrixXcd M=param.get_as_operator_check("apply_Operator_right", i, 1);

        bool apply_before=param.get_as_bool("apply_Operator_right", false, i, 2);
        add_MultitimeOp(time, Eigen::MatrixXcd::Identity(M.rows(), M.cols()), M, apply_before);
      }
    }

    //pulses with general operators
    if(param.is_specified("add_Pulse")){
      std::vector<std::vector<std::string> >  svv=param.get("add_Pulse");
      for(size_t r=0; r<svv.size(); r++){
        if(svv[r].size()<1)continue;
       
        std::pair<ComplexFunctionPtr, Eigen::MatrixXcd> res=Pulse_Selector(svv[r]);
        set_dim(res.second.rows(), "at 'add_Pulse'");
        add_Pulse(res.first, res.second);
      }
    } 
    if(param.is_specified("add_forward_Pulse")){
      std::vector<std::vector<std::string> >  svv=param.get("add_forward_Pulse");
      for(size_t r=0; r<svv.size(); r++){
        if(svv[r].size()<1)continue;
       
        std::pair<ComplexFunctionPtr, Eigen::MatrixXcd> res=Pulse_Selector(svv[r]);
        set_dim(res.second.rows(), "at 'add_forward_Pulse'");
        add_forward_Pulse(res.first, res.second);
      }
    }
    if(param.is_specified("add_backward_Pulse")){
      std::vector<std::vector<std::string> >  svv=param.get("add_backward_Pulse");
      for(size_t r=0; r<svv.size(); r++){
        if(svv[r].size()<1)continue;
       
        std::pair<ComplexFunctionPtr, Eigen::MatrixXcd> res=Pulse_Selector(svv[r]);
        set_dim(res.second.rows(), "at 'add_backward_Pulse'");
        add_backward_Pulse(res.first, res.second);
      }
    } 
    print_pulses(param);
  }

}//namespace
