#include "DynamicalMap.hpp"
#include "DummyException.hpp"
#include "BinaryReader.hpp"
#include "InfluenceFunctional_Vector.hpp" 
#include "Simulation_QUAPI.hpp"

namespace ACE{

const Eigen::MatrixXcd &DynamicalMap::get(int i)const{ 
  if(i<E.size()){
    return E[i];
  }
  if(Pade.size()!=0){
    std::cerr<<"DynamicalMap::get with Pade NOT IMPLEMENTED YET!"<<std::endl;
    throw DummyException();
  }
  std::cerr<<"DynamicalMap::get: Out of bounds ("<<i<<"/"<<E.size()<<")!"<<std::endl;
  throw DummyException();
}
Eigen::MatrixXcd DynamicalMap::invert(const Eigen::MatrixXcd &M, double regularize, std::ostream *os){
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::VectorXcd inv_sigma(svd.singularValues().rows());
    if(os!=NULL){
      *os<<svd.singularValues().transpose()<<std::endl;
    }
    for(int j=0; j<inv_sigma.rows(); j++){
      if(regularize>0.){
        double s=svd.singularValues()(j);
        inv_sigma(j)=s/(s*s+regularize*regularize);
      }else{
        inv_sigma(j)=1./svd.singularValues()(j);
      }
    }
    return (svd.matrixV()*(inv_sigma.asDiagonal())*svd.matrixU().adjoint()).eval();
}
Eigen::MatrixXcd DynamicalMap::get_dE(int i, double regularize, std::ostream *os)const{ 
  if(i<E.size()){
    if(i<1){
      if(os!=NULL){
        invert(get(0), regularize, os);//just for output
      }
      return get(0);
    }
//    Eigen::JacobiSVD<Eigen::MatrixXcd> svd(get(i-1), Eigen::ComputeFullU | Eigen::ComputeFullV);
//    Eigen::VectorXcd inv_sigma(svd.singularValues().rows());
//    for(int j=0; j<inv_sigma.rows(); j++){
//      inv_sigma(j)=1./svd.singularValues()(j);
//    }
//    return get(i)*svd.matrixV()*(inv_sigma.asDiagonal())*svd.matrixU().adjoint();
    return get(i)*invert(get(i-1), regularize, os);
  }
/*
  if(Pade.size()!=0){
    std::cerr<<"DynamicalMap::get with Pade NOT IMPLEMENTED YET!"<<std::endl;
    throw DummyException();
  }
*/
  std::cerr<<"DynamicalMap::get: Out of bounds ("<<i<<"/"<<E.size()<<")!"<<std::endl;
  throw DummyException();
}

Eigen::MatrixXcd DynamicalMap::get_dE_ref(int i, double eps, double tref)const{ 
  if(i<1){
    return get(0);
  }
  Eigen::JacobiSVD<Eigen::MatrixXcd> svd(get(i-1), Eigen::ComputeFullU | Eigen::ComputeFullV);
  int fulldim=svd.singularValues().rows();
  int newdim=0;
  for(int j=0; j<fulldim; j++){
    if(svd.singularValues()(j)>=eps)newdim++;
    else break;
  }
  if(newdim<1){
    std::cerr<<"DynamicalMap::get_dE_ref: newdim<1!"<<std::endl;
    throw DummyException();
  }
  std::cout<<"t="<<i*dt<<" fulldim="<<fulldim<<" newdim="<<newdim<<std::endl;

  Eigen::VectorXcd inv_sigma(newdim);
  for(int j=0; j<newdim; j++){
    inv_sigma(j)=1./svd.singularValues()(j);
  }
  Eigen::MatrixXcd Vtilde=svd.matrixV().block(0,0,fulldim,newdim);
  Eigen::MatrixXcd Utilde=svd.matrixU().block(0,0,fulldim,newdim);

  Eigen::MatrixXcd result=get(i)*Vtilde*inv_sigma.asDiagonal()*Utilde.adjoint();
  if(newdim==fulldim){
    return result;
  }

  int i_ref=round((tref-ta)/dt);
  if(i_ref<0||i_ref>=i-1){
    std::cerr<<"DynamicalMap::get_dE_ref: i_ref<0||i_ref>=i-1!"<<std::endl;
    throw DummyException();
  }
  std::cout<<"i_ref="<<i_ref<<std::endl;
  Eigen::MatrixXcd Uhat=svd.matrixU().block(0,newdim,fulldim,fulldim-newdim);
  Eigen::MatrixXcd Vhat=svd.matrixV().block(0,newdim,fulldim,fulldim-newdim);
  //result+=Vhat*Vhat.adjoint()*get_dE(i_ref)*Uhat*Uhat.adjoint();
  result+=get_dE(i_ref)*Uhat*Uhat.adjoint();
  return result;
}

void DynamicalMap::calculate(Propagator &prop, ProcessTensorForwardList &PT, Simulation_PT &sim, const TimeGrid &tgrid, int verbosity){
   
  ta=tgrid.ta;
  dt=tgrid.dt;

  int N=prop.get_dim(); int NL=N*N;
  if(N<2){
    std::cerr<<"DynamicalMap::calculate: N<2!"<<std::endl;
    throw DummyException();
  }
  if(tgrid.n_tot<1){
    std::cerr<<"DynamicalMap::calculate: tgrid.n_tot<2!"<<std::endl;
    throw DummyException();
  }

  //Initialize E
  {std::vector<Eigen::MatrixXcd> tmp(tgrid.n_tot, Eigen::MatrixXcd::Zero(NL,NL));E.swap(tmp);}
 
  //Calculate E
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      if(verbosity>0){
        std::cout<<"("<<i<<","<<j<<")/("<<N<<","<<N<<")"<<std::endl;
      }
      Eigen::MatrixXcd initial_rho=Eigen::MatrixXcd::Zero(N,N);
      initial_rho(i,j)=1;

      OutputPrinter extractor; 
      extractor.do_extract=true; extractor.start_extract=0; 
      extractor.full_densmat=true;
      sim.run_std(prop, PT, initial_rho, tgrid, extractor);

      for(int l=0; l<tgrid.n_tot; l++){
        E[l].col(i*N+j)=extractor.rho_t[l+1];  //Note: rho_t[0]=Identity
      }
    }
  }
}
void DynamicalMap::calculate_TEMPO(Parameters &param){
  TimeGrid                 tgrid(param);
  FreePropagator           prop(param);
  DiagBB                   diagBB(param,"Boson");
  InfluenceFunctional_Vector IF(tgrid.n_tot, tgrid.dt, diagBB);
  Eigen::MatrixXcd initial_rho=InitialState(param);
  bool use_symmetric_Trotter=param.get_as_bool("use_symmetric_Trotter", true);
  int verbosity=param.get_as_int("verbosity",1);
  bool silent=false;
  Simulation_TEMPO sim;
  sim.setup_output(param);
  sim.compressor=RankCompressor_Selector(param);

  ta=tgrid.ta;
  dt=tgrid.dt;

  int N=prop.get_dim(); int NL=N*N;
  if(N<2){
    std::cerr<<"DynamicalMap::calculate: N<2!"<<std::endl;
    throw DummyException();
  }
  if(tgrid.n_tot<1){
    std::cerr<<"DynamicalMap::calculate: tgrid.n_tot<2!"<<std::endl;
    throw DummyException();
  }

  //Initialize E
  {std::vector<Eigen::MatrixXcd> tmp(tgrid.n_tot, Eigen::MatrixXcd::Zero(NL,NL));E.swap(tmp);}
 
  //Calculate E
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      if(verbosity>0){
        std::cout<<"("<<i<<","<<j<<")/("<<N<<","<<N<<")"<<std::endl;
      }
      Eigen::MatrixXcd initial_rho=Eigen::MatrixXcd::Zero(N,N);
      initial_rho(i,j)=1;
  
      sim.output_Op.ops=std::vector<Eigen::MatrixXcd>(N*N, Eigen::MatrixXcd::Zero(N,N));
      for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
          sim.output_Op.ops[i*N+j](j*N+i)=1.;
        }
      }  
      sim.run(prop, IF, tgrid.ta, tgrid.dt, tgrid.get_t_tot(), initial_rho, silent, use_symmetric_Trotter);

      for(int l=0; l<tgrid.n_tot; l++){
        for(int i1=0; i1<N; i1++){
          for(int j1=0; j1<N; j1++){
            E[l](i1*N+j1, i*N+j)=sim.results[l+1].second[i1*N+j1];
          }
        } 
      }
    }
  }
}
void DynamicalMap::calculate(Parameters &param){
  TimeGrid                 tgrid(param);
  ProcessTensorForwardList PT(param);
  FreePropagator           prop(param);
  Simulation_PT            sim(param);
 
  calculate(prop, PT, sim, tgrid, 1);
}
void DynamicalMap::calculate_stepwise(Parameters &param){
  TimeGrid                 tgrid(param);
  FreePropagator           prop(param);
  Simulation_PT            sim(param);
 
  ta=tgrid.ta;
  dt=tgrid.dt;
  int n_tot=tgrid.n_tot;
  {std::vector<Eigen::MatrixXcd> E2(n_tot); E.swap(E2);}
  for(int n=1; n<=n_tot; n++){
    Parameters param2=param;
    param2.override_param("te", tgrid.get_t(n));
    ProcessTensorForwardList PT(param2);
    TimeGrid tgrid2(param2);

    DynamicalMap other; 
    other.calculate(prop, PT, sim, tgrid2, 1);
    E[n-1]=other.E[n-1];
  }
}

void DynamicalMap::make_physical_dE(Eigen::MatrixXcd & dE){
    //check if general hermitian density matrices are mapped to something hermitian:
    int DL=dE.rows(); int D=sqrt(DL);
    for(int i=0; i<D; i++){
      Eigen::MatrixXcd rho = L_Vector_to_H_Matrix( dE.col(i*D+i) );
      Eigen::MatrixXcd phys = 0.5*(rho+rho.adjoint()); 
      dE.col(i*D+i) = H_Matrix_to_L_Vector( phys/phys.trace() );
      dE.col(i*D+i) = H_Matrix_to_L_Vector( phys );
    }
    for(int i=0; i<D; i++){
      for(int j=i+1; j<D; j++){
        Eigen::MatrixXcd rho1 = L_Vector_to_H_Matrix( dE.col(i*D+j) );
        Eigen::MatrixXcd rho2 = L_Vector_to_H_Matrix( dE.col(j*D+i) );
        Eigen::MatrixXcd phys1 = 0.5*(rho1+rho2.adjoint());
        Eigen::MatrixXcd phys2 = phys1.adjoint(); //0.5*(rho2+rho1.adjoint());
//        dE.col(i*D+j) = H_Matrix_to_L_Vector( phys1 );
//        dE.col(j*D+i) = H_Matrix_to_L_Vector( phys2 );
        dE.col(i*D+j) = H_Matrix_to_L_Vector( phys1 - phys1.trace()/(double)D*Eigen::MatrixXcd::Identity(D,D) );
        dE.col(j*D+i) = H_Matrix_to_L_Vector( phys2 - phys2.trace()/(double)D*Eigen::MatrixXcd::Identity(D,D) );
      }
    }
}

void DynamicalMap::make_physical(){
 std::vector<Eigen::MatrixXcd> E2(E.size());
   for(size_t l=0; l<E.size(); l++){
    Eigen::MatrixXcd dE=get_dE(l);
    make_physical_dE(dE);
    if(l<1){
      E2[0]=dE;
    }else{
      E2[l]=E2[l-1]*dE;
    }
  } 
  E.swap(E2);
}
void DynamicalMap::propagate(const TimeGrid &tgrid, const Eigen::MatrixXcd &initial_state,  OutputPrinter &printer)const{
  Eigen::VectorXcd rho=H_Matrix_to_L_Vector(initial_state);
  
  int n_tot=tgrid.n_tot;
  if((int)E.size()<n_tot){n_tot=E.size();}
  
  if(initial_state.rows()*initial_state.cols()!=get(0).cols()){
    std::cerr<<"initial_state.rows()="<<initial_state.rows()<<std::endl;
    std::cerr<<"versus get(0).cols="<<get(0).cols()<<std::endl;
    throw DummyException();
  }

  printer.print(0, tgrid.ta, rho);
  for(int n=0; n<n_tot; n++){
    if(get(n).cols()!=rho.rows()){
      std::cerr<<"DynamicalMap::propagate: E.cols()!=rho.rows()!"<<std::endl;
      throw DummyException();
    } 
    Eigen::VectorXcd rho2 = get(n)*rho;
    printer.print(n+1, tgrid.get_t(n+1), rho2);
  }
}

void DynamicalMap::propagate_dE(const TimeGrid &tgrid, const Eigen::MatrixXcd &initial_state,  OutputPrinter &printer)const{
  Eigen::VectorXcd rho=H_Matrix_to_L_Vector(initial_state);
  
  int n_tot=tgrid.n_tot;
  if((int)E.size()<n_tot){n_tot=E.size();}

  printer.print(0, tgrid.ta, rho);
  for(int n=0; n<n_tot; n++){
    if(get(n).cols()!=rho.rows()){
      std::cerr<<"DynamicalMap::propagate: dE.cols()!=rho.rows()!"<<std::endl;
      throw DummyException();
    } 
    Eigen::VectorXcd rho2 = get_dE(n)*rho;
    rho=rho2;
    printer.print(n+1, tgrid.get_t(n+1), rho);
  }
}

Eigen::MatrixXcd DynamicalMap::get_Liouvillian(int i, double regularize)const{ 
  if(i<1){
    return (get(0)).log()/dt;
  }
  return (get_dE(i,regularize).log()+get_dE(i-1,regularize).log())/(2.*dt);
}


LindbladMasterEquation DynamicalMap::get_LindbladMasterEquation_Hall(int l)const{
  LindbladMasterEquation LME;
  LME.set_from_Liouvillian_Hall(get_Liouvillian(l), 0,0);
  return LME;
}
LindbladMasterEquation DynamicalMap::get_LindbladMasterEquation(int l)const{
  LindbladMasterEquation LME;
  LME.set_from_Liouvillian(get_Liouvillian(l), 0,0);
  return LME;
}

void DynamicalMap::print_Lindblad_Hall(const std::string & prefix)const{
  std::ofstream ofs(prefix.c_str());
  for(size_t l=0; l<E.size(); l++){
    LindbladMasterEquation LME=get_LindbladMasterEquation_Hall(l);
    ofs<<ta+l*dt;
    for(size_t k=0; k<LME.L.size(); k++){
      ofs<<" "<<LME.L[k].first;
    }
    ofs<<std::endl;
  } 
}
void DynamicalMap::print_Lindblad(const std::string & prefix)const{
  std::ofstream ofs(prefix.c_str());
  for(size_t l=0; l<E.size(); l++){
    LindbladMasterEquation LME=get_LindbladMasterEquation(l);
    ofs<<ta+l*dt;
    for(size_t k=0; k<LME.L.size(); k++){
      ofs<<" "<<LME.L[k].first;
    }
    ofs<<std::endl;
  } 
}

void DynamicalMap::print_normdiff(std::ostream &os)const{
  for(size_t l=0; l<E.size(); l++){
    Eigen::MatrixXcd dE=get_dE(l);
    Eigen::MatrixXcd dEref;
    if(l<1){
      dEref=Eigen::MatrixXcd::Identity(dE.rows(), dE.cols());
    }else{
      dEref=get_dE(l-1);
    }
    os << ta+l*dt <<" "<< (dE-dEref).norm()/dt << std::endl;
  }
}
void DynamicalMap::print_TT_normdiff(std::ostream &os)const{
  std::vector<Eigen::MatrixXcd> TT=get_TT();
  for(size_t l=0; l<TT.size(); l++){
    Eigen::MatrixXcd TTref;
    if(l<1){
      TTref=Eigen::MatrixXcd::Identity(TT[l].rows(), TT[l].cols());
    }else{
      TTref=TT[l-1];
    }
    os << ta+l*dt <<" "<< (TT[l]-TTref).norm()/dt << std::endl;
  }
}
void DynamicalMap::print_TT_norm(std::ostream &os)const{
  std::vector<Eigen::MatrixXcd> TT=get_TT();
  for(size_t l=0; l<TT.size(); l++){
    os << ta+l*dt <<" "<< TT[l].norm()/dt << std::endl;
  }
}
void DynamicalMap::print_eigenvalues(std::ostream &os, double regularize)const{
  for(size_t l=0; l<E.size(); l++){
//std::cout<<"print_eigenvalues: "<<l<<"/"<<E.size()<<std::endl;
    Eigen::MatrixXcd Liou=get_Liouvillian(l, regularize);
    int DL=Liou.rows();
    Eigen::ComplexSchur<Eigen::MatrixXcd> schur(Liou);
    Eigen::VectorXcd lambda=schur.matrixT().diagonal();
    std::vector<std::complex<double> > evals(lambda.rows());
    for(int i=0; i<DL; i++){
      evals[i]=lambda(i);
    }
    std::sort(evals.begin(), evals.end(), 
      [](const std::complex<double> &a, const std::complex<double> &b){
      return a.real()<b.real(); } );

    os<<ta+l*dt;
    for(int i=0; i<DL; i++){
      os<<" "<<evals[i].real()<<" "<<evals[i].imag();
    }
    os<<std::endl;
  }
}
void DynamicalMap::print_dE_eigenvalues(std::ostream &os)const{
  for(size_t l=0; l<E.size(); l++){
    Eigen::MatrixXcd dE=get_dE(l);
    int DL=dE.rows();
    Eigen::ComplexSchur<Eigen::MatrixXcd> schur(dE);
    Eigen::VectorXcd lambda=schur.matrixT().diagonal();
    std::vector<std::complex<double> > evals(lambda.rows());
    for(int i=0; i<DL; i++){
      evals[i]=lambda(i);
    }
    std::sort(evals.begin(), evals.end(), 
      [](const std::complex<double> &a, const std::complex<double> &b){
      return a.real()<b.real(); } );

    os<<ta+l*dt;
    for(int i=0; i<DL; i++){
      os<<" "<<evals[i].real()<<" "<<evals[i].imag();
    }
    os<<std::endl;
  }
}
void DynamicalMap::print_E_eigenvalues(std::ostream &os)const{
  for(size_t l=0; l<E.size(); l++){
    Eigen::MatrixXcd dE=get(l);
    int DL=dE.rows();
    Eigen::ComplexSchur<Eigen::MatrixXcd> schur(dE);
    Eigen::VectorXcd lambda=schur.matrixT().diagonal();
    std::vector<std::complex<double> > evals(lambda.rows());
    for(int i=0; i<DL; i++){
      evals[i]=lambda(i);
    }
    std::sort(evals.begin(), evals.end(), 
      [](const std::complex<double> &a, const std::complex<double> &b){
      return a.real()<b.real(); } );

    os<<ta+l*dt;
    for(int i=0; i<DL; i++){
      os<<" "<<evals[i].real()<<" "<<evals[i].imag();
    }
    os<<std::endl;
  }
}
void DynamicalMap::print_eigenvalues_ref(std::ostream &os, double eps, double tref)const{
  for(size_t l=0; l<E.size(); l++){
    Eigen::MatrixXcd Liou=get_dE_ref(l,eps,tref).log()/dt;
    int DL=Liou.rows();
    Eigen::ComplexSchur<Eigen::MatrixXcd> schur(Liou);
    Eigen::VectorXcd lambda=schur.matrixT().diagonal();
    std::vector<std::complex<double> > evals(lambda.rows());
    for(int i=0; i<DL; i++){
      evals[i]=lambda(i);
    }
    std::sort(evals.begin(), evals.end(), 
      [](const std::complex<double> &a, const std::complex<double> &b){
      return a.real()<b.real(); } );

    os<<ta+l*dt;
    for(int i=0; i<DL; i++){
      os<<" "<<evals[i].real()<<" "<<evals[i].imag();
    }
    os<<std::endl;
  }
}
void DynamicalMap::print_dE_eigenvalues_ref(std::ostream &os, double eps, double tref)const{
  for(size_t l=0; l<E.size(); l++){
    Eigen::MatrixXcd dE=get_dE_ref(l,eps,tref);
    int DL=dE.rows();
    Eigen::ComplexSchur<Eigen::MatrixXcd> schur(dE);
    Eigen::VectorXcd lambda=schur.matrixT().diagonal();
    std::vector<std::complex<double> > evals(lambda.rows());
    for(int i=0; i<DL; i++){
      evals[i]=lambda(i);
    }
    std::sort(evals.begin(), evals.end(), 
      [](const std::complex<double> &a, const std::complex<double> &b){
      return a.real()<b.real(); } );

    os<<ta+l*dt;
    for(int i=0; i<DL; i++){
      os<<" "<<evals[i].real()<<" "<<evals[i].imag();
    }
    os<<std::endl;
  }
}
void DynamicalMap::print_fullE_eigenvalues(std::ostream &os)const{
  for(size_t l=0; l<E.size(); l++){
    Eigen::MatrixXcd Liou=(get(l)).log()/((l+1)*dt);
    int DL=Liou.rows();
    Eigen::ComplexSchur<Eigen::MatrixXcd> schur(Liou);
    Eigen::VectorXcd lambda=schur.matrixT().diagonal();
    std::vector<std::complex<double> > evals(lambda.rows());
    for(int i=0; i<DL; i++){
      evals[i]=lambda(i);
    }
    std::sort(evals.begin(), evals.end(), 
      [](const std::complex<double> &a, const std::complex<double> &b){
      return a.real()<b.real(); } );

    os<<ta+l*dt;
    for(int i=0; i<DL; i++){
      os<<" "<<evals[i].real()<<" "<<evals[i].imag();
    }
    os<<std::endl;
  }
}
void DynamicalMap::print_dE_singularvalues(std::ostream &os)const{
  for(size_t l=0; l<E.size(); l++){
    os<<ta+l*dt<<" ";
    invert(get_dE(l),0,&os);
  }
}
void DynamicalMap::print_singularvalues(std::ostream &os)const{
  for(size_t l=0; l<E.size(); l++){
    os<<ta+l*dt<<" ";
    invert(get(l),0,&os);
  }
}

Eigen::MatrixXcd DynamicalMap::get_average_dE(double tstart, double tend)const{
  int nstart=round((tstart-ta)/dt);
  int nend=round((tend-ta)/dt);
  if(nstart<0)nstart=0;
  if(nend>(int)E.size())nend=E.size();
  if(nstart>=nend){
    std::cerr<<"DynamicalMap::print_average_eigenvalues: nstart>=nend!"<<std::endl;
    throw DummyException();
  }

  int DL=get(0).rows();
  Eigen::MatrixXcd dE = Eigen::MatrixXcd::Zero(DL,DL);
  for(int l=nstart; l<nend; l++){
    dE += get_dE(l) / ((double)(nend-nstart));
  }
  return dE;
}

void DynamicalMap::Richardson_combine(const DynamicalMap &other){
  int step_ratio = round(dt/other.dt);
  if(step_ratio<2){
    std::cerr<<"DynamicalMap::Richardson_combine step_ratio<2!"<<std::endl;
  }

  if(other.E.size()<E.size()*step_ratio){
    std::cerr<<"DynamicalMap::Richardson_combine other.E.size()<E.size()*step_ratio!"<<std::endl;
    throw DummyException();
  }
  int DL=get(0).rows();
  std::vector<Eigen::MatrixXcd> dE2(E.size(), Eigen::MatrixXcd::Identity(DL,DL));
  for(size_t l=0; l<E.size(); l++){
//    Eigen::MatrixXcd E2=(pow(step_ratio,2)*other.get(l*step_ratio)-get(l)) /
//                        (pow(step_ratio,2)-1.);
//    E[l]=E2;

    for(int i=0; i<step_ratio; i++){
      dE2[l]=other.get_dE(l*step_ratio+i)*dE2[l];
    }
    dE2[l]=(pow(step_ratio,2)*dE2[l]-get_dE(l))/(pow(step_ratio,2)-1.);
  } 
  for(size_t l=0; l<E.size(); l++){
    Eigen::MatrixXcd Eref=Eigen::MatrixXcd::Identity(DL,DL);
    if(l>0)Eref=E[l-1];
    E[l]=dE2[l]*Eref;
  }
}

void DynamicalMap::print_average_eigenvalues(std::ostream &os, double tstart, double tend)const{

  Eigen::MatrixXcd Liou = get_average_dE(tstart,tend).log()/dt;
  int DL=Liou.rows();
  Eigen::ComplexSchur<Eigen::MatrixXcd> schur(Liou);
  Eigen::VectorXcd lambda=schur.matrixT().diagonal();
  std::vector<std::complex<double> > evals(lambda.rows());
  for(int i=0; i<DL; i++){
    evals[i]=lambda(i);
  }
  std::sort(evals.begin(), evals.end(), 
      [](const std::complex<double> &a, const std::complex<double> &b){
      return a.real()<b.real(); } );

  os<<tstart;
  for(int i=0; i<DL; i++){
    os<<" "<<evals[i].real()<<" "<<evals[i].imag();
  }
  os<<std::endl;
}

void DynamicalMap::analyze_physicality(std::ostream &os)const{
  for(size_t l=0; l<E.size(); l++){
    os<<"E["<<l<<"] (t="<<ta+l*dt<<"):"<<std::endl;
    Eigen::MatrixXcd dE=get_dE(l);
    //check if general hermitian density matrices are mapped to something hermitian:
    int DL=dE.rows(); int D=sqrt(DL);
    for(int i=0; i<D; i++){
      Eigen::VectorXcd rho_in=Eigen::VectorXcd::Zero(DL);
      rho_in(i*D+i)=1.;
      Eigen::MatrixXcd rho_out=L_Vector_to_H_Matrix(dE*rho_in);
      os<<"|rho_out-rho_out.adjoint()|="<<(rho_out-rho_out.adjoint()).norm()<<std::endl;
      os<<"rho_out.trace()-1="<<rho_out.trace()-1.<<std::endl;
    }
    for(int i=0; i<D; i++){
      for(int j=i+1; j<D; j++){
        Eigen::VectorXcd rho_in=Eigen::VectorXcd::Zero(DL);
        rho_in(i*D+i)=0.5;
        rho_in(j*D+j)=0.5;
        rho_in(i*D+j)=0.5;
        rho_in(j*D+i)=0.5;
        Eigen::MatrixXcd rho_out=L_Vector_to_H_Matrix(dE*rho_in);
        os<<"|rho_out-rho_out.adjoint()|="<<(rho_out-rho_out.adjoint()).norm()<<std::endl;
        os<<"rho_out.trace()-1="<<rho_out.trace()-1.<<std::endl;

        rho_in(i*D+j)=std::complex<double>(0.,-0.5);
        rho_in(j*D+i)=std::complex<double>(0.,0.5);
        rho_out=L_Vector_to_H_Matrix(dE*rho_in);
        os<<"|rho_out-rho_out.adjoint()|="<<(rho_out-rho_out.adjoint()).norm()<<std::endl;
        os<<"rho_out.trace()-1="<<rho_out.trace()-1.<<std::endl;
      }
    }
  }
}

std::vector<Eigen::MatrixXcd> DynamicalMap::get_TT(int length)const{
  if(length<1)length=E.size();
  std::vector<Eigen::MatrixXcd> TT(length);
  for(int l=0; l<length; l++){
    TT[l]=E[l];
    for(int k=0; k<l; k++){
      TT[l]-=TT[l-k-1]*E[k];
    }
  }
  return TT;
}
void DynamicalMap::extrapolate_TT(double te2, int length){
  std::vector<Eigen::MatrixXcd> TT=get_TT(length);
  length=TT.size();
  int DL=TT[0].rows(); int D=sqrt(DL);
    
  int n_new=round( (te2-ta)/dt );
  {std::vector<Eigen::MatrixXcd> E2(n_new, Eigen::MatrixXcd::Zero(DL,DL));E.swap(E2);}
  for(int i=0; i<D; i++){
    for(int j=0; j<D; j++){
      std::vector<Eigen::VectorXcd> rho(length, Eigen::VectorXcd::Zero(DL));
      rho[0](i*D+j)=1;
      for(int n=0; n<n_new; n++){
        Eigen::VectorXcd rho_next=Eigen::VectorXcd::Zero(DL);
        for(int l=0; l<length; l++){
          rho_next+=TT[l]*rho[l]; 
        } 
        for(int l=length-1; l>0; l--){
          rho[l]=rho[l-1];
        }
        rho[0]=rho_next;
        E[n].col(i*D+j)=rho_next;
      }           
    }
  }
}

void DynamicalMap::calculate_Pade(const std::vector<int> &steps){
  int dim=get(0).rows();
  int S=steps.size();


  Eigen::MatrixXcd id=Eigen::MatrixXcd::Identity(dim,dim);
  Eigen::MatrixXcd bigMat=Eigen::MatrixXcd::Zero(S*dim, 3*dim);
  Eigen::MatrixXcd bigB=Eigen::MatrixXcd::Zero(S*dim, dim);
  for(int s=0; s<S; s++){
    Eigen::MatrixXcd dET=get_dE(steps[s]).transpose();  
    bigMat.block(s*dim,  0,    dim,dim)=id;
    bigMat.block(s*dim,  dim,  dim,dim)=steps[s]*id;
    bigMat.block(s*dim,  2*dim,dim,dim)=-steps[s]*dET;

    bigB.block(s*dim,   0, dim, dim)=dET;
  }
  
  Eigen::JacobiSVD<Eigen::MatrixXcd> svd( bigMat , Eigen::ComputeFullU | Eigen::ComputeFullV);
  std::cout<<"Pade: singular values: "<<svd.singularValues().transpose()<<std::endl;

  Eigen::MatrixXcd bigX = svd.solve(bigB);

  Pade = std::vector<Eigen::MatrixXcd>(3, Eigen::MatrixXcd::Zero(dim,dim));
  for(int d=0; d<dim; d++){
    Pade[0] = bigX.block(0,    0,dim,dim).transpose();
    Pade[1] = bigX.block(dim,  0,dim,dim).transpose();
    Pade[2] = bigX.block(2*dim,0,dim,dim).transpose();
  }
}
void DynamicalMap::calculate_Pade2(const std::vector<int> &steps){
  int dim=get(0).rows();
  int S=steps.size();

  Eigen::MatrixXcd id=Eigen::MatrixXcd::Identity(dim,dim);
  Eigen::MatrixXcd bigMat=Eigen::MatrixXcd::Zero(S*dim, 5*dim);
  Eigen::MatrixXcd bigB=Eigen::MatrixXcd::Zero(S*dim, dim);
  for(int s=0; s<S; s++){
    int n=steps[s];
    Eigen::MatrixXcd dET=get_dE(n).transpose();  
    bigMat.block(s*dim,  0,    dim,dim)=id;
    bigMat.block(s*dim,  dim,  dim,dim)=n*id;
    bigMat.block(s*dim,  2*dim,dim,dim)=n*n*id;
    bigMat.block(s*dim,  3*dim,dim,dim)=-n*dET;
    bigMat.block(s*dim,  4*dim,dim,dim)=-n*n*dET;

    bigB.block(s*dim,   0, dim, dim)=dET;
  }
  
  Eigen::JacobiSVD<Eigen::MatrixXcd> svd( bigMat , Eigen::ComputeFullU | Eigen::ComputeFullV);
  std::cout<<"Pade: singular values: "<<svd.singularValues().transpose()<<std::endl;

  Eigen::MatrixXcd bigX = svd.solve(bigB);

  Pade = std::vector<Eigen::MatrixXcd>(5, Eigen::MatrixXcd::Zero(dim,dim));
  for(int d=0; d<dim; d++){
    Pade[0] = bigX.block(0,    0,dim,dim).transpose();
    Pade[1] = bigX.block(dim,  0,dim,dim).transpose();
    Pade[2] = bigX.block(2*dim,0,dim,dim).transpose();
    Pade[3] = bigX.block(3*dim,0,dim,dim).transpose();
    Pade[4] = bigX.block(4*dim,0,dim,dim).transpose();
  }
  
}
void DynamicalMap::Pade_extrapolate(int add_steps){
std::cout<<"DynamicalMap: NOT IMPLEMENTED YET!"<<std::endl;
  throw DummyException();
}

void DynamicalMap::read(const std::string &fname){
 try{
  std::ifstream ifs(fname.c_str());
  if(!file_exists(fname)){ 
    std::cerr<<"DynamicalMap::read: File '"<<fname<<"' cannot be read!"<<std::endl;
    throw DummyException();
  }
  std::string str=binary_read_fixedSizeString(ifs, 4, "DynamicalMap");
  if(str!="DM00"){
    std::cerr<<"File \""<<fname<<"\" does not seem to contain a DynamicalMap!"<<std::endl;
    throw DummyException();
  }
  ta=binary_read<double>(ifs, "DynamicalMap");
  dt=binary_read<double>(ifs, "DynamicalMap");

  size_t sz=binary_read<size_t>(ifs, "DynamicalMap");
  {std::vector<Eigen::MatrixXcd> tmp(sz); E.swap(tmp);}
  for(size_t s=0; s<E.size(); s++){
    E[s]=binary_read_EigenMatrixXcd(ifs, "DynamicalMap");
  }

  sz=binary_read<size_t>(ifs, "DynamicalMap");
  {std::vector<Eigen::MatrixXcd> tmp(sz); Pade.swap(tmp);}
  for(size_t s=0; s<Pade.size(); s++){
    Pade[s]=binary_read_EigenMatrixXcd(ifs, "DynamicalMap");
  }

 }catch(std::exception &e){
   std::cerr<<"called from DynamicalMap::read("<<fname<<")"<<std::endl;
   throw e;
 }
}

void DynamicalMap::write(const std::string &fname)const{
  std::ofstream ofs(fname.c_str());
  binary_write_fixedSizeString(ofs, 4, "DM00");
  binary_write<double>(ofs, ta);
  binary_write<double>(ofs, dt);

  binary_write<size_t>(ofs, E.size());
  for(size_t s=0; s<E.size(); s++){
    binary_write_EigenMatrixXcd(ofs, E[s]);
  }

  binary_write<size_t>(ofs, Pade.size());
  for(size_t s=0; s<Pade.size(); s++){
    binary_write_EigenMatrixXcd(ofs, Pade[s]);
  }
}

}//namespace
