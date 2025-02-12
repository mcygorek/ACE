#include "Largest_EV.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include "DummyException.hpp"

namespace ACE{

std::complex<double> Largest_EV_Arnoldi(Eigen::VectorXcd &vec, int m, 
      double epsilon,
      std::function<Eigen::VectorXcd(const Eigen::VectorXcd &v)> Afunc, 
      int verbosity, bool reortho){

  if(verbosity>1)std::cout<<"m="<<m<<" epsilon="<<epsilon<<std::endl;
  Eigen::MatrixXcd H=Eigen::MatrixXcd::Zero(m, m);
  Eigen::MatrixXcd Q=Eigen::MatrixXcd::Zero(vec.size(),m);
  Q.col(0)=vec; Q.col(0).normalize();
  for(int k=1; k<m; k++){
//    Q.col(k-1).normalize();
    Q.col(k)=Afunc(Q.col(k-1));
    for(int j=0; j<k; j++){
      H(j,k-1)=Q.col(j).dot(Q.col(k));
      Q.col(k)-=H(j,k-1)*Q.col(j);
    }
    if(reortho){
      for(int j=0; j<k; j++){
        std::complex<double> overlap=Q.col(j).dot(Q.col(k));
        Q.col(k)-=overlap*Q.col(j);
        H(j,k-1)+=overlap;
      }
    }

    H(k,k-1)=Q.col(k).norm(); 
    Q.col(k)/=H(k,k-1);
 
    if(H(k,k-1).real()<epsilon){
      if(verbosity>0){
        std::cout<<"Largest_EV_Arnoldi converged to "<<H(k,k-1).real()<<"<"<<epsilon<<" at iteration "<<k<<"/"<<m<<std::endl;
      }
      m=k;
      Eigen::MatrixXcd tmp = H.block(0,0,m,m); H.swap(tmp);
      break;
    }
    if(verbosity>1){
      std::cout<<"Arnoldi k="<<k<<" H(k,k-1)="<<H(k,k-1)<<std::endl;
    }

  }

  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solveH(H);
//  int imax=0; for(int i=1; i<m; i++){if(std::abs(solveH.eigenvalues()(i))>std::abs(solveH.eigenvalues()(imax))){imax=i;}}
    int imax=0; for(int i=1; i<m; i++){if(std::real(solveH.eigenvalues()(i))>std::real(solveH.eigenvalues()(imax))){imax=i;}}

  vec=Eigen::VectorXcd::Zero(vec.size());
  for(int i=0; i<m; i++){
    vec+=solveH.eigenvectors()(i, imax)*Q.col(i);
  }
//  return vec.dot(Afunc(vec));  //Ritz value
  return solveH.eigenvalues()(imax);
}


std::complex<double> Largest_EV_Arnoldi_BLAS(Eigen::VectorXcd &vec, int m, 
      double epsilon,
      std::function<Eigen::VectorXcd(const Eigen::VectorXcd &v)> Afunc, 
      int verbosity){

  int dim=vec.size();
  if(dim<1){ 
    std::cerr<<"Largest_EV_Arnoldi_BLAS called with empty vec!"<<std::endl;
    throw DummyException();
  }
  if(dim==1){
    vec(0)=1;
    std::complex<double> max_eval=Afunc(vec)(0);
    if(!std::isfinite(std::real(max_eval)) || !std::isfinite(std::imag(max_eval))){
      std::cerr<<"Arnoldi: max_eval="<<max_eval<<std::endl;
      throw DummyException();
    }
    return max_eval;
  }
  if(m>dim){
    m=dim;
    //Direct diagonalization:
    Eigen::MatrixXcd Afull(dim,dim);
    Eigen::VectorXcd v(dim);
    for(int i=0; i<dim; i++){
      if(i==0){
        v=Eigen::VectorXcd::Zero(dim);
      }else{
        v(i-1)=0.;
      }
      v(i)=1;
      Afull.col(i)=Afunc(v);
    }

    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solveH(Afull);
    int imax=0; for(int i=1; i<m; i++){if(std::real(solveH.eigenvalues()(i))>std::real(solveH.eigenvalues()(imax))){imax=i;}}
    vec=solveH.eigenvectors().col(imax);
    return solveH.eigenvalues()(imax);
  }

if(verbosity>0){std::cout<<"Arnoldi: dim="<<dim<<" m="<<m<<std::endl;}
  Eigen::MatrixXcd H=Eigen::MatrixXcd::Zero(m, m);
  Eigen::MatrixXcd Q=Eigen::MatrixXcd::Zero(dim,m);
  Q.col(0)=vec; Q.col(0).normalize();
  Eigen::VectorXcd r=Afunc(Q.col(0));
  Eigen::VectorXcd c(dim);
  H(0,0)=Q.col(0).dot(r);
  r-=H(0,0)*Q.col(0);
  for(int j=0; j<m-1; j++){
    H(j+1, j)=r.norm();

    if(H(j+1,j).real()<epsilon){
      if(verbosity>0){
        std::cout<<"Largest_EV_Arnoldi_BLAS converged to "<<H(j+1,j).real()<<"<"<<epsilon<<" at iteration "<<j+1<<"/"<<m<<std::endl;
      }
      m=j+1;
      Eigen::MatrixXcd tmp = H.block(0,0,m,m); H.swap(tmp);
      break;
    }
    if(verbosity>1){
      std::cout<<"Arnoldi j="<<j<<" H(j+1,j)="<<H(j+1,j)<<std::endl;
    }

    Q.col(j+1)=r/H(j+1, j);
    r=Afunc(Q.col(j+1));
    H.block(0, j+1, j+2, 1)=Q.block(0,0,dim,j+2).adjoint()*r;
    r-=Q.block(0,0,dim,j+2)*H.block(0,j+1, j+2, 1);

    c.head(j+2)=Q.block(0,0,dim,j+2).adjoint()*r; 
    r-=Q.block(0,0,dim,j+2)*c.head(j+2);
    H.block(0,j+1, j+2, 1)+=c.head(j+2);
  }

  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solveH(H);
  int imax=0; for(int i=1; i<m; i++){if(std::real(solveH.eigenvalues()(i))>std::real(solveH.eigenvalues()(imax))){imax=i;}}
  vec=Eigen::VectorXcd::Zero(vec.size());
  for(int i=0; i<m; i++){
    vec+=solveH.eigenvectors()(i, imax)*Q.col(i);
  }
  std::complex<double> max_eval=solveH.eigenvalues()(imax);
/*
  std::cout<<std::boolalpha<<std::isfinite(std::real(max_eval))<<" "<<std::isfinite(std::imag(max_eval))<<std::endl;
  if(!std::isfinite(std::real(max_eval)) || !std::isfinite(std::imag(max_eval))){
    std::cerr<<"H="<<H<<std::endl;   
    std::cerr<<"vec="<<vec.transpose()<<std::endl;
    std::cerr<<"Arnoldi: max_eval!=max_eval"<<std::endl;
    throw DummyException();
  }
std::cout<<"Arnoldi: max_eval="<<max_eval<<std::endl;
*/
  return max_eval;
}

std::complex<double> Largest_EV_KrylovSchur(Eigen::VectorXcd &vec, 
      int maxiter, int m, int k, double epsilon,
      std::function<Eigen::VectorXcd(const Eigen::VectorXcd &v)> Afunc, 
      int verbosity){

  if(k>=m){
    std::cerr<<"Largest_EV_KrylovSchur: k>=m!"<<std::endl;
    throw DummyException();
  }
  if(k<1){
    std::cerr<<"Largest_EV_KrylovSchur: k<1!"<<std::endl;
    throw DummyException();
  }

  int dim=vec.size();
  int m_orig=m;
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solveH;
//  Eigen::ComplexSchur<Eigen::MatrixXcd> Schur;
  std::vector<std::pair<std::complex<double>,int> > srt(m);

  Eigen::MatrixXcd H=Eigen::MatrixXcd::Zero(m, m);
  Eigen::MatrixXcd Q=Eigen::MatrixXcd::Zero(dim,m);
  Q.col(0)=vec; Q.col(0).normalize();
  Eigen::VectorXcd r(dim);
  Eigen::VectorXcd c(dim);
  r=Afunc(Q.col(0));
  H(0,0)=Q.col(0).dot(r);
  r-=H(0,0)*Q.col(0);
  int start=0;
  bool stopped=false;
  for(int outeriter=0; outeriter<((maxiter+m-1)/m); outeriter++){
    //Arnoldi to fill up to size m:
    for(int j=start; j<m-1; j++){
      H(j+1, j)=r.norm();

      if(H(j+1,j).real()<epsilon){
        if(verbosity>0){
          std::cout<<"Largest_EV_KrylovSchur converged to "<<H(j+1,j).real()<<"<"<<epsilon<<" at iteration "<<j+m*outeriter+1<<"/"<<maxiter<<std::endl;
        }
        m=j+1;
        Eigen::MatrixXcd tmp = H.block(0,0,m,m); H.swap(tmp);
        stopped=true;
        break;
      }
      if(verbosity>1){
        std::cout<<"KrylovSchur: outeriter="<<outeriter<<" j="<<j<<" H(j+1,j)="<<H(j+1,j)<<std::endl;
      }

      Q.col(j+1)=r/H(j+1, j);
      r=Afunc(Q.col(j+1));
      H.block(0, j+1, j+2, 1)=Q.block(0,0,dim,j+2).adjoint()*r;
      r-=Q.block(0,0,dim,j+2)*H.block(0,j+1, j+2, 1);

      c.head(j+2)=Q.block(0,0,dim,j+2).adjoint()*r; 
      r-=Q.block(0,0,dim,j+2)*c.head(j+2);
      H.block(0,j+1, j+2, 1)+=c.head(j+2);
    }

    //Schur form and keep the k dominant pairs: For now, just diagonalize
//    Schur.compute(H);
    solveH.compute(H);
    srt.resize(m);
    for(int i=0; i<m; i++){srt[i]=std::make_pair(solveH.eigenvalues()(i), i);}
    std::sort(srt.begin(),srt.end(), [](const std::pair<std::complex<double>,int> &a, const std::pair<std::complex<double>,int> &b){ return std::real(a.first)>std::real(b.first); });
   
    if(verbosity>1){
      std::cout<<"Largest "<<k<<" eigenvalue approximations:";
      for(int i=0; i<k; i++){ std::cout<<" "<<srt[i].first; } 
      std::cout<<std::endl;
      for(int i=0; i<k; i++){ std::cout<<" "<<srt[i].second; } 
      std::cout<<std::endl;
    }
    Eigen::MatrixXcd newQ(m, k);
    for(int i=0; i<k; i++){
      newQ.col(i)=solveH.eigenvectors().col(srt[i].second); 
      for(int j=0; j<i; j++){
        newQ.col(i)-=(newQ.col(j).dot(newQ.col(i)))*newQ.col(j);
      }
      newQ.col(i).normalize();
    }
    Eigen::MatrixXcd tmp=newQ.adjoint()*H*newQ;
//      H(i,i)=srt[i].first;
    H=Eigen::MatrixXcd::Zero(m,m);
    H.block(0,0,k,k)=tmp;
//    H=newQ.adjoint()*H*newQ;
//    if(verbosity>1){std::cout<<"Q2^*Q2="<<std::endl<<newQ.adjoint()*newQ<<std::endl;}
//    if(verbosity>1){std::cout<<"Q2Q2^*="<<std::endl<<newQ*newQ.adjoint()<<std::endl;}
//    if(verbosity>1){std::cout<<"H="<<std::endl<<H.block(0,0,k,k)<<std::endl;}
    tmp=Q*newQ;
    Q=Eigen::MatrixXcd::Zero(dim,m);
    Q.block(0,0,dim,k)=tmp;

    start=k-1;
  }
  Eigen::VectorXcd test=Afunc(Q.col(0));
  std::cout<<"srt[0].first="<<srt[0].first<<" H(0,0)="<<H(0,0)<<" Q.col(0).dot(test)="<<Q.col(0).dot(test)<<std::endl;
  vec=Q.col(0);
  return srt[0].first;
}

std::complex<double> Largest_EV_Arnoldi_restart(Eigen::VectorXcd &vec, 
      int maxiter, double epsilon, 
      std::function<Eigen::VectorXcd(const Eigen::VectorXcd &v)> Afunc, 
      int verbosity, bool reortho, int m){

  int m_orig=m;
  if(m>maxiter){m=maxiter;}
  int restarts=maxiter/m;
  if(verbosity>1)std::cout<<"Arnoldi m="<<m<<" maxiter="<<maxiter<<" epsilon="<<epsilon<<std::endl;
  Eigen::MatrixXcd Q=Eigen::MatrixXcd::Zero(vec.size(),m);
  std::complex<double> max_eval=0;
  double diff=0;

  for(int run_=0; run_<restarts; run_++){
    if(maxiter<(run_+1)*m_orig){ m=maxiter-run_*m_orig; }
    Eigen::MatrixXcd H=Eigen::MatrixXcd::Zero(m, m);
    bool stopped=false;
    Q.col(0)=vec; Q.col(0).normalize();
    for(int k=1; k<m; k++){
      Q.col(k)=Afunc(Q.col(k-1));
      for(int j=0; j<k; j++){
        H(j,k-1)=Q.col(j).dot(Q.col(k));
        Q.col(k)-=H(j,k-1)*Q.col(j);
      }
      if(reortho){
        for(int j=0; j<k; j++){
          std::complex<double> overlap=Q.col(j).dot(Q.col(k));
          Q.col(k)-=overlap*Q.col(j);
          H(j,k-1)+=overlap;
        }
      }
      H(k,k-1)=Q.col(k).norm(); 
      Q.col(k)/=H(k,k-1);
 
      if(H(k,k-1).real()<epsilon){
        if(verbosity>0){
          std::cout<<"Largest_EV_Arnoldi converged to "<<H(k,k-1).real()<<"<"<<epsilon<<" at iteration "<<k+run_*m_orig<<"/"<<maxiter<<std::endl;
        }
        m=k;
        Eigen::MatrixXcd tmp = H.block(0,0,m,m); H.swap(tmp);
        stopped=true;
        break;
      }
      if(verbosity>1){
        std::cout<<"Arnoldi k="<<k<<" H(k,k-1)="<<H(k,k-1)<<std::endl;
      }
    }

    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solveH(H);
 //   int imax=0; for(int i=1; i<m; i++){if(std::abs(solveH.eigenvalues()(i))>std::abs(solveH.eigenvalues()(imax))){imax=i;}}
    int imax=0; for(int i=1; i<m; i++){if(std::real(solveH.eigenvalues()(i))>std::real(solveH.eigenvalues()(imax))){imax=i;}}

    Eigen::VectorXcd tmp=Eigen::VectorXcd::Zero(vec.size());
    for(int i=0; i<m; i++){
      tmp+=solveH.eigenvectors()(i, imax)*Q.col(i);
    }
    vec=Afunc(tmp);
    max_eval=tmp.dot(vec);
//    max_eval = solveH.eigenvalues()(imax);
    diff=(vec-max_eval*tmp).norm();
    vec.normalize();
    if(stopped){break;}

    if(diff<epsilon){
      if(verbosity>0){
          std::cout<<"Largest_EV_Arnoldi converged in run "<<run_<<"/"<<restarts<<std::endl;
      }
      break;
    }
  }
  if(verbosity>0){
    std::cout<<"max_eval="<<max_eval<<" diff="<<diff<<std::endl;
  }
  return max_eval;
}


}//namespace
