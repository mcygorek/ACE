#include <Eigen/Eigenvalues>
#include "LindbladMasterEquation.hpp"
#include "Constants.hpp"
#include "Reader.hpp"
#include "DummyException.hpp"
#include "HermitianLiouvilleBasis.hpp"
#include "LiouvilleTools.hpp"

namespace ACE{
int LindbladMasterEquation::get_dim()const{
  return H.rows();
}
Eigen::MatrixXcd LindbladMasterEquation::construct_Liouvillian()const{
    constexpr double hbar=hbar_in_meV_ps;
    int D=get_dim();
    Eigen::MatrixXcd Lrec=Eigen::MatrixXcd::Zero(D*D,D*D);
    for(int mu1=0; mu1<D; mu1++){
      for(int mu0=0; mu0<D; mu0++){
        for(int nu=0; nu<D; nu++){
          Lrec(mu1*D+nu, mu0*D+nu)+=std::complex<double>(0,-1./hbar)*H(mu1,mu0);
        }
      }
    }
    for(int nu1=0; nu1<D; nu1++){
      for(int nu0=0; nu0<D; nu0++){
        for(int mu=0; mu<D; mu++){
          Lrec(mu*D+nu1, mu*D+nu0)+=std::complex<double>(0,1./hbar)*std::conj(H(nu1,nu0));
        }
      }
    }
    for(size_t j=0; j<L.size(); j++){
      for(int mu1=0; mu1<D; mu1++){
        for(int nu1=0; nu1<D; nu1++){
          for(int mu0=0; mu0<D; mu0++){
            for(int nu0=0; nu0<D; nu0++){
              Lrec(mu1*D+nu1,mu0*D+nu0)+= 
L[j].first*L[j].second(mu1,mu0)*std::conj(L[j].second(nu1,nu0));
            }
          }
        }
      }
      Eigen::MatrixXcd LdL=L[j].first*L[j].second.adjoint()*L[j].second;
      
      for(int mu1=0; mu1<D; mu1++){
        for(int mu0=0; mu0<D; mu0++){
          for(int nu=0; nu<D; nu++){
            Lrec(mu1*D+nu,mu0*D+nu)+=-0.5*LdL(mu1*D+mu0);
            Lrec(nu*D+mu1,nu*D+mu0)+=-0.5*std::conj(LdL(mu1*D+mu0));
          }
        }
      }
    }
    return Lrec;
}

void LindbladMasterEquation::set_from_Liouvillian_Hall(const Eigen::MatrixXcd &Liou, double threshold, int verbosity){

  constexpr double hbar=hbar_in_meV_ps;
  int D=sqrt(Liou.rows()); int DL=D*D;
  if(Liou.rows()!=DL || Liou.rows()!=DL){
    std::cerr<<"Liouvillian.rows()!=DL || Liouvillian.rows()!=DL"<<std::endl;
    throw DummyException();
  }

  //Matrix containing basis of operators as columns:
  Eigen::MatrixXcd G = HermitianLiouvilleBasis(D).get_Matrix();
  //as list of matrices on Hilbert space:
  std::vector<Eigen::MatrixXcd> GS(DL);
  for(int m=0; m<DL; m++){
    GS[m]=L_Vector_to_H_Matrix( G.col(m) );
  }

  //Hall et al. PRA 89, 042120, Eq (A7):
  Eigen::MatrixXcd deph_mat = Eigen::MatrixXcd::Zero(DL-1, DL-1);
  for(int i=1; i<DL; i++){
    for(int j=1; j<DL; j++){
      for(int m=0; m<DL; m++){
        deph_mat(i-1,j-1) +=
          (GS[m]*GS[i]*L_Vector_to_H_Matrix(Liou*G.col(m))*GS[j]).trace();
      }
    }
  }

//  std::cout<<"(deph_mat-deph_mat.adjoint()).norm()="<<(deph_mat-deph_mat.adjoint()).norm()<<std::endl;
  Eigen::MatrixXcd B = 0.5*(deph_mat+deph_mat.adjoint());
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(B);
//std::cout<<"Eigenvalues: "<<solver.eigenvalues().transpose()<<std::endl;
  if(verbosity>0)std::cout<<"Eigenvalues: "<<solver.eigenvalues().transpose()<<std::endl;
  
  {std::vector<std::pair<double,Eigen::MatrixXcd> > tmp(DL-1, 
               std::pair<double,Eigen::MatrixXcd>(0,
               Eigen::MatrixXcd::Zero(D,D)));         L.swap(tmp);}
  for(int k=0; k<DL-1; k++){
    L[k].first=solver.eigenvalues()(k);
    for(int j=0; j<DL-1; j++){
      L[k].second+=solver.eigenvectors()(j,k)*GS[j+1];
    }
    //normalize to 2-norm (instead of Frobenius norm):
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd(L[k].second);
    double norm=svd.singularValues()(0);
    //But only we don't divide by zero
    if(norm>1e-9){
      L[k].second/=norm;
      L[k].first*=norm*norm;
    }
  }

  std::sort(L.begin(), L.end(), 
     [](const std::pair<double,Eigen::MatrixXcd>  &a, 
        const std::pair<double,Eigen::MatrixXcd>  &b){
        return a.first>b.first;});

  //get H from the remainder, assuming Tr(H)=0
  H = Eigen::MatrixXcd::Zero(D,D);
  Eigen::MatrixXcd LiouH = Liou - construct_Liouvillian();
  for(int mu1=0; mu1<D; mu1++){
    for(int mu0=0; mu0<D; mu0++){
      for(int nu=0; nu<D; nu++){
        H(mu1,mu0)+=LiouH(mu1*D+nu, mu0*D+nu);
      }
    }
  }
  H*=std::complex<double>(0,hbar)/((double)D); 
  
}

void LindbladMasterEquation::set_from_Liouvillian(const Eigen::MatrixXcd &Liou, double threshold, int verbosity){

  constexpr double hbar=hbar_in_meV_ps;
  int D=sqrt(Liou.rows()); int DL=D*D;
  if(Liou.rows()!=DL || Liou.rows()!=DL){
    std::cerr<<"Liouvillian.rows()!=DL || Liouvillian.rows()!=DL"<<std::endl;
    throw DummyException();
  }

  std::complex<double> totaltrace=0;
  for(int mu=0; mu<D; mu++){
    for(int nu=0; nu<D; nu++){
      totaltrace+=Liou(mu*D+nu, mu*D+nu);
    }
  }
  if(verbosity>0)std::cout<<"Total trace: "<<totaltrace<<std::endl;
  //extract Hamiltonian
  Eigen::MatrixXcd Hfw=Eigen::MatrixXcd::Zero(D,D);
  for(int mu1=0; mu1<D; mu1++){
    for(int mu0=0; mu0<D; mu0++){
      for(int nu=0; nu<D; nu++){
        Hfw(mu1,mu0)+=Liou(mu1*D+nu, mu0*D+nu)/((double)D);
      }
    }
    Hfw(mu1,mu1)-=0.5*totaltrace/((double)D*D);
  }
  Hfw*=std::complex<double>(0,hbar);

  Eigen::MatrixXcd Hbw=Eigen::MatrixXcd::Zero(D,D);
  for(int nu1=0; nu1<D; nu1++){
    for(int nu0=0; nu0<D; nu0++){
      for(int mu=0; mu<D; mu++){
        Hbw(nu1,nu0)+=std::conj(Liou(mu*D+nu1, mu*D+nu0))/((double)D);
      }
    }
    Hbw(nu1,nu1)-=0.5*std::conj(totaltrace)/((double)D*D);
  }
  Hbw*=std::complex<double>(0,hbar);

//    if(verbosity>0)std::cout<<"Hfw:"<<std::endl<<Hfw<<std::endl<<std::endl;
//    if(verbosity>0)std::cout<<"Hbw:"<<std::endl<<Hbw<<std::endl<<std::endl;

  Eigen::MatrixXcd Htilde=(Hfw+Hbw)/2.;
  double norm_fwbw=(Hfw-Hbw).norm();
  if(verbosity>0)std::cout<<"|Hfw-Hbw|="<<norm_fwbw<<std::endl;
  H=(Htilde+Htilde.adjoint())/2.;
  
  //extract collapse operators
  Eigen::MatrixXcd Ltilde=Liou;
  for(int mu1=0; mu1<D; mu1++){
    for(int mu0=0; mu0<D; mu0++){
      for(int nu=0; nu<D; nu++){
        Ltilde(mu1*D+nu, mu0*D+nu)-=std::complex<double>(0,-1./hbar)*Hfw(mu1,mu0);
      }
    }
  }
  for(int nu1=0; nu1<D; nu1++){
    for(int nu0=0; nu0<D; nu0++){
      for(int mu=0; mu<D; mu++){
        Ltilde(mu*D+nu1, mu*D+nu0)-=std::complex<double>(0,1./hbar)*std::conj(Hbw(nu1,nu0));
      }
    }
  }

  Eigen::MatrixXcd A(D*D,D*D);
  for(int mu1=0; mu1<D; mu1++){
    for(int nu1=0; nu1<D; nu1++){
      for(int mu0=0; mu0<D; mu0++){
        for(int nu0=0; nu0<D; nu0++){
          A(mu1*D+mu0,nu1*D+nu0)=Ltilde(mu1*D+nu1,mu0*D+nu0);
        }
      }
    }
  }
   
  double norm_A_antiH=(A-A.adjoint()).norm();
    
  if(verbosity>0)std::cout<<"|A-A^+|="<<norm_A_antiH<<std::endl;
  Eigen::MatrixXcd B=((A+A.adjoint())/2.).eval();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(B);
  if(verbosity>0)std::cout<<"Eigenvalues: "<<solver.eigenvalues().transpose()<<std::endl;
 

  {std::vector<std::pair<double,Eigen::MatrixXcd> > tmp(DL, 
               std::pair<double,Eigen::MatrixXcd>(0,
               Eigen::MatrixXcd::Zero(D,D)));         L.swap(tmp);}
  for(int k=0; k<DL; k++){
    L[k].first=solver.eigenvalues()(k);
    for(int mu1=0; mu1<D; mu1++){
      for(int mu0=0; mu0<D; mu0++){
        L[k].second(mu1,mu0)+=solver.eigenvectors()(mu1*D+mu0,k);
      }
    }
    //normalize to 2-norm (instead of Frobenius norm):
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd(L[k].second);
    double norm=svd.singularValues()(0);
    //But only we don't divide by zero
    if(norm>1e-9){
      L[k].second/=norm;
      L[k].first*=norm*norm;
    }
  }

  std::sort(L.begin(), L.end(), 
     [](const std::pair<double,Eigen::MatrixXcd>  &a, 
        const std::pair<double,Eigen::MatrixXcd>  &b){
        return a.first>b.first;});

/*
  L.clear();
  for(int j=0; j<D*D; j++){
    if(threshold >0. && std::abs(solver.eigenvalues()(j))<threshold){
      continue;
    }
    Eigen::MatrixXcd thisL(D,D);
    int mu1max=0; int mu0max=0;
    for(int mu1=0; mu1<D; mu1++){
      for(int mu0=0; mu0<D; mu0++){
        thisL(mu1,mu0)=solver.eigenvectors()(mu1*D+mu0,j);
        if(std::abs(thisL(mu1,mu0))>std::abs(thisL(mu1max,mu0max))){
          mu1max=mu1; mu0max=mu0;
        }
      }
    }
    //phase of largest element to 1:
    {std::complex<double> phase=std::abs(thisL(mu1max,mu0max))/thisL(mu1max,mu0max);
    thisL*=phase;}
     
    L.push_back(std::make_pair(solver.eigenvalues()(j), thisL));
  }
*/
}

void LindbladMasterEquation::print(std::ostream &os, double epsilon)const{
    os<<"H:"<<std::endl<<round_Matrix(H, epsilon)<<std::endl<<std::endl;
    for(size_t j=0; j<L.size(); j++){
      os<<"gamma["<<j<<"]="<<L[j].first<<std::endl;
      os<<" L["<<j<<"]:"<<std::endl<<round_Matrix(L[j].second,epsilon)<<std::endl<<std::endl;
    }
}

void LindbladMasterEquation::print_param(const std::string &fname, const TimeGrid &tgrid, bool print_Markov, double epsilon)const{
  std::ofstream ofs(fname.c_str());
  ofs<<"ta "<<tgrid.ta<<std::endl;
  ofs<<"dt "<<tgrid.dt<<std::endl;
  ofs<<"te "<<tgrid.n_tot*tgrid.dt<<std::endl;
  ofs<<"add_Hamiltonian "<<Matrix_as_parameter(H ,epsilon)<<std::endl;
  for(size_t j=0; j<L.size(); j++){
    if(!print_Markov || L[j].first>0.){
      ofs<<"add_Lindblad "<<L[j].first<<" "<<Matrix_as_parameter(L[j].second, epsilon)<<std::endl;
    }
  }
}

}//namespace


