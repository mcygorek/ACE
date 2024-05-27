#include "Trafo_Chain.hpp"
#include "MPS.hpp"
#include "ModePropagatorGenerator.hpp"
#include "Tensor.hpp"
#include "BinaryReader.hpp"
#include "otimes.hpp"

namespace ACE{

Eigen::MatrixXd pseudoinverse(const Eigen::MatrixXd &M){
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(M, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::VectorXd sv=svd.singularValues();
  int maxk=sv.size();
  Eigen::VectorXd sv_inv(maxk);
  for(int i=0; i<maxk; i++){
    sv_inv(i)=(fabs(sv(i))>1e-30) ? (1./sv(i)) : 0.;
  }
  return svd.matrixV().block(0,0,M.cols(), maxk) * sv_inv.head(maxk).asDiagonal() * svd.matrixU().adjoint().block(0,0,maxk, M.rows());
}

  void Trafo_Chain::print_info(std::ostream &ofs){
    ofs<<"Trafo_Chain: size="<<T.size()<<":";
    for(size_t i=0; i<T.size(); i++){
      ofs<<" ("<<T[i].rows()<<","<<T[i].cols()<<")";
    }
    ofs<<std::endl;
    ofs<<"Trafo_Chain (inverse): size="<<Tinv.size()<<":";
    for(size_t i=0; i<Tinv.size(); i++){
      ofs<<" ("<<Tinv[i].rows()<<","<<Tinv[i].cols()<<")";
    }
    ofs<<std::endl;
  }

  void Trafo_Chain::add_low_to_high(const Eigen::MatrixXd &R){
    T.push_back(R);
    
    Eigen::MatrixXd Rinv=R.adjoint();
    for(int i=0; i<Rinv.cols(); i++){
      double norm2=Rinv.col(i).dot(Rinv.col(i));
      Rinv.col(i)/=norm2;
    }
    Tinv.push_back(Rinv);
  }


  void Trafo_Chain::add_high_to_low(const Eigen::MatrixXd &L){
    if(L.rows()!=T.back().rows()){
      std::cerr<<"Trafo_Chain: L.rows()!=T.back().rows() ("<<L.rows()<<" vs. "<<T.back().rows()<<")!"<<std::endl;
      exit(1);
    }
  
    Eigen::MatrixXd Linv=L.adjoint();
    for(int i=0; i<Linv.rows(); i++){
      double norm2=Linv.row(i).dot(Linv.row(i));
      Linv.row(i)/=norm2;
    }

    Eigen::MatrixXd tmp=Linv*T.back(); 
    T.back()=tmp;
    
    tmp=Tinv.back()*L;
    Tinv.back()=tmp;   
  }


  std::vector<int> Trafo_Chain::get_dims()const{
    std::vector<int> dims(T.size());
    int lastdim=1;
    for(size_t i=0; i<T.size(); i++){
      if(i>0)lastdim=T[i-1].rows();
      dims[i]=T[i].cols()/lastdim;
    }
    return dims;
  }

  void Trafo_Chain::check_compatible(const Trafo_Chain &other)const{
    if(other.size()!=size()){
      std::cerr<<"Trafo_Chain::check_compatible: other.size()!=size()!"<<std::endl;
      exit(1);
    }
    std::vector<int> dims=get_dims();
    std::vector<int> dims2=other.get_dims();
    for(size_t i=0; i<dims.size(); i++){
      if(dims[i]!=dims2[i]){
        std::cerr<<"Trafo_Chain::check_compatible: dims["<<i<<"]!=dims2["<<i<<"] ("<<dims[i]<<" vs. "<<dims2[i]<<")!"<<std::endl;
        exit(1);
      }
    }
  }
  
  Eigen::MatrixXd Trafo_Chain::overlap_matrix(const Trafo_Chain &other, bool skiplastmultiply)const{
    Eigen::MatrixXd v(1,1); v(0,0)=1;
//std::cout<<"START:: overlap_matrix"<<std::endl;
    std::vector<int> dims=get_dims(); 
    if(T.size()!=dims.size() || Tinv.size() !=dims.size()){
      std::cerr<<"Trafo_Chain::overlap_matrix: T.size()!=dims.size() || Tinv.size() !=dims.size()"<<std::endl;
      exit(1);
    }
    if(other.T.size()!=dims.size() || other.Tinv.size() !=dims.size()){
      std::cerr<<"Trafo_Chain::overlap_matrix: other.T.size()!=dims.size() || other.Tinv.size() !=dims.size()"<<std::endl;
      exit(1);
    }


    Eigen::MatrixXd v2;
    for(int o=0; o<dims.size(); o++){
      if(v.rows()*dims[o]!=T[o].cols()){
        std::cerr<<"Trafo_Chain::overlap_matrix: v.rows()*dims["<<o<<"]!=T["<<o<<"].cols()"<<std::endl;
        exit(1);
      }
      if(v.cols()*dims[o]!=other.Tinv[o].rows()){
        std::cerr<<"Trafo_Chain::overlap_matrix: v.cols()*dims["<<o<<"]!=other.Tinv["<<o<<"].cols()"<<std::endl;
        exit(1);
      }
      v2=Eigen::MatrixXd::Zero(T[o].cols(), other.Tinv[o].rows());
      for(size_t i=0; i<dims[o]; i++){
        for(int d1=0; d1<v.rows(); d1++){
          for(int d2=0; d2<v.cols(); d2++){
            v2(d1*dims[o]+i, d2*dims[o]+i) += v(d1, d2);
          }
        }
      }
      v=T[o]*v2*other.Tinv[o];
    }
//std::cout<<"END:: overlap_matrix"<<std::endl;
    if(skiplastmultiply)return v2;
    return v;
  }

  double Trafo_Chain::overlap(const Trafo_Chain &other)const{
    Eigen::MatrixXd v=overlap_matrix(other);
//std::cout<<"TEST: got overlap matrix!"<<std::endl;
    double res=0;
    for(int k1=0; k1<lastdim(); k1++){
      for(int k2=0; k2<other.lastdim(); k2++){
        res+=v(k1,k2)*v(k1,k2);
      }
    }
    return res;
  }

  void Trafo_Chain::combine(const Trafo_Chain &other, int max_add){
    check_compatible(other);
    if(T.size()<1)return;
    if(max_add==0)return;
    std::vector<int> dims=get_dims(); 

    std::vector<Eigen::MatrixXd> R(T.size());
    std::vector<Eigen::MatrixXd> Rinv(Tinv.size());
    for(size_t o=0; o<T.size(); o++){
      if(o==0){
        Rinv[o].resize(Tinv[o].rows(), Tinv[o].cols()+other.Tinv[o].cols());
        Rinv[o].block(0,0,Tinv[o].rows(),Tinv[o].cols())=Tinv[o];
        Rinv[o].block(0,Tinv[o].cols(),other.Tinv[o].rows(),other.Tinv[o].cols())=other.Tinv[o];
        R[o].resize(T[o].rows()+other.T[o].rows(), T[o].cols());
        R[o].block(0,0,T[o].rows(),T[o].cols())=T[o];
        R[o].block(T[o].rows(),0,other.T[o].rows(),other.T[o].cols())=other.T[o];
      }else{
        Rinv[o]=Eigen::MatrixXd::Zero(Tinv[o].rows()+other.Tinv[o].rows(),
                                   Tinv[o].cols()+other.Tinv[o].cols());
        Rinv[o].block(0,0,Tinv[o].rows(),Tinv[o].cols())=Tinv[o];
        Rinv[o].block(Tinv[o].rows(),Tinv[o].cols(),other.Tinv[o].rows(),other.Tinv[o].cols())=other.Tinv[o];
        R[o]=Eigen::MatrixXd::Zero(T[o].rows()+other.T[o].rows(), 
                                      T[o].cols()+other.T[o].cols());
        R[o].block(0,0,T[o].rows(),T[o].cols())=T[o];
        R[o].block(T[o].rows(),T[o].cols(),other.T[o].rows(),other.T[o].cols())=other.T[o];
      }
    }
    if(max_add>0 && max_add<other.T.back().rows()){
      Eigen::MatrixXd tmp;
      tmp=R.back().block(0,0,T.back().rows()+max_add,T.back().cols()+other.T.back().cols());
      R.back()=tmp;
      tmp=Rinv.back().block(0,0,Tinv.back().rows()+other.Tinv.back().rows(),Tinv.back().cols()+max_add);
      Rinv.back()=tmp;
    }
    T=R;
    Tinv=Rinv;
  }

  void Trafo_Chain::orthogonalize_single(int i, int j){
    Eigen::MatrixXd OM=overlap_matrix(*this); // v=T[o]*v2*other.Tinv[o];

//std::cout<<"ORTHO_SINGLE: "<<i<<", "<<j<<": "<<OM(i,j)<<std::endl;
    if(i==j){
//      if(fabs(OM(j,j))<swap_thr){
//std::cout<<"ORTHO_SINGLE: "<<i<<", "<<j<<": "<<OM(i,j)<<std::endl;
//        Tinv.back().col(j)*=0.;
//      }else{
        Tinv.back().col(j)/=OM(j,j);
//      }
    }else{
//      T.back().row(j)-=OM(i,j)*T.back().row(i);
//      Tinv.back().col(j)-=OM(j,i)*Tinv.back().col(i);
      T.back().row(j)-=T.back().row(i)*OM(j,i);
      Tinv.back().col(j)-=OM(i,j)*Tinv.back().col(i);
    }
  }

  void Trafo_Chain::orthogonalize(){
    if(T.size()<1)return;
    for(int j=0; j<T.back().rows(); j++){
      for(int i=0; i<=j; i++){
        orthogonalize_single(i,j);
      }
    }
  }


  //Perform SVD on all but the last T:
  void Trafo_Chain::compress_all_but_last_T(double epsilon){ 
    std::vector<int> dims=get_dims();
    for(size_t o=0; o<T.size()-1; o++){
      Eigen::JacobiSVD<Eigen::MatrixXd> svd( T[o] , Eigen::ComputeThinU | Eigen::ComputeThinV);
      Eigen::VectorXd sval=svd.singularValues();
      int newdim=1;
      for(int i=1; i<sval.rows(); i++){
        if(sval(i)>epsilon*sval(0))newdim++;
        else break;
      }
      T[o]=svd.matrixV().adjoint().block(0,0,newdim, T[o].cols());
      Eigen::MatrixXd pass_on=svd.matrixU().block(0,0,svd.matrixU().rows(),newdim)*sval.head(newdim).asDiagonal();
      Eigen::MatrixXd Id=Eigen::MatrixXd::Identity(dims[o+1],dims[o+1]);
      T[o+1]=T[o+1]*otimes_real(pass_on, Id);
    }
  }

  //Perform SVD on all but the last Tinv:
  void Trafo_Chain::compress_all_but_last_Tinv(double epsilon){ 
    std::vector<int> dims=get_dims();
    for(size_t o=0; o<Tinv.size()-1; o++){
      Eigen::JacobiSVD<Eigen::MatrixXd> svd( Tinv[o] , Eigen::ComputeThinU | Eigen::ComputeThinV);
      Eigen::VectorXd sval=svd.singularValues();
      int newdim=1;
      for(int i=1; i<sval.rows(); i++){
        if(sval(i)>epsilon*sval(0))newdim++;
        else break;
      }
      Tinv[o]=svd.matrixU().block(0,0,Tinv[o].rows(),newdim);
      Eigen::MatrixXd pass_on=sval.head(newdim).asDiagonal()*svd.matrixV().adjoint().block(0,0,newdim,svd.matrixV().rows());
      Eigen::MatrixXd Id=Eigen::MatrixXd::Identity(dims[o+1],dims[o+1]);
      Tinv[o+1]=otimes_real(pass_on, Id)*Tinv[o+1];
    }
  }
  
  void Trafo_Chain::compress_weight_to_Tinv(double epsilon){  //make orthogonal. Residual SVD to last Tinv.
std::cout<<"COMPRESS:: START (compress_weight_to_Tinv)"<<std::endl;
    std::vector<int> dims=get_dims(); 
    if(T.size()<1)return;

    compress_all_but_last_T(epsilon);
    compress_all_but_last_Tinv(epsilon);
    
    {
      Eigen::JacobiSVD<Eigen::MatrixXd> svd( T.back() , Eigen::ComputeThinU | Eigen::ComputeThinV);
      Eigen::VectorXd sval=svd.singularValues();
      int newdim=1;
      for(int i=1; i<sval.rows(); i++){
        if(sval(i)>epsilon*sval(0))newdim++;
        else break;
      }
      T.back()=svd.matrixV().adjoint().block(0,0,newdim, T.back().cols());
      Tinv.back()=Tinv.back()*svd.matrixU().block(0,0,svd.matrixU().rows(),newdim)*sval.head(newdim).asDiagonal();
    }
    {
      Eigen::JacobiSVD<Eigen::MatrixXd> svd( Tinv.back() , Eigen::ComputeThinU | Eigen::ComputeThinV);
      Eigen::VectorXd sval=svd.singularValues();
      std::cout<<"COMPRESS:: SVDs: "<<sval.transpose()<<std::endl;
      double sum=0; for(int i=0; i<sval.rows(); i++)sum+=sval(i);
      std::cout<<"sum of SVDs: "<<sum<<std::endl;
    }

std::cout<<"COMPRESS:: END"<<std::endl;
  }

  void Trafo_Chain::compress_weight_sym(double epsilon){  
std::cout<<"COMPRESS:: START (compress_weight_sym)"<<std::endl;
    std::vector<int> dims=get_dims(); 
    if(T.size()<1)return;

    compress_all_but_last_T(epsilon);
    compress_all_but_last_Tinv(epsilon);

    Eigen::MatrixXd A=Tinv.back()*T.back();
    Eigen::JacobiSVD<Eigen::MatrixXd> svd( A , Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd sval=svd.singularValues();
    int newdim=1;
    for(int i=1; i<sval.rows(); i++){
      if(sval(i)>epsilon*sval(0))newdim++;
      else break;
    }
    Tinv.back()=svd.matrixU().block(0,0,svd.matrixU().rows(),newdim);
    T.back()=svd.matrixV().adjoint().block(0,0,newdim, svd.matrixV().adjoint().cols());

    Eigen::VectorXd sqrtsval(newdim);
    for(int i=0; i<newdim; i++){
      sqrtsval(i)=sqrt(sval(i));
    }
    T.back()=sqrtsval.asDiagonal()*T.back();
    Tinv.back()=Tinv.back()*sqrtsval.asDiagonal();
    
    std::cout<<"COMPRESS:: SVDs: "<<sval.transpose()<<std::endl;
    double sum=0; for(int i=0; i<sval.rows(); i++)sum+=sval(i);
    std::cout<<"sum of SVDs: "<<sum<<std::endl;

//    for(int o=0; o<T.size(); o++)Tinv[o]=T[o].adjoint();
std::cout<<"COMPRESS:: END"<<std::endl;
  }

  void Trafo_Chain::SVD_sweep_T_to_Tinv(double epsilon){ 
    std::vector<int> dims=get_dims();
    if(T.size()<1)return;
    for(size_t o=0; o<T.size(); o++){
      Eigen::JacobiSVD<Eigen::MatrixXd> svd( T[o] , Eigen::ComputeThinU | Eigen::ComputeThinV);
      Eigen::VectorXd sval=svd.singularValues();
      int newdim=1;
      for(int i=1; i<sval.rows(); i++){
        if(sval(i)>epsilon*sval(0))newdim++;
        else break;
      }
      if(o<T.size()-1){
        T[o]=svd.matrixV().adjoint().block(0,0,newdim, T[o].cols());
        Eigen::MatrixXd pass_on=svd.matrixU().block(0,0,svd.matrixU().rows(),newdim)*sval.head(newdim).asDiagonal();
        Eigen::MatrixXd Id=Eigen::MatrixXd::Identity(dims[o+1],dims[o+1]);
        T[o+1]=T[o+1]*otimes_real(pass_on, Id);
      }else{
        T.back()=svd.matrixV().adjoint().block(0,0,newdim, T.back().cols());
        Tinv.back()=Tinv.back()*svd.matrixU().block(0,0, T.back().rows(), newdim)*sval.head(newdim).asDiagonal();
      }
    }
    for(int o=Tinv.size(); o>=0; o--){
      Eigen::JacobiSVD<Eigen::MatrixXd> svd( Tinv[o] , Eigen::ComputeFullU | Eigen::ComputeFullV);
      Eigen::VectorXd sval=svd.singularValues();
      int newdim=1;
      for(int i=1; i<sval.rows(); i++){
        if(sval(i)>epsilon*sval(0))newdim++;
        else break;
      }
      if(o>0){
        int facdim=Tinv[o].rows()/dims[o];
        Tinv[o]=svd.matrixV().adjoint().block(0,0,newdim, Tinv[o].cols());

std::cerr<<"PROBLEM: How to un-expand row dimension of U?"<<std::endl; exit(1);

        Eigen::MatrixXd pass_on=svd.matrixU().block(0,0,svd.matrixU().rows(),newdim)*sval.head(newdim).asDiagonal();
        Eigen::MatrixXd rearrange(dims[o], dims[o]);
        Eigen::JacobiSVD<Eigen::MatrixXd> svd( rearrange , Eigen::ComputeFullU | Eigen::ComputeFullV);

        Eigen::MatrixXd Id=Eigen::MatrixXd::Identity(facdim,facdim);
        Tinv[o-1]=Tinv[o]*otimes_real(pass_on, Id);
      }else{
        Tinv[0]=svd.matrixV().adjoint().block(0,0, Tinv[0].rows(), Tinv[0].cols());
std::cout<<"SVD_sweep_T_to_Tinv("<<epsilon<<"): singular values: "<<sval.transpose()<<std::endl;
//        Tinv.back()=Tinv.back()*svd.matrixU().block(0,0, T.back().rows(), newdim)*sval.head(newdim).asDiagonal();
      }
    }
  }

  //Perform SVD on all but the last Tinv:
  void Trafo_Chain::compress(double epsilon){
    compress_weight_to_Tinv(epsilon);
//    compress_weight_sym(epsilon); 
    
//    Eigen::MatrixXd OM=overlap_matrix(*this);
//std::cout<<"Overlaps: max_diff_from_ortho="<<max_diff_from_ortho(OM)<<std::endl;
//    print_diff_from_ortho(OM,1e-8);
  }

  void Trafo_Chain::SVD_orthogonalize(double epsilon){
    Eigen::MatrixXd OM=overlap_matrix(*this);
    Eigen::JacobiSVD<Eigen::MatrixXd> svd( OM , Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd sval=svd.singularValues();
    int newdim=1;
    for(int i=1; i<sval.rows(); i++){
      if(sval(i)>epsilon*sval(0))newdim++;
      else break;
    }
std::cout<<"SVDs of overlap: "<<sval.transpose()<<std::endl;
    Eigen::VectorXd svalinv(newdim);
    for(int i=0; i<newdim; i++){
      svalinv(i)=1./sval(i);
    }
 
    T.back()=svd.matrixU().adjoint().block(0,0,newdim, T.back().rows())*T.back();
    Tinv.back()=Tinv.back()*svd.matrixV().block(0,0,Tinv.back().cols(),newdim)
                           *svalinv.asDiagonal();
}
 
  void Trafo_Chain::Eigen_orthogonalize(double epsilon){
    Eigen::MatrixXd OM=overlap_matrix(*this);
    Eigen::EigenSolver<Eigen::MatrixXd> solver( OM );
    Eigen::VectorXcd sval=solver.eigenvalues();
    int newdim=1;
    for(int i=1; i<sval.rows(); i++){
      if(abs(sval(i))>epsilon*abs(sval(0)))newdim++;
      else break;
    }
std::cout<<"eigenvalues of overlap: "<<sval.transpose()<<std::endl;
    Eigen::VectorXd svalinv(newdim);
    for(int i=0; i<newdim; i++){
      svalinv(i)=1./sqrt(abs(sval(i)));
    }
    Eigen::MatrixXd V=solver.pseudoEigenvectors();
    T.back()=svalinv.asDiagonal()
            *V.inverse().block(0,0,newdim, T.back().rows())*T.back();
    Tinv.back()=Tinv.back()*V.block(0,0,Tinv.back().cols(),newdim)
                           *svalinv.asDiagonal();
}

  void Trafo_Chain::add_ortho(Trafo_Chain other, double epsilon, double epsilon2, int max_add){
    check_compatible(other);
std::cout<<"ADD_ORTHO:: START"<<std::endl;

    combine(other, max_add);

    Eigen_orthogonalize(epsilon);
    if(epsilon2>0){
      compress(epsilon2);
      Eigen_orthogonalize(epsilon);
    }
//    orthogonalize();

std::cout<<"ADD_ORTHO:: END"<<std::endl;
  } 

  void Trafo_Chain::add_ortho2(Trafo_Chain other, double epsilon, int max_add){
    check_compatible(other);
std::cout<<"ADD_ORTHO2:: START"<<std::endl;

    //Rearrange "other" so that components with least overlap are first
    {
      Eigen::MatrixXd OM=overlap_matrix(other); //v=T[o]*v2*other.Tinv[o];
      Eigen::JacobiSVD<Eigen::MatrixXd> svd( OM , Eigen::ComputeFullU | Eigen::ComputeFullV);

      Eigen::MatrixXd tmp(other.Tinv.back().cols(), other.Tinv.back().cols());
      for(int j=0; j<other.Tinv.back().cols(); j++){
        for(int k=0; k<other.Tinv.back().cols(); k++){
          tmp(j,k)=svd.matrixV()(j,other.Tinv.back().cols()-1-k);
        }
      }
      other.Tinv.back()=other.Tinv.back()*tmp;
      other.T.back()=tmp.adjoint()*other.T.back();
    }
    //Orthogonalize "other" wrt. *this; redefine "max_add" if residual norm smaller than epsilon

    //combine 
    combine(other, max_add);

std::cout<<"ADD_ORTHO2:: END"<<std::endl;
  } 


  void Trafo_Chain::combine_Tinv_T(double epsilon){
    if(T.size()<1)return;
    std::vector<int> dims=get_dims(); 

    std::vector<Eigen::MatrixXd> R(T.size());
    for(size_t o=0; o<T.size(); o++){
      if(o==0){
        R[o].resize(T[o].rows()+Tinv[o].cols(), T[o].cols());
        R[o].block(0,0,T[o].rows(),T[o].cols())=T[o];
        R[o].block(T[o].rows(),0,Tinv[o].cols(),Tinv[o].rows())=Tinv[o].adjoint();
      }else{
        R[o]=Eigen::MatrixXd::Zero(T[o].rows()+Tinv[o].cols(), 
                                   T[o].cols()+Tinv[o].rows());
        R[o].block(0,0,T[o].rows(),T[o].cols())=T[o];
        R[o].block(T[o].rows(),T[o].cols(),Tinv[o].cols(),Tinv[o].rows())=Tinv[o].adjoint();
      }
    }

    T=R;

    for(size_t o=0; o<T.size(); o++){
      Eigen::JacobiSVD<Eigen::MatrixXd> svd( T[o] , Eigen::ComputeThinU | Eigen::ComputeThinV);
      Eigen::VectorXd sval=svd.singularValues();
      int newdim=1;
      for(int i=1; i<sval.rows(); i++){
        if(sval(i)>epsilon*sval(0))newdim++;
        else break;
      }
std::cout<<"o="<<o<<" svals: "<<sval.transpose()<<std::endl;
      T[o]=svd.matrixV().adjoint().block(0,0,newdim, T[o].cols());
      if(o<T.size()-1){
        Eigen::MatrixXd pass_on=svd.matrixU().block(0,0,svd.matrixU().rows(),newdim)*sval.head(newdim).asDiagonal();
        Eigen::MatrixXd Id=Eigen::MatrixXd::Identity(dims[o+1],dims[o+1]);
        T[o+1]=T[o+1]*otimes_real(pass_on, Id);
      }
    
      Tinv[o]=T[o].adjoint();
    }
  }

  void Trafo_Chain::read(const std::string &fname){
    T.clear();
    std::ifstream ifs(fname.c_str());
    if(!ifs.good()){
      std::cerr<<"Cannot read Trafo_Chain file '"<<fname<<"'!"<<std::endl;
      exit(1);
    }

    if(binary_read_int(ifs, "Trafo_Chain")!=9090){
      std::cerr<<"Cannot interpret '"<<fname<<"' as Trafo_Chain file!"<<std::endl;
      exit(1);
    }

    T.resize(binary_read_int(ifs, "Trafo_Chain"));
    for(size_t i=0; i<T.size(); i++){
      T[i]=binary_read_EigenMatrixXd(ifs, "Trafo_Chain");
    }

    Tinv.resize(binary_read_int(ifs, "Trafo_Chain"));
    for(size_t i=0; i<Tinv.size(); i++){
      Tinv[i]=binary_read_EigenMatrixXd(ifs, "Trafo_Chain");
    }
  }

  void Trafo_Chain::write(const std::string &fname){
    std::ofstream ofs(fname.c_str());

    binary_write_int(ofs, 9090); //let this be our magic number

    binary_write_int(ofs, T.size());
    for(size_t i=0; i<T.size(); i++){
      binary_write_EigenMatrixXd(ofs, T[i]);
    }
    binary_write_int(ofs, Tinv.size());
    for(size_t i=0; i<Tinv.size(); i++){
      binary_write_EigenMatrixXd(ofs, Tinv[i]);
    }
  }

}//namespace
