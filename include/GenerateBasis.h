#ifndef ACE_GENERATE_BASIS_DEFINED_H
#define ACE_GENERATE_BASIS_DEFINED_H

namespace ACE{

Eigen::MatrixXcd generate_basis(std::vector<Eigen::VectorXcd> v, double epsilon=1e-12){

  if(v.size()<1){
    std::cerr<<"gernerate_basis: need at least one vector!"<<std::endl;
    exit(1);
  }
  int N=v[0].rows();
  if(N<1){std::cerr<<"generate_basis: N<1!"<<std::endl; exit(1);}
  for(int i=1; i<(int)v.size(); i++){
    if(v[i].rows()!=N){std::cerr<<"generate_basis: v[i].rows()!=N!"<<std::endl; exit(1);}
  }
  
  Eigen::MatrixXcd U=Eigen::MatrixXcd::Identity(N,N);
  for(int i=0; i<(int)v.size() && i<N; i++){
    U.col(i)=v[i];
  }
  

 int counter=0;
 for(int loop=0; loop<3; loop++){

  for(int i=0; i<N; i++){
    for(int j=0; j<i; j++){
      std::complex<double> ov=U.col(j).dot(U.col(i));
      U.col(i)-=ov*U.col(j);
    } 

    double norm=U.col(i).norm();
//std::cout<<"v["<<i<<"].norm(): "<<norm<<std::endl;
    if(norm<epsilon){
      U.col(i)=Eigen::VectorXcd::Zero(N); 
      U.col(i)((counter++)%N)=1.;
      i--; continue;
    }else{
      U.col(i).normalize();
    }
  }
 }

  //check again:
  double d=max_diff_from_ortho(U);
  if(d>1e-15){ 
    std::cerr<<"generate_basis: max_diff_from_ortho(U)>1e-15  ("<<d<<")!"<<std::endl;
    print_diff_from_ortho(U, 1e-15);
    exit(1);
  }
  return U;
}
Eigen::MatrixXcd generate_basis(const Eigen::VectorXcd &v){
  return generate_basis(std::vector<Eigen::VectorXcd>(1,v));
}
Eigen::MatrixXcd generate_basis(const Eigen::VectorXcd &v, 
                                const Eigen::MatrixXcd &M){
  std::vector<Eigen::VectorXcd> list(1,v);
  for(int i=0; i<M.cols(); i++)list.push_back(M.col(i));
  return generate_basis(list);
}

}//namespace
#endif
