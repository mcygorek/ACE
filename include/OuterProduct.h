#ifndef OUTER_PRODUCT_DEFINED_H
#define OUTER_PRODUCT_DEFINED_H

#include <Eigen/Core>
#include <vector>

int loop_pow(int x, int n){
  int r=1;
  for(int i=0; i<n; i++)r*=x;
  return r;
}

Eigen::MatrixXcd OuterProduct(const Eigen::MatrixXcd &o1, 
                              const Eigen::MatrixXcd &o2){

  Eigen::MatrixXcd A(o1.rows()*o2.rows(), o1.cols()*o2.cols());
  for(int r1=0; r1<o1.rows(); r1++){
    for(int c1=0; c1<o1.cols(); c1++){
      for(int r2=0; r2<o2.rows(); r2++){
        for(int c2=0; c2<o2.cols(); c2++){
          A(r1*o2.rows()+r2, c1*o2.cols()+c2) = o1(r1, c1) * o2(r2, c2);
        }
      }
    }
  }
  return A;
}
Eigen::MatrixXcd OuterProduct(const Eigen::MatrixXcd &o1, 
                              const Eigen::MatrixXcd &o2,
                              const Eigen::MatrixXcd &o3){

  return OuterProduct(OuterProduct(o1,o2),o3);
}

Eigen::MatrixXcd ExpandSingleOp(int dim_before, const Eigen::MatrixXcd &m, int dim_after){
  if(dim_before<1 && dim_after<1)return m;
  if(dim_before<1)return OuterProduct(m, Eigen::MatrixXcd::Identity(dim_after,dim_after));
  if(dim_after<1)return OuterProduct(Eigen::MatrixXcd::Identity(dim_before,dim_before),m);
  return OuterProduct(Eigen::MatrixXcd::Identity(dim_before,dim_before), m,
                      Eigen::MatrixXcd::Identity(dim_after,dim_after));
}


Eigen::VectorXcd VectorOuterProduct(const Eigen::VectorXcd &o1, 
                              const Eigen::VectorXcd &o2){

  Eigen::VectorXcd A(o1.size()*o2.size());
  for(int r1=0; r1<o1.size(); r1++){
    for(int r2=0; r2<o2.size(); r2++){
      A(r1*o2.size()+r2) = o1(r1) * o2(r2);
    }
  }
  return A;
}

Eigen::VectorXcd H_Matrix_to_L_Vector(const Eigen::MatrixXcd &M){
  Eigen::VectorXcd v(M.rows()*M.cols());
  for(int i=0; i<M.rows(); i++){
    for(int j=0; j<M.cols(); j++){
      v(i*M.cols()+j)=M(i,j);
    } 
  }
  return v;
}

#endif
