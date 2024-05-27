#pragma once
#ifndef ACE_CHECK_MATRIX_DEFINED_H
#define ACE_CHECK_MATRIX_DEFINED_H

#include <Eigen/Dense>
	
namespace ACE{

extern void check_matrix_square(const Eigen::MatrixXcd & A, const std::string & name);

inline int get_dim_check_square(const Eigen::MatrixXcd & A, const std::string & name){
  check_matrix_square(A, name);
  return A.rows();
}

extern void check_matrix_rows_min(const Eigen::MatrixXcd & A, int dim, const std::string & name);

extern void check_matrix_cols_min(const Eigen::MatrixXcd & A, int dim, const std::string & name);

//both dimensions must be at least dim
inline void check_matrix_dims_min(const Eigen::MatrixXcd & A, int dim, const std::string & name){
  check_matrix_rows_min(A, dim, name);
  check_matrix_cols_min(A, dim, name);
}

//check if matrix is square and at least of dimension dim
inline void check_matrix_square_min(const Eigen::MatrixXcd & A, int dim, const std::string & name){
  check_matrix_square(A, name);
  check_matrix_rows_min(A, dim, name);
}

inline int get_dim_check_square_min(const Eigen::MatrixXcd & A, int dim, const std::string & name){
  check_matrix_square(A, name);
  check_matrix_rows_min(A, dim, name);
  return A.rows();
}

extern void check_matrix_rows_eq(const Eigen::MatrixXcd & A, int dim, const std::string & name);

extern void check_matrix_cols_eq(const Eigen::MatrixXcd & A, int dim, const std::string & name);

//check if matrix is square and dimension equals dim
inline void check_matrix_square_eq(const Eigen::MatrixXcd & A, int dim, const std::string & name){
  check_matrix_square(A, name);
  check_matrix_rows_eq(A, dim, name);
}


//vector tests
extern int get_vector_dim_min(const Eigen::VectorXcd & A, int dim, const std::string & name);

inline void check_vector_dim_min(const Eigen::VectorXcd & A, int dim, const std::string & name){
  get_vector_dim_min(A, dim, name);
}

extern void check_vector_dim_eq(const Eigen::VectorXcd & A, int dim, const std::string & name);


//check if dimension is a square of an integer number
int get_dim_sqrt(int i, const std::string & name);

inline void check_dim_is_square(int i, const std::string & name){
  get_dim_sqrt(i, name);
}


//generally: check bounds:
extern void check_bounds(int i, int dim, const std::string & name);

extern void check_at_least(int i, int cmp, const std::string & name);

}//namespace
#endif
