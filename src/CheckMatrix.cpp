#include "PCH.hpp"
#include "CheckMatrix.hpp"
#include <stdexcept>
#include <sstream>
#include <Eigen/Dense>
	
namespace ACE{

void check_matrix_square(const Eigen::MatrixXcd & A, const std::string & name){
  if(A.rows()!=A.cols()){
    std::stringstream ss;
    ss<<"Matrix "<<name<<" is not square (Dimensions: "<<A.rows()<<", "<<A.cols()<<")!"<<std::endl;

    throw std::runtime_error(ss.str()); 
  }
}

void check_matrix_rows_min(const Eigen::MatrixXcd & A, int dim, const std::string & name){
  if(A.rows()<dim){
    std::stringstream ss;
    ss<<"Matrix "<<name<<" has not enough rows ("<<A.rows()<<" while >= "<<dim<<" are required)!"<<std::endl;

    throw std::runtime_error(ss.str()); 
  }
}

void check_matrix_cols_min(const Eigen::MatrixXcd & A, int dim, const std::string & name){
  if(A.cols()<dim){
    std::stringstream ss;
    ss<<"Matrix "<<name<<" has not enough columns ("<<A.cols()<<" while >= "<<dim<<" are required)!"<<std::endl;

    throw std::runtime_error(ss.str()); 
  }
}

void check_matrix_rows_eq(const Eigen::MatrixXcd & A, int dim, const std::string & name){
  if(A.rows()!=dim){
    std::stringstream ss;
    ss<<"Matrix "<<name<<" has not the right amount of rows ("<<A.rows()<<" while = "<<dim<<" are required)!"<<std::endl;

    throw std::runtime_error(ss.str()); 
  }
}

void check_matrix_cols_eq(const Eigen::MatrixXcd & A, int dim, const std::string & name){
  if(A.cols()!=dim){
    std::stringstream ss;
    ss<<"Matrix "<<name<<" has not the right amount of columns ("<<A.cols()<<" while = "<<dim<<" are required)!"<<std::endl;

    throw std::runtime_error(ss.str()); 
  }
}


//vector tests
int get_vector_dim_min(const Eigen::VectorXcd & A, int dim, const std::string & name){
  if(A.size()<dim){
    std::stringstream ss;
    ss<<"Vector "<<name<<" is not long enough ("<<A.size()<<" while >= "<<dim<<" are required)!"<<std::endl;

    throw std::runtime_error(ss.str()); 
  }
  return A.size();
}

void check_vector_dim_eq(const Eigen::VectorXcd & A, int dim, const std::string & name){
  if(A.size()!=dim){
    std::stringstream ss;
    ss<<"Vector "<<name<<" has not the right length ("<<A.size()<<" while = "<<dim<<" are required)!"<<std::endl;

    throw std::runtime_error(ss.str()); 
  }
}


//check if dimension is a square of an integer number
int get_dim_sqrt(int i, const std::string & name){
  int j=sqrt(i);
  if( j*j !=i ){
    std::stringstream ss;
    ss<<"Dimension "<<name<<"="<<i<<" is not a square of on integer!"<<std::endl;

    throw std::runtime_error(ss.str());
  }
  return j;
}

//generally: check bounds:
void check_bounds(int i, int dim, const std::string & name){
  if( i<0 ){
    std::stringstream ss;
    ss<<"Out of bounds: Index "<<name<<"="<<i<<" must not be smaller than 0!"<<std::endl;

    throw std::runtime_error(ss.str());
  }else if( i>=dim ){
    std::stringstream ss;
    ss<<"Out of bounds: Index "<<name<<"="<<i<<" must be smaller than "<<dim<<"!"<<std::endl;

    throw std::runtime_error(ss.str());
  }
}

void check_at_least(int i, int cmp, const std::string & name){
  if(i<cmp){
    std::stringstream ss;
    ss<<"Variable: '"<<name<<"'="<<i<<" but it should be at least "<<cmp<<"!"<<std::endl;
    throw std::runtime_error(ss.str());
  }
}

}//namespace
