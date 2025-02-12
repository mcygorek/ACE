#include "MPS_Matrix.hpp"
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include "CheckMatrix.hpp"

namespace ACE{


template <typename T>
  void MPS_Matrix_ScalarType<T>::resize(int dim_i_, int dim_d1_, int dim_d2_){
    dim_i=dim_i_;
    dim_d1=dim_d1_;
    dim_d2=dim_d2_;
    deallocate();
    allocate();
  }
  
template <typename T>
  void MPS_Matrix_ScalarType<T>::fill(const T &c){
    for(int i=0; i<dim_i*dim_d1*dim_d2; i++)mem[i]=c;
  }

template <typename T>
  void MPS_Matrix_ScalarType<T>::resize_fill_one(int dim_i_, int dim_d1_, int dim_d2_){
    resize(dim_i_, dim_d1_, dim_d2_);
    fill(1.);
  }
template <typename T>
  void MPS_Matrix_ScalarType<T>::set_zero(){
    //for(int i=0; i<dim_i*dim_d1*dim_d2; i++)mem[i]=0.;
    memset(mem, 0., sizeof(T)*dim_i*dim_d1*dim_d2);
  }

template <typename T>
  double MPS_Matrix_ScalarType<T>::max_element_abs()const{
    double a=0;
    int end=dim_d1*dim_d2*dim_i;
    for(int x=0; x<end; x++){
      double this_a=abs(mem[x]);
      if(this_a>a)a=this_a;
      if(!std::isfinite(this_a))return this_a;
    }
    return a; 
  }

 template<>
  double MPS_Matrix_ScalarType<double>::norm2_d1(int d1)const{
    double n2=0;
    for(int i=0; i<dim_i; i++){
      for(int d2=0; d2<dim_d2; d2++){
        double tmp=mem[(d1*dim_i+i)*dim_d2+d2];
        n2+=tmp*tmp;
      }
    }
    return n2;  
  }
template<>
  double MPS_Matrix_ScalarType<std::complex<double> >::norm2_d1(int d1)const{
    double n2=0;
    for(int i=0; i<dim_i; i++){
      for(int d2=0; d2<dim_d2; d2++){
        std::complex<double> tmp=mem[(d1*dim_i+i)*dim_d2+d2];
        n2+=tmp.real()*tmp.real()+tmp.imag()*tmp.imag();
      }
    }
    return n2;  
  }

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
                      MPS_Matrix_ScalarType<T>::get_Matrix_d1i_d2()const{

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> M(dim_d1*dim_i, dim_d2);
  for(int d1=0; d1<dim_d1; d1++){
    for(int i=0; i<dim_i; i++){
      for(int d2=0; d2<dim_d2; d2++){
        M(d1*dim_i+i, d2) = mem[(d1*dim_i+i)*dim_d2+d2];
      }
    }
  } 
  return M;
}
template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
                      MPS_Matrix_ScalarType<T>::get_Matrix_d1_id2()const{

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> M(dim_d1, dim_i*dim_d2);
  for(int d1=0; d1<dim_d1; d1++){
    for(int i=0; i<dim_i; i++){
      for(int d2=0; d2<dim_d2; d2++){
        M(d1, i*dim_d2+d2) = mem[(d1*dim_i+i)*dim_d2+d2];
      }
    }
  } 
  return M;
}
template <typename T> void MPS_Matrix_ScalarType<T>::set_from_Matrix_d1i_d2(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &M, int dimi ){

  resize( dimi, M.rows()/dimi, M.cols() );
  for(int d1=0; d1<dim_d1; d1++){
    for(int i=0; i<dim_i; i++){
      for(int d2=0; d2<dim_d2; d2++){
        mem[(d1*dim_i+i)*dim_d2+d2] = M(d1*dim_i+i, d2);
      }
    }
  } 
}
template <typename T> void MPS_Matrix_ScalarType<T>::set_from_Matrix_d1_id2(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &M, int dimi ){

  resize( dimi, M.rows(), M.cols()/dimi );
  for(int d1=0; d1<dim_d1; d1++){
    for(int i=0; i<dim_i; i++){
      for(int d2=0; d2<dim_d2; d2++){
        mem[(d1*dim_i+i)*dim_d2+d2] = M(d1, i*dim_d2+d2);
      }
    }
  } 
}


template <typename T>
  void MPS_Matrix_ScalarType<T>::print_dims(std::ostream &os)const{
    os<<dim_i<<" "<<dim_d1<<" "<<dim_d2;
  }

template <typename T>
  void MPS_Matrix_ScalarType<T>::print_HR(const std::string &fname, double thr)const{
    std::ofstream ofs(fname.c_str());
    ofs<<"Dimensions: "<<dim_i<<" "<<dim_d1<<" "<<dim_d2<<std::endl;
    for(int i=0; i<dim_i; i++){
      for(int d1=0; d1<dim_d1; d1++){
        for(int d2=0; d2<dim_d2; d2++){
          if(thr<=0. || abs(operator()(i,d1,d2))>thr){
            ofs<<i<<" "<<d1<<" "<<d2<<" "<<operator()(i,d1,d2)<<std::endl;
          }
        }
      }
    }
  }

template <typename T>
  void MPS_Matrix_ScalarType<T>::copy(const MPS_Matrix_ScalarType<T> &other){
    resize(other.dim_i, other.dim_d1, other.dim_d2);
    for(int i=0; i<dim_i*dim_d1*dim_d2; i++){
      mem[i]=other.mem[i];
    }
  }

template <typename T>
  void MPS_Matrix_ScalarType<T>::swap(MPS_Matrix_ScalarType<T> &other){
    int dim_i_=dim_i;
    int dim_d1_=dim_d1;
    int dim_d2_=dim_d2;
    T *mem_=mem;
    dim_i=other.dim_i;
    dim_d1=other.dim_d1;
    dim_d2=other.dim_d2;
    mem=other.mem;
    other.dim_i=dim_i_;
    other.dim_d1=dim_d1_;
    other.dim_d2=dim_d2_;
    other.mem=mem_;
  }

template <typename T>
  void MPS_Matrix_ScalarType<T>::inner_multiply_left(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & M){
    if(dim_d1!=M.cols()){
      std::cerr<<"MPS_Matrix::inner_multiply_left: dim_d1!=M.cols() ("<<dim_d1<<" vs. "<<M.cols()<<")!"<<std::endl;
      exit(1);
    }
    int new_d1=M.rows();
    MPS_Matrix_ScalarType<T> A(dim_i, new_d1, dim_d2); 
/*
    for(int i=0; i<dim_i; i++){
      for(int r=0; r<new_d1; r++){
        for(int d2=0; d2<dim_d2; d2++){
          A(i,r,d2)=0.;
          for(int d1=0; d1<dim_d1; d1++){
            A(i,r,d2)+=M(r,d1)*operator()(i,d1,d2);
          }
        }
      }
    }
*/
    A.set_zero();
    for(int i=0; i<dim_i; i++){
      Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 0, Eigen::OuterStride<> >(A.mem+i*dim_d2,new_d1,dim_d2,Eigen::OuterStride<>(dim_i*dim_d2)).noalias() += M * \
      Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 0, Eigen::OuterStride<> >(mem+i*dim_d2,dim_d1,dim_d2,Eigen::OuterStride<>(dim_i*dim_d2));
    }

    swap(A);
  }

template <typename T>
  void MPS_Matrix_ScalarType<T>::inner_multiply_right(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & M){
    if(dim_d2!=M.rows()){
      std::cerr<<"MPS_Matrix::inner_multiply_right: dim_d2!=M.rows() ("<<dim_d2<<" vs. "<<M.rows()<<")!"<<std::endl;
      exit(1);
    }

    int new_d2=M.cols();
    MPS_Matrix_ScalarType<T> A(dim_i, dim_d1, new_d2); 
/*
    for(int i=0; i<dim_i; i++){
      for(int d1=0; d1<dim_d1; d1++){
        for(int c=0; c<new_d2; c++){
          A(i,d1,c)=0.;
          for(int d2=0; d2<dim_d2; d2++){
            A(i,d1,c)+=operator()(i,d1,d2)*M(d2,c);
          }
        }
      }
    }
*/
    A.set_zero();
    for(int i=0; i<dim_i; i++){
      Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 0, Eigen::OuterStride<> >(A.mem+i*new_d2,dim_d1,new_d2,Eigen::OuterStride<>(dim_i*new_d2)).noalias() += \
      Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 0, Eigen::OuterStride<> >(mem+i*dim_d2,dim_d1,dim_d2,Eigen::OuterStride<>(dim_i*dim_d2))  *  M;
    }

    swap(A);
  }

template <typename T>
  void MPS_Matrix_ScalarType<T>::contract_front(){
    if(dim_d1<2)return;
    MPS_Matrix_ScalarType<T> A(dim_i, 1, dim_d2);
    A.set_zero();
    for(int i=0; i<dim_i; i++){
      for(int d1=0; d1<dim_d1; d1++){
        for(int d2=0; d2<dim_d2; d2++){
          A(i,0,d2)+=operator()(i,d1,d2);
        }
      }
    }
    swap(A);
  }

template <typename T>
  void MPS_Matrix_ScalarType<T>::contract_back(){
    if(dim_d2<2)return;
    MPS_Matrix_ScalarType<T> A(dim_i, dim_d1, 1);
    A.set_zero();
    for(int i=0; i<dim_i; i++){
      for(int d1=0; d1<dim_d1; d1++){
        for(int d2=0; d2<dim_d2; d2++){
          A(i,d1,0)+=operator()(i,d1,d2);
        }
      }
    }
    swap(A);
  }

template <typename T>
  void MPS_Matrix_ScalarType<T>::multiply(const MPS_Matrix_ScalarType<T> &other){
    if(other.dim_i!=dim_i){
      std::cerr<<"MPS_Matrix::multiply: other.dim_i!=dim_i!"<<std::endl;
      exit(1);
    }

    MPS_Matrix_ScalarType<T> A(dim_i, dim_d1*other.dim_d1, dim_d2*other.dim_d2);
#ifdef PRINT_MPS_MATRIX_MULTIPLY_DIMS
std::cout<<"MPS_Matrix::multiply: dims: "<<A.dim_i<<" "<<A.dim_d1<<" "<<A.dim_d2<<" <- "<<dim_d1<<"*"<<other.dim_d1<<", "<<dim_d2<<"*"<<other.dim_d2<<std::endl;
#endif
    //A.set_zero();
    for(int i=0; i<dim_i; i++){
      for(int d1=0; d1<dim_d1; d1++){
        for(int d2=0; d2<dim_d2; d2++){
          for(int od1=0; od1<other.dim_d1; od1++){
            for(int od2=0; od2<other.dim_d2; od2++){
              A(i, d1*other.dim_d1 + od1, d2*other.dim_d2 + od2)=
                 operator()(i,d1,d2) * other(i,od1,od2);
            }
          }
        }
      }
    }
    swap(A);
  }

template <typename T>
 void MPS_Matrix_ScalarType<T>::read_binary(std::istream &ifs){
    ifs.read((char*)&dim_i, sizeof(int));
    ifs.read((char*)&dim_d1, sizeof(int));
    ifs.read((char*)&dim_d2, sizeof(int));
    resize(dim_i, dim_d1, dim_d2);
    
    ifs.read((char*)mem, sizeof(T)*dim_i*dim_d1*dim_d2);
}
template <typename T>
  void MPS_Matrix_ScalarType<T>::write_binary(std::ostream &ofs)const{
    ofs.write((char*)&dim_i, sizeof(int));
    ofs.write((char*)&dim_d1, sizeof(int));
    ofs.write((char*)&dim_d2, sizeof(int));
    ofs.write((char*)mem, sizeof(T)*dim_i*dim_d1*dim_d2);
}


template <typename T>
std::ostream &operator<<(std::ostream &os, const MPS_Matrix_ScalarType<T> &a){
  a.print_dims(os);
  return os;
}

template class MPS_Matrix_ScalarType<double>;
template class MPS_Matrix_ScalarType<std::complex<double>>;
}//namespace

