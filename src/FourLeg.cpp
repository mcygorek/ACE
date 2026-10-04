#include "FourLeg.hpp"
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include "CheckMatrix.hpp"
#include "DummyException.hpp"

namespace ACE{


template <typename T>
  void FourLeg_ScalarType<T>::resize(int dim_d1_, int dim_d2_, int dim_d3_, int dim_d4_){
    dim_d1=dim_d1_;
    dim_d2=dim_d2_;
    dim_d3=dim_d3_;
    dim_d4=dim_d4_;
    deallocate();
    allocate();
  }
  
template <typename T>
  void FourLeg_ScalarType<T>::fill(const T &c){
    for(int i=0; i<dim_d1*dim_d2*dim_d3*dim_d4; i++)mem[i]=c;
  }

template <typename T>
  void FourLeg_ScalarType<T>::resize_fill_one(int dim_d1_, int dim_d2_, int dim_d3_, int dim_d4_){
    resize(dim_d1_, dim_d2_, dim_d3_, dim_d4_);
    fill(1.);
  }
template <typename T>
  void FourLeg_ScalarType<T>::set_zero(){
    //for(int i=0; i<dim_i*dim_d1*dim_d2; i++)mem[i]=0.;
    memset(mem, 0., sizeof(T)*dim_d1*dim_d2*dim_d3*dim_d4);
  }

template <typename T>
  double FourLeg_ScalarType<T>::max_element_abs()const{
    double a=0;
    int end=dim_d1*dim_d2*dim_d3*dim_d4;
    for(int x=0; x<end; x++){
      double this_a=abs(mem[x]);
      if(this_a>a)a=this_a;
      if(!std::isfinite(this_a))return this_a;
    }
    return a; 
  }

/*
template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
                      FourLeg_ScalarType<T>::get_Matrix_d1i_d2()const{

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
                      FourLeg_ScalarType<T>::get_Matrix_d1_id2()const{

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
template <typename T> void FourLeg_ScalarType<T>::set_from_Matrix_d1i_d2(
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
template <typename T> void FourLeg_ScalarType<T>::set_from_Matrix_d1_id2(
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
*/


template <typename T>
  void FourLeg_ScalarType<T>::print_dims(std::ostream &os)const{
    os<<dim_d1<<" "<<dim_d2<<" "<<dim_d3<<" "<<dim_d4;
  }

template <typename T>
  void FourLeg_ScalarType<T>::print_HR(const std::string &fname, double thr)const{
    std::ofstream ofs(fname.c_str());
    ofs<<"Dimensions: "<<dim_d1<<" "<<dim_d2<<" "<<dim_d3<<" "<<dim_d4<<std::endl;
    for(int d1=0; d1<dim_d1; d1++){
      for(int d2=0; d2<dim_d2; d2++){
        for(int d3=0; d3<dim_d3; d3++){
          for(int d4=0; d4<dim_d4; d4++){
            if(thr<=0. || abs(operator()(d1,d2,d3,d4))>thr){
              ofs<<d1<<" "<<d2<<" "<<d3<<" "<<d4<<" "<<operator()(d1,d2,d3,d4)<<std::endl;
            }
          }
        }
      }
    }
  }

template <typename T>
  void FourLeg_ScalarType<T>::copy(const FourLeg_ScalarType<T> &other){
    resize(other.dim_d1, other.dim_d2, other.dim_d3, other.dim_d4);
    for(int i=0; i<dim_d1*dim_d2*dim_d3*dim_d4; i++){
      mem[i]=other.mem[i];
    }
  }

template <typename T>
  void FourLeg_ScalarType<T>::resize_copy(int new_d1, int new_d2, int new_d3, int new_d4, const FourLeg_ScalarType<T> &other){
    resize(new_d1, new_d2, new_d3, new_d4);
    set_zero();
    int min_d1=new_d1; if(other.dim_d1<min_d1){min_d1=other.dim_d1;}
    int min_d2=new_d2; if(other.dim_d2<min_d2){min_d2=other.dim_d2;}
    int min_d3=new_d3; if(other.dim_d3<min_d3){min_d3=other.dim_d3;}
    int min_d4=new_d4; if(other.dim_d4<min_d4){min_d4=other.dim_d4;}
    for(int d1=0; d1<min_d1; d1++){
      for(int d2=0; d2<min_d2; d2++){
        for(int d3=0; d3<min_d3; d3++){
          for(int d4=0; d4<min_d4; d4++){
            operator()(d1,d2,d3,d4)=other(d1,d2,d3,d4);
          }
        }
      } 
    }
  }


template <typename T>
  void FourLeg_ScalarType<T>::swap(FourLeg_ScalarType<T> &other){
    int dim_d1_=dim_d1;
    int dim_d2_=dim_d2;
    int dim_d3_=dim_d3;
    int dim_d4_=dim_d4;
    T *mem_=mem;
    dim_d1=other.dim_d1;
    dim_d2=other.dim_d2;
    dim_d3=other.dim_d3;
    dim_d4=other.dim_d4;
    mem=other.mem;
    other.dim_d1=dim_d1_;
    other.dim_d2=dim_d2_;
    other.dim_d3=dim_d3_;
    other.dim_d4=dim_d4_;
    other.mem=mem_;
  }


template <typename T>
  void FourLeg_ScalarType<T>::multiply_d1(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & M){
    if(dim_d1!=M.rows()){
      std::cerr<<"FourLeg::multiply_d1: dim_d1!=M.rows() ("<<dim_d1<<" vs. "<<M.rows()<<")!"<<std::endl;
      exit(1);
    }

    int new_d1=M.cols();
    FourLeg_ScalarType<T> A(new_d1, dim_d2, dim_d3, dim_d4); 
    A.set_zero();
    for(int d2=0; d2<dim_d2; d2++){
      for(int d3=0; d3<dim_d3; d3++){
        for(int d4=0; d4<dim_d4; d4++){
          for(int c=0; c<new_d1; c++){
            for(int d1=0; d1<dim_d1; d1++){
              A(c,d2,d3,d4)+=operator()(d1,d2,d3,d4)*M(d1,c);
            }
          }
        }
      }
    }
    swap(A);
  }


template <typename T>
  void FourLeg_ScalarType<T>::multiply_d2(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & M){
    if(dim_d2!=M.rows()){
      std::cerr<<"FourLeg::multiply_d2: dim_d2!=M.rows() ("<<dim_d2<<" vs. "<<M.rows()<<")!"<<std::endl;
      exit(1);
    }

    int new_d2=M.cols();
    FourLeg_ScalarType<T> A(dim_d1, new_d2, dim_d3, dim_d4); 
    A.set_zero();
    for(int d1=0; d1<dim_d1; d1++){
      for(int d3=0; d3<dim_d3; d3++){
        for(int d4=0; d4<dim_d4; d4++){
          for(int c=0; c<new_d2; c++){
            for(int d2=0; d2<dim_d2; d2++){
              A(d1,c,d3,d4)+=operator()(d1,d2,d3,d4)*M(d2,c);
            }
          }
        }
      }
    }
    swap(A);
  }


template <typename T>
  void FourLeg_ScalarType<T>::multiply_d3(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & M){
    if(dim_d3!=M.rows()){
      std::cerr<<"FourLeg::multiply_d3: dim_d3!=M.rows() ("<<dim_d3<<" vs. "<<M.rows()<<")!"<<std::endl;
      exit(1);
    }

    int new_d3=M.cols();
    FourLeg_ScalarType<T> A(dim_d1, dim_d2, new_d3, dim_d4); 
    A.set_zero();
    for(int d1=0; d1<dim_d1; d1++){
      for(int d2=0; d2<dim_d2; d2++){
        for(int d4=0; d4<dim_d4; d4++){
          for(int c=0; c<new_d3; c++){
            for(int d3=0; d3<dim_d3; d3++){
              A(d1,d2,c,d4)+=operator()(d1,d2,d3,d4)*M(d3,c);
            }
          }
        }
      }
    }
    swap(A);
  }


template <typename T>
  void FourLeg_ScalarType<T>::multiply_d4(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & M){
    if(dim_d4!=M.rows()){
      std::cerr<<"FourLeg::multiply_d4: dim_d4!=M.rows() ("<<dim_d4<<" vs. "<<M.rows()<<")!"<<std::endl;
      exit(1);
    }

    using RowMatrixMap = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >;
    int new_d4=M.cols();
    FourLeg_ScalarType<T> A(dim_d1, dim_d2, dim_d3, new_d4); 
/*
    A.set_zero();
    for(int d1=0; d1<dim_d1; d1++){
      for(int d2=0; d2<dim_d2; d2++){
        for(int d3=0; d3<dim_d3; d3++){
          for(int c=0; c<new_d4; c++){
            for(int d4=0; d4<dim_d4; d4++){
              A(d1,d2,d3,c)+=operator()(d1,d2,d3,d4)*M(d4,c);
    } } } } }
*/
    RowMatrixMap(A.mem, dim_d1*dim_d2*dim_d3, new_d4) = 
      RowMatrixMap(mem, dim_d1*dim_d2*dim_d3, dim_d4) * M;
    swap(A);
  }


template <typename T>
  void FourLeg_ScalarType<T>::join_d3_d1(const FourLeg_ScalarType<T> &other){

  if(dim_d3!=other.dim_d1){
    std::cerr<<"FourLeg::join_d3_d1: dim_d3!=other.dim_d1 ("<<dim_d3<<" vs. "<<other.dim_d1<<")!"<<std::endl;
    exit(1);
  }

  FourLeg_ScalarType<T> tmp(dim_d1, dim_d2*other.dim_d2, other.dim_d3, dim_d4*other.dim_d4);
  tmp.set_zero();
  for(int d1=0; d1<dim_d1; d1++){
   for(int d2=0; d2<dim_d2; d2++){
    for(int od2=0; od2<other.dim_d2; od2++){
     for(int d4=0; d4<dim_d4; d4++){
      for(int od4=0; od4<other.dim_d4; od4++){
       for(int od3=0; od3<other.dim_d3; od3++){
        for(int d3=0; d3<dim_d3; d3++){
//          tmp.mem[((((d1*dim_d2+d2)*other.dim_d2+od2)*other.dim_d3+od3)*dim_d4+d4)*other.dim_d4+od4] += operator()(d1, d2, d3, d4) * other(d3, od2, od3, od4);
          tmp(d1, d2*other.dim_d2+od2, od3, d4*other.dim_d4+od4) += 
                 operator()(d1, d2, d3, d4) * other(d3, od2, od3, od4);
        } 
       }
      }
     }
    }
   }
  }
  tmp.swap(*this);  
} 
template <typename T>
  void FourLeg_ScalarType<T>::join_d3_d1_multiply_d2(const FourLeg_ScalarType<T> &other, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & P){
  // assuming dim(d1,d3) << dim(d2,d4) and *this.dim_d4<other.dim_d4, it is
  // best to first contract P with *this, and then the result with other
  // Note: multiply defined over first index (row) of P
  if(P.rows() != dim_d2*other.dim_d2){
    std::cerr<<"FourLeg::join_d3_d1_multiply_d2: P.rows() != dim_d2*other.dim_d2 ("<<P.cols()<<" vs. "<<dim_d2<<"*"<<other.dim_d2<<")!"<<std::endl;
    throw DummyException();
  }
  // contract (P *this)
  using RowMatrixMap = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >;
  using OuterStrideMap = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 0, Eigen::OuterStride<Eigen::Dynamic> >;
  using cInnerStrideMap = Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 0, Eigen::InnerStride<Eigen::Dynamic> >;
  using DynamicStrideMap = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> >;


  FourLeg_ScalarType<T> tmp(dim_d1, P.cols()*other.dim_d2, dim_d3, dim_d4);
  tmp.set_zero();
/*
  for(int d1=0; d1<dim_d1; d1++){
    for(int c=0; c<P.cols(); c++){
      for(int od2=0; od2<other.dim_d2; od2++){
        for(int d3=0; d3<dim_d3; d3++){
          for(int d4=0; d4<dim_d4; d4++){
            for(int d2=0; d2<dim_d2; d2++){
              tmp(d1, c*other.dim_d2+od2, d3, d4)+=P(d2*other.dim_d2+od2, c)*operator()(d1,d2,d3,d4);
  } } } } } }
*/
  
  for(int d1=0; d1<dim_d1; d1++){
    for(int od2=0; od2<other.dim_d2; od2++){
      OuterStrideMap(tmp.mem+(d1*tmp.dim_d2+od2)*dim_d3*dim_d4, P.cols(), dim_d3*dim_d4, Eigen::OuterStride<Eigen::Dynamic>(other.dim_d2*dim_d3*dim_d4)) =
        cInnerStrideMap(P.data()+od2, P.cols(), dim_d2, Eigen::InnerStride<Eigen::Dynamic>(other.dim_d2)) *
        RowMatrixMap(mem+d1*dim_d2*dim_d3*dim_d4, dim_d2, dim_d3*dim_d4);
  } } 
   

  if(dim_d3!=other.dim_d1){
    std::cerr<<"FourLeg::join_d3_d1_multiply_d2: dim_d3!=other.dim_d1 ("<<dim_d3<<" vs. "<<other.dim_d1<<")!"<<std::endl;
    throw DummyException();
  }
  int old_dim_d4=dim_d4;
  resize(tmp.dim_d1, P.cols(), other.dim_d3, old_dim_d4*other.dim_d4);
  set_zero();
/*
  for(int d1=0; d1<dim_d1; d1++){
    for(int c=0; c<P.cols(); c++){
      for(int od2=0; od2<other.dim_d2; od2++){
        for(int d4=0; d4<old_dim_d4; d4++){
          for(int od4=0; od4<other.dim_d4; od4++){
            for(int od3=0; od3<other.dim_d3; od3++){
              for(int od1=0; od1<other.dim_d1; od1++){
//                operator()(d1, c, od3, d4*other.dim_d4+od4) += 
//tmp(d1, c*other.dim_d2+od2, od1, d4) * other(od1, od2, od3, od4);
      mem[(d1*dim_d2+c)*dim_d3*dim_d4+ (od3*dim_d4+d4*other.dim_d4) +od4]+= 
tmp.mem[(d1*dim_d2+c)*other.dim_d2*dim_d3*old_dim_d4+ (od2*dim_d3*old_dim_d4+d4) +od1*old_dim_d4] * 
other.mem[od1*other.dim_d2*other.dim_d3*other.dim_d4+ (od2*other.dim_d3+od3)*other.dim_d4 +od4];
  } } } } } } }
*/
  for(int od2=0; od2<other.dim_d2; od2++){
    for(int d4=0; d4<old_dim_d4; d4++){
      for(int od3=0; od3<other.dim_d3; od3++){
        OuterStrideMap(mem+od3*dim_d4+d4*other.dim_d4, dim_d1*P.cols(), other.dim_d4, Eigen::OuterStride<Eigen::Dynamic>(dim_d3*dim_d4)) +=

          DynamicStrideMap(tmp.mem+od2*tmp.dim_d3*tmp.dim_d4+d4, dim_d1*P.cols(), other.dim_d1, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(other.dim_d2*other.dim_d1*old_dim_d4, old_dim_d4)) * 

          OuterStrideMap(other.mem+(od2*other.dim_d3+od3)*other.dim_d4, other.dim_d1, other.dim_d4, Eigen::OuterStride<Eigen::Dynamic>(other.dim_d2*other.dim_d3*other.dim_d4));
  } } } 

} 


template <typename T>
  MPS_Matrix_ScalarType<T> FourLeg_ScalarType<T>::reduce_fix_d1(int d1) const{
    if(d1>=dim_d1){
      std::cerr<<"FourLeg::reduce_fix_d1: d1>=dim_d1 ("<<d1<<" vs. "<<dim_d1<<")!"<<std::endl;
      exit(1);
    }
  MPS_Matrix_ScalarType<T> M(dim_d3, dim_d2, dim_d4);
  for(int d2=0; d2<dim_d2; d2++){
    for(int d3=0; d3<dim_d3; d3++){
      for(int d4=0; d4<dim_d4; d4++){
        M(d3,d2,d4)=operator()(d1,d2,d3,d4);
      }
    }
  }
  return M;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
                      FourLeg_ScalarType<T>::get_Matrix_d1_d2d3d4()const{

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> M(dim_d1, dim_d2*dim_d3*dim_d4);
/*
  for(int d1=0; d1<dim_d1; d1++){
    for(int d2=0; d2<dim_d2; d2++){
      for(int d3=0; d3<dim_d3; d3++){
        for(int d4=0; d4<dim_d4; d4++){
          M(d1,(d2*dim_d3+d3)*dim_d4+d4) = operator()(d1,d2,d3,d4);
  } } } } 
*/
  using RowMatrixMap = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >;
  M = RowMatrixMap(mem, dim_d1, dim_d2*dim_d3*dim_d4);
  return M;
}
template <typename T> void FourLeg_ScalarType<T>::set_from_Matrix_d1_d2d3d4(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &M, int dim2, int dim3 ){

  resize( M.rows(), dim2, dim3, M.cols()/(dim2*dim3) );
/*
  for(int d1=0; d1<dim_d1; d1++){
    for(int d2=0; d2<dim_d2; d2++){
      for(int d3=0; d3<dim_d3; d3++){
        for(int d4=0; d4<dim_d4; d4++){
          operator()(d1,d2,d3,d4) = M(d1,(d2*dim_d3+d3)*dim_d4+d4);
  } } } } 
*/
  using RowMatrixMap = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >;
  RowMatrixMap(mem, dim_d1, dim_d2*dim_d3*dim_d4) = M;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
                      FourLeg_ScalarType<T>::get_Matrix_d1d2d3_d4()const{

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> M(dim_d1*dim_d2*dim_d3, dim_d4);
/*  
  for(int d1=0; d1<dim_d1; d1++){
    for(int d2=0; d2<dim_d2; d2++){
      for(int d3=0; d3<dim_d3; d3++){
        for(int d4=0; d4<dim_d4; d4++){
          M(((d1*dim_d2)+d2)*dim_d3+d3, d4) = operator()(d1,d2,d3,d4);
  } } } }
*/
  using RowMatrixMap = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >;
  M = RowMatrixMap(mem, dim_d1*dim_d2*dim_d3, dim_d4);
  return M;
}
template <typename T> void FourLeg_ScalarType<T>::set_from_Matrix_d1d2d3_d4(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &M, int dim2, int dim3 ){

  resize( M.rows()/(dim2*dim3), dim2, dim3, M.cols() );
/*
  for(int d1=0; d1<dim_d1; d1++){
    for(int d2=0; d2<dim_d2; d2++){
      for(int d3=0; d3<dim_d3; d3++){
        for(int d4=0; d4<dim_d4; d4++){
          operator()(d1,d2,d3,d4) = M((d1*dim_d2+d2)*dim_d3+d3, d4);
  } } } }
*/  
 using RowMatrixMap = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >;
 RowMatrixMap(mem, dim_d1*dim_d2*dim_d3, dim_d4) = M;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
                      FourLeg_ScalarType<T>::get_Matrix_d2_d1d3d4()const{

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> M(dim_d2,dim_d1*dim_d3*dim_d4);
/*
  for(int d1=0; d1<dim_d1; d1++){
    for(int d2=0; d2<dim_d2; d2++){
      for(int d3=0; d3<dim_d3; d3++){
        for(int d4=0; d4<dim_d4; d4++){
          M(d2,((d1*dim_d3)+d3)*dim_d4+d4) = operator()(d1,d2,d3,d4);
  } } } } 
*/
  using RowMatrixMap = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >;
  using ColMatrixMap = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> >;

  for(int d1=0; d1<dim_d1; d1++){
    ColMatrixMap(M.data()+d1*dim_d2*dim_d3*dim_d4, dim_d2, dim_d3*dim_d4) =
      RowMatrixMap(mem+d1*dim_d2*dim_d3*dim_d4, dim_d2, dim_d3*dim_d4);
  }

  return M;
}
template <typename T> void FourLeg_ScalarType<T>::set_from_Matrix_d2_d1d3d4(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &M, int dim3, int dim4 ){

  resize( M.cols()/(dim3*dim4), M.rows(), dim3, dim4);
/*
  for(int d1=0; d1<dim_d1; d1++){
    for(int d2=0; d2<dim_d2; d2++){
      for(int d3=0; d3<dim_d3; d3++){
        for(int d4=0; d4<dim_d4; d4++){
          operator()(d1,d2,d3,d4) = M(d2, (d1*dim_d3+d3)*dim_d4+d4);
  } } } }  
*/
  using RowMatrixMap = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >;
  using cColMatrixMap = Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> >;

  for(int d1=0; d1<dim_d1; d1++){
    RowMatrixMap(mem+d1*dim_d2*dim_d3*dim_d4, dim_d2, dim_d3*dim_d4) = 
      cColMatrixMap(M.data()+d1*dim_d2*dim_d3*dim_d4, dim_d2, dim_d3*dim_d4);
  }
}

template <typename T> 
  Eigen::VectorXcd FourLeg_ScalarType<T>::get_reduced_d3(const Eigen::VectorXcd &cap_d1, const Eigen::VectorXcd &cap_d2, const Eigen::VectorXcd &cap_d4)const{

  if(cap_d1.rows()!=dim_d1){
    std::cerr<<"FourLeg::get_reduced_d3: cap_d1.rows()!=dim_d1!"<<std::endl;
    throw DummyException();
  }
  if(cap_d2.rows()!=dim_d2){
    std::cerr<<"FourLeg::get_reduced_d3: cap_d2.rows()!=dim_d2!"<<std::endl;
    throw DummyException();
  }
  if(cap_d4.rows()!=dim_d4){
    std::cerr<<"FourLeg::get_reduced_d3: cap_d4.rows()!=dim_d4!"<<std::endl;
    throw DummyException();
  }
  Eigen::VectorXcd rho=Eigen::VectorXcd::Zero(dim_d3);
  for(int d1=0; d1<dim_d1; d1++){
    for(int d2=0; d2<dim_d2; d2++){
      for(int d3=0; d3<dim_d3; d3++){
        for(int d4=0; d4<dim_d4; d4++){
          rho(d3)+=cap_d1(d1)*cap_d2(d2)*cap_d4(d4)*operator()(d1,d2,d3,d4);
        }
      }
    }
  }
  return rho;
}


template <typename T>
 void FourLeg_ScalarType<T>::read_binary(std::istream &ifs){
    ifs.read((char*)&dim_d1, sizeof(int));
    ifs.read((char*)&dim_d2, sizeof(int));
    ifs.read((char*)&dim_d3, sizeof(int));
    ifs.read((char*)&dim_d4, sizeof(int));
    resize(dim_d1, dim_d2, dim_d3, dim_d4);
    
    ifs.read((char*)mem, sizeof(T)*dim_d1*dim_d2*dim_d3*dim_d4);
}
template <typename T>
  void FourLeg_ScalarType<T>::write_binary(std::ostream &ofs)const{
    ofs.write((char*)&dim_d1, sizeof(int));
    ofs.write((char*)&dim_d2, sizeof(int));
    ofs.write((char*)&dim_d3, sizeof(int));
    ofs.write((char*)&dim_d4, sizeof(int));
    ofs.write((char*)mem, sizeof(T)*dim_d1*dim_d2*dim_d3*dim_d4);
}

template <typename T>
  FourLeg_ScalarType<T> FourLeg_ScalarType<T>::Identity_d1_d3(int d){
    FourLeg_ScalarType<T> ret(d, 1, d, 1);
    ret.fill(0.);
    for(int d1=0; d1<d; d1++){
      ret(d1, 0, d1, 0)=1.;
    }
    return ret;
}
template <typename T>
  FourLeg_ScalarType<T> FourLeg_ScalarType<T>::Operator_d1_d3(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & M){
  FourLeg_ScalarType<T> ret(M.cols(), 1, M.rows(), 1);
  for(int j=0; j<M.rows(); j++){
    for(int i=0; i<M.cols(); i++){
      ret(i, 0, j, 0) = M(j, i);
    }
  }
  return ret;
}
template <typename T>
  FourLeg_ScalarType<T> FourLeg_ScalarType<T>::Operator_forward_d1_d3(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & M){
  FourLeg_ScalarType<T> ret(M.cols()*M.cols(), 1, M.rows()*M.cols(), 1);
  ret.fill(0.);
  for(int d=0; d<M.rows(); d++){
    for(int j=0; j<M.rows(); j++){
      for(int i=0; i<M.cols(); i++){
//multiplication col->row corresponds to d1->d3
        ret(i*M.cols()+d, 0, j*M.cols()+d, 0) = M(j, i);
      }
    }
  }
  return ret;
}
template <typename T>
  FourLeg_ScalarType<T> FourLeg_ScalarType<T>::Operator_backward_d1_d3(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & M){
  FourLeg_ScalarType<T> ret(M.rows()*M.rows(), 1, M.rows()*M.cols(), 1);
  ret.fill(0.);
  for(int d=0; d<M.rows(); d++){
    for(int i=0; i<M.rows(); i++){
      for(int j=0; j<M.cols(); j++){
//"backward" means multiplication from the right -> d1=row, d3=col
        ret(d*M.cols()+i, 0, d*M.rows()+j, 0) = M(i, j);
      }
    }
  }
  return ret;
}
template <typename T>
  FourLeg_ScalarType<T> FourLeg_ScalarType<T>::Operator_forward_backward_d1_d3(
           const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & op1,
           const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & op2){
  FourLeg_ScalarType<T> ret(op1.cols()*op2.rows(), 1, op1.rows()*op2.cols(), 1);
  ret.fill(0.);
  for(int j1=0; j1<op1.rows(); j1++){
    for(int i1=0; i1<op1.cols(); i1++){
      for(int j2=0; j2<op2.rows(); j2++){
        for(int i2=0; i2<op2.cols(); i2++){
          ret(i1*op2.rows()+i2, 0, j1*op2.cols()+j2, 0) = op1(j1,i1)*op2(i2,j2);
        }
      }
    }
  }
  return ret;
}


template <typename T>
std::ostream &operator<<(std::ostream &os, const FourLeg_ScalarType<T> &a){
  a.print_dims(os);
  return os;
}


template class FourLeg_ScalarType<double>;
template class FourLeg_ScalarType<std::complex<double>>;

FourLeg MPS_Matrix_to_FourLeg(const MPS_Matrix &M, int dim_d2){
    int dim_d4=M.dim_d1/dim_d2;
    if(dim_d4*dim_d2!=M.dim_d1){
      std::cerr<<"MPS_Matrix_to_Fourleg: dim_d4*dim_d2!=M.dim_d1!"<<std::endl;
      exit(1);
    }
    FourLeg f(M.dim_d2, dim_d2, M.dim_i, dim_d4);

    /*for(int d1=0; d1<f.dim_d1; d1++){
      for(int d2=0; d2<f.dim_d2; d2++){
        for(int d3=0; d3<f.dim_d3; d3++){
          for(int d4=0; d4<f.dim_d4; d4++){
            f(d1, d2, d3, d4)=M(d3, d2*dim_d4+d4, d1);
    } } } } */
    for(int d1=0; d1<f.dim_d1; d1++){
      for(int d3=0; d3<f.dim_d3; d3++){
        Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 0, Eigen::OuterStride<Eigen::Dynamic> >(f.mem+(d1*f.dim_d2*f.dim_d3+d3)*f.dim_d4, f.dim_d2, f.dim_d4, Eigen::OuterStride<Eigen::Dynamic>(f.dim_d3*f.dim_d4))
        =
        Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 0, Eigen::InnerStride<Eigen::Dynamic> >(M.mem+d3*M.dim_d2+d1, f.dim_d2, f.dim_d4, Eigen::InnerStride<Eigen::Dynamic>(M.dim_i*M.dim_d2)); 
    } }
    return f;
}
MPS_Matrix FourLeg_to_MPS_Matrix(const FourLeg &f){
    MPS_Matrix M(f.dim_d3, f.dim_d2*f.dim_d4, f.dim_d1);

/*    for(int d1=0; d1<f.dim_d1; d1++){
      for(int d2=0; d2<f.dim_d2; d2++){
        for(int d3=0; d3<f.dim_d3; d3++){
          for(int d4=0; d4<f.dim_d4; d4++){
            M(d3, d2*f.dim_d4+d4, d1)=f(d1, d2, d3, d4);
    } } } } */
    for(int d1=0; d1<f.dim_d1; d1++){
      for(int d3=0; d3<f.dim_d3; d3++){
        Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 0, Eigen::InnerStride<Eigen::Dynamic> >(M.mem+d3*M.dim_d2+d1, f.dim_d2, f.dim_d4, Eigen::InnerStride<Eigen::Dynamic>(M.dim_i*M.dim_d2))        = 
        Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 0, Eigen::OuterStride<Eigen::Dynamic> >(f.mem+(d1*f.dim_d2*f.dim_d3+d3)*f.dim_d4, f.dim_d2, f.dim_d4, Eigen::OuterStride<Eigen::Dynamic>(f.dim_d3*f.dim_d4));
    } }
    return M;
}

}//namespace
