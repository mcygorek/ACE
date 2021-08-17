#ifndef ACE_SPACE_DEFINED_H
#define ACE_SPACE_DEFINED_H

/** Structure to deal with product spaces:
On the one hand, the total system is a product of system and environment parts.
On the other hand, the Liouville space is a product of Hilbert spaces.
Then again, we could have systems (e.g., two-level system + cavity) that
are themselves product spaces.

This class tracks product spaces, mostly to perform consistency checks.
The reference will be Hilbert spaces. 
*/

class Space{
private:
  std::vector<int> dims;
  int totaldim;

public:
  
  int get_dim()const{ return totaldim; } //total dimension
  bool is_empty()const{ return dims.size()<1; }

  void update_total(){ //updates value of totaldim 
    totaldim=1;
    for(std::vector<int>::iterator it=dims.begin(); it!=dims.end; ++it){
      if(*it<1){
        std::cerr<<"Space::update_total: dimension<1!"<<std::endl;
        exit(1);
      }
      totaldim*=*it;
    }
  }

  Space otimes(const Space & sp2) const{
    Space res(dims);
    res.dims.insert(res.dims.end(), sp2.dims.begin(), sp2.dims.end());
    res.update_total();
    return res;
  }
  Space operator*(const Space & sp2) const{
    return otimes(sp2);
  }

  Space(){totaldim=1;}
  Space(int n){
    dims=std::vector<int>(1,n);
    update_total();
  }
  Space(const std::vector<int> &dims_) : dims(dims_) {
    update_total();
  }
  virtual ~Space(){}
};

bool are_compatible(const Space & sp1, const Space & sp2){
  return sp1.get_dim()==sp2.get_dim();
}


//Hilbert space matrix:
class HMatrix{
private:
  Eigen::MatrixXcd EMatrix;
  Space            sp;
  
public: 
//  Eigen::MatrixXcd & get_EMatrix(){ return EMatrix; }
  const Eigen::MatrixXcd & get_EMatrix() const{ return EMatrix; }

  int get_dim()const{ return sp.get_dim(); }


  std::complex<double> & operator()(int i, int j){
    if(i<0 || i>=get_dim() || j<0||j>=get_dim()){
      std::cerr<<"HMatrix: out of bounds: ("<<i<<", "<<j<<") vs. dim "<<get_dim()<<"!"<<std::endl;
      exit(1);
    }
    return Ematrix(i,j);
  }
  const std::complex<double> & operator()(int i, int j) const{
    if(i<0 || i>=get_dim() || j<0||j>=get_dim()){
      std::cerr<<"HMatrix: out of bounds: ("<<i<<", "<<j<<") vs. dim "<<get_dim()<<"!"<<std::endl;
      exit(1);
    }
    return Ematrix(i,j);
  }  


//TODO:   operator*, otimes, ... 
// same for LMatrix: additionally: expand MPSMatrix



  HMatrix(const Space & sp_=Space() ) : sp(sp_) {
    EMatrix=Eigen::MatrixXcd::Zero(sp.get_dim());
  }
  HMatrix(const Eigen::MatrixXcd & EMatrix_, const Space & sp_=Space() ) 
    : Ematrix(EMatrix_), sp(sp_) {
 
    if(EMatrix.rows()!=EMatrix.cols()){
      std::cerr<<"Constructing HMatrix: EMatrix.rows()!=EMatrix.cols()!"<<std::endl;
      exit(1);
    }
    if(sp.is_empty()){
      sp=Space(EMatrix.rows());
    }else{
      if(EMatrix.rows()!=sp.get_dim()){
        std::cerr<<"Constructing HMatrix: EMatrix.rows()!=sp.get_dim()!"<<std::endl; 
        exit(1);
      }
    }
  } 
  virtual ~HMatrix();
};




#endif 
