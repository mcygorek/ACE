#ifndef ABSTRACT_PROPAGATOR_DEFINED_H
#define ABSTRACT_PROPAGATOR_DEFINED_H

class Propagator{
public:
  /// Propagator in Liouville space. To be updated in each time step:
  Eigen::MatrixXcd M;


  virtual int get_dim()const=0;

  virtual bool is_time_independent() const=0;

  virtual void update(double t, double dt)=0;


  virtual void check_dimensions()const{}

  virtual Eigen::MatrixXcd propagate_noupdate(const Eigen::MatrixXcd &rho){
    int dim=rho.rows();
    Eigen::MatrixXcd tmp=Eigen::MatrixXcd::Zero(dim,dim);
    for(int i=0; i<dim; i++){
      for(int j=0; j<dim; j++){
        for(int k=0; k<dim; k++){
          for(int l=0; l<dim; l++){
            tmp(i,j)+=M(i*dim+j, k*dim+l)*rho(k,l);
          }
        }
      }
    }
    return tmp;
  }
  virtual Eigen::MatrixXcd propagate(const Eigen::MatrixXcd &rho, double t, double dt){
    update(t, dt);
    return propagate_noupdate(rho);
  }
  virtual Eigen::MatrixXcd get_Htot(double t)const{
    return Eigen::MatrixXcd::Zero(get_dim(), get_dim());
  }

  virtual ~Propagator(){};
};

#endif
