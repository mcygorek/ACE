#ifndef ABSTRACT_PROPAGATOR_DEFINED_H
#define ABSTRACT_PROPAGATOR_DEFINED_H

#include <Eigen/Dense>

namespace ACE{

class Propagator{
public:
  /// Propagator in Liouville space. To be updated in each time step:
  Eigen::MatrixXcd M;

  virtual int get_dim()const=0;

  virtual bool is_time_independent() const=0;

  virtual void update(double t, double dt)=0;


  inline virtual void check_dimensions()const{}

  virtual Eigen::MatrixXcd propagate_noupdate(const Eigen::MatrixXcd &rho);
  
  virtual Eigen::MatrixXcd propagate(const Eigen::MatrixXcd &rho, double t, double dt);
  
  virtual Eigen::MatrixXcd get_Htot(double t)const;

  inline virtual ~Propagator(){};
};

}//namespace
#endif
