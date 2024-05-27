#ifndef ACE_RUNGEKUTTA_DEFINED_H
#define ACE_RUNGEKUTTA_DEFINED_H
#include <Eigen/Dense>

/**
Basic Runge Kutta solver. Requires a function that produces a Vector, which
is interpreted as the derivative. Changes initial value to result in place.
*/

namespace ACE{

template <typename VecType>
class EigenRungeKutta_Function{
public:
  virtual VecType f_RK(const VecType &in, double t) const = 0;
};

template <typename VecType>
class EigenRungeKutta{
public:
  typedef VecType (*f_func)(const VecType &in, double t);

  void do_step(VecType &in, f_func f, double t, double dt){
    VecType k1, k2, k3, k4;
    k1 = f(in, t);
    k2 = f(in + 0.5*dt*k1, t+0.5*dt);
    k3 = f(in + 0.5*dt*k2, t+0.5*dt);
    k4 = f(in + dt*k3, t+dt);
    in += (k1+2.*k2+2.*k3+k4)*dt/6.;
  }
  void do_step(VecType &in, EigenRungeKutta_Function<VecType> &f, double t, double dt){
    VecType k1, k2, k3, k4;
    k1 = f.f_RK(in, t);
    k2 = f.f_RK(in + 0.5*dt*k1, t+0.5*dt);
    k3 = f.f_RK(in + 0.5*dt*k2, t+0.5*dt);
    k4 = f.f_RK(in + dt*k3, t+dt);
    in += (k1+2.*k2+2.*k3+k4)*dt/6.;
  }



  EigenRungeKutta(VecType &in, f_func f, double t, double dt){
    do_step(in, f, t, dt);
  }
  EigenRungeKutta(VecType &in, EigenRungeKutta_Function<VecType> &f, double t, double dt){
    do_step(in, f, t, dt);
  }

};


}//namespace

#endif
