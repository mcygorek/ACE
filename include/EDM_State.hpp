#ifndef ACE_EDM_STATE_DEFINED_H
#define ACE_EDM_STATE_DEFINED_H

#include "EDM_Index.hpp"
#include "EDM_Filter.hpp"
#include <Eigen/Core>
#include <memory>
#include <complex>
#include <unordered_map>
#include <map>
#include <fstream>
#include <iostream>

namespace ACE {

class EDM_State{
//Describes current state via CoeffList and site dimensions.
//As trace closures also change, these are stored here as well

public:
//  typedef std::map<EDM_Index, std::complex<double> >::iterator Iterator;
//  typedef std::map<EDM_Index, std::complex<double> >::const_iterator cIterator;
//  typedef std::unordered_map<EDM_Index, std::complex<double> >::iterator Iterator;
//  typedef std::unordered_map<EDM_Index, std::complex<double> >::const_iterator cIterator;


  EDM_Index Ldim; //List of dimensions for each site in the network

//  std::map<EDM_Index, std::complex<double> > coeffs; //coefficients of states with non-zero contribution
//  std::map<EDM_Index, std::complex<double> > candidates; //coefficients of states with non-zero contribution
  std::unordered_map<EDM_Index, std::complex<double> > coeffs; //coefficients of states with non-zero contribution
  std::unordered_map<EDM_Index, std::complex<double> > candidates; //coefficients of states with non-zero contribution


  std::vector<Eigen::VectorXcd> closures; //<-to trace out degrees of freedom

  std::vector<Eigen::VectorXd> relevance; //<-e.g. SVs for PT-MPO in forwardNF
 
  inline void clear(){
    Ldim.clear();
    coeffs.clear();
    candidates.clear();
    closures.clear();
    relevance.clear();
  }
  inline void swap(EDM_State & other){
    Ldim.swap(other.Ldim);
    coeffs.swap(other.coeffs);
    candidates.swap(other.candidates);
    closures.swap(other.closures);
    relevance.swap(other.relevance);
  }
  inline void copy_empty(const EDM_State & other){
    coeffs.clear();
    candidates.clear();
    Ldim=other.Ldim;
    closures=other.closures;
    relevance=other.relevance;
  }
  inline bool empty(){
    return coeffs.empty() && candidates.empty();
  }

  void add_to(const std::pair<EDM_Index, std::complex<double> > & P);
  inline void add_to(const EDM_Index & I, const std::complex<double> &val){
    add_to(std::make_pair(I,val));
  }
  

  inline void add_filtered(const EDM_Index & I, const std::complex<double> &val, EDM_Filter &filter){
    if(filter.increment_passes(val)){
      add_to(I,val);
    }
  }

  double get_relevance(const EDM_Index &I)const;
  void reduce(const EDM_Filter &filter);

  void print(std::ostream &os=std::cout)const;
  void print_closures(std::ostream &os=std::cout)const;

  void write(const std::string &fname)const;

  void set_single(const Eigen::VectorXcd & rho, double thr=0);
  void set_from_product(const std::vector<Eigen::VectorXcd> & rhos,double thr=0);
  void set_from_product(const Eigen::VectorXcd & rho0, const Eigen::VectorXcd & rho1, double thr=0);

  EDM_State trace_over(int r)const;
  Eigen::VectorXcd get_reduced(int r)const;

  void check_coeffs_consistent()const;
  void check_closures_consistent()const;

  EDM_State(){}
};
}
#endif
