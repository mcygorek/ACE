#pragma once
#ifndef PARAMETERS_DEFINED_H
#define PARAMETERS_DEFINED_H

#include <map>
#include <set>
#include <cmath>
#include <iosfwd>
#include <vector>

//#include "Eigen_fwd.hpp"
#include <Eigen/Core>
//#include "Reader.hpp"
#include "Printable.hpp"

namespace ACE{

typedef std::vector<std::vector<std::string> > Parameters_Entry;

class Parameters: public Printable{
public:
  std::map<std::string, Parameters_Entry> map;
  std::set<std::string> requested;
  bool register_requested;

  typedef std::map<std::string, Parameters_Entry>::iterator Iterator;
  typedef std::map<std::string, Parameters_Entry>::const_iterator cIterator;
  
  const Parameters_Entry dummy; //<-Empty entry for function "get"

  inline bool is_empty()const{return map.size()<1;}

  //append string to Parameters_Entry of key; creates key if it doesn't exist
  void add_to(const std::string & key, const std::vector<std::string> &arg);
    
  inline void add_to(const std::string & key, const std::string &arg){
    add_to(key, std::vector<std::string>(1, arg));
  }  
  inline void add_to(const std::string & key, const double &d){
    std::stringstream ss; ss<<d; 
    add_to(key, std::vector<std::string>(1, ss.str()));
  }
  inline void erase(const std::string & key){
    map.erase(key);
  } 
 
  void override_param(const std::string & key, const std::vector<std::string> &arg);
  void override_param(const std::string & key, const std::string &arg);
  void override_param(const std::string & key, double d);

  void add_if_not_specified(const std::string & key, const std::vector<std::string> &arg);
    
  inline void add_if_not_specified(const std::string & key, const std::string &arg){
    add_if_not_specified(key, std::vector<std::string>(1, arg));
  }
 
  inline void set_requested(const std::string &key);
  
  const Parameters_Entry &get(const std::string &key);
  

  inline bool is_specified(const std::string &key)const{
    return map.find(key)!=map.end();
  }
  void complain_if_not_specified(const std::string &key)const;
  void complain_if_row_shorter(const std::string &key, int n, int row=0, const std::string & context="")const;
//  void complain_if_rows_shorter(const std::string &key, int n)const;
  
  void complain_if_conflict(const std::string &key1, const std::string &key2)const;

  int get_nr_rows(const std::string &key)const;
  
  std::vector<std::string> get_row(const std::string &key, int row);
  
  int get_nr_cols(const std::string &key, int row)const;
 
  std::vector<std::string> get_all_strings(const std::string &key);
   
  inline std::vector<std::string> get_all_strings_check(const std::string &key){
    complain_if_not_specified(key);
    return get_all_strings(key); 
  }
  
  std::vector<double> get_all_double(const std::string &key);
  
  std::vector<size_t> get_all_size_t(const std::string &key);

  std::string get_as_single_string(const std::string &key, int row=0);

  std::vector<double> get_row_doubles(const std::string &key, int row=0, int min=0);

  std::string get_as_string(const std::string &key, const std::string &def="", int row=0, int col=0);

  std::string get_as_string_check(const std::string &key, int row=0, int col=0);

  Eigen::MatrixXcd get_as_operator(const std::string &key, Eigen::MatrixXcd def=Eigen::MatrixXcd::Zero(1,1), int row=0, int col=0);

  inline Eigen::MatrixXcd get_as_operator_check(const std::string &key, int row=0, int col=0){
    complain_if_not_specified(key);
    return get_as_operator(key, Eigen::MatrixXcd::Zero(1,1), row, col);
  }


  double get_as_double(const std::string &key, double def=0., int row=0, int col=0);
  inline double get_as_double_check(const std::string &key, int row=0, int col=0){
    complain_if_not_specified(key);
    return get_as_double(key, 0, row, col);
  }


  inline int get_as_int(const std::string &key, int def=0., int row=0, int col=0){
    return round(get_as_double(key, def, row, col));
  }
  inline int get_as_int_check(const std::string &key, int row=0, int col=0){
    complain_if_not_specified(key);
    return get_as_int(key, 0, row, col);
  }

  inline size_t get_as_size_t(const std::string &key, int def=0, int row=0, int col=0){
    int i=get_as_int(key, def, row, col);
    if(i<0){ std::cerr<<"Parameter '"<<key<<"' must not be smaller than zero!"<<std::endl; exit(1); }
    return (size_t)i;
  }
  inline size_t get_as_size_t_check(const std::string &key, int row=0, int col=0){
    complain_if_not_specified(key);
    return get_as_size_t(key, 0, row, col);
  }


  bool get_as_bool(const std::string &key, bool def=false, int row=0, int col=0);
  inline bool get_as_bool_check(const std::string &key, int row=0, int col=0){
    complain_if_not_specified(key);
    return get_as_bool(key, false, row, col);
  }


  void add_from_stringvec(const std::string &key, const std::vector<std::string> &toks);

  void add_from_line(const std::string &line);

  void add_from_file(const std::string &fname);
  
  void setup(int args, char** argv, bool read_first_as_driver=true);
  
  std::vector<std::string> get_lines()const;

  virtual std::ostream & print(std::ostream & ofs=std::cout)const;
  
  void print(const std::string &fname)const;
  
  void add_from_prefix(const std::string prefix, Parameters &other);
 
  inline void clear(){
    map.clear();
    requested.clear();
  } 
  Parameters(const std::string &fname){
    register_requested=false;
    add_from_file(fname);
  }
  Parameters(const Parameters &other){
    map=other.map;
  }
  Parameters(int args, char** argv, bool reg_req=false, bool read_first_as_driver=true){
    setup(args, argv, read_first_as_driver);
    register_requested=false;
  }
  Parameters(){ 
    register_requested=false;
  }
  
};

}//namespace
#endif
