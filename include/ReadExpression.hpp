#pragma once
#ifndef READEXPRESSION_DEFINED_H
#define READEXPRESSION_DEFINED_H
#include <complex>
#include <Eigen/Dense>
#include <vector>
#include <iostream>

/* Class to turn a string into a valid complex number */

namespace ACE{

class ExpressionOperand{
public:
  enum TYPE { None, Copy, Value, Add, Subtract, Multiply, Divide, Otimes,
              Sqrt, Exp, Ln, Sin, Cos, Tan} type;
 
  Eigen::MatrixXcd value;
  std::vector<ExpressionOperand> operands;


  void print_tree(std::ostream &os=std::cout)const;

  void set_scalar(std::complex<double> c);
 
  inline bool is_scalar()const{ 
    return ((value.rows()==1) && (value.cols()==1));
  }

  void complain_if_not_scalar()const;
  
  inline std::complex<double> get_scalar()const{
    complain_if_not_scalar();
    return value(0,0);
  }

  int needed_operands()const;
  std::string get_name()const;
  bool is_complete()const;

  Eigen::MatrixXcd eval();
  
  std::complex<double> eval_scalar();

  inline void add_operand(const ExpressionOperand &op){
    operands.push_back(op);
  }

  void add_last_operand(const ExpressionOperand &op);
  
  inline void add_last_operand(const Eigen::MatrixXcd & c){
    add_last_operand(ExpressionOperand(c));
  }
  inline void add_last_operand(std::complex<double> c){
    add_last_operand(ExpressionOperand(c));
  }
  inline void add_last_operand(double c){
    add_last_operand(std::complex<double>(c,0.));
  }


  ExpressionOperand() : type(None), value(Eigen::MatrixXcd::Zero(1,1)){}
  ExpressionOperand(TYPE tp) : type(tp), value(Eigen::MatrixXcd::Zero(1,1)){}
  ExpressionOperand(std::complex<double> c) : type(Value), value(Eigen::MatrixXcd::Zero(1,1)){
    set_scalar(c);
  }
  ExpressionOperand(double c) : type(Value), value(Eigen::MatrixXcd::Zero(1,1)){
    set_scalar(c);
  }
  ExpressionOperand(const Eigen::MatrixXcd &M) : type(Value){
    value=M;
  }

};

class ReadExpression {
public:
  Eigen::MatrixXcd number;

  inline bool is_scalar()const{ 
    return number.rows()==1 && number.cols()==1;
  }
  void complain_if_not_scalar()const;

  static int predefined_operator_dim(std::stringstream &ss, const std::string &str);

  Eigen::MatrixXcd read(const std::string &str, bool force_print=false);
  
  inline operator Eigen::MatrixXcd () const { 
    return number; 
  }
  inline operator std::complex<double> () const { 
    complain_if_not_scalar();
    return number(0,0);
  }
 
  ReadExpression(const std::string &str="", bool force_print=false) : number(Eigen::MatrixXcd::Zero(1,1)){
    read(str, force_print); 
  }

};

}//namespace
#endif
