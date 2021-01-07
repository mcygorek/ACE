#ifndef MATH_READER_DEFINED_H
#define MATH_READER_DEFINED_H

#include "Reader.h"
#include <Eigen/Core>
#include <fstream>
#include <complex>
#include <cmath>

/** Parse mathematical expressions such as
  1 + 3 * 7

1.: Identify types
  1 -> number
  + -> operator ADD
  3 -> number
  * -> operator MULT
  7 -> number


2.: Evaluate in correct order

*/

namespace MathReader{
  enum MathOpType {NONE, NUMBER, UNARY, BINARY_MID, PARENTHESIS};

  class MathOp{
  public:
    MathOpType type;
    std::string name;

    int get_nr_args()const{
      switch(type){
        case NONE:
        case NUMBER:
        case PARENTHESIS:
          return 0;
        case UNARY:
          return 1;
        case BINARY_MID:
          return 2;
      }
      return 0;
    }
    void check_nr_args(const std::vector<Eigen::MatrixXcd> &args)const{
      if(args.size()!=get_nr_args()){
        std::cerr<<"MathReader: wrong number of arguments!"<<std::endl; 
        exit(1);
      }
    }
    
    virtual Eigen::MatrixXcd evaluate(const std::vector<Eigen::MatrixXcd> &args)const{
      return Eigen::MatrixXcd();
    }
    MathOp(const std::string &str="none"){
      type=NONE;
      name=str;
    }
  };
  class MathOp_Number: public MathOp{
  public:
    Eigen::MatrixXcd value;
     
    virtual Eigen::MatrixXcd evaluate(const std::vector<Eigen::MatrixXcd> &args)const{
      return value;
    }
    void setup(const Eigen::MatrixXcd &mat=Eigen::MatrixXcd()){
      value=mat;
//      std::stringstream ss; ss<<c.real();
//      if(fabs(c.imag())>1e-20)ss<<"+i_"<<c.imag();
      std::stringstream ss; ss<<mat;
      name=ss.str();
      type=NUMBER;
    }
    MathOp_Number(std::complex<double> c=0.){
      value=Eigen::MatrixXcd(1,1); value(0,0)=c;
      setup(value);
    }
    MathOp_Number(const Eigen::MatrixXcd &mat){
      setup(mat);
    }

  };
  class MathOp_Plus: public MathOp{
  public:
    virtual Eigen::MatrixXcd evaluate(const std::vector<Eigen::MatrixXcd> &args)const{
      if(args.size()<2){std::cerr<<"MathOp + needs 2 arguments!"<<std::endl;exit(1);}
      return args[0]+args[1];
    }
    MathOp_Plus() : MathOp("+") {type=BINARY_MID;}
  };
  class MathOp_Minus: public MathOp{
  public:
    virtual Eigen::MatrixXcd evaluate(const std::vector<Eigen::MatrixXcd> &args)const{
      if(args.size()<2){std::cerr<<"MathOp - needs 2 arguments!"<<std::endl;exit(1);}
      return args[0]-args[1];
    }
    MathOp_Minus() : MathOp("-") {type=BINARY_MID;}
  };
  class MathOp_Times: public MathOp{
  public:
    virtual Eigen::MatrixXcd evaluate(const std::vector<Eigen::MatrixXcd> &args)const{
      if(args.size()<2){std::cerr<<"MathOp * needs 2 arguments!"<<std::endl;exit(1);}
      if(args[0].rows()==1 && args[0].cols()==1){
        return args[0](0,0)*args[1];
      }
      if(args[1].rows()==1 && args[1].cols()==1){
        return args[0]*args[1](0,0);
      }
      return args[0]*args[1];
    }
    MathOp_Times() : MathOp("*") {type=BINARY_MID;}
  };
  class MathOp_Divide: public MathOp{
  public:
    virtual Eigen::MatrixXcd evaluate(const std::vector<Eigen::MatrixXcd> &args)const{
      if(args.size()<2){std::cerr<<"MathOp * needs 2 arguments!"<<std::endl;exit(1);}
      if(args[1].rows()==1 && args[1].cols()==1){
        return args[0]/args[1](0,0);
      }
      std::cerr<<"Error parsing divison: second argument not scalar!"<<std::endl;
      exit(1);
    }
    MathOp_Divide() : MathOp("/") {type=BINARY_MID;}
  };

  class MathOp_Parenthesis: public MathOp{
  public:
    MathOp_Parenthesis(const std::string &nm){
      type=PARENTHESIS;
      name=nm;
      std::cerr<<"Parentheses not implemented yes!"<<std::endl;
      exit(1);
    }
  };



  class IdentifiedOp{
  public:
    std::vector<MathOp*> list;
   
    ~IdentifiedOp(){
      for(size_t i=0; i<list.size(); i++)delete list[i];
    }
  };

  std::ostream &operator<<(std::ostream &os, const IdentifiedOp &id){
    if(id.list.size()<1)return os;
    os<<id.list[0]->name;
    for(size_t i=1; i<id.list.size(); i++)os<<" "<<id.list[i]->name;
    return os;
  }

  //split up some tokens, e.g. 1*2 -> 1 * 2
  void preprocess(std::vector<std::string> &toks){
    for(size_t i=0; i<toks.size(); i++){
      if(toks[i].length()>1)for(size_t j=0; j<toks[i].length(); j++){
        if(toks[i][j]=='+'||toks[i][j]=='*'||toks[i][j]=='/'||
           toks[i][j]=='('||toks[i][j]==')'){
          std::string s1=""; if(j>0)s1=toks[i].substr(0,j);
          std::string s2(1,toks[i][j]);
          std::string s3=""; if(j<toks[i].length()-1)s3=toks[i].substr(j+1);
          toks[i]=s2;
//std::cout<<"j: "<<j<<" s1: '"<<s1<<"' s2: '"<<s2<<"' s3: '"<<s3<<"'"<<std::endl;
          if(s3!="")toks.insert(toks.begin()+i+1,s3);
          if(s1!="")toks.insert(toks.begin()+i,s1);
        }
      }
    }
  }

  //attribute meaning to the individual tokens
  IdentifiedOp identify(std::vector<std::string> toks){
    preprocess(toks);
    IdentifiedOp results;
    for(size_t i=0; i<toks.size(); i++){
      if(Reader::isDouble(toks[i])){
        results.list.push_back(new MathOp_Number(Reader::readDouble(toks[i])));
      }else if(toks[i]=="+"){
        results.list.push_back(new MathOp_Plus());
      }else if(toks[i]=="-"){
        results.list.push_back(new MathOp_Minus());
      }else if(toks[i]=="*"){
        results.list.push_back(new MathOp_Times());
      }else if(toks[i]=="/"){
        results.list.push_back(new MathOp_Divide());
      }else if(toks[i]=="("){
        results.list.push_back(new MathOp_Parenthesis("("));
      }else if(toks[i]==")"){
        results.list.push_back(new MathOp_Parenthesis(")"));
      }else{
        std::cerr<<"Cannot parse '"<<toks[i]<<"'!"<<std::endl;
        exit(1);
//        results.list.push_back(new MathOp(toks[i]));
      }
    }

    return results;
  }
  
  void partial_eval_binary_mid(IdentifiedOp & id, int i){
          if(id.list[i-1]->type!=NUMBER){
            std::cerr<<"id.list[i-1]->type!=NUMBER"<<std::endl;
            exit(1);
          }
          if(id.list[i+1]->type!=NUMBER){
            std::cerr<<"id.list[i+1]->type!=NUMBER"<<std::endl;
            exit(1);
          }
        
          std::vector<Eigen::MatrixXcd> MV(2);
          MV[0]=id.list[i-1]->evaluate(std::vector<Eigen::MatrixXcd>());
          MV[1]=id.list[i+1]->evaluate(std::vector<Eigen::MatrixXcd>());
          MathOp *mop=new MathOp_Number(id.list[i]->evaluate(MV));
 
          delete id.list[i-1];
          delete id.list[i];
          delete id.list[i+1];
          id.list.erase(id.list.begin()+i-1, id.list.begin()+i+1);
          id.list[i-1]=mop;
  }
 
  void check_parentheses(const std::vector<std::string> &toks){
    int nr_open=0;
    for(size_t i=0; i<toks.size(); i++){
      for(size_t j=0; j<toks[i].length(); j++){
        if(toks[i][j]=='(')nr_open++;
        if(toks[i][j]==')')nr_open--;
        if(nr_open<0){
          std::cerr<<"Check parentheses: ')' found without '('!"<<std::endl;
          exit(1);
        }
      }
    }
    if(nr_open!=0){
      std::cerr<<"Check parentheses: "<<nr_open<<" parentheses not closed'!"<<std::endl;
      exit(1);
    }
  }
  Eigen::MatrixXcd evaluate(const std::vector<std::string> &toks){
    check_parentheses(toks);
    IdentifiedOp id=identify(toks);
    while(true){
      std::cout<<id<<std::endl;
      if(id.list.size()<1){
        std::cerr<<"MathReader::evaluate: Empty!"<<std::endl;
        exit(1);
      }else if(id.list.size()==1){
        if(id.list[0]->get_nr_args()==0){
          return id.list[0]->evaluate(std::vector<Eigen::MatrixXcd>());
        }else{
           std::cerr<<"id.list.get_nr_args()!=0!"<<std::endl;
           exit(1);
        }
      }
      if(id.list[0]->type==BINARY_MID||id.list.back()->type==BINARY_MID){
        std::cerr<<"Error parsing: operator "<<id.list.back()->name<<" at the end of expression!"<<std::endl;
        exit(1);
      }

      bool wcont=false;
      for(int i=1; i<(int)id.list.size()-1; i++){
        if(id.list[i]->name=="*"||id.list[i]->name=="/"){
         partial_eval_binary_mid(id, i);
         wcont=true;
        }
      }
      if(wcont)continue;

      for(int i=1; i<(int)id.list.size()-1; i++){
        if(id.list[i]->name=="+"||id.list[i]->name=="-"){
         partial_eval_binary_mid(id, i);
         wcont=true;
        }
      }
      if(wcont)continue;

/*
      for(int i=0; i<(int)id.list.size(); i++){
        if(id.list[i]->name=="("){
          int nr_to_replace=-1;
          std::vector<string> inside;
          int nr_open=1;
          for(int j=i+1; (int)id.list.size(); j++){
            if(id.list[i]->name==")"){
              nr_open--;
              if(nr_open<1){
                nr_to_replace=j-i+1;
                break;
              }
            }else if(id.list[i]->name=="("){
              nr_open++;
            }
            inside.push_back(tok
          }
        }
      }
*/      

std::cout<<"MathReader: mismatch in number of operators!"<<std::endl; exit(1); 
      return Eigen::MatrixXcd();
    }
  }

};

#endif
