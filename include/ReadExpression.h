#ifndef READEXPRESSION_DEFINED_H
#define READEXPRESSION_DEFINED_H
#include "Reader.h"
#include <complex>
#include <Eigen/Dense>
#include "OuterProduct.h"

/* Class to turn a string into a valid complex number */


class ExpressionOperand{
public:
  enum TYPE { None, Copy, Value, Add, Subtract, Multiply, Divide, Otimes,
              Sqrt, Exp, Sin, Cos, Tan} type;
 
  Eigen::MatrixXcd value;
  std::vector<ExpressionOperand> operands;

  void set_scalar(std::complex<double> c){
    value=Eigen::MatrixXcd::Zero(1,1);
    value(0,0)=c;
  }
  bool is_scalar()const{ 
    return ((value.rows()==1) && (value.cols()==1));
  }
  void complain_if_not_scalar()const{
    if(!is_scalar()){
      std::cerr<<"Expression: cannot convert Expression of dimension ("<<value.rows()<<", "<<value.cols()<<") to scalar!"<<std::endl;
      exit(1);
    }
  }
  std::complex<double> get_scalar()const{
    complain_if_not_scalar();
    return value(0,0);
  }

  int needed_operands()const{
    switch(type){
      case Value: return 0; //access stored variable
      case Copy: return 1; 
      case Add: return 2; 
      case Subtract: return 2; 
      case Multiply: return 2; 
      case Divide: return 2; 
      case Otimes: return 2; 
      case Sqrt: return 1; 
      case Exp: return 1; 
      case Sin: return 1; 
      case Cos: return 1; 
      case Tan: return 1; 
      default: return 0; //E.g.: None
    }
  }
  std::string get_name()const{
    switch(type){
      case Value: return "Value";
      case Copy: return "Copy";
      case Add: return "Add";
      case Subtract: return "Subtract";
      case Multiply: return "Multiply";
      case Divide: return "Divide";
      case Otimes: return "Otimes";
      case Sqrt: return "Sqrt";
      case Exp: return "Exp";
      case Sin: return "Sin";
      case Cos: return "Cos";
      case Tan: return "Tan";
      default: return "None";
    }
  }

  bool is_complete()const{
    if(needed_operands()==0) return true;
    if( (int)operands.size()!=needed_operands() ) return false;
    for(size_t i=0; i<operands.size(); i++){
      if(!operands[i].is_complete())return false;
    }
    return true;
  }

  Eigen::MatrixXcd eval(){
 
    if(type==None){ std::cerr<<"Trying to evaluation expression 'None'!"<<std::endl; exit(1);}

    if(type==Value)return value;
 
    if(!is_complete()){
      std::cerr<<"Cannot evaluate expression '"<<get_name()<<"': Not enough operands!"<<std::endl;
      exit(1);
    }
    for(size_t i=0; i<operands.size(); i++){
      operands[i].eval(); //All required operands are now of type "Value"
    }
   
    if(type==Copy){
      value=operands[0].value;
      type=Value;
      return value;
    }else if(type==Add){
      value=operands[0].value+operands[1].value;
      type=Value;
      return value;
    }else if(type==Subtract){
      value=operands[0].value-operands[1].value;
      type=Value;
      return value;
    }else if(type==Multiply){
      type=Value;
      if(operands[0].is_scalar() && operands[1].is_scalar()){
        set_scalar(operands[0].get_scalar()*operands[1].get_scalar());
      }else if(operands[0].is_scalar()){
        value=operands[0].get_scalar()*operands[1].value;
      }else if(operands[1].is_scalar()){
        value=operands[0].value*operands[1].get_scalar();
      }else{
        value=operands[0].value*operands[1].value;
      }
      return value;
    }else if(type==Divide){
      std::complex<double> c=operands[1].get_scalar();
      value=operands[0].value*(1./c);
      type=Value;
      return value;
    }else if(type==Otimes){
      value=OuterProduct(operands[0].value, operands[1].value);
      type=Value;
      return value;
    }else if(type==Sqrt){
      set_scalar(std::sqrt(operands[0].get_scalar()));
      type=Value;
      return value;
    }else if(type==Exp){
      set_scalar(std::exp(operands[0].get_scalar()));
      type=Value;
      return value;
    }else if(type==Sin){
      set_scalar(std::sin(operands[0].get_scalar()));
      type=Value;
      return value;
    }else if(type==Cos){
      set_scalar(std::cos(operands[0].get_scalar()));
      type=Value;
      return value;
    }else if(type==Tan){
      set_scalar(std::tan(operands[0].get_scalar()));
      type=Value;
      return value;

    }else{
      std::cerr<<"eval() not implemented for expression '"<<get_name()<<"'!"<<std::endl;
      exit(1);
    }
  }
  std::complex<double> eval_scalar(){
    value=eval();
    if(!is_scalar()){
      std::cerr<<"Cannot cast Matrix to scalar!"<<std::endl;
      exit(1);
    }
    return value(0,0);
  }

  void add_operand(const ExpressionOperand &op){
    operands.push_back(op);
  }
  void add_last_operand(const ExpressionOperand &op){
    for(size_t i=0; i<operands.size(); i++){
      if(!operands[i].is_complete()){
        operands[i].add_last_operand(op);
        return;
      }
    }
    if( (int)operands.size()==needed_operands() ){
      std::cerr<<"Too many operands in expression!"<<std::endl;
      exit(1);
    }
    add_operand(op);
  }
  void add_last_operand(const Eigen::MatrixXcd & c){
    add_last_operand(ExpressionOperand(c));
  }
  void add_last_operand(std::complex<double> c){
    add_last_operand(ExpressionOperand(c));
  }


  ExpressionOperand() : type(None), value(Eigen::MatrixXcd::Zero(1,1)){}
  ExpressionOperand(TYPE tp) : type(tp), value(Eigen::MatrixXcd::Zero(1,1)){}
  ExpressionOperand(std::complex<double> c) : type(Value), value(Eigen::MatrixXcd::Zero(1,1)){
    set_scalar(c);
  }
  ExpressionOperand(const Eigen::MatrixXcd &M) : type(Value){
    value=M;
  }

};

class ReadExpression {
public:
  Eigen::MatrixXcd number;

  bool is_scalar()const{ 
    return number.rows()==1 && number.cols()==1;
  }
  void complain_if_not_scalar()const{
    if(!is_scalar()){
      std::cerr<<"ReadExpression: cannot convert Expression '"<<number<<"' to scalar!"<<std::endl;
      exit(1);
    }
  }
 
  Eigen::MatrixXcd read(const std::string &str){

#ifdef DEBUG_EXPRESSIONS
    std::cout<<"ReadExpression called with argument '"<<str<<"'"<<std::endl;
#endif
    number=Eigen::MatrixXcd::Zero(1,1);

    double d=0.;
    std::stringstream ss(str);
 
    ExpressionOperand op(ExpressionOperand::Copy);

    while(!ss.eof()){
      ss>> std::ws;  // eat up any leading white spaces
      
      int peek=ss.peek(); //peek character

      if(peek=='#')break; //end line if comment found
      
      if(peek=='('){
        std::string stmp=Reader::find_matching_brace_block_complain(ss);
#ifdef DEBUG_EXPRESSIONS
        std::cout<<"Term within parentheses: '"<<stmp<<"'"<<std::endl;
#endif
        op.add_last_operand((Eigen::MatrixXcd)ReadExpression(stmp));
        continue;
      }     
      if(peek=='{'){
        std::string stmp=Reader::find_matching_brace_block_complain(ss,'{','}');
#ifdef DEBUG_EXPRESSIONS
        std::cout<<"Term within curly braces: '"<<stmp<<"'"<<std::endl;
#endif
        op.add_last_operand((Eigen::MatrixXcd)ReadExpression(stmp));
        continue;
      }     

      if(peek=='+' || peek=='-' || peek=='*' || peek=='/'){
#ifdef DEBUG_EXPRESSIONS
        std::cout<<"Found operator: "<<(char)peek<<std::endl;
#endif
        //If op.is_complete(), then '+' or '-' are operators. If not, they are part of a number:
        if(!op.is_complete()){
          if(peek=='*'||peek=='/'){ 
            std::cerr<<"'*' or '/' after incomplete expression!"<<std::endl;
            exit(1);
          }
          //else: go further in the loop and try to read a number.
        }else{
          if(peek=='+'){
            ExpressionOperand optmp=op; 
            op=ExpressionOperand(ExpressionOperand::Add);
            op.add_operand(optmp);
          }else if(peek=='-'){
            ExpressionOperand optmp=op; 
            op=ExpressionOperand(ExpressionOperand::Subtract);
            op.add_operand(optmp);
          }else{
            if(op.operands.size()<1){
              std::cerr<<"Error prioritizing expressions!"<<std::endl;
              exit(1);
            }
            ExpressionOperand op2(ExpressionOperand::Multiply);
            if(peek=='/')op2.type=ExpressionOperand::Divide;
            op2.add_operand(op.operands.back());
            op.operands.back()=op2;
          }
          ss.seekg(((int)ss.tellg())+1,std::ios_base::beg);
          continue;
        }
      }
      //process "otimes"
      if(Reader::is_next_in_stream(ss, "otimes")){
        if(!op.is_complete()){
          std::cerr<<"'otimes' after incomplete expression!"<<std::endl;
          exit(1);
        }else{
          if(op.operands.size()<1){
            std::cerr<<"Error prioritizing expressions!"<<std::endl;
            exit(1);
          }
          ExpressionOperand op2(ExpressionOperand::Otimes);
          op2.add_operand(op.operands.back());
          op.operands.back()=op2;
        }
        continue;
      }
  
      //try to identify basic matrices
      if(Reader::is_next_in_stream(ss, "Id_")){
        int dim=0;
        ss>>dim;
        if(ss.fail() || dim<1){
          std::cerr<<"No valid dimension number after underscore at 'Id_'!"<<std::endl;
          exit(1);
        }
        op.add_last_operand(Eigen::MatrixXcd::Identity(dim,dim));
        continue;
      }
      if(peek=='|'){
        const std::string failstr="Cannot read matrix basis element |i><j|_n!";
        char c; ss>>c;
//std::cout<<"TEST: "<<c<<std::endl;
        int i; ss>>i;
//std::cout<<"TEST: "<<i<<std::endl;
        if(ss.fail() || i<0){std::cerr<<failstr<<std::endl; exit(1);}
        ss>>std::ws>>c;
//std::cout<<"TEST: "<<c<<std::endl;
        if(ss.fail() || c!='>'){std::cerr<<failstr<<std::endl; exit(1);}
        ss>>std::ws>>c;
//std::cout<<"TEST: "<<c<<std::endl;
        if(ss.fail() || c!='<'){std::cerr<<failstr<<std::endl; exit(1);}
        int j; ss>>j;
//std::cout<<"TEST: "<<j<<std::endl;
        if(ss.fail() || j<0){std::cerr<<failstr<<std::endl; exit(1);}
        ss>>std::ws>>c;
//std::cout<<"TEST: "<<c<<std::endl;
        if(ss.fail() || c!='|'){std::cerr<<failstr<<std::endl; exit(1);}
        ss>>std::ws>>c;
//std::cout<<"TEST: "<<c<<std::endl;
        if(ss.fail() || c!='_'){std::cerr<<failstr<<std::endl; exit(1);}
        int dim; ss>>dim;
//std::cout<<"TEST: "<<dim<<std::endl;
        if(ss.fail() || dim<1){std::cerr<<failstr<<std::endl; exit(1);}

        if(i>=dim||j>=dim){
          std::cerr<<"Error reading matrix in expression: i>=dim||j>=dim!"<<std::endl; 
          exit(1);
        }

        Eigen::MatrixXcd mat=Eigen::MatrixXcd::Zero(dim,dim);
        mat(i,j)=1;
        op.add_last_operand(mat);
        continue;
      } 

      ss>>d;   //try to read a number 
      if(!ss.fail()){   // a number has been found
#ifdef DEBUG_EXPRESSIONS
        std::cout<<"Found operand: number: "<<d<<std::endl;
#endif
        //useful for internal purposes: interpret 2_i-3 as "2-3*i"
        double im=0;
        int peek=ss.peek();
        if(!ss.eof() && peek=='_'){
          int got=ss.get(); got=ss.get();
          if(!ss.eof() && (char)got=='i'){
            ss>>im; 
            if(ss.fail()){
              std::cerr<<"Can't read imaginary part after '_i'!"<<std::endl;
              exit(1);
            }
          }else{
            std::cerr<<"Can't interpret underscore without '_i'!"<<std::endl;
            exit(1);
          }
        }
        op.add_last_operand(std::complex<double>(d,im));
        
        continue;
      }
      ss.clear();
      //no number:

      int spos=ss.tellg(); //remember position in steam
      std::string stmp; ss>>stmp;
      if(stmp=="")break;

#ifdef DEBUG_EXPRESSIONS
      std::cout<<"stmp: '"<<stmp<<"'"<<std::endl;
#endif

      if(stmp[0]=='i'){
        if(stmp.length()<2||stmp[1]=='+'||stmp[1]=='-'||stmp[1]=='*'||stmp[1]=='/'){
#ifdef DEBUG_EXPRESSIONS
        std::cout<<"imaginary unit found: "<<stmp[0]<<std::endl;
#endif
          op.add_last_operand(ExpressionOperand(std::complex<double>(0.,1.)));
          ss.seekg(spos+1,std::ios_base::beg);
          continue;
        }
      }
      if(stmp.length()>=2 && stmp.substr(0,2)=="pi"){
        op.add_last_operand(ExpressionOperand(M_PI));
        ss.seekg(spos+2,std::ios_base::beg);
        continue;
      }else if(stmp.length()>=2 && stmp.substr(0,2)=="kB"){
        op.add_last_operand(ExpressionOperand(8.6173303e-2));
        ss.seekg(spos+2,std::ios_base::beg);
        continue;
      }else if(stmp.length()>=4 && stmp.substr(0,4)=="hbar"){
        op.add_last_operand(ExpressionOperand(0.6582119569));
        ss.seekg(spos+4,std::ios_base::beg);
        continue;
      }

      if(stmp.length()>=4 && stmp.substr(0,4)=="sqrt"){
        op.add_last_operand(ExpressionOperand(ExpressionOperand::Sqrt));
        ss.seekg(spos+4,std::ios_base::beg);
        continue;
      }else if(stmp.length()>=3 && stmp.substr(0,3)=="exp"){
        op.add_last_operand(ExpressionOperand(ExpressionOperand::Exp));
        ss.seekg(spos+3,std::ios_base::beg);
        continue;
      }else if(stmp.length()>=3 && stmp.substr(0,3)=="sin"){
        op.add_last_operand(ExpressionOperand(ExpressionOperand::Sin));
        ss.seekg(spos+3,std::ios_base::beg);
        continue;
      }else if(stmp.length()>=3 && stmp.substr(0,3)=="cos"){
        op.add_last_operand(ExpressionOperand(ExpressionOperand::Cos));
        ss.seekg(spos+3,std::ios_base::beg);
        continue;
      }else if(stmp.length()>=3 && stmp.substr(0,3)=="tan"){
        op.add_last_operand(ExpressionOperand(ExpressionOperand::Tan));
        ss.seekg(spos+3,std::ios_base::beg);
        continue;
      }

      std::cerr<<"'"<<stmp<<"': "<<"No Expression!"<<std::endl; exit(1);
    }

    number=op.eval();
    return number;
  }
  operator Eigen::MatrixXcd () const { 
    return number; 
  }
  operator std::complex<double> () const { 
    complain_if_not_scalar();
    return number(0,0);
  }
 
  ReadExpression(const std::string &str="") : number(Eigen::MatrixXcd::Zero(1,1)){
    read(str); 
  }

};

#endif
