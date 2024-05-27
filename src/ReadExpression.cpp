#include "PCH.hpp"
#include "ReadExpression.hpp"
#include "ReaderBasics.hpp"
//#include <complex>
//#include <Eigen/Dense>
#include "otimes.hpp"
#include "Operators_Boson.hpp"
#include <iostream>
//#include <vector>


namespace ACE{

 
void ExpressionOperand::set_scalar(std::complex<double> c){
    value=Eigen::MatrixXcd::Zero(1,1);
    value(0,0)=c;
}
void ExpressionOperand::complain_if_not_scalar()const{
    if(!is_scalar()){
      std::cerr<<"Expression: cannot convert Expression of dimension ("<<value.rows()<<", "<<value.cols()<<") to scalar!"<<std::endl;
      exit(1);
    }
  }

int ExpressionOperand::needed_operands()const{
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

std::string ExpressionOperand::get_name()const{
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

bool ExpressionOperand::is_complete()const{
    if(needed_operands()==0) return true;
    if( (int)operands.size()!=needed_operands() ) return false;
    for(size_t i=0; i<operands.size(); i++){
      if(!operands[i].is_complete())return false;
    }
    return true;
}

Eigen::MatrixXcd ExpressionOperand::eval(){
 
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
      value=otimes(operands[0].value, operands[1].value);
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
void ExpressionOperand::print_tree(std::ostream &os)const{
  os<<get_name();
  if(type==Value){
    os<<"_"<<value.rows();
  }
  if(operands.size()<1){
    os<<" ";
    return;
  }
  os<<"(";
  for(size_t i=0; i<operands.size(); i++)operands[i].print_tree(os);
  os<<") ";
}
std::complex<double> ExpressionOperand::eval_scalar(){
    value=eval();
    if(!is_scalar()){
      std::cerr<<"Cannot cast Matrix to scalar!"<<std::endl;
      exit(1);
    }
    return value(0,0);
}

void ExpressionOperand::add_last_operand(const ExpressionOperand &op){
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



void ReadExpression::complain_if_not_scalar()const{
    if(!is_scalar()){
      std::cerr<<"ReadExpression: cannot convert Expression '"<<number<<"' to scalar!"<<std::endl;
      exit(1);
    }
}


int ReadExpression::predefined_operator_dim(std::stringstream &ss, const std::string &str){
    int dim=0;
    std::string with_us=std::string(str+"_");
    if(is_next_in_stream(ss, with_us)){
      ss>>dim;
      if(ss.fail() || dim<1){
        std::cerr<<"No valid dimension number after underscore at '"<<str<<"_'!"<<std::endl;   
        exit(1);
      }
    }
    return dim;
}
 
Eigen::MatrixXcd ReadExpression::read(const std::string &str, bool force_print){

#ifdef DEBUG_EXPRESSIONS
    std::cout<<"ReadExpression called with argument '"<<str<<"'"<<std::endl;
#endif
    number=Eigen::MatrixXcd::Zero(1,1);

    double d=0.;
    std::stringstream ss(str);
 
    ExpressionOperand op(ExpressionOperand::Copy);

    while(!ss.eof()){
if(force_print){op.print_tree(); std::cout<<std::endl;}
      ss>> std::ws;  // eat up any leading white spaces
      
      int peek=ss.peek(); //peek character

      if(peek=='#')break; //end line if comment found
      
      if(peek=='('){
        std::string stmp=find_matching_brace_block_complain(ss);
#ifdef DEBUG_EXPRESSIONS
        std::cout<<"Term within parentheses: '"<<stmp<<"'"<<std::endl;
#endif
        op.add_last_operand((Eigen::MatrixXcd)ReadExpression(stmp,force_print));
        continue;
      }     
      if(peek=='{'){
        std::string stmp=find_matching_brace_block_complain(ss,'{','}');
#ifdef DEBUG_EXPRESSIONS
        std::cout<<"Term within curly braces: '"<<stmp<<"'"<<std::endl;
#endif
        op.add_last_operand((Eigen::MatrixXcd)ReadExpression(stmp,force_print));
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
          if(peek=='-'){ //deal with cases such as "{-hbar}"
            op.add_operand(ExpressionOperand::Multiply);
            op.add_last_operand(-1.);
            ss.seekg(((int)ss.tellg())+1,std::ios_base::beg);
            continue;
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
      if(is_next_in_stream(ss, "otimes")){
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
      int id_dim=predefined_operator_dim(ss,"Id");
      if(id_dim){
        op.add_last_operand(Eigen::MatrixXcd::Identity(id_dim,id_dim));
        continue;
      }

      id_dim=predefined_operator_dim(ss,"bdagger");
      if(id_dim){
        op.add_last_operand(Operators_Boson::adagger(id_dim));
        continue;
      }

      id_dim=predefined_operator_dim(ss,"b");
      if(id_dim){
        op.add_last_operand(Operators_Boson::a(id_dim));
        continue;
      }

      id_dim=predefined_operator_dim(ss,"n");
      if(id_dim){
        op.add_last_operand(Operators_Boson::n(id_dim));
        continue;
      }

      if(is_next_in_stream(ss, "sigma_x")){
        op.add_last_operand(sigma_x());
        continue;
      }
      if(is_next_in_stream(ss, "sigma_y")){
        op.add_last_operand(sigma_y());
        continue;
      }
      if(is_next_in_stream(ss, "sigma_z")){
        op.add_last_operand(sigma_z());
        continue;
      }
      if(is_next_in_stream(ss, "sigma_plus")){
        op.add_last_operand(sigma_plus());
        continue;
      }
      if(is_next_in_stream(ss, "sigma_minus")){
        op.add_last_operand(sigma_minus());
        continue;
      }

      if(peek=='|'){
        const std::string failstr="Cannot read matrix basis element |i><j|_n!";
        char c; ss>>c;
        int i; ss>>i;
        if(ss.fail() || i<0){std::cerr<<failstr<<std::endl; exit(1);}
        ss>>std::ws>>c;
        if(ss.fail() || c!='>'){std::cerr<<failstr<<std::endl; exit(1);}
        ss>>std::ws>>c;
        if(ss.fail() || c!='<'){std::cerr<<failstr<<std::endl; exit(1);}
        int j; ss>>j;
        if(ss.fail() || j<0){std::cerr<<failstr<<std::endl; exit(1);}
        ss>>std::ws>>c;
        if(ss.fail() || c!='|'){std::cerr<<failstr<<std::endl; exit(1);}
        ss>>std::ws>>c;
        if(ss.fail() || c!='_'){std::cerr<<failstr<<std::endl; exit(1);}
        int dim; ss>>dim;
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

      //try to read a number 
      std::stringstream ss2(ss.str());
      ss2.seekg(ss.tellg());
      ss2>>d;   
      if(!ss2.fail()){   // a number has been found
        ss>>d;
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
      ss2.clear();
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
      }else if(stmp.length()>=2 && stmp.substr(0,2)=="wn"){
        op.add_last_operand(ExpressionOperand(0.188364997625644));
        ss.seekg(spos+2,std::ios_base::beg);
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


}//namespace
