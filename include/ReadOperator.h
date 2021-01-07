#ifndef READOPERATOR_DEFINED_H
#define READOPERATOR_DEFINED_H
#include "ReadExpression.h"
#include <Eigen/Dense>

/**
Class to read an Operator from file:

We would like to be able to specify, e.g., an initial state or a Hamiltonian
from a human readable file. It is assumed that all matrix elements are zero
except the ones specified directly via expressions like |0><1|.
Furthermore, we want to introduce a few appreviations, like "Id" for the 
identity matrix. 

Consider the expression:

3:  0.5 * |0><0| + 0.5 * |1><1|  

should set up a 3x3 matrix with 0.5 in the first two diagonal elements.

The parsing will be pretty tough!
Steps: 
1. coarse tokening: we need to split up the expressions into different 
   sub-expression. First, the term before ":" has to be separated. Then,
   we decompose into segments belonging to the individual basis elements |i><j|

->  { 3 },  {{0.5 },{|0><0|}}, {{0.5},{|1><1}}

2. the prefactor before each basis element (in the above example: {0.5}) needs
   to be evaluated as an expression for a general complex number. This task
   can be passed on to a more general routine that knows what to do with
   more complex representations such as " 1/sqrt(2) + i/sqrt(2) "...    

3. Once elementary operators are implemented, we can include compositions of 
   elementary operators to shorten the notation.

*/

class ReadOperator{
public:
  Eigen::MatrixXcd M;
 

  std::complex<double> eval_prefac(const std::string &str){        
    std::complex<double> prefac=1.;
    std::string stmp=Reader::trim(str);
    if(stmp!=""){
      if(stmp[stmp.length()-1]=='*'){
        stmp=stmp.substr(0, stmp.length()-1);
      }
      if(stmp[0]=='+'){
        stmp=stmp.substr(1, stmp.length());
      }else if(stmp[0]=='-'){
        if(stmp.length()<2)stmp="-1";
        else stmp="0 "+stmp;
      }
      if(stmp!="")prefac=ReadExpression(stmp);
    }
    return prefac;
  }

  Eigen::MatrixXcd read_within_braces(const std::string &str){
    int dim=0;
  
    {
      std::stringstream ss(str);
      if(Reader::is_next_in_stream(ss,"product")){
        std::string stmp=Reader::find_matching_brace_block_complain(ss,'{','}');
        M=ReadOperator("{"+stmp+"}");
        while(Reader::find_matching_brace_block(ss, stmp, '{', '}')){
          Eigen::MatrixXcd M2=OuterProduct(M, ReadOperator("{"+stmp+"}"));
          M=M2;
        }
        return M;
      }
    }

    int colon_pos=str.find_first_of(':');
    if(colon_pos==std::string::npos || colon_pos==0){
      std::cerr<<"Error reading operator: no ':' found!"<<std::endl;
      exit(1);
    }

    dim=Reader::readSizeT(str.substr(0,colon_pos));
    M=Eigen::MatrixXcd::Zero(dim,dim);
     
    const std::string BSERROR="Error reading operator: cannot identify basis element!";

    std::stringstream ss(str.substr(colon_pos+1));
    std::stringstream ss2;
    char c;
    bool final_was_basis=false;

    while(c=ss.get(), !ss.eof()){
      if(c=='I' && Reader::is_next_in_stream(ss, "d")){
        std::cout<<"Id found"<<std::endl;
        M += eval_prefac(ss2.str()) * Eigen::MatrixXcd::Identity(dim, dim); 
        ss2.str("");
        final_was_basis=true;
        continue;
      }
      if(c=='|'){
        int i,j; 
        ss>>i; 
        if(ss.fail()){std::cerr<<BSERROR<<std::endl;exit(1);}
//        std::cout<<"i: "<<i<<std::endl;
        ss>>c;
        if(c!='>'){std::cerr<<BSERROR<<std::endl;exit(1);}
        ss>>c;
        if(c!='<'){std::cerr<<BSERROR<<std::endl;exit(1);}
        ss>>j; 
        if(ss.fail()){std::cerr<<BSERROR<<std::endl;exit(1);}
//        std::cout<<"j: "<<j<<std::endl;
        ss>>c;
        if(c!='|'){std::cerr<<BSERROR<<std::endl;exit(1);}

        if(i>=dim||j>=dim||i<0||j<0){
          std::cerr<<"Error reading operator: out of bounds!"<<std::endl;
          exit(1);
        }

        M(i,j)+=eval_prefac(ss2.str());       
        ss2.str("");
        final_was_basis=true;
        continue;
      }else{
        ss2<<c;
        if(final_was_basis )final_was_basis=false;
      }
    } 
     
    if( (!final_was_basis) && Reader::trim(ss2.str())!=""){
      std::cerr<<"Operator expression must end in basis element!"<<std::endl;
      std::cerr<<"'"<<Reader::trim(ss2.str())<<"'"<<std::endl;
      exit(1);
    }
 
    return M;
  } 

  Eigen::MatrixXcd read(const std::string &str){
    std::stringstream ss(str);
    std::string str2=Reader::find_matching_brace_block_complain(ss,'{','}');
    return read_within_braces(str2);
  }

  operator Eigen::MatrixXcd ()const{return M;}

  ReadOperator(const std::string &str){
    read(str);
  }
};


#endif
