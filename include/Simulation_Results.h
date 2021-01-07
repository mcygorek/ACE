#ifndef SIMULATION_RESULTS_DEFINED_H
#define SIMULATION_RESULTS_DEFINED_H

#include <Eigen/Core>
#include <Eigen/Eigenvalues> 
#include "Operators.h"
#include "FT_Parameters.h"

class Output_Ops{
public:
  std::vector<Eigen::MatrixXcd> ops;
  std::vector<int> proj;

  size_t size()const{return ops.size();}
  Eigen::MatrixXcd & operator[](size_t i){return ops[i];}
  const Eigen::MatrixXcd & operator[](size_t i)const {return ops[i];}

  void add(const Eigen::MatrixXcd &op){
    if(op.rows()!=op.cols()){
      std::cerr<<"Output_Ops::add: op.rows()!=op.cols()!"<<std::endl;
      exit(1);
    }
    if(ops.size()>0 && ops[0].rows()!=op.rows()){
      std::cerr<<"Output_Ops::add: ops[0].rows()!=op.rows()!"<<std::endl;
      exit(1);
    }
    ops.push_back(op);
  }
  void add_proj(int i){
    if(ops.size()>0 && i>=ops[0].rows()){
      std::cerr<<"Output_Ops::add_proj: i>=ops[0].rows()!"<<std::endl;
      exit(1);
    }
    proj.push_back(i);
  }

  void setup(Parameters &param){
    Operators2x2 op;
    if(param.is_specified("add_Output")){
      for(int i=0; i<param.get_nr_rows("add_Output"); i++){
        std::string str=param.get_as_single_string("add_Output", i);
std::cout<<"add_Output: "<<str<<std::endl;
        add(ReadExpression(str)); 
      }
    }else if(param.is_specified("MTLS")){
      int MTLS=param.get_as_size_t("MTLS",0);
      for(size_t i=0; i<MTLS; i++){
        add( ExpandSingleOp(loop_pow(2, i), op.ketbra(1,1), 
                                        loop_pow(2, MTLS-1-i)) );
      }
    }else{
      add( op.ketbra(1,1) );
      add( op.ketbra(0,0) );
      add( op.ketbra(0,1) );
    }
  }
  Output_Ops(){}
};


typedef std::pair<double, std::vector<std::complex<double> > > Simulation_Results_Entry;


class Simulation_Results{
public: 
  std::vector<Simulation_Results_Entry> list;

  size_t size()const{return list.size();}
  void clear(){return list.clear();}
  void resize(size_t s){return list.resize(s);}
  Simulation_Results_Entry & operator[](size_t i){return list[i];}
  const Simulation_Results_Entry & operator[](size_t i)const {return list[i];}
  void push_back(const Simulation_Results_Entry & r){list.push_back(r);}
  Simulation_Results_Entry & back(){return list.back();}
  const Simulation_Results_Entry & back()const{return list.back();}



  void print(std::ostream &ofs)const{
    for(size_t i=0; i<list.size(); i++){
      ofs<<list[i].first;
      for(size_t j=0; j<list[i].second.size(); j++){
        ofs<<" "<<list[i].second[j].real();
        ofs<<" "<<list[i].second[j].imag();
      }
      ofs<<std::endl;
    }
  }
  void print(const std::string &fname)const{
    std::ofstream ofs(fname.c_str());
    print(ofs);
  }
  void read(const std::string &fname){
    std::ifstream ifs(fname.c_str());
    std::vector<std::string> toks;
    list.clear();
    while(Reader::getRelevantLineTokens(ifs, toks)){
      if(toks.size()%2==0){ 
        std::cerr<<"Error reading Simulation_Results '"<<fname<<"': even number of columns!"<<std::endl;
        exit(1);
      }
      if(toks.size()<1){
        std::cerr<<"Error reading Simulation_Results '"<<fname<<"': no column!"<<std::endl;
        exit(1);
      }
      Simulation_Results_Entry entry;
      entry.first=Reader::readDouble(toks[0]," Simulation_results: time");
      for(int i=0; (2*i+2)<toks.size(); i++){
        std::complex<double> c;
        c.real(Reader::readDouble(toks[2*i+1]," Simulation_results: real"));
        c.imag(Reader::readDouble(toks[2*i+2]," Simulation_results: imag"));
        entry.second.push_back(c);
      }
      list.push_back(entry);
    }
  }

  void set(int step, double t, const Output_Ops & output_Op, const Eigen::MatrixXcd &rho, const Eigen::MatrixXcd *Hamil=NULL){
    if(step>=list.size())list.resize(step+1);
    list[step].first=t;

    if(output_Op.size()<1){
      list[step].second.resize(rho.rows()*rho.cols());
      for(int i=0; i<rho.rows(); i++){
        for(int j=0; j<rho.cols(); j++){
          list[step].second[i*rho.cols()+j]=rho(i,j);
        }
      }
    }else{
      list[step].second.resize(output_Op.size());
      for(size_t o=0; o<output_Op.size(); o++){
        std::complex<double> res=0;
        for(int i=0; i<rho.rows(); i++){
          for(int j=0; j<rho.cols(); j++){
            res+=rho(i,j)*output_Op[o](j,i);
          }
        }
        list[step].second[o]=res;
      }
    }  

    if(Hamil!=NULL && output_Op.proj.size()>0){
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(*Hamil);
      int dim=solver.eigenvalues().size();
      if(dim!=rho.rows()){
        std::cerr<<"Error calculating projection: dim!=rho.rows()!"<<std::endl;
        exit(1);
      }
      for(size_t i=0; i<output_Op.proj.size(); i++){
        if(output_Op.proj[i]>=dim){
          std::cerr<<"Error calculating projection: output_Op.proj[i]>=dim!"<<std::endl;
          exit(1);
        }
        Eigen::VectorXcd EV=solver.eigenvectors().col(output_Op.proj[i]);
        std::complex<double> c=EV.adjoint()*rho*EV;
        list[step].second.push_back(c);
      }
    }
  }

  void add_back(int step, const std::complex<double> &c){
    if(step>=list.size()){
      std::cerr<<"Simulation_Results::add_back: step>=list.size()!"<<std::endl;
      exit(1);
    }
    list[step].second.push_back(c);
  }
  void add_back_nan(int step){
    if(step>=list.size()){
      std::cerr<<"Simulation_Results::add_back: step>=list.size()!"<<std::endl;
      exit(1);
    }
    list[step].second.push_back(std::complex<double>(1./0., 1./0.));
  }

  Simulation_Results(){}
  Simulation_Results(const std::string &str){
    read(str);
  }
};


#endif
