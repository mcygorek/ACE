#include "Parameters.hpp"
#include "PCH.hpp"
//#include <map>
//#include <set>
#include <cmath>
#include "Reader.hpp"
#include "ReadExpression.hpp"
#include "DummyException.hpp"
//#include "Printable.hpp"

namespace ACE{

  //append string to Parameters_Entry of key; creates key if it doesn't exist
  void Parameters::add_to(const std::string & key, const std::vector<std::string> &arg){
    Iterator it=map.find(key);
    if(it==map.end()){
      map.insert(std::make_pair(key,Parameters_Entry(1,arg)));
    }else{
      it->second.push_back(arg);
    }
  }
  
  void Parameters::add_if_not_specified(const std::string & key, const std::vector<std::string> &arg){
    Iterator it=map.find(key);
    if(it==map.end()){
      map.insert(std::make_pair(key,Parameters_Entry(1,arg)));
    }
  }  
 
  void Parameters::set_requested(const std::string &key){
    if(requested.find(key)==requested.end()){
      requested.insert(key);
    }
  }
  const Parameters_Entry & Parameters::get(const std::string &key){
    Iterator it=map.find(key);
    if(it==map.end()){
      return dummy;
    }else{
      set_requested(key);
      return it->second;
    }
  }

  void Parameters::complain_if_not_specified(const std::string &key)const{
    if(!is_specified(key)){
      std::cerr<<"Please specify parameter '"<<key<<"'!"<<std::endl;
      throw DummyException();
    }
  }
  void Parameters::complain_if_row_shorter(const std::string &key, int n, int row, const std::string & context)const{
    if(get_nr_rows(key)<row+1 || get_nr_cols(key,row)<n){
      std::cerr<<"Please specify "<<n<<" arguments for parameter '"<<key<<"'!"<<std::endl;
      if(context!="")std::cerr<<"Context: "<<context<<std::endl;
      throw DummyException();
    }
  }
  void Parameters::complain_if_conflict(const std::string &key1, const std::string &key2)const{
    if(is_specified(key1) && is_specified(key2)){
      std::cerr<<"Please do not specify two conflicting parameters '"<<key1<<"' and '"<<key2<<"'!"<<std::endl;
      throw DummyException();
    }
  }
  void Parameters::override_param(const std::string & key, const std::vector<std::string> &arg){
    erase(key);
    add_to(key, arg);
  }
  void Parameters::override_param(const std::string & key, const std::string &arg){
    override_param(key, std::vector<std::string>(1, arg));
  }
  void Parameters::override_param(const std::string & key, double d){
    override_param(key, std::vector<std::string>(1, double_to_string(d)));
  }

  int Parameters::get_nr_rows(const std::string &key)const{
    cIterator it=map.find(key);
    if(it==map.end()){
      return 0;
    }else{
      return it->second.size();
    }
  }

  std::vector<std::string> Parameters::get_row(const std::string &key, int row){
    Parameters_Entry pe=get(key);
    if(row>=(int)pe.size()||row<0){
      std::cerr<<"Parameters: Error: get_row: row>=get(\""<<key<<"\").size() || row<0!"<<std::endl;
      throw DummyException();
    }
    return pe[row];
  }

  int Parameters::get_nr_cols(const std::string &key, int row)const{
    cIterator it=map.find(key);
    int nr_rows;
    if(it==map.end()){
      return 0;
    }else{
      nr_rows=it->second.size();
    }
    if(row>=nr_rows){
      std::cerr<<"Parameters: Error: get_nr_cols: row>=nr_rows!"<<std::endl;
      throw DummyException();
    }
    return it->second[row].size();
  }
 
  std::vector<std::string> Parameters::get_all_strings(const std::string &key){
    std::vector<std::string> sv;
    Parameters_Entry got=get(key);
    for(size_t i=0; i<got.size(); i++){
      for(size_t j=0; j<got[i].size(); j++){
        sv.push_back(got[i][j]);
      }
    }
    return sv;
  } 
  
  std::vector<double> Parameters::get_all_double(const std::string &key){
    std::vector<std::string> sv=get_all_strings(key);
    std::vector<double> dv;
    for(size_t i=0; i<sv.size(); i++){
      dv.push_back(readDouble(sv[i], key));
    }
    return dv;
  }

  std::vector<size_t> Parameters::get_all_size_t(const std::string &key){
    std::vector<std::string> sv=get_all_strings(key);
    std::vector<size_t> iv;
    for(size_t i=0; i<sv.size(); i++){
      iv.push_back(readSizeT(sv[i],key));
    }
    return iv;
  }

  std::string Parameters::get_as_single_string(const std::string &key, int row){
    Parameters_Entry pe=get(key);
    if(row>=(int)pe.size()||row<0){
      std::cerr<<"Parameters: Error: get_as_single_string: row>=get(\""<<key<<"\").size() || row<0!"<<std::endl;
      throw DummyException();
    }
    std::stringstream ss;
    for(size_t i=0; i<pe[row].size(); i++){
      if(i>0)ss<<" ";
      ss<<pe[row][i];
    }
    return ss.str();
  }
  
  std::vector<double> Parameters::get_row_doubles(const std::string &key, int row, int min){
    std::vector<std::string> row_s=get_row(key, row);
    std::vector<double> ret;
    double d;
    for(size_t i=0; i<row_s.size(); i++){
     if(!canReadDouble(row_s[i], d))break;
     ret.push_back(d);
    }
    if((int)ret.size()<min){
      std::cerr<<"Reading parameter '"<<key<<"': "<<ret.size()<<" doubles found where "<<min<<" are required!"<<std::endl;
      throw DummyException();
    }
    return ret;
  } 


  std::string Parameters::get_as_string(const std::string &key, const std::string &def, int row, int col){
   set_requested(key);
    Iterator it=map.find(key);
    if(it==map.end())return def;
    if(row>=(int)it->second.size()||row<0)return def;
    if(col>=(int)it->second[row].size()||col<0)return def;
    return it->second[row][col];
  }

  std::string Parameters::get_as_string_check(const std::string &key, int row, int col){
    set_requested(key);
    Iterator it=map.find(key);
    if(it==map.end()){
      std::cerr<<"Parameter '"<<key<<"' not specified!"<<std::endl;
      throw DummyException();
    }else{
      if(row>=(int)it->second.size()||row<0){
        std::cerr<<"Parameters::get: '"<<key<<"' row>=it->second.size()||row<0!"<<std::endl;
        throw DummyException();
      }
      if(col>=(int)it->second[row].size()||col<0){
        std::cerr<<"Parameters::get: '"<<key<<"' col>=it->second[row].size()||col<0!"<<std::endl;
        throw DummyException();
      }
      return it->second[row][col];
    }
  }

  Eigen::MatrixXcd Parameters::get_as_operator(const std::string &key, Eigen::MatrixXcd def, int row, int col){

    std::string str=get_as_string(key, "", row, col);
#ifdef DEBUG_EXPRESSIONS
    std::cout<<"get_as_operator: key: '"<<key<<"' string: '"<<str<<"'"<<std::endl;
#endif
    if(str=="")return def;
    return ReadExpression(str); //readDouble(str, key);
  }

  double Parameters::get_as_double(const std::string &key, double def, int row, int col){
    std::string str=get_as_string(key, "", row, col);
    if(str=="")return def;
    return readDouble(str, key);
  }

  bool Parameters::get_as_bool(const std::string &key, bool def, int row, int col){
    std::string str=get_as_string(key, "", row, col);
    if(str=="")return def;
    if(str=="true"||str=="TRUE"||str=="True")return true;
    if(str=="false"||str=="FALSE"||str=="False")return false;
    std::cerr<<"Parameters::get_as_bool: Cannot interpret '"<<str<<"' as Boolean value!"<<std::endl;
    throw DummyException();
    return def;
  }

  void Parameters::add_from_stringvec(const std::string &key, const std::vector<std::string> &toks){
    //take argument in { } as a single parameter, (also save {})
    std::vector<std::string> svec;
    for(int i=0; i<(int)toks.size(); i++){
      if(toks[i].size()<1){
        std::cerr<<"Parameters::add_from_stringvec: toks[i].size()<1!"<<std::endl;
        throw DummyException();
      }
      if(toks[i][0]!='{'){svec.push_back(toks[i]); continue;}

      int nr_open=1;
      std::stringstream ss;
      for(; i<(int)toks.size(); i++){
        for(size_t j=0; j<toks[i].length(); j++){
          if(toks[i][j]=='{')nr_open++;
          if(toks[i][j]=='}'){nr_open--;}//if(nr_open<=1) break; <- let this be handled when accessing the parameter
          ss<<toks[i][j];
        }
        if(nr_open==1)break;
        ss<<' ';
      }
      if(nr_open!=1){
        std::cerr<<"Error: Parameters: can't find closing '}'!"<<std::endl;
        throw DummyException();
      }
      svec.push_back(ss.str());
    }

    add_to(key, svec);
  }
  void Parameters::add_from_line(const std::string &line){
      std::vector<std::string> toks;
      tokenize(line, toks);
      if(toks.size()<1)return;
      std::string key=toks[0];
      toks.erase(toks.begin());
      add_from_stringvec(key, toks);
  }
  void Parameters::add_from_file(const std::string &fname){
    std::ifstream ifs(fname.c_str());
    if(!ifs.good()){
      std::cerr<<"Cannot open driver file '"<<fname<<"'!"<<std::endl;
      throw DummyException();
    }
    std::string line;
    while(getRelevantLine(ifs,line)){
      add_from_line(line);
    }
  }

  void Parameters::setup(int args, char** argv, bool read_first_as_driver){
    std::string key="";
    std::vector<std::string> svec;
    for(int i=1; i<args; i++){
      //check if argument is a key (starts with "-" and is not a number)
      std::string arg(argv[i]);

      if(arg.length()<2 || arg[0]!='-' || (arg[1]>='0' && arg[1] <= '9') 
         || arg[1]=='.'){
        if(key==""){
          if(read_first_as_driver){
            key=std::string("driver");
          }else{
            //TODO: complain?
          }
        }   
        svec.push_back(arg);
      }else{
        if(svec.size()>0){
          add_from_stringvec(key, svec);
          svec.clear();
        }
        key=arg.substr(1);
      }
    }
    if(svec.size()>0)add_from_stringvec(key,svec);
   
    Iterator it=map.find("driver");
    if(it!=map.end()){
      for(size_t i=0; i<it->second.size(); i++){
        for(size_t j=0; j<it->second[i].size(); j++){
          add_from_file(it->second[i][j]);
        }
      }
    } 
  }
  std::vector<std::string> Parameters::get_lines()const{
    std::vector<std::string> list;
    for(cIterator it=map.begin(); it!=map.end(); ++it){
      if(it->first=="driver")continue;
      if(it->first=="print_param")continue;
      for(size_t row=0; row<it->second.size(); row++){
        std::stringstream ss;
        ss<<it->first;
        for(size_t col=0; col<it->second[row].size(); col++){
          ss<<" "<<it->second[row][col];
        } 
        list.push_back(ss.str());
      }
    }
    return list;
  }
  std::ostream & Parameters::print(std::ostream & ofs)const{
    for(cIterator it=map.begin(); it!=map.end(); ++it){
      if(it->first=="driver")continue;
      if(it->first=="print_param")continue;
      for(size_t row=0; row<it->second.size(); row++){
        ofs<<it->first;
        for(size_t col=0; col<it->second[row].size(); col++){
          ofs<<" "<<it->second[row][col];
        }
        ofs<<std::endl;
      }
    }
    return ofs;
  }
  void Parameters::print(const std::string &fname)const{
    std::ofstream ofs(fname.c_str());
    print(ofs);
  }
  
  void Parameters::add_from_prefix(const std::string prefix, Parameters &other){
    std::string prefix_=std::string(prefix+"_");

    for(cIterator it=other.map.begin(); it!=other.map.end(); ++it){
      size_t pos=it->first.find(prefix_);
      if(pos==std::string::npos)continue;
 
      std::string key=it->first.substr(prefix_.length());
     
      Iterator it2=map.find(key);
      if(it2==map.end()){
        map.insert(std::make_pair(key,it->second));
      }else{
        for(size_t i=0; i<it->second.size(); i++){
          it2->second.push_back(it->second[i]);
        }
      }
    }
  }
  

}//namespace
