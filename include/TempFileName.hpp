#ifndef ACE_TEMP_FILE_NAME_DEFINED_H
#define ACE_TEMP_FILE_NAME_DEFINED_H

#include <string>
#include <iostream>
#include "ReaderBasics.hpp"

namespace ACE{
/* choose temporary file name. Remove file with this name when object is destroyed." */
class TempFileName{
public:
  std::string fname;
  bool noremove;
  
  static std::string randomString(int length);
  
  inline const std::string & get()const{ return fname; }
  inline const char* c_str()const{ return fname.c_str(); }
  inline operator const std::string & () const {return fname;}
  inline void set_explicit(const std::string & fname_){ fname = fname_; }
  TempFileName & set_noremove(bool nr);
 
  void swap(TempFileName &other); 
  void remove();
  void initialize(const std::string &str="");

  inline TempFileName(const std::string &str=""){
    initialize(str);
  }
  ~TempFileName(){
//std::cout<<"TempFileName destructor for '"<<get()<<"' called."<<std::endl;
//print_file_exists(get());
    remove();
  }
};
}//namespace
#endif
