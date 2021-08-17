#ifndef ACE_PRINTABLE_DEFINED_H_
#define ACE_PRINTABLE_DEFINED_H_

#include <iostream>

class Printable{
public:
  virtual std::ostream & print(std::ostream & os=std::cout)const=0;

  friend std::ostream& operator<<(std::ostream& os, const Printable& prn);
};

std::ostream& operator<<(std::ostream& os, const Printable& prn){
  return prn.print(os);
}
#endif
