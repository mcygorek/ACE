#include "PCH.hpp"
#include <iostream>
#include "Printable.hpp"

namespace ACE{

std::ostream& operator<<(std::ostream& os, const Printable& prn){
  return prn.print(os);
}

}//namespace
