#include "DummyException.hpp"

namespace ACE{

const char * DummyException::what () const throw (){
  return "ACE exception";
}

};
