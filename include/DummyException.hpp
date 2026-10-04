#ifndef ACE_DUMMY_EXCEPTION_DEFINED_H
#define ACE_DUMMY_EXCEPTION_DEFINED_H

#include "PCH.hpp"

namespace ACE{
/* Dummy exception: 
   Reason for exception is printed already before throwing. 
   Needed to clean up file streams
*/
class DummyException: public std::exception{
public:   
  const char * what () const throw (); 
};
}//namespace
#endif
