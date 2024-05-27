#ifndef ACE_BUFFERED_ELEMENT_DEFINED_H
#define ACE_BUFFERED_ELEMENT_DEFINED_H
#include <fstream>

namespace ACE{
class BufferedElement{
public:
//  virtual void BufferedElement_read(std::istream &ifs)=0;
//  virtual void BufferedElement_write(std::ostream &ofs)const=0;
  virtual void read_binary(std::istream &is)=0;//, const std::string &context="");
  virtual void write_binary(std::ostream &os)const=0;
  virtual ~BufferedElement()=0;
};

//Example: BufferedInt:
class BufferedInt: public BufferedElement{
public:
  int value;
  virtual void read_binary(std::istream &is);
  virtual void write_binary(std::ostream &os)const;
  BufferedInt(int i=0): value(i){};
  virtual ~BufferedInt();
};

}//namespace
#endif
