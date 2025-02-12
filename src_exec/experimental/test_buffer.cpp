#include "BufferedContainer.hpp"
#include "Parameters.hpp"
#include "DummyException.hpp"

using namespace ACE;
/*
class BufferedInt: public BufferedElement{
public:
  int value;

  virtual void read_binary(std::istream &is){
    std::string line;
    std::getline(is, line);
    std::istringstream iss(line);
    if(!(iss>>value)){ 
      std::cerr<<"Can't read BufferedInt!"<<std::endl;
      throw DummyException();
    }
  }
  virtual void write_binary(std::ostream &os)const{
    os<<value<<std::endl;
  }
};
*/

void print_all(BufferedContainer<BufferedInt> &container){
  for(int i=0; i<container.n_tot; i++){
    std::cout<<"container.get("<<i<<").value=";
    std::cout<<container.get(i).value<<std::endl; 
  }
}

int main(int args, char ** argv){
  Parameters param(args, argv);

  std::cout<<"** Creating empty container **"<<std::endl;
  BufferedContainer<BufferedInt> container;
  container.print_info(); std::cout<<std::endl;


  std::cout<<std::endl;
  std::cout<<"** Appending numbers **"<<std::endl;
  container.push_back(BufferedInt(12));
  container.push_back(BufferedInt(42));
  container.push_back(BufferedInt(-17));
  container.push_back(BufferedInt(101));
  container.push_back(BufferedInt(404));
  container.print_info(); std::cout<<std::endl;
  print_all(container);
 
 
  std::cout<<std::endl;
  std::cout<<"** Resizing **"<<std::endl;
  container.resize(7);
  container.print_info(); std::cout<<std::endl;
  print_all(container);


  std::cout<<std::endl;
  std::cout<<"** Copy to writing container **"<<std::endl;
 {
  BufferedContainer<BufferedInt> container_write;
  container_write.initialize("test.dat", 3);
  container_write.copy_content(container);
  container_write.print_info(); std::cout<<std::endl;
  print_all(container_write);
  std::cout<<"written"<<std::endl;
 }


  std::cout<<"** Reading container **"<<std::endl;
  BufferedContainer<BufferedInt> container_read("test.dat");
  container_read.print_info(); std::cout<<std::endl;
  print_all(container_read);
 
  std::cout<<"container.get(2)="<<container_read.get(2).value<<std::endl;

  return 0;
}
