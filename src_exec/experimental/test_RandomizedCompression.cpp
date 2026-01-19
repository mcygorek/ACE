#include "ACE.hpp"
#include "ProcessTensorBuffer.hpp"
#include "RandomizedCompression.hpp"

using namespace ACE;


int main(int args, char ** argv){
  Parameters param(args, argv, true);
  try{
    std::vector<std::shared_ptr<ModePropagatorGenerator> > mpgs = \
                                                MPG_Selector(param);
    if(mpgs.size()<=0 || !mpgs[0]){
      std::cerr<<"mpgs.size()<=0 || !mpgs[0]!"<<std::endl;
      throw DummyException();
    }
    std::cout<<"mpgs[0]->get_N_modes()="<<mpgs[0]->get_N_modes()<<std::endl;

    TimeGrid tgrid(param);
    TruncationLayout trunc(param);
    int random_chi=param.get_as_double_check("random_chi");
    double dict_zero=param.get_as_double("dict_zero",-1.);
    std::string write_PT=param.get_as_string_check("write_PT");
    int buffer_blocksize=param.get_as_int("buffer_blocksize",-1);

    if(mpgs[0]->get_N_modes()>1){
      ProcessTensorBuffer PTB;
      PTB.set_new_file(write_PT, buffer_blocksize);
      PTB.set_from_ModePropagator(*mpgs[0]->get_ModePropagator(0).get(), tgrid, dict_zero);
 
      ProcessTensorBuffer PTB2;
      PTB2.set_new_temporary(PTB);
      PTB2.set_from_ModePropagator(*mpgs[0]->get_ModePropagator(1).get(), tgrid, dict_zero);
     
std::cout<<"BEFORE"<<std::endl;
      Randomized_Combine(PTB, PTB2, random_chi, trunc.get_base());
std::cout<<"AFTER"<<std::endl;
    }


  }catch (DummyException &e){
    return 1;
  }

  return 0;
}
