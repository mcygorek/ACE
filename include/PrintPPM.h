#ifndef PRINT_PPM_DEFINED_H_
#define PRINT_PPM_DEFINED_H_

#include "Reader.hpp"

namespace ACE{

class PPM_Canvas{
  int width, height;
  double *mem;

public:

  int get_width()const{return width;}
  int get_height()const{return height;}
 
  double & operator()(int i, int j){
    if(i<0||i>=width || j<0||j>=height){
      std::cerr<<"PPM_Canvas::("<<i<<","<<j<<") out of bounds ("<<width<<","<<height<<")!"<<std::endl;
      exit(1);
    }
    return mem[j*width+i];
  } 
  double operator()(int i, int j) const{
    if(i<0||i>=width || j<0||j>=height){
      std::cerr<<"PPM_Canvas::("<<i<<","<<j<<") out of bounds ("<<width<<","<<height<<")!"<<std::endl;
      exit(1);
    }
    return mem[j*width+i];
  } 
  
  void fill(double val){
    for(int i=0; i<width*height; i++)mem[i]=val;
  }
  void allocate(){
    if(width<0 || height<0){
      std::cerr<<"PPM_Canvas::allocate: dth<0 || height<0!"<<std::endl;
    }
    mem=new double[width*height];
  }
  void deallocate(){
    delete [] mem;
  }
  void reallocate(){
    deallocate();
    allocate(); 
  }
  void resize(int w, int h){
    width=w; height=h; reallocate();
    fill(0.);
  }


  void print_pgm(const std::string &fname, double maxval=1., double minval=0.)const{
  
    if(width<1 || height<1){
      std::cerr<<"print_pgm: width<1 || height<1 !"<<std::endl;
      exit(1);
    }

    std::ofstream ofs(fname.c_str());
    ofs<<"P2"<<std::endl;
    ofs<<width<<" "<<height<<std::endl<<"255"<<std::endl;
    for(int r=0; r<width; r++){
      for(int c=0; c<height; c++){
        int x=(255.*(mem[r*height+c]-minval)/(maxval-minval)+0.5);
        if(x>255)x=255; 
        if(x<0)x=0;
        ofs<<x<<" ";
      }
      ofs<<std::endl;
    }
  }

  PPM_Canvas(int w=0, int h=0): width(w), height(h){
    allocate();
  }
};



}//namespace
#endif
