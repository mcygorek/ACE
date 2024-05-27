#ifndef SMART_PTR_DEFINED_H
#define SMART_PTR_DEFINED_H

#include <vector>
#include <iostream>


/** ATTENTION: This is an experimental structure to implement a basic smart 
    pointer without requiring C++11.

    It should work for managing functions, but it may cause issues with other
    types. It will cause memory leaks if "T" is a type whose destructor is 
    not _virtual_. 

*/
namespace ACE{

template <typename T> class Smart_Ptr{
public:
  size_t *counter;
  T *obj;

public:

  T* operator->(){
    return obj;
  }
  const T* operator->()const{
    return obj;
  }
  operator T*() {
    return obj;
  } 
  const  T & ref()const { return *obj;}
  T & ref() { return *obj;}

  size_t get_count()const{ 
    if(counter==NULL)return 0;
    else return *counter; 
  }

  
  void dealloc(){

    if(counter==NULL){
     //do nothing if not initialized (counter==0)   
    }else if(*counter==1){
      delete counter;
      delete obj;
      counter=NULL;
    }else if(*counter>1){
      --(*counter);
//      counter=NULL;  //<-This particular instance shall have no access
    }

  }
/*  void copy_obj(const T &o){
    counter=new size_t(1);
    obj=new T(o);
  }
  Smart_Ptr<T> & operator=(const T &o){
    dealloc();
    copy_obj(o);
    return *this;
  }
  Smart_Ptr(const T &o){
    copy_obj(o);
  }
*/
  template<typename T2> void copy(const Smart_Ptr<T2> &p){
    if(p.counter==NULL){
      std::cerr<<"Smart_Ptr: Error: copy uninstantiated Smart_Ptr!"<<std::endl;
      exit(1);
    }
    counter=p.counter;
    ++(*counter);
    obj=(T*)p.obj;
  }


  Smart_Ptr<T> & operator=(const Smart_Ptr<T> &p){ 

#ifdef DEBUG_SMARTPTR
std::cout<<"SMART POINTER ASSIGNMENT OPERATOR (from other Smart_Ptr) CALLED!"<<std::endl;
std::cout<<"BEFORE: get_count(): "<<p.get_count()<<std::endl;
#endif

    dealloc();
    copy(p);

#ifdef DEBUG_SMARTPTR
std::cout<<"NOW: get_count(): "<<get_count()<<" other: "<<p.get_count()<<std::endl;
#endif

    return *this;
  }

  template<typename T2> Smart_Ptr(const Smart_Ptr<T2> &p){

#ifdef DEBUG_SMARTPTR
std::cout<<"SMART POINTER COPY CONSTRUCTOR CALLED!"<<std::endl;
std::cout<<"BEFORE: get_count(): "<<p.get_count()<<std::endl;
#endif
    
    copy(p);

#ifdef DEBUG_SMARTPTR
std::cout<<"NOW: get_count(): "<<get_count()<<" other: "<<p.get_count()<<std::endl;
#endif
  }
 
  Smart_Ptr(const Smart_Ptr<T> &p){

#ifdef DEBUG_SMARTPTR
std::cout<<"SMART POINTER COPY CONSTRUCTOR CALLED!"<<std::endl;
std::cout<<"BEFORE: get_count(): "<<p.get_count()<<std::endl;
#endif
    
    copy(p);

#ifdef DEBUG_SMARTPTR
std::cout<<"NOW: get_count(): "<<get_count()<<" other: "<<p.get_count()<<std::endl;
#endif
  }

  Smart_Ptr(){
    obj=NULL;
    counter=NULL;
#ifdef DEBUG_SMARTPTR
std::cout<<"SMART POINTER DEFAULT CONSTRUCTOR CALLED! get_count(): "<<get_count()<<std::endl;
#endif
  }

  Smart_Ptr<T> & operator=(T *p){
    dealloc();
    obj=p;
    counter=new size_t(1);

#ifdef DEBUG_SMARTPTR
std::cout<<"SMART POINTER ASSIGNMENT FROM POINTER CALLED! get_count(): "<<get_count()<<std::endl;
#endif
    return *this;
  }
  Smart_Ptr(T *t){
    obj=t;
    counter=new size_t(1);
#ifdef DEBUG_SMARTPTR
std::cout<<"SMART POINTER CONSTRUCTOR FROM POINTER CALLED! get_count(): "<<get_count()<<std::endl;
#endif
  }

  virtual ~Smart_Ptr(){
#ifdef DEBUG_SMARTPTR
std::cout<<"SMART POINTER DESTRUCTOR CALLED WITH get_count(): "<<get_count()<<std::endl;
#endif
    dealloc();
  }
};

}//namespace
#endif
