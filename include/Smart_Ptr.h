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

template <typename T> class Smart_Ptr{
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
  void copy(const Smart_Ptr<T> &p){
    counter=p.counter;
    ++(*counter);
    obj=p.obj;
  }
  Smart_Ptr<T> & operator=(const Smart_Ptr<T> &p){ 
    dealloc();
    copy(p);
    return *this;
  }
  Smart_Ptr(const Smart_Ptr<T> &p){
    copy(p);
  }
  Smart_Ptr(){
    obj=NULL;
    counter=NULL;
  }
  Smart_Ptr<T> & operator=(T *p){
    dealloc();
    obj=p;
    counter=new size_t(1);
    return *this;
  }
  Smart_Ptr(T *t){
    obj=t;
    counter=new size_t(1);
  }
  virtual ~Smart_Ptr(){
    dealloc();
  }
};


#endif
