#ifndef ACE_INDEXED_PT_DEFINED_H
#define ACE_INDEXED_PT_DEFINED_H

/** Purpose: 
Store indexed PT where multiple matrices can be identical and have to be stored
only once. Used for log(t) contraction of PTs with finite memory using 
divide-and-conquer algorithm. 
At its core is essentially an implementation of shared_ptr, but also the DnC 
algorithm has to be part of this class


Optional (TODO): implement efficient IO: keep recently touched elements in buffer; preload prospective next PT elements, in particular in the active region
*/

namespace ACE{

class IndexedPT_Element{
private:
  int count;
public:
  ProcessTensorElement *e;

// Functions: 
// - destruct only when no other element points to it
// - when modified, other indexed elements that share the same phyiscal "e" should be left untouched => let go of ownership and create new one.
// This implies that we need a "get_rw" and a "get_ro" function, where the former triggers the physical copy.
// - implement "copy" function that only copies ownership
// - implement binary read/write functions.
  
};

//TODO: inherit from "ProcessTensorForward" for easy integration with rest of the code
class IndexedPT{
public:
  std::vector<IndexedPT_Element> list;

};


}//namespace
#endif
