#ifndef ACE_PT_SUM_DIFF_DEFINED_H
#define ACE_PT_SUM_DIFF_DEFINED_H

#include "ProcessTensorElement.hpp"

namespace ACE{

/*
Implements primivites for additive combination and splitting of PTs.
Functions operate on a connection involving 4 ProcessTensorElements

-E11-E12-       -|         |-
  |   |    <=>   |E'11-E'12|
-E21-E22-       -|         |-

Functions will use same memory for E11<> E'11 and E12<>E'12
*/

//returns false if inner bond E11-E12 has dimension < 2
namespace AddPT{
bool split_alternating(ProcessTensorElement &E11, ProcessTensorElement &E12,
                       ProcessTensorElement &E21, ProcessTensorElement &E22);

void add(ProcessTensorElement &E11,
         ProcessTensorElement &E21);

void add_head(ProcessTensorElement &E11,
                ProcessTensorElement &E21);

void add_tail(ProcessTensorElement &E11,
                ProcessTensorElement &E21);



}//namespace
}//namespace
#endif
