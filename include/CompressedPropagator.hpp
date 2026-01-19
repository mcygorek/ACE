#ifndef ACE_COMPRESSED_PROPAGATOR_DEFINED_H
#define ACE_COMPRESSED_PROPAGATOR_DEFINED_H

#include "CompressionTree.hpp"
#include "ProcessTensorElement.hpp"
#include "ModePropagatorGenerator.hpp"

namespace ACE{

namespace CompressedPropagator{

typedef std::pair<ProcessTensorElement, Eigen::VectorXcd> PTE_rhoE0;

void compress(PTE_rhoE0 &e_rhoE0,
                       std::shared_ptr<CompressionTree> TTree,
                       std::shared_ptr<CompressionTree> TTree_inv);

PTE_rhoE0 calculate_single(
                       std::shared_ptr<CompressionTree> TTree,
                       std::shared_ptr<CompressionTree> TTree_inv,
//                       std::shared_ptr<ModePropagatorGenerator> mpg, int mode,
                       std::shared_ptr<ModePropagator> mpp,
                       double ta, double dt, double dict_zero);

PTE_rhoE0 calculate_ACE_sequential(
                       std::shared_ptr<CompressionTree> TTree, 
                       std::shared_ptr<CompressionTree> TTree_inv,
                       std::shared_ptr<ModePropagatorGenerator> mpg, 
                       int mode, double ta, double dt, double dict_zero);


}
}
#endif
