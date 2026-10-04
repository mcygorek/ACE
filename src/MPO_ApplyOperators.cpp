#include "MPO_op.hpp"
#include "Parameters.hpp"
#include "CCTEMPO.hpp"
#include "otimes.hpp"
#include "MPO_ApplyOperators.hpp"

namespace ACE{


void MPO_ApplyOperators::apply(int step, Comb & comb, const TruncatedSVD &trunc, int verbosity){
  for(size_t i=0; i<list.size(); i++){
    if(list[i].first==step){
      std::cout<<"MPO_ApplyOperator at step="<<step<<std::endl;
      if(comb.size()!=list[i].second.size()){
        std::cerr<<"MPO_ApplyOperators::apply: comb.size()!=list[i].second.size()!"<<std::endl;
        throw DummyException();
      }
      comb.join_d3_d1_compress_d2_d4(list[i].second, trunc, verbosity );
      comb.sweep_d4_d2(trunc);
    }
  }
 
}

void MPO_ApplyOperators::setup(Parameters & param){
  list.clear();
  int N_sites=CCTEMPO::get_N_sites(param);

  TimeGrid tgrid(param);
  double threshold=param.get_as_double("apply_Operator_threshold", 1e-14);
  TruncatedSVD trunc(threshold);

  std::vector<std::vector<std::string> > entry = param.get("apply_Operator_sum_left");
  for(size_t i=0; i<entry.size(); i++){
    if(entry[i].size()<N_sites+1){
      std::cerr<<"apply_Operator_left TIME {OP_1} {OP_2} .. {OP_N_sites}!"<<std::endl;
      throw DummyException();
    }
    double t=readDouble(entry[i][0], "apply_Operator_sum_left: TIME");
    std::vector<Eigen::MatrixXcd> mat_list(N_sites);
    for(int l=0; l<N_sites; l++){
      Eigen::MatrixXcd mat=ReadExpression(entry[i][1+l]);
      mat_list[l] = otimes(mat, Eigen::MatrixXcd::Identity(mat.rows(), mat.cols()));
    }
    list.push_back(std::make_pair<int, MPO>(
          tgrid.get_closest_n(t), 
          MPO_op::OneBody_Sum(mat_list)));
    list.back().second.sweep_d2_d4(trunc,1);
    list.back().second.sweep_d4_d2(trunc,1);
/*
    std::cout<<"TEST:"<<std::endl;
    std::cout<<MPO_op::evaluate(list.back().second,{0,0},{0,0})<<std::endl;
    std::cout<<MPO_op::evaluate(list.back().second,{1,0},{0,0})<<std::endl;
    std::cout<<MPO_op::evaluate(list.back().second,{0,0},{1,0})<<std::endl;
    std::cout<<MPO_op::evaluate(list.back().second,{1,1},{0,0})<<std::endl;
    std::cout<<MPO_op::evaluate(list.back().second,{1,0},{1,0})<<std::endl;
    std::cout<<mat_list[0]<<std::endl;
    std::cout<<mat_list[1]<<std::endl;
*/
  }
   
  entry = param.get("apply_Operator_sum_right");
  for(size_t i=0; i<entry.size(); i++){
    if(entry[i].size()<N_sites+1){
      std::cerr<<"apply_Operator_left TIME {OP_1} {OP_2} .. {OP_N_sites}!"<<std::endl;
      throw DummyException();
    }
    double t=readDouble(entry[i][0], "apply_Operator_sum_right: TIME");
    std::vector<Eigen::MatrixXcd> mat_list(N_sites);
    for(int l=0; l<N_sites; l++){
      Eigen::MatrixXcd mat=ReadExpression(entry[i][1+l]);
      mat_list[l] = otimes(Eigen::MatrixXcd::Identity(mat.rows(), mat.cols()), mat.transpose());
    }
    list.push_back(std::make_pair<int, MPO>(
          tgrid.get_closest_n(t), 
          MPO_op::OneBody_Sum(mat_list)));
    list.back().second.sweep_d2_d4(trunc,1);
    list.back().second.sweep_d4_d2(trunc,1);
  }

}

}//namespace
