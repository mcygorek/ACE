#ifndef ACE_ENVIRONMENT_PROJECTOR_SCHUR
#define ACE_ENVIRONMENT_PROJECTOR_SCHUR

#include "EnvironmentProjector.h"

namespace ACE{

class EnvironmentProjector_Schur: public EnvironmentProjector{
public:
  RankCompressor_SVD compr;
  double keep_nonoverlapping;


  struct relevantVector{
    Eigen::VectorXcd v;   //Schur Vector orthonormalized wrt. selected vectors
    double row_weight;    //sqrt(sum (row of T) )
    double init_overlap; //overlap (modsqr) with initial bath state
    double norm; //overlap (modsqr) with selected relevant vectors
  };

  double calculate_relevance(const relevantVector &rV){
    return rV.row_weight * (rV.init_overlap + keep_nonoverlapping) * rV.norm;
  }

  virtual Eigen::MatrixXcd get_Q(const Eigen::MatrixXcd &Mdis, 
                                 const Eigen::VectorXcd &bath_init){ 

     int ML=get_vector_dim_min(bath_init, 1, "EnvironmentProjector_Schur::get_Q: bath_init");
     int NL=Mdis.rows()/ML;
     check_matrix_square_eq(Mdis, NL*ML, "EnvironmentProjector_Schur::get_Q: Mdis");
     check_at_least(NL, 4, "EnvironmentProjector_Schur::get_Q: NL");
     

     Eigen::MatrixXcd Q=Eigen::MatrixXcd::Identity(ML,ML);
     if(compr.has_effect()){
       std::vector<Eigen::VectorXcd> selected;
//       std::vector<relevantVector> relVec(NL*NL*ML);

       for(int a1=0; a1<NL; a1++){
         for(int a2=0; a2<NL; a2++){
           //From environment part of propagator with fixed a1,a2:
           Eigen::MatrixXcd tmp(ML,ML);
           for(int d1=0; d1<ML; d1++){
             for(int d2=0; d2<ML; d2++){
               tmp(d1,d2)=Mdis(a1*ML+d1, a2*ML+d2);
             } 
           }
           //Calculate Schur form
           Eigen::ComplexSchur<Eigen::MatrixXcd> schur(tmp);
           Eigen::MatrixXcd T=schur.matrixT();
           Eigen::MatrixXcd U=schur.matrixU();

//std::cout<<"TEST! "<<a1<<" "<<a2<<std::endl;
#ifdef PRINT_SCHUR_T
std::stringstream ss; ss<<"SchurT_"<<a1<<"_"<<a2<<".dat";
std::ofstream ofs(ss.str().c_str());
ofs<<std::setprecision(3);
for(int i=0; i<T.rows(); i++){
  for(int j=0; j<T.cols(); j++){
    ofs<<std::abs(T(i,j))<<" ";
  }
  ofs<<std::endl;
}
#endif
//           std::cout<<"'relevance' for a1="<<a1<<" and a2="<<a2<<":";
           for(int d=0; d<ML; d++){
             relevantVector rV;
              
             rV.v=U.col(d);
             rV.row_weight=0.;
             for(int d2=d; d2<ML; d2++)rV.row_weight+=abs(T(d,d2));
             rV.row_weight=sqrt(rV.row_weight);
             std::complex<double> ov=rV.v.dot(bath_init);
             rV.init_overlap=(std::conj(ov)*ov).real();

             //remove overlap with selected vectors
             rV.norm=1.;
             for(size_t l=0; l<selected.size(); l++){
               std::complex<double> ov=selected[l].dot(rV.v);
               rV.v-=ov*selected[l];
               double norm_=rV.v.norm();
               rV.norm*=norm_;
               rV.v/=norm_;
             }
             double relevance=calculate_relevance(rV);

 //            std::cout<<" "<<relevance;

             if(calculate_relevance(rV)>compr.threshold){
               selected.push_back(rV.v);
             }
//             relVec[(a1*NL+a2)*ML+d]=rV;
           }
//           std::cout<<std::endl;  

         } 
       }
       std::cout<<"Selected: "<<selected.size()<<"/"<<ML<<std::endl;
       Q=Eigen::MatrixXcd(ML, selected.size());
       for(size_t l=0; l<selected.size(); l++){
         Q.col(l)=selected[l];
       }
       std::cout<<"max_diff_from_ortho(Q): "<<max_diff_from_ortho(Q)<<std::endl;
     }
     return Q;
  }
  virtual Eigen::MatrixXcd get_inv_Q(const Eigen::MatrixXcd &Q) {
    return Q.adjoint();
  }
  void setup(Parameters &param){
    compr.setup(param);
    keep_nonoverlapping=param.get_as_double("keep_nonoverlapping");
    if(keep_nonoverlapping<0){
      std::stringstream ss; ss<<"EnvironmentProjector_Schur::setup: keep_nonoverlapping must be non-negative!"<<std::endl;
      throw(std::runtime_error(ss.str()));
    }
  }
  EnvironmentProjector_Schur(Parameters &param){
    setup(param);
  }
  virtual ~EnvironmentProjector_Schur(){}
};

}//namespace

#endif
