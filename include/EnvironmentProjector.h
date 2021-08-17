#ifndef ACE_ENVIRONMENT_PROJECTOR_H
#define ACE_ENVIRONMENT_PROJECTOR_H

#include "InfluenceFunctional_OD.h"
#include <iomanip>


MPS_Matrix Disentangled_MatrixXcd_to_MPS_Matrix(
                            const Eigen::MatrixXcd & Mdis, int N){
  int NL=N*N;
  check_matrix_square_min(Mdis, NL, "Disentangled_MatrixXcd_to_MPS_Matrix");
  int ML=Mdis.rows()/NL;
  check_matrix_rows_eq(Mdis, ML*NL, "Disentangled_MatrixXcd_to_MPS_Matrix");

  MPS_Matrix A(NL*NL, ML, ML);
  for(int a1=0; a1<NL; a1++){
    for(int a2=0; a2<NL; a2++){
      for(int d1=0; d1<ML; d1++){
        for(int d2=0; d2<ML; d2++){
          A(a1*NL+a2, d1, d2)=Mdis(a1*ML+d2, a2*ML+d1);
        }
      }
    }
  }
  return A;
}

InfluenceFunctional_OD IF_from_Disentangled_MatrixXcd(
        const IF_TimeGrid &tgrid, const Eigen::MatrixXcd & Mdis, 
        const Eigen::VectorXcd & bath_init_vector, 
        const Eigen::VectorXcd & bath_end_vector){

  int ML=get_vector_dim_min(bath_init_vector, 1, "IF_from_Disentangled_MatrixXcd: init");
  check_vector_dim_eq(bath_end_vector, ML, "IF_from_Disentangled_MatrixXcd: end vs. init");
//  int M=get_dim_sqrt(ML, "IF_from_Disentangled_MatrixXcd: init");

  check_matrix_square_min(Mdis, ML, "IF_from_Disentangled_MatrixXcd");
  int NL=Mdis.rows()/ML;
  int N=get_dim_sqrt(NL, "IF_from_Disentangled_MatrixXcd");
  check_matrix_rows_eq(Mdis, ML*NL, "IF_from_Disentangled_MatrixXcd");


  InfluenceFunctional_OD IF(tgrid, N);

  MPS_Matrix A=Disentangled_MatrixXcd_to_MPS_Matrix(Mdis, N);

  for(size_t n=0; n<IF.a.size(); n++){
    IF.a[n]=A;
  }

  IF.reduce_first(bath_init_vector);
  IF.reduce_last(bath_end_vector);

  
  IF.dict.set_default(N);
//  IF.print_dims();
  IF.check_consistency();
  IF.calculate_closures();

  return IF;
}

InfluenceFunctional_OD IF_from_Disentangled_MatrixXcd(
        const IF_TimeGrid &tgrid, const Eigen::MatrixXcd & Mdis, 
        const Eigen::VectorXcd & bath_init_vector){
  
  int M=sqrt(bath_init_vector.rows());
  Eigen::VectorXcd bath_end_vector=
      H_Matrix_to_L_Vector(Eigen::MatrixXcd::Identity(M,M));

  return IF_from_Disentangled_MatrixXcd(tgrid, Mdis, bath_init_vector, bath_end_vector);
}


class EnvironmentProjector{
public:
  //Q works on environment part of disentangled propagator
  virtual Eigen::MatrixXcd get_Q(const Eigen::MatrixXcd &Mdis, 
                                 const Eigen::VectorXcd &bath_init)=0;

  //We also need inverse of Q; if Q is hermitian, adjoint is enough
  virtual Eigen::MatrixXcd get_inv_Q(const Eigen::MatrixXcd &Q)=0;



  //Generates an InfluenceFunctional_OD from the defined implementation of Q
  InfluenceFunctional_OD get_IF_from_MPG(ModePropagatorGenerator &mpg, IF_TimeGrid tgrid){

    int N = mpg.get_N();
    int NL = N*N; 

    Eigen::VectorXcd bath_init = Eigen::VectorXcd::Ones(1);
    Eigen::VectorXcd bath_end = Eigen::VectorXcd::Ones(1);
    Eigen::MatrixXcd Mdis = Eigen::MatrixXcd::Identity(NL,NL);
 
    for(int k=0; k<mpg.get_N_modes(); k++){ 
      std::cout<<"Mode: "<<k<<"/"<<mpg.get_N_modes()<<std::endl;
    
      ModePropagatorPtr mp=mpg.getModePropagator(k);
      mp->update(tgrid.ta, tgrid.dt/2.);

      Eigen::VectorXcd new_bath_init=H_Matrix_to_L_Vector(mp->get_bath_init());
      int ML = new_bath_init.size();
      int M = get_dim_sqrt(ML, "get_IF_test: ML");
    
      bath_init = Vector_otimes(bath_init, new_bath_init);
      bath_end = Vector_otimes(bath_end,
                   H_Matrix_to_L_Vector(Eigen::MatrixXcd::Identity(M,M)));

      //Recall: Disentangle: 
      //from (H_S otimes H_E) otimes (H_S otimes H_E)
      //to   (H_S otimes H_S) otimes (H_E otimes H_E)
      //
      //Here: Expand new propagator for HE2 by inserting an identity matrix 
      // corresponding to HE1 in the middle. 
      //Mdis is expandend by a identity for HE2 at the end:

      Eigen::MatrixXcd Mdis2 = ExpandMatrix( 
             Disentangle_Propagator(mp->M, N), NL, Mdis.rows()/NL );
     
      Mdis = Mdis2 * ExpandMatrix(Mdis, Mdis.rows(), ML) * Mdis2;

 
      Eigen::MatrixXcd Q = get_Q(Mdis, bath_init);
      Eigen::MatrixXcd Qinv = get_inv_Q(Q);
    
      //Recall: Q acts on HE1+HE2 only:
      Mdis = otimes(Eigen::MatrixXcd::Identity(NL,NL),Qinv) 
              * Mdis * otimes(Eigen::MatrixXcd::Identity(NL,NL),Q);

      bath_init = Qinv * bath_init;
      bath_end = Q.transpose() * bath_end;
    }

    return IF_from_Disentangled_MatrixXcd(tgrid, Mdis, bath_init, bath_end);
  }

  virtual ~EnvironmentProjector(){}
};


#endif

