#include "ACE.hpp"
#include "PT_infinite.hpp"
#include "ProcessTensorBuffer.hpp"
#include "ProcessTensorRepeat.hpp"
#include "Largest_EV.hpp"
#include "Timings.hpp"
#include <slepceps.h>

using namespace ACE;

struct SLEPC_MatMult_Context{
  const Eigen::VectorXcd *LambdaB;
  const Eigen::MatrixXcd *GLG2;
};

PetscErrorCode myMatMult(Mat A, Vec vr, Vec vl){  //vl = A*vr
  PetscInt rhigh, rlow;
  const Eigen::MatrixXcd *M;
  PetscCall(MatShellGetContext(A, &M));
  PetscCall(VecGetOwnershipRange(vr,&rlow,&rhigh));
  PetscCall(VecZeroEntries(vl));
  PetscCall(PetscBarrier((PetscObject)vl));
  const std::complex<double>* ar;
  PetscCall(VecGetArrayRead(vr, &ar));

  Eigen::VectorXcd L=(*M).block(0,rlow,M->rows(), rhigh-rlow)*Eigen::Map<const Eigen::VectorXcd>(ar, rhigh-rlow);

  PetscCall(VecRestoreArrayRead(vr, &ar));
  std::vector<PetscInt> indices(L.size());
  for(int i=0; i<L.size(); i++){
    indices[i]=i;
  }
  PetscCall(VecSetValues(vl, L.size(), &indices[0], &L(0), ADD_VALUES));
  PetscCall(VecAssemblyBegin(vl)); 
  PetscCall(VecAssemblyEnd(vl)); 
  return PETSC_SUCCESS;
}

PetscErrorCode SLEPC_PT_iTEBD_R_times(Mat A, Vec vr, Vec vl){  //vl = A*vr
  PetscInt rank;
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));

  const SLEPC_MatMult_Context *ctx;
  PetscCall(MatShellGetContext(A, &ctx));

  const Eigen::VectorXcd & LambdaB=*(ctx->LambdaB); 
  const Eigen::MatrixXcd & GLG2=*(ctx->GLG2); 

  PetscInt rhigh, rlow;
  PetscCall(VecGetOwnershipRange(vr,&rlow,&rhigh));
  PetscCall(VecZeroEntries(vl));
  PetscCall(PetscBarrier((PetscObject)vl));
  const std::complex<double>* ar;
  PetscCall(VecGetArrayRead(vr, &ar));

////  Eigen::VectorXcd VR=Eigen::Map<const Eigen::VectorXcd>(ar, rhigh-rlow);
//  Eigen::MatrixXcd M=PT_iTEBD_calc_R(LambdaB, GLG2);
//  Eigen::VectorXcd VL=M.block(0,rlow,M.rows(), rhigh-rlow)*Eigen::Map<const Eigen::VectorXcd>(ar, rhigh-rlow);

  int dim_d=LambdaB.size();
  int dim_i=GLG2.rows()/dim_d;
  //limited first dimension of intermediate matrix
  int d1low=rlow/dim_d; int d1high=(rhigh-1)/dim_d+1; 
  int dfirst=rlow-d1low*dim_d; int dlast=rhigh-(d1high-1)*dim_d;

if(rhigh-rlow>0){
  typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXcdRM;
  MatrixXcdRM VR=MatrixXcdRM::Zero((d1high-d1low),dim_d);
  Eigen::Map<Eigen::VectorXcd>(&VR(0,dfirst),rhigh-rlow)=Eigen::Map<const Eigen::VectorXcd>(ar, rhigh-rlow);

  for(int d1=d1low; d1<d1high; d1++){
    for(int d2=0; d2<dim_d; d2++){
      if(d1==d1low && d2<dfirst){  d2=dfirst-1; continue;  }
      if(d1==d1high-1 && d2>=dlast){ break;}
      VR((d1-d1low),d2)*=LambdaB(d1)*std::conj(LambdaB(d2));
    }
  }

  MatrixXcdRM VL=MatrixXcdRM::Zero(dim_d,dim_d);
  for(int i=0; i<dim_i; i++){
   for(int j=0; j<dim_i; j++){
      VL.noalias()+=  \
          Eigen::Map<const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>, 0, Eigen::OuterStride<> > (&GLG2(i*dim_d,j*dim_d), dim_d, dim_d, Eigen::OuterStride<>(dim_i*dim_d)).block(0,d1low,dim_d,d1high-d1low) * VR * \
          Eigen::Map<const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 0, Eigen::OuterStride<> > (&GLG2(i*dim_d,j*dim_d), dim_d, dim_d, Eigen::OuterStride<>(dim_i*dim_d)).conjugate();

    }
  }

  PetscCall(VecRestoreArrayRead(vr, &ar));
  std::vector<PetscInt> indices(dim_d*dim_d);
  for(int i=0; i<dim_d*dim_d; i++){
    indices[i]=i;
  }
  PetscCall(VecSetValues(vl, dim_d*dim_d, &indices[0], &VL(0,0), ADD_VALUES));
 }
  PetscCall(VecAssemblyBegin(vl)); 
  PetscCall(VecAssemblyEnd(vl)); 

/*
  Eigen::VectorXcd VL=Eigen::VectorXcd::Zero(dim_d*dim_d);
  for(int i=0; i<dim_i; i++){
   for(int j=0; j<dim_i; j++){
      Eigen::VectorXcd tmp=Eigen::VectorXcd::Zero((d1high-d1low)*dim_d);
      for(int d1=d1low; d1<d1high; d1++){
        for(int d2=0; d2<dim_d; d2++){
          for(int d=0; d<dim_d; d++){
            if(d1==d1low && d<dfirst){  d=dfirst-1; continue;  }
            if(d1==d1high-1 && d>=dlast){ break;}
            tmp((d1-d1low)*dim_d+d2)+=std::conj(GLG2(i*dim_d+d2,j*dim_d+d)*LambdaB(d))*ar[d1*dim_d+d-rlow];
          }
        }
      }
      for(int d1=0; d1<dim_d; d1++){
        for(int d2=0; d2<dim_d; d2++){
          for(int d=d1low; d<d1high; d++){
            VL(d1*dim_d+d2)+=GLG2(i*dim_d+d1, j*dim_d+d)*LambdaB(d)*tmp((d-d1low)*dim_d+d2);
          }
        }
      }
    }
  }

  PetscCall(VecRestoreArrayRead(vr, &ar));
  std::vector<PetscInt> indices(VL.size());
  for(int i=0; i<VL.size(); i++){
    indices[i]=i;
  }
  PetscCall(VecSetValues(vl, VL.size(), &indices[0], &VL(0), ADD_VALUES));
  PetscCall(VecAssemblyBegin(vl)); 
  PetscCall(VecAssemblyEnd(vl)); 
*/
  return PETSC_SUCCESS;
}

PetscErrorCode SLEPC_PT_iTEBD_LT_times(Mat A, Vec vr, Vec vl){  //vl = A*vr
  PetscInt rank;
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));

  const SLEPC_MatMult_Context *ctx;
  PetscCall(MatShellGetContext(A, &ctx));

  const Eigen::VectorXcd & LambdaB=*(ctx->LambdaB); 
  const Eigen::MatrixXcd & GLG2=*(ctx->GLG2); 

  PetscInt rhigh, rlow;
  PetscCall(VecGetOwnershipRange(vr,&rlow,&rhigh));
  PetscCall(VecZeroEntries(vl));
  PetscCall(PetscBarrier((PetscObject)vl));
  const std::complex<double>* ar;
  PetscCall(VecGetArrayRead(vr, &ar));

////  Eigen::VectorXcd VR=Eigen::Map<const Eigen::VectorXcd>(ar, rhigh-rlow);
//  Eigen::MatrixXcd M=PT_iTEBD_calc_R(LambdaB, GLG2);
//  Eigen::VectorXcd VL=M.block(0,rlow,M.rows(), rhigh-rlow)*Eigen::Map<const Eigen::VectorXcd>(ar, rhigh-rlow);

  int dim_d=LambdaB.size();
  int dim_i=GLG2.rows()/dim_d;
  //limited first dimension of intermediate matrix
  int d1low=rlow/dim_d; int d1high=(rhigh-1)/dim_d+1; 
  int dfirst=rlow-d1low*dim_d; int dlast=rhigh-(d1high-1)*dim_d;

if(rhigh-rlow>0){
  typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXcdRM;
  MatrixXcdRM VR=MatrixXcdRM::Zero((d1high-d1low),dim_d);
  Eigen::Map<Eigen::VectorXcd>(&VR(0,dfirst),rhigh-rlow)=Eigen::Map<const Eigen::VectorXcd>(ar, rhigh-rlow);

  for(int d1=d1low; d1<d1high; d1++){
    for(int d2=0; d2<dim_d; d2++){
      if(d1==d1low && d2<dfirst){  d2=dfirst-1; continue;  }
      if(d1==d1high-1 && d2>=dlast){ break;}
      VR((d1-d1low),d2)*=LambdaB(d1)*std::conj(LambdaB(d2));
    }
  }

  MatrixXcdRM VL=MatrixXcdRM::Zero(dim_d,dim_d);
  for(int i=0; i<dim_i; i++){
   for(int j=0; j<dim_i; j++){
      VL.noalias()+=  \
          Eigen::Map<const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 0, Eigen::OuterStride<> > (&GLG2(i*dim_d,j*dim_d), dim_d, dim_d, Eigen::OuterStride<>(dim_i*dim_d)).block(0,d1low,dim_d,d1high-d1low) * VR * \
          Eigen::Map<const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>, 0, Eigen::OuterStride<> > (&GLG2(i*dim_d,j*dim_d), dim_d, dim_d, Eigen::OuterStride<>(dim_i*dim_d)).conjugate();

    }
  }

  PetscCall(VecRestoreArrayRead(vr, &ar));
  std::vector<PetscInt> indices(dim_d*dim_d);
  for(int i=0; i<dim_d*dim_d; i++){
    indices[i]=i;
  }
  PetscCall(VecSetValues(vl, dim_d*dim_d, &indices[0], &VL(0,0), ADD_VALUES));
 }
  PetscCall(VecAssemblyBegin(vl)); 
  PetscCall(VecAssemblyEnd(vl)); 


  //Eigen::VectorXcd VR=Eigen::Map<const Eigen::VectorXcd>(ar, rhigh-rlow);
/*
  Eigen::VectorXcd VL=Eigen::VectorXcd::Zero(dim_d*dim_d);
  for(int i=0; i<dim_i; i++){
   for(int j=0; j<dim_i; j++){
      Eigen::VectorXcd tmp=Eigen::VectorXcd::Zero((d1high-d1low)*dim_d);
      for(int d1=d1low; d1<d1high; d1++){
        for(int d2=0; d2<dim_d; d2++){
          for(int d=0; d<dim_d; d++){
            if(d1==d1low && d<dfirst){  d=dfirst-1; continue;  }
            if(d1==d1high-1 && d>=dlast){ break;}
            tmp((d1-d1low)*dim_d+d2)+=std::conj(GLG2(i*dim_d+d,j*dim_d+d2)*LambdaB(d))*ar[d1*dim_d+d-rlow];
          }
        }
      }
      for(int d1=0; d1<dim_d; d1++){
        for(int d2=0; d2<dim_d; d2++){
          for(int d=d1low; d<d1high; d++){
            VL(d1*dim_d+d2)+=GLG2(i*dim_d+d, j*dim_d+d1)*LambdaB(d)*tmp((d-d1low)*dim_d+d2);
          }
        }
      }
    }
  }

  PetscCall(VecRestoreArrayRead(vr, &ar));
  std::vector<PetscInt> indices(VL.size());
  for(int i=0; i<VL.size(); i++){
    indices[i]=i;
  }
  PetscCall(VecSetValues(vl, VL.size(), &indices[0], &VL(0), ADD_VALUES));
  PetscCall(VecAssemblyBegin(vl)); 
  PetscCall(VecAssemblyEnd(vl)); 
*/
  return PETSC_SUCCESS;
}

PetscErrorCode Apply_PETSC_Shell_to_Eigen_Vector(Eigen::VectorXcd &V, PetscVoidFn *fct, void *ctx){
  int dim=V.size();

  Vec tmp_vec, tmp_vec2, tmp_vec3;
  Mat tmp_A;
  PetscCall(MatCreateShell(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,dim,dim,(void*)ctx,&tmp_A));
  PetscCall(MatShellSetOperation(tmp_A,MATOP_MULT,(PetscVoidFn *)fct));     
  PetscCall(MatCreateVecs(tmp_A,NULL,&tmp_vec));
  PetscCall(MatCreateVecs(tmp_A,NULL,&tmp_vec2));

  VecScatter sf;
  VecScatterCreateToAll(tmp_vec2, &sf, &tmp_vec3);
  PetscCall(VecZeroEntries(tmp_vec));
  for(int d=0; d<dim; d++){ 
    PetscCall(VecSetValue(tmp_vec,d,V(d),INSERT_VALUES));
  }
  PetscCall(VecAssemblyBegin(tmp_vec)); 
  PetscCall(VecAssemblyEnd(tmp_vec)); 

  PetscErrorCode (*fct2)(Mat, Vec, Vec)=(PetscErrorCode (*)(Mat, Vec, Vec))fct;
  PetscCall(fct2(tmp_A, tmp_vec, tmp_vec2));
  //Collect local eigenvector to global xr->yr 
  VecScatterBegin(sf, tmp_vec2, tmp_vec3, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(sf, tmp_vec2, tmp_vec3, INSERT_VALUES, SCATTER_FORWARD);
  //copy into Eigen structure
  const std::complex<double>* ar;
  PetscCall(VecGetArrayRead(tmp_vec3, &ar));
  V=Eigen::Map<const Eigen::VectorXcd>(ar, dim);
  PetscCall(VecRestoreArrayRead(tmp_vec3, &ar));
  //Destroy
  PetscCall(VecScatterDestroy(&sf));
  PetscCall(VecDestroy(&tmp_vec));
  PetscCall(VecDestroy(&tmp_vec2));
  PetscCall(VecDestroy(&tmp_vec3));
  PetscCall(MatDestroy(&tmp_A));

  return PETSC_SUCCESS;
} 

PetscErrorCode Eigen_Mat_from_PETSC_Shell(Eigen::MatrixXcd &R2, int dim, PetscVoidFn *fct, void *ctx){

  R2=Eigen::MatrixXcd::Zero(dim,dim);
  Vec tmp_vec, tmp_vec2, tmp_vec3;
  Mat tmp_A;
  PetscCall(MatCreateShell(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,dim,dim,(void*)ctx,&tmp_A));
  PetscCall(MatShellSetOperation(tmp_A,MATOP_MULT,(PetscVoidFn *)fct));     
  PetscCall(MatCreateVecs(tmp_A,NULL,&tmp_vec));
  PetscCall(MatCreateVecs(tmp_A,NULL,&tmp_vec2));

  VecScatter sf;
  VecScatterCreateToAll(tmp_vec2, &sf, &tmp_vec3);
  for(int d=0; d<dim; d++){ 
    PetscCall(VecZeroEntries(tmp_vec));
    PetscCall(VecSetValue(tmp_vec,d,1.,INSERT_VALUES));
    PetscCall(VecAssemblyBegin(tmp_vec)); 
    PetscCall(VecAssemblyEnd(tmp_vec)); 
    PetscErrorCode (*fct2)(Mat, Vec, Vec)=(PetscErrorCode (*)(Mat, Vec, Vec))fct;
    PetscCall(fct2(tmp_A, tmp_vec, tmp_vec2));
    //Collect local eigenvector to global xr->yr 
    VecScatterBegin(sf, tmp_vec2, tmp_vec3, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(sf, tmp_vec2, tmp_vec3, INSERT_VALUES, SCATTER_FORWARD);
    //copy into Eigen structure
    const std::complex<double>* ar;
    PetscCall(VecGetArrayRead(tmp_vec3, &ar));
    R2.col(d)=Eigen::Map<const Eigen::VectorXcd>(ar, dim);
    PetscCall(VecRestoreArrayRead(tmp_vec3, &ar));
  }
  //Destroy
  PetscCall(VecScatterDestroy(&sf));
  PetscCall(VecDestroy(&tmp_vec));
  PetscCall(VecDestroy(&tmp_vec2));
  PetscCall(VecDestroy(&tmp_vec3));
  PetscCall(MatDestroy(&tmp_A));

  return PETSC_SUCCESS;
}

std::complex<double> Largest_EV_SLEPC(Eigen::VectorXcd &vec, PetscVoidFn *fct, void *ctx, int argc, char **argv) {

  Mat            A;           
  EPS            eps;         
  EPSType        type;
  PetscReal      error,tol,re,im;
  PetscScalar    kr,ki;
  Vec            xr,xi,yr;
  PetscInt       i,Istart,Iend,nev,maxit,its,nconv;
  PetscInt       n=vec.size();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Compute the operator matrix that defines the eigensystem, Ax=kx
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  int dim=vec.size();
  PetscCall(MatCreateShell(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,dim,dim,ctx,&A));
  PetscCall(MatShellSetOperation(A,MATOP_MULT,fct));     


  PetscCall(MatCreateVecs(A,NULL,&xr));
  PetscCall(MatCreateVecs(A,NULL,&xi));

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //            Create the eigensolver and set various options
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  
  // Create eigensolver context
  PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));

  // Set operators. In this case, it is a standard eigenvalue problem
  PetscCall(EPSSetOperators(eps,A,NULL));
  PetscCall(EPSSetProblemType(eps,EPS_NHEP));
  PetscCall(EPSSetWhichEigenpairs(eps,EPS_LARGEST_REAL));
  
  // Set solver parameters at runtime
  PetscCall(EPSSetFromOptions(eps));

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //                  Solve the eigensystem
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  PetscCall(EPSSolve(eps));
  // Optional: Get some information from the solver and display it
  PetscCall(EPSGetIterationNumber(eps,&its));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %" PetscInt_FMT "\n",its));
  PetscCall(EPSGetType(eps,&type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type));
  PetscCall(EPSGetDimensions(eps,&nev,NULL,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %" PetscInt_FMT "\n",nev));
  PetscCall(EPSGetTolerances(eps,&tol,&maxit));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%" PetscInt_FMT "\n",(double)tol,maxit));

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //                Display solution and clean up
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  
  // Get number of converged approximate eigenpairs
  PetscCall(EPSGetConverged(eps,&nconv));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %" PetscInt_FMT "\n\n",nconv));


  std::complex<double> res(0.,0.);
  if (nconv>0) {
    // Display eigenvalues and relative errors
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
         "           k          ||Ax-kx||/||kx||\n"
         "   ----------------- ------------------\n"));

    for (i=0;i<nconv;i++) {
     // Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and
     // ki (imaginary part)
      PetscCall(EPSGetEigenpair(eps,i,&kr,&ki,xr,xi));
      // Compute the relative error associated to each eigenpair
      PetscCall(EPSComputeError(eps,i,EPS_ERROR_RELATIVE,&error));

#if defined(PETSC_USE_COMPLEX)
      re = PetscRealPart(kr);
      im = PetscImaginaryPart(kr);
#else
      re = kr;
      im = ki;
#endif
      if (im!=0.0) PetscCall(PetscPrintf(PETSC_COMM_WORLD," %9f%+9fi %12g\n",(double)re,(double)im,(double)error));
      else PetscCall(PetscPrintf(PETSC_COMM_WORLD,"   %12f       %12g\n",(double)re,(double)error));

      if(i==0){
        res=std::complex<double>(re, im);

        //Collect local eigenvector to global xr->yr 
        VecScatter sf;
        VecScatterCreateToAll(xr, &sf, &yr);
        VecScatterBegin(sf, xr, yr, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(sf, xr, yr, INSERT_VALUES, SCATTER_FORWARD);

        //copy into Eigen structure
        const std::complex<double>* ar;
        PetscCall(VecGetArrayRead(yr, &ar));
        vec=Eigen::Map<const Eigen::VectorXcd>(ar, vec.size());
        PetscCall(VecRestoreArrayRead(yr, &ar));
 
        //Destroy
        VecScatterDestroy(&sf);
        VecDestroy(&yr);
      }
    }
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\n"));
  }

  // Free work space
  PetscCall(EPSDestroy(&eps));
  PetscCall(MatDestroy(&A));
  PetscCall(VecDestroy(&xr));
  PetscCall(VecDestroy(&xi));
  return res;
}


PT_iTEBD_X_Result SLEPC_PT_iTEBD_X(const Eigen::VectorXcd & LambdaB, const Eigen::MatrixXcd &GLG2, int args, char **argv){

  int dim_d=LambdaB.size();
  int dim_i=GLG2.rows()/dim_d;

  Eigen::MatrixXcd R;

  bool use_full_matrices=false;
#ifdef USE_FULL_MATRICES 
  use_full_matrices=true;
#endif

   //Largest eigenvalue -> matrix VR:
  Eigen::MatrixXcd VR=Eigen::MatrixXcd::Identity(dim_d, dim_d)/sqrt((double)dim_d);

  std::complex<double> max_eval=0;
//Old exact:
/*
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solveR(R);
  int imax=0; for(int i=1; i<dim_d*dim_d; i++){if(std::real(solveR.eigenvalues()(i))>std::real(solveR.eigenvalues()(imax))){ imax=i; }}

  for(int d1=0; d1<dim_d; d1++){ 
    for(int d2=0; d2<dim_d; d2++){
      VR(d1,d2)=solveR.eigenvectors()(d1*dim_d+d2,imax);
    }
  } 
  max_eval=solveR.eigenvalues()(imax);
*/

//SLEPC:
  Eigen::VectorXcd vec(dim_d*dim_d);
  if(use_full_matrices){
    R=PT_iTEBD_calc_R(LambdaB, GLG2);
    max_eval=Largest_EV_SLEPC(vec, (PetscVoidFn *)myMatMult, &R, args, argv);
  }else{
    SLEPC_MatMult_Context ctx; ctx.LambdaB=&LambdaB; ctx.GLG2=&GLG2;
    max_eval=Largest_EV_SLEPC(vec, (PetscVoidFn *)SLEPC_PT_iTEBD_R_times, &ctx, args, argv);
  }

//CHECK: Construct R from new context
//  Eigen::MatrixXcd R2;
//  Eigen_Mat_from_PETSC_Shell(R2, dim_d*dim_d, (PetscVoidFn *)SLEPC_PT_iTEBD_R_times, (void*)&ctx);
//  std::cout<<"|R2-R|="<<(R2-R).norm()<<std::endl;


  //Check quality of eigenvector
  if(use_full_matrices){
    Eigen::VectorXcd vec2=R*vec;
    std::complex<double> val2=vec2.norm();
    std::cout<<"|vec2/val2-vec|="<<(vec2/val2-vec).norm()<<" val2="<<val2<<std::endl;
  }

  for(int d1=0; d1<dim_d; d1++){ 
    for(int d2=0; d2<dim_d; d2++){
      VR(d1,d2)=vec(d1*dim_d+d2);
    }
  }
//  std::cout<<"1-std::abs(vec.dot(exactEV))="<<1-std::abs(vec.dot(solveR.eigenvectors().col(imax)))<<std::endl;


  //Decompose VR:
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solve(VR);
  Eigen::VectorXcd sqrt_D(dim_d);
  for(int d=0; d<dim_d; d++){  
    std::complex<double> eval=solve.eigenvalues()(d);
    if(std::abs(eval)<1e-15){
      std::cerr<<"PT_iTEBD_X: not invertible: "<<solve.eigenvalues().transpose()<<std::endl;
      throw DummyException();
    }
    sqrt_D(d)=sqrt(eval);
  }
  return {solve.eigenvectors(), sqrt_D, max_eval};
}

PT_iTEBD_X_Result SLEPC_PT_iTEBD_Y(const Eigen::VectorXcd & LambdaB, const Eigen::MatrixXcd &GLG2, int args, char **argv){

  int dim_d=LambdaB.size();
  int dim_i=GLG2.rows()/dim_d;
  Eigen::MatrixXcd L,LT;

  bool use_full_matrices=false;
#ifdef USE_FULL_MATRICES 
  use_full_matrices=true;
#endif

   //Largest eigenvalue -> matrix VL:
  Eigen::MatrixXcd VL=Eigen::MatrixXcd::Identity(dim_d, dim_d)/sqrt((double)dim_d);

  std::complex<double> max_eval=0;

  if(use_full_matrices){
    L=PT_iTEBD_calc_L(LambdaB, GLG2);
  }

//Old exact:
/*
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solveL(L.transpose());
  int imax=0; for(int i=1; i<dim_d*dim_d; i++){if(std::real(solveL.eigenvalues()(i))>std::real(solveL.eigenvalues()(imax))){ imax=i; }}

  for(int d1=0; d1<dim_d; d1++){ 
    for(int d2=0; d2<dim_d; d2++){
      VL(d1,d2)=solveL.eigenvectors()(d1*dim_d+d2,imax);
    }
  } 
  max_eval=solveL.eigenvalues()(imax);
*/
//SLEPC:
  Eigen::VectorXcd vec(dim_d*dim_d);

  if(use_full_matrices){
    LT =L.transpose();
    max_eval=Largest_EV_SLEPC(vec, (PetscVoidFn *)myMatMult, &LT, args, argv);
  }else{
    SLEPC_MatMult_Context ctx; ctx.LambdaB=&LambdaB; ctx.GLG2=&GLG2;
    max_eval=Largest_EV_SLEPC(vec, (PetscVoidFn *)SLEPC_PT_iTEBD_LT_times, &ctx, args, argv);
  }

//CHECK: Construct LT from new context
//  Eigen::MatrixXcd LT2;
//  Eigen_Mat_from_PETSC_Shell(LT2, dim_d*dim_d, (PetscVoidFn *)SLEPC_PT_iTEBD_LT_times, (void*)&ctx);
//  std::cout<<"|LT2-LT|="<<(LT2-L.transpose()).norm()<<std::endl;

  if(use_full_matrices){
    Eigen::VectorXcd vec2=L.transpose()*vec;
    std::complex<double> val2=vec2.norm();
    std::cout<<"|vec2/val2-vec|="<<(vec2/val2-vec).norm()<<" val2="<<val2<<std::endl;
  }

  for(int d1=0; d1<dim_d; d1++){ 
    for(int d2=0; d2<dim_d; d2++){
      VL(d1,d2)=vec(d1*dim_d+d2);
    }
  }
//  std::cout<<"1-std::abs(vec.dot(exactEV))="<<1-std::abs(vec.dot(solveR.eigenvectors().col(imax)))<<std::endl;

  //Decompose VL:
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solve(VL);
  Eigen::VectorXcd sqrt_D(dim_d);
  for(int d=0; d<dim_d; d++){  
    std::complex<double> eval=solve.eigenvalues()(d);
    if(std::abs(eval)<1e-15){
      std::cerr<<"PT_iTEBD_Y: not invertible: "<<solve.eigenvalues().transpose()<<std::endl;
      throw DummyException();
    }
    sqrt_D(d)=sqrt(eval);
  }
  return {solve.eigenvectors(), sqrt_D, max_eval};
}

void SLEPC_PT_iTEBD_step(Eigen::VectorXcd & LambdaA, Eigen::VectorXcd &LambdaB,
                   MPS_Matrix & GammaA, MPS_Matrix & GammaB, 
                   Eigen::MatrixXcd expS, const TruncatedSVD &trunc,
                   const infinite_normalize_specs & specs,
                   int args, char** argv){
  //see appendix of Orus & Vidal, PRB 78, 155117 (2008)
  //Check dimensions
  if(LambdaA.size()!=GammaA.dim_d2){
    std::cerr<<"LambdaA.size()!=GammaA.dim_d2  ("<<LambdaA.size()<<" vs. "<<GammaA.dim_d2<<")!"<<std::endl;
    throw DummyException();
  }
  if(LambdaB.size()!=GammaB.dim_d2){
    std::cerr<<"LambdaB.size()!=GammaB.dim_d2  ("<<LambdaB.size()<<" vs. "<<GammaB.dim_d2<<")!"<<std::endl;
    throw DummyException();
  }
  if(LambdaA.size()!=GammaB.dim_d1){
    std::cerr<<"LambdaA.size()!=GammaB.dim_d1  ("<<LambdaA.size()<<" vs. "<<GammaB.dim_d1<<")!"<<std::endl;
    throw DummyException();
  }
  if(LambdaB.size()!=GammaA.dim_d1){
    std::cerr<<"LambdaB.size()!=GammaA.dim_d1  ("<<LambdaB.size()<<" vs. "<<GammaA.dim_d1<<")!"<<std::endl;
    throw DummyException();
  }
  int dim_i=expS.rows();
  if(expS.cols()!=dim_i){
    std::cerr<<"expS.cols()!=dim_i  ("<<expS.cols()<<" vs. "<<dim_i<<")!"<<std::endl;
    throw DummyException();
  }
  if(GammaA.dim_i!=dim_i){
    std::cerr<<"GammaA.dim_i!=dim_i  ("<<GammaA.dim_i<<" vs. "<<dim_i<<")!"<<std::endl;
    throw DummyException();
  }
  if(GammaB.dim_i!=dim_i){
    std::cerr<<"GammaB.dim_i!=dim_i  ("<<GammaB.dim_i<<" vs. "<<dim_i<<")!"<<std::endl;
    throw DummyException();
  }

  //Construct fragments:
  //Start with progagated Gamma_A LambdaA GammaB
  Eigen::MatrixXcd GLG2=Eigen::MatrixXcd::Zero(dim_i*GammaA.dim_d1, dim_i*GammaB.dim_d2);
  for(int i1=0; i1<dim_i; i1++){
    for(int i2=0; i2<dim_i; i2++){
      Eigen::Map<Eigen::MatrixXcd, 0, Eigen::OuterStride<> > (&GLG2(i2*GammaA.dim_d1,i1*GammaB.dim_d2), GammaA.dim_d1, GammaB.dim_d2, Eigen::OuterStride<>(dim_i*GammaA.dim_d1)) = expS(i1, i2) * \
GammaA.get_Matrix_d1_d2(i1) * LambdaA.asDiagonal() * GammaB.get_Matrix_d1_d2(i2);
/*
      for(int d1=0; d1<GammaA.dim_d1; d1++){
        for(int d2=0; d2<GammaB.dim_d2; d2++){
          for(int d3=0; d3<GammaA.dim_d2; d3++){
            GLG2(i2*GammaA.dim_d1+d1, i1*GammaB.dim_d2+d2) += \
         expS(i1, i2) * GammaA(i1, d1, d3) * LambdaA(d3) * GammaB(i2,d3,d2);
          }
        }
      } 
*/
    }
  }

  //Normalize: 
  Eigen::MatrixXcd X, Xinv;
  { 
    time_point norm_start=now();
    PT_iTEBD_X_Result Xs=SLEPC_PT_iTEBD_X(LambdaB, GLG2, args, argv);
    std::cout<<"Right canonicalize: runtime="<<time_diff(now()-norm_start)<<"ms"<<std::endl;
    X=std::get<0>(Xs)*std::get<1>(Xs).asDiagonal();
    Xinv=std::get<1>(Xs).cwiseInverse().asDiagonal()*std::get<0>(Xs).adjoint();

  }

  Eigen::MatrixXcd YT, YTinv;
  {
    time_point norm_start=now();
    PT_iTEBD_X_Result Ys=SLEPC_PT_iTEBD_Y(LambdaB, GLG2, args, argv);
    std::cout<<"Left canonicalize: runtime="<<time_diff(now()-norm_start)<<"ms"<<std::endl;
    YT=(std::get<0>(Ys)*std::get<1>(Ys).asDiagonal()).transpose();
    YTinv=(std::get<1>(Ys).cwiseInverse().asDiagonal()*std::get<0>(Ys).adjoint()).transpose();


  }
 
  time_point SVD_start=now();
  TruncatedSVD_Return ret=trunc.compress(YT*LambdaB.asDiagonal()*X);
  std::cout<<"SVD: runtime="<<time_diff(now()-SVD_start)<<"ms"<<std::endl;

  int new_d=ret.sigma.size();
  LambdaB=ret.sigma;

  //Build Sigma matrix
  Eigen::MatrixXcd lVXi=ret.sigma.asDiagonal()*ret.Vdagger*Xinv;
  Eigen::MatrixXcd YTiUl=YTinv*ret.U*ret.sigma.asDiagonal();

  Eigen::MatrixXcd Sigma_tmp=Eigen::MatrixXcd::Zero(dim_i*new_d, dim_i*GammaB.dim_d2);
  for(int i1=0; i1<dim_i; i1++){
    for(int i2=0; i2<dim_i; i2++){
      for(int d1=0; d1<new_d; d1++){
        for(int d=0; d<GammaA.dim_d1; d++){
          for(int d2=0; d2<GammaB.dim_d2; d2++){
           Sigma_tmp(i1*new_d+d1, i2*GammaB.dim_d2+d2) += \
             lVXi(d1, d)*GLG2(i1*GammaA.dim_d1+d, i2*GammaB.dim_d2+d2);
          }
        }
      } 
    }
  }
  Eigen::MatrixXcd Sigma=Eigen::MatrixXcd::Zero(dim_i*new_d, dim_i*new_d);
  for(int i1=0; i1<dim_i; i1++){
    for(int i2=0; i2<dim_i; i2++){
      for(int d1=0; d1<new_d; d1++){
        for(int d=0; d<GammaB.dim_d2; d++){
          for(int d2=0; d2<new_d; d2++){
           Sigma(i1*new_d+d1, i2*new_d+d2) += \
             Sigma_tmp(i1*new_d+d1, i2*GammaB.dim_d2+d)*YTiUl(d, d2);
          }
        }
      } 
    }
  }

  //Split Sigma matrix   
  TruncatedSVD_Return ret2=trunc.compress(Sigma);
  int new_d2=ret2.sigma.size();
  LambdaA=ret2.sigma;

  GammaA.resize(dim_i, new_d, new_d2);
  for(int i=0; i<dim_i; i++){
    for(int d1=0; d1<new_d; d1++){
      for(int d2=0; d2<new_d2; d2++){
        GammaA(i,d1,d2)=ret2.U(i*new_d+d1, d2)/LambdaB(d1);
      }
    }
  }
  GammaB.resize(dim_i, new_d2, new_d);
  for(int i=0; i<dim_i; i++){
    for(int d1=0; d1<new_d2; d1++){
      for(int d2=0; d2<new_d; d2++){
        GammaB(i,d1,d2)=ret2.Vdagger(d1, i*new_d+d2)/LambdaB(d2);
      }
    }
  }

  //Balance weight between LambdaA and LambdaB:
  std::complex<double> LambdaMean=sqrt(LambdaA(0)*LambdaB(0));
  LambdaA*=LambdaMean/LambdaA(0);
  LambdaB*=LambdaMean/LambdaB(0);
}


int main(int args, char** argv){

  bool is_root=false;
  PetscFunctionBeginUser;
  PetscCall(SlepcInitialize(&args,&argv,NULL,NULL));
  {PetscInt rank;
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
  if(rank==0)is_root=true;}

  Parameters param(args, argv, true);
  std::string Gaussian_prefix=param.get_as_string("Gaussian_prefix", "Boson");
  DiagBB diagBB(param,Gaussian_prefix);

  param.complain_if_not_specified("write_PT");  

  double dict_zero=param.get_as_double("dict_zero",0);
  TruncationLayout trunc_layout(param);

  if(!diagBB.is_set_up()){
    std::cerr<<"DiagBB not set up!"<<std::endl;
    throw DummyException();
  }
  int N=diagBB.get_dim();  
  int NL=N*N;

  TimeGrid tgrid(param);
  int n_mem=tgrid.n_mem;
  if(tgrid.n_mem<2){
    std::cerr<<"tgrid.n_mem<2!"<<std::endl;
    throw DummyException();
  }
 
  infinite_normalize_specs specs;
  specs.iter=param.get_as_int("infinite_normalize_iter", 0); // 0: brute-force normalize, >0 power iteration, <0 Arnoldi
  specs.dont_normalize=param.get_as_bool("infinite_normalize_dont",false);
  specs.use_iter=param.get_as_bool("infinite_normalize_use_iter", specs.iter!=0);
  specs.use_Arnoldi=param.get_as_bool("infinite_normalize_Arnoldi", specs.use_iter);
  specs.eps=param.get_as_double("infinite_normalize_eps", 1e-15);

  Eigen::VectorXcd LambdaA(1); LambdaA(0)=1;
  Eigen::VectorXcd LambdaB(1); LambdaB(0)=1;
  MPS_Matrix GammaA(NL+1, 1,1); GammaA.fill(1.);
  MPS_Matrix GammaB(NL+1, 1,1); GammaB.fill(1.);

  for(int n=n_mem-1; n>0; n--){ 
    TruncatedSVD trunc=trunc_layout.get_base_line(n_mem-1-n, n_mem-1);
    trunc.keep=-1.;
    std::cout<<"n="<<n<<"/"<<n_mem<<" thr="<<trunc.threshold<<" dim_A="<<LambdaA.size()<<" dim_B="<<LambdaB.size()<<std::endl;
    Eigen::MatrixXcd expS=Eigen::MatrixXcd::Ones(NL+1,NL+1);
    expS.block(0,0,NL,NL)=diagBB.calculate_expS(n,tgrid.dt);
    if(n%2==0){
      SLEPC_PT_iTEBD_step(LambdaA, LambdaB, GammaA, GammaB, expS, trunc, specs, args, argv);
    }else{
      SLEPC_PT_iTEBD_step(LambdaB, LambdaA, GammaB, GammaA, expS, trunc, specs, args, argv);
    }
  }

  //First time step;
  std::cout<<"n="<<0<<"/"<<n_mem<<" dim_A="<<LambdaA.size()<<" dim_B="<<LambdaB.size()<<std::endl;
  MPS_Matrix f(NL+1, GammaA.dim_d1, GammaB.dim_d2); 
  f.set_zero();
  Eigen::MatrixXcd expS=Eigen::MatrixXcd::Ones(NL+1,NL+1);
  expS.block(0,0,NL,NL)=diagBB.calculate_expS(0,tgrid.dt);
  for(int i=0; i<NL+1; i++){
    f.get_Matrix_d1_d2(i) = expS(i,i) * \
                           GammaA.get_Matrix_d1_d2(i)*LambdaA.asDiagonal()*\
                           GammaB.get_Matrix_d1_d2(i)*LambdaB.asDiagonal();
  }

  Eigen::MatrixXcd frag0 = f.get_Matrix_d1_d2(NL);

  Eigen::VectorXcd vr, vl;
  { Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver(frag0);
//    std::cout<<"Right eigenvalues: "<<solver.eigenvalues().transpose()<<std::endl;
    int dim_d=LambdaB.size();
    std::vector<std::pair<int, std::complex<double> > > srt(dim_d);
    for(int d=0; d<dim_d; d++){
      srt[d].first=d; srt[d].second=solver.eigenvalues()(d);
    }
    std::sort(srt.begin(), srt.end(), [] (const std::pair<int, std::complex<double> > & a, const std::pair<int, std::complex<double> > &b){ return std::abs(a.second)>std::abs(b.second);});
    std::cout<<"Largest eigenvalue: "<<srt[0].second<<" (nr="<<srt[0].first<<")"<<std::endl;
    vr=solver.eigenvectors().col(srt[0].first);
  }

  { Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver(frag0.transpose());
//    std::cout<<"Left eigenvalues: "<<solver.eigenvalues().transpose()<<std::endl;
    int dim_d=LambdaB.size();
    std::vector<std::pair<int, std::complex<double> > > srt(dim_d);
    for(int d=0; d<dim_d; d++){
      srt[d].first=d; srt[d].second=solver.eigenvalues()(d);
    }
    std::sort(srt.begin(), srt.end(), [] (const std::pair<int, std::complex<double> > & a, const std::pair<int, std::complex<double> > &b){ return std::abs(a.second)>std::abs(b.second);});
    std::cout<<"Largest eigenvalue: "<<srt[0].second<<" (nr="<<srt[0].first<<")"<<std::endl;
    vl=solver.eigenvectors().col(srt[0].first).transpose();
  }


  //write:
  std::string write_PT=param.get_as_string("write_PT");
  int buffer_blocksize=param.get_as_int("buffer_blocksize",-1);

 if(is_root){
  if(param.get_as_bool("write_as_PTB",false)){
    std::shared_ptr<ProcessTensorForward> result(new ProcessTensorBuffer());
    ProcessTensorBuffer *PTB = dynamic_cast<ProcessTensorBuffer*>(result.get());
    if(write_PT!=""){
      PTB->set_new_file(write_PT, buffer_blocksize);
    }
    PTB->resize(tgrid.n_tot);
    for(int n=0; n<tgrid.n_tot; n++){
      ProcessTensorElement &e = PTB->get(n);
      e.clear();
      e.accessor.dict.set_default_diag(N);
      e.M.resize(NL,f.dim_d1,f.dim_d2);
      for(int i=0; i<NL; i++){
        e.M.get_Matrix_d1_d2(i)=f.get_Matrix_d1_d2(i);
      }
    } 
    std::complex<double> prod=vl.transpose()*vr;
    vr/=prod;
    std::cout<<"vl.transpose()*vr="<<vl.transpose()*vr<<std::endl;
    std::cout<<"vl.transpose()*frag0*vr="<<vl.transpose()*frag0*vr<<std::endl;
    PTB->get(0).M.inner_multiply_left(vl.transpose());
    PTB->get(tgrid.n_tot-1).M.inner_multiply_right(vr);
    PTB->expand_DiagBB(diagBB, dict_zero);
    PTB->calculate_closures();
    
  }else{ //write as ProcessTensorRepeat
    std::shared_ptr<ProcessTensorForward> result(new ProcessTensorRepeat());
    ProcessTensorRepeat *PTR = dynamic_cast<ProcessTensorRepeat*>(result.get());
    PTR->set_specs(write_PT, buffer_blocksize);

    ProcessTensorBuffer PTB;
    PTB.resize(3);
    for(int n=0; n<3; n++){
      ProcessTensorElement &e = PTB.get(n);
      e.clear();
      e.accessor.dict.set_default_diag(N);
      e.M.resize(NL,f.dim_d1,f.dim_d2);
      for(int i=0; i<NL; i++){
        e.M.get_Matrix_d1_d2(i)=f.get_Matrix_d1_d2(i);
      }
    } 
    std::complex<double> prod=vl.transpose()*vr;
    vr/=prod;
    std::cout<<"vl.transpose()*vr="<<vl.transpose()*vr<<std::endl;
    std::cout<<"vl.transpose()*frag0*vr="<<vl.transpose()*frag0*vr<<std::endl;
    PTB.get(0).M.inner_multiply_left(vl.transpose());
    PTB.get(2).M.inner_multiply_right(vr);
    PTB.expand_DiagBB(diagBB, dict_zero);
    PTB.calculate_closures();
  
    PTR->initial.resize(1);
    PTR->initial.get(0)=PTB.get(0);
    PTR->repeated.resize(1);
    PTR->repeated.get(0)=PTB.get(1);
  }
 }

  PetscCall(SlepcFinalize());
  return 0;
}

