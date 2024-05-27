#include "InfluenceFunctional_OD.hpp"
#include "BlockCombine_Basics.hpp"

namespace ACE{

void BlockCombine_FixedOrder(InfluenceFunctional_OD & IF1, const InfluenceFunctional_OD & IF2, double threshold, bool print_timesteps){

  if(IF1.a.size() != IF2.a.size()){
    std::cerr<<"BlockCombine: IF1.a.size() != IF2.a.size()!"<<std::endl;
    exit(1);
  }
  if(IF1.a.size() < 1 ){
    return;
  }
  if(IF1.dict.get_NL()!=IF2.dict.get_NL()){
    std::cerr<<"BlockCombine: IF1.dict.get_NL()!=IF2.dict.get_NL()!"<<std::endl;
    exit(1);
  }
  int NL=IF1.dict.get_NL();

  IF_OD_Dictionary newdict(IF1.dict);
  newdict.join(IF2.dict);
  std::vector<std::vector<int> > newrev=newdict.get_reverse_beta();

 
  std::vector<std::pair<int, int>> k_list, k_list_last;
  Eigen::MatrixXcd R1 = Eigen::MatrixXcd::Identity(IF1.a[0].dim_d1, IF1.a[0].dim_d1);
  Eigen::MatrixXcd R2 = Eigen::MatrixXcd::Identity(IF2.a[0].dim_d1, IF2.a[0].dim_d1);

  for(int k1=0; k1<IF1.a[0].dim_d1; k1++){
    for(int k2=0; k2<IF2.a[0].dim_d1; k2++){
      k_list_last.push_back(std::make_pair(k1,k2));
    }
  }

  for(int n=0; n<(int)IF1.a.size(); n++){ 
    if(n>0){
      if(print_timesteps){
        if(n==1){std::cout<<"Sweep forward: "<<n<<std::flush;}
        else{
          std::stringstream ss_last; ss_last<<n-1;
          for(int i=ss_last.str().length(); i>0; i--)std::cout<<'\b';
          std::cout<<n<<std::flush;
        }
        if(n==IF1.a.size()-1)std::cout<<std::endl;
      }
    }

    MPS_Matrix M1 = IF1.a[n];
    MPS_Matrix M2 = IF2.a[n];
    M1.inner_multiply_left(R1);
    M2.inner_multiply_left(R2);
    Eigen::MatrixXcd A1 = M1.get_Matrix_d1i_d2();
    Eigen::MatrixXcd A2 = M2.get_Matrix_d1i_d2();
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd1(A1, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd2(A2, Eigen::ComputeThinU | Eigen::ComputeThinV);
    
    double svd_max=svd1.singularValues()(0)*svd2.singularValues()(0);
    k_list.clear();
    for(int k1=0; k1<svd1.singularValues().size(); k1++){
      for(int k2=0; k2<svd2.singularValues().size(); k2++){
        if(svd1.singularValues()(k1)*svd2.singularValues()(k2)>=threshold*svd_max){
          k_list.push_back(std::make_pair(k1,k2));
        }else{
          break;
        }
      }
    }
    
    int dim_i_orig=IF1.a[n].dim_i;
    IF1.a[n].resize(newrev.size(), k_list_last.size(), k_list.size());
    IF1.a[n].set_zero();
    for(int d=0; d<(int)k_list.size(); d++){
      for(int i=0; i<NL; i++){
        for(int j=0; j<NL; j++){
          int i_ind=IF1.dict.beta[i*NL+j];         if(i_ind<0)continue;

          for(int k=0; k<NL; k++){
            int i_ind_new=newdict.beta[i*NL+k];  if(i_ind_new<0)continue;
            if(newrev[i_ind_new][0]!=i*NL+k)continue;

            int i_ind2=IF2.dict.beta[j*NL+k];  if(i_ind2<0)continue;
              
            for(int d_=0; d_<(int)k_list_last.size(); d_++){
              IF1.a[n](i_ind_new, d_, d) += 
  svd1.matrixU()(k_list_last[d_].first*dim_i_orig+i_ind, k_list[d].first) 
* svd2.matrixU()(k_list_last[d_].second*IF2.a[n].dim_i+i_ind2, k_list[d].second);
            }
          }
        } 
      }
    }
  
    k_list_last=k_list;
    R1 = svd1.singularValues().head(svd1.matrixV().cols()).asDiagonal()
       * svd1.matrixV().adjoint();
    R2 = svd2.singularValues().head(svd2.matrixV().cols()).asDiagonal()
       * svd2.matrixV().adjoint();

  }

  {//last: close 
    Eigen::MatrixXcd R=Eigen::MatrixXcd::Zero(k_list_last.size(), R1.cols()*R2.cols() );
  
    for(int k=0; k<(int)k_list.size(); k++){ 
      for(int c1=0; c1<R1.cols(); c1++){
        for(int c2=0; c2<R2.cols(); c2++){
          R(k, c1*R2.cols()+c2) = R1(k_list[k].first, c1) 
                                * R2(k_list[k].second, c2);
        }
      }
    }
    IF1.a.back().inner_multiply_right(R);
  }
  
  IF1.dict=newdict;
  IF1.env_ops.clear();
  IF1.calculate_closures();
}

void BlockCombine_AlternateOrder(InfluenceFunctional_OD & IF1, const InfluenceFunctional_OD & IF2, double threshold, bool print_timesteps){

  if(IF1.a.size() != IF2.a.size()){
    std::cerr<<"BlockCombine: IF1.a.size() != IF2.a.size()!"<<std::endl;
    exit(1);
  }
  if(IF1.a.size() < 1 ){
    return;
  }
  if(IF1.dict.get_NL()!=IF2.dict.get_NL()){
    std::cerr<<"BlockCombine: IF1.dict.get_NL()!=IF2.dict.get_NL()!"<<std::endl;
    exit(1);
  }
  int NL=IF1.dict.get_NL();

  IF_OD_Dictionary newdict(IF1.dict);
  newdict.join(IF2.dict);
  std::vector<std::vector<int> > newrev=newdict.get_reverse_beta();

 
  std::vector<std::pair<int, int>> k_list, k_list_last;
  Eigen::MatrixXcd R1 = Eigen::MatrixXcd::Identity(IF1.a[0].dim_d1, IF1.a[0].dim_d1);
  Eigen::MatrixXcd R2 = Eigen::MatrixXcd::Identity(IF2.a[0].dim_d1, IF2.a[0].dim_d1);

  for(int k1=0; k1<IF1.a[0].dim_d1; k1++){
    for(int k2=0; k2<IF2.a[0].dim_d1; k2++){
      k_list_last.push_back(std::make_pair(k1,k2));
    }
  }

  for(int n=0; n<(int)IF1.a.size(); n++){ 
    if(n>0){
      if(print_timesteps){
        if(n==1){std::cout<<"Sweep forward: "<<n<<std::flush;}
        else{
          std::stringstream ss_last; ss_last<<n-1;
          for(int i=ss_last.str().length(); i>0; i--)std::cout<<'\b';
          std::cout<<n<<std::flush;
        }
        if(n==IF1.a.size()-1)std::cout<<std::endl;
      }
    }

    MPS_Matrix M1 = IF1.a[n];
    MPS_Matrix M2 = IF2.a[n];
    M1.inner_multiply_left(R1);
    M2.inner_multiply_left(R2);
    Eigen::MatrixXcd A1 = M1.get_Matrix_d1i_d2();
    Eigen::MatrixXcd A2 = M2.get_Matrix_d1i_d2();
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd1(A1, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd2(A2, Eigen::ComputeThinU | Eigen::ComputeThinV);
    
    double svd_max=svd1.singularValues()(0)*svd2.singularValues()(0);
    k_list.clear();
    for(int k1=0; k1<svd1.singularValues().size(); k1++){
      for(int k2=0; k2<svd2.singularValues().size(); k2++){
        if(svd1.singularValues()(k1)*svd2.singularValues()(k2)>=threshold*svd_max){
          k_list.push_back(std::make_pair(k1,k2));
        }else{
          break;
        }
      }
    }
    
    int dim_i_orig=IF1.a[n].dim_i;
    IF1.a[n].resize(newrev.size(), k_list_last.size(), k_list.size());
    IF1.a[n].set_zero();

    if(n%2==0){
      for(int d=0; d<(int)k_list.size(); d++){
        for(int i=0; i<NL; i++){
          for(int j=0; j<NL; j++){
            int i_ind=IF1.dict.beta[i*NL+j];         if(i_ind<0)continue;

            for(int k=0; k<NL; k++){
              int i_ind_new=newdict.beta[i*NL+k];  if(i_ind_new<0)continue;
              if(newrev[i_ind_new][0]!=i*NL+k)continue;

              int i_ind2=IF2.dict.beta[j*NL+k];  if(i_ind2<0)continue;
               
              for(int d_=0; d_<(int)k_list_last.size(); d_++){
                IF1.a[n](i_ind_new, d_, d) += 
  svd1.matrixU()(k_list_last[d_].first*dim_i_orig+i_ind, k_list[d].first) 
* svd2.matrixU()(k_list_last[d_].second*IF2.a[n].dim_i+i_ind2, k_list[d].second);
              }
            }
          }
        }
      }
    }else{
      for(int d=0; d<(int)k_list.size(); d++){
        for(int j=0; j<NL; j++){
          for(int k=0; k<NL; k++){
            int i_ind=IF1.dict.beta[j*NL+k];         if(i_ind<0)continue;

            for(int i=0; i<NL; i++){
              int i_ind_new=newdict.beta[i*NL+k];  if(i_ind_new<0)continue;
              if(newrev[i_ind_new][0]!=i*NL+k)continue;

              int i_ind2=IF2.dict.beta[i*NL+j];  if(i_ind2<0)continue;
               
              for(int d_=0; d_<(int)k_list_last.size(); d_++){
                IF1.a[n](i_ind_new, d_, d) += 
  svd1.matrixU()(k_list_last[d_].first*dim_i_orig+i_ind, k_list[d].first) 
* svd2.matrixU()(k_list_last[d_].second*IF2.a[n].dim_i+i_ind2, k_list[d].second);
              }
            }
          }
        } 
      }
    }
  
    k_list_last=k_list;
    R1 = svd1.singularValues().head(svd1.matrixV().cols()).asDiagonal()
       * svd1.matrixV().adjoint();
    R2 = svd2.singularValues().head(svd2.matrixV().cols()).asDiagonal()
       * svd2.matrixV().adjoint();

  }

  {//last: close 
    Eigen::MatrixXcd R=Eigen::MatrixXcd::Zero(k_list_last.size(), R1.cols()*R2.cols() );
  
    for(int k=0; k<(int)k_list.size(); k++){ 
      for(int c1=0; c1<R1.cols(); c1++){
        for(int c2=0; c2<R2.cols(); c2++){
          R(k, c1*R2.cols()+c2) = R1(k_list[k].first, c1) 
                                * R2(k_list[k].second, c2);
        }
      }
    }
    IF1.a.back().inner_multiply_right(R);
  }  
  IF1.dict=newdict;
  IF1.env_ops.clear();
  IF1.calculate_closures();
}

void BlockCombineCompress_NF(
  InfluenceFunctional_OD & IF1, const std::vector<Eigen::VectorXd> &NF1, 
  const InfluenceFunctional_OD & IF2, const std::vector<Eigen::VectorXd> &NF2,
  double threshold, bool print_timesteps, int *maxdim_intermediate){

  if(IF1.a.size() != IF2.a.size()){
    std::cerr<<"BlockCombine: IF1.a.size() != IF2.a.size()!"<<std::endl;
    exit(1);
  }
  if(IF1.a.size() != NF1.size()){
    std::cerr<<"BlockCombine: IF1.a.size() != NF1.size()!"<<std::endl;
    exit(1);
  }
  if(IF2.a.size() != NF2.size()){
    std::cerr<<"BlockCombine: IF2.a.size() != NF2.size()!"<<std::endl;
    exit(1);
  }
  if(IF1.a.size() < 1 ){
    return;
  }
  if(IF1.dict.get_NL()!=IF2.dict.get_NL()){
    std::cerr<<"BlockCombine: IF1.dict.get_NL()!=IF2.dict.get_NL()!"<<std::endl;
    exit(1);
  }
  int NL=IF1.dict.get_NL();

  IF_OD_Dictionary newdict(IF1.dict);
  newdict.join(IF2.dict);
  std::vector<std::vector<int> > newrev=newdict.get_reverse_beta();

 
  std::vector<std::pair<int, int>> k_list, k_list_last;

  for(int k1=0; k1<IF1.a[0].dim_d1; k1++){
    for(int k2=0; k2<IF2.a[0].dim_d1; k2++){
      k_list_last.push_back(std::make_pair(k1,k2));
    }
  }
  Eigen::MatrixXcd R = Eigen::MatrixXcd::Identity(k_list_last.size(), k_list_last.size());

  for(int n=0; n<(int)IF1.a.size(); n++){ 
    if(n>0){
      if(print_timesteps){
        if(n==1){std::cout<<"Sweep forward: "<<n<<std::flush;}
        else{
          std::stringstream ss_last; ss_last<<n-1;
          for(int i=ss_last.str().length(); i>0; i--)std::cout<<'\b';
          std::cout<<n<<std::flush;
        }
        if(n==IF1.a.size()-1)std::cout<<std::endl;
      }
    }

    if(NF1[n].size()!=IF1.a[n].dim_d2){
      std::cerr<<"NF1["<<n<<"].size()!=IF1.a["<<n<<"].dim_d2 ("<<NF1[n].size()<<" vs. "<<IF1.a[n].dim_d2<<")!"<<std::endl;
      exit(1);
    }
    if(NF2[n].size()!=IF2.a[n].dim_d2){
      std::cerr<<"NF2["<<n<<"].size()!=IF2.a["<<n<<"].dim_d2 ("<<NF2[n].size()<<" vs. "<<IF2.a[n].dim_d2<<")!"<<std::endl;
      exit(1);
    }

    double svd_max=NF1[n](0)*NF2[n](0);
    k_list.clear();
    for(int k1=0; k1<NF1[n].size(); k1++){
      for(int k2=0; k2<NF2[n].size(); k2++){
        if(NF1[n](k1)*NF2[n](k2)>=threshold*svd_max){
          k_list.push_back(std::make_pair(k1,k2));
        }else{
          break;
        }
      }
    }
    if(maxdim_intermediate!=NULL){
      if(*maxdim_intermediate<k_list.size())*maxdim_intermediate=k_list.size();
    }

    MPS_Matrix M(newrev.size(), k_list_last.size(), k_list.size());
    M.set_zero();

    if(n%2==0){
      for(int d=0; d<(int)k_list.size(); d++){
        for(int i=0; i<NL; i++){
          for(int j=0; j<NL; j++){
            int i_ind=IF1.dict.beta[i*NL+j];         if(i_ind<0)continue;

            for(int k=0; k<NL; k++){
              int i_ind_new=newdict.beta[i*NL+k];  if(i_ind_new<0)continue;
              if(newrev[i_ind_new][0]!=i*NL+k)continue;

              int i_ind2=IF2.dict.beta[j*NL+k];  if(i_ind2<0)continue;
               
              for(int d_=0; d_<(int)k_list_last.size(); d_++){
                M(i_ind_new, d_, d) += 
                   IF1.a[n](i_ind, k_list_last[d_].first, k_list[d].first)
                 * IF2.a[n](i_ind2, k_list_last[d_].second, k_list[d].second);
              }
            }
          }
        }
      }
    }else{
      for(int d=0; d<(int)k_list.size(); d++){
        for(int j=0; j<NL; j++){
          for(int k=0; k<NL; k++){
            int i_ind=IF1.dict.beta[j*NL+k];         if(i_ind<0)continue;

            for(int i=0; i<NL; i++){
              int i_ind_new=newdict.beta[i*NL+k];  if(i_ind_new<0)continue;
              if(newrev[i_ind_new][0]!=i*NL+k)continue;

              int i_ind2=IF2.dict.beta[i*NL+j];  if(i_ind2<0)continue;
               
              for(int d_=0; d_<(int)k_list_last.size(); d_++){
                M(i_ind_new, d_, d) += 
                   IF1.a[n](i_ind, k_list_last[d_].first, k_list[d].first)
                 * IF2.a[n](i_ind2, k_list_last[d_].second, k_list[d].second);
              }
            }
          }
        } 
      }
    }

    
    if(n<IF1.a.size()-1){//Compress:
      M.inner_multiply_left(R);
      Eigen::MatrixXcd A = M.get_Matrix_d1i_d2();
      Eigen::JacobiSVD<Eigen::MatrixXcd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
      int maxdim = svd.singularValues().size();
      for(int i=1; i<svd.singularValues().size(); i++){
        if(svd.singularValues()(i)<threshold*svd.singularValues()(0)){
          maxdim=i;
          break;
        }
      }
      IF1.a[n].set_from_Matrix_d1i_d2(svd.matrixU().block(0,0,svd.matrixU().rows(), maxdim), M.dim_i);

      R = svd.singularValues().head(maxdim).asDiagonal()
        * svd.matrixV().block(0,0,svd.matrixV().rows(), maxdim).adjoint();

    }else{ //last step:
      M.inner_multiply_left(R);
      IF1.a[n].swap(M);
    }
    k_list_last=k_list;
  }
  
  IF1.dict=newdict;
  IF1.env_ops.clear();
  IF1.calculate_closures();
  
}

void BlockCombineCompressBack_NF(
  InfluenceFunctional_OD & IF1, const std::vector<Eigen::VectorXd> &NF1, 
  const InfluenceFunctional_OD & IF2, const std::vector<Eigen::VectorXd> &NF2,
  double threshold, bool print_timesteps, int *maxdim_intermediate,
  bool noSVD){
 
  if(maxdim_intermediate!=NULL){*maxdim_intermediate=0;}
  if(IF1.a.size() != IF2.a.size()){
    std::cerr<<"BlockCombine: IF1.a.size() != IF2.a.size()!"<<std::endl;
    exit(1);
  }
  if(IF1.a.size() != NF1.size()){
    std::cerr<<"BlockCombine: IF1.a.size() != NF1.size()!"<<std::endl;
    exit(1);
  }
  if(IF2.a.size() != NF2.size()){
    std::cerr<<"BlockCombine: IF2.a.size() != NF2.size()!"<<std::endl;
    exit(1);
  }
  if(IF1.a.size() < 1 ){
    return;
  }
  if(IF1.dict.get_NL()!=IF2.dict.get_NL()){
    std::cerr<<"BlockCombine: IF1.dict.get_NL()!=IF2.dict.get_NL()!"<<std::endl;
    exit(1);
  }
  int NL=IF1.dict.get_NL();

  IF_OD_Dictionary newdict(IF1.dict);
  newdict.join(IF2.dict);
  std::vector<std::vector<int> > newrev=newdict.get_reverse_beta();

 
  int n_max=IF1.a.size();
  std::vector<std::vector<Eigen::VectorXcd> > new_ops(n_max);

  std::vector<std::pair<int, int>> k_list, k_list_last;
  for(int k1=0; k1<IF1.a[n_max-1].dim_d2; k1++){
    for(int k2=0; k2<IF2.a[n_max-1].dim_d2; k2++){
      k_list_last.push_back(std::make_pair(k1,k2));
    }
  }
  Eigen::MatrixXcd L = Eigen::MatrixXcd::Identity(k_list_last.size(), k_list_last.size());
  Eigen::MatrixXcd Linv = Eigen::MatrixXcd::Identity(k_list_last.size(), k_list_last.size());

  for(int n=n_max-1; n>=0; n--){ 
    if(print_timesteps){
      if(n==n_max-1){std::cout<<"Sweep backward: "<<n<<std::flush;}
      if(n==0){
        std::cout<<std::endl;
      }else{
        std::stringstream ss_last; ss_last<<n+1;
        for(int i=ss_last.str().length(); i>0; i--)std::cout<<"\b \b";
        std::cout<<n<<std::flush;
      }
    }

    if(NF1[n].size()!=IF1.a[n].dim_d2){
      std::cerr<<"NF1["<<n<<"].size()!=IF1.a["<<n<<"].dim_d2 ("<<NF1[n].size()<<" vs. "<<IF1.a[n].dim_d2<<")!"<<std::endl;
      exit(1);
    }
    if(NF2[n].size()!=IF2.a[n].dim_d2){
      std::cerr<<"NF2["<<n<<"].size()!=IF2.a["<<n<<"].dim_d2 ("<<NF2[n].size()<<" vs. "<<IF2.a[n].dim_d2<<")!"<<std::endl;
      exit(1);
    }

    k_list.clear();
    if(n==0){
      for(int k1=0; k1<IF1.a[0].dim_d1; k1++){
        for(int k2=0; k2<IF2.a[0].dim_d1; k2++){
          k_list.push_back(std::make_pair(k1,k2));
        }
      }
    }else{
      double svd_max=NF1[n-1](0)*NF2[n-1](0);
      for(int k1=0; k1<NF1[n-1].size(); k1++){
        for(int k2=0; k2<NF2[n-1].size(); k2++){
          if(NF1[n-1](k1)*NF2[n-1](k2)>=threshold*svd_max){
            k_list.push_back(std::make_pair(k1,k2));
          }else{
            break;
          }
        }
      }
    }
    if(maxdim_intermediate!=NULL){
      if(*maxdim_intermediate<k_list.size())*maxdim_intermediate=k_list.size();
    }

    MPS_Matrix M(newrev.size(), k_list.size(), k_list_last.size());
    M.set_zero();

    if(n%2==0){
      for(int d=0; d<(int)k_list.size(); d++){
        for(int i=0; i<NL; i++){
          for(int j=0; j<NL; j++){
            int i_ind=IF1.dict.beta[i*NL+j];         if(i_ind<0)continue;

            for(int k=0; k<NL; k++){
              int i_ind_new=newdict.beta[i*NL+k];  if(i_ind_new<0)continue;
              if(newrev[i_ind_new][0]!=i*NL+k)continue;

              int i_ind2=IF2.dict.beta[j*NL+k];  if(i_ind2<0)continue;
               
              for(int d_=0; d_<(int)k_list_last.size(); d_++){
                M(i_ind_new, d, d_) += 
                   IF1.a[n](i_ind, k_list[d].first, k_list_last[d_].first)
                 * IF2.a[n](i_ind2, k_list[d].second, k_list_last[d_].second);
              }
            }
          }
        }
      }
    }else{
      for(int d=0; d<(int)k_list.size(); d++){
        for(int j=0; j<NL; j++){
          for(int k=0; k<NL; k++){
            int i_ind=IF1.dict.beta[j*NL+k];         if(i_ind<0)continue;

            for(int i=0; i<NL; i++){
              int i_ind_new=newdict.beta[i*NL+k];  if(i_ind_new<0)continue;
              if(newrev[i_ind_new][0]!=i*NL+k)continue;

              int i_ind2=IF2.dict.beta[i*NL+j];  if(i_ind2<0)continue;
               
              for(int d_=0; d_<(int)k_list_last.size(); d_++){
                M(i_ind_new, d, d_) += 
                   IF1.a[n](i_ind, k_list[d].first, k_list_last[d_].first)
                 * IF2.a[n](i_ind2, k_list[d].second, k_list_last[d_].second);
              }
            }
          }
        } 
      }
    }

    //environment observables:
//std::cout<<"IF1.env_ops.size()="<<IF1.env_ops.size();
//std::cout<<" IF2.env_ops.size()="<<IF2.env_ops.size()<<std::endl;
    if(IF1.env_ops.size()>0 || IF2.env_ops.size()>0){
      if(IF2.env_ops.size()<1 || IF2.env_ops[n].size()<1){
        new_ops[n]=std::vector<Eigen::VectorXcd>(IF1.env_ops[n].size(),Eigen::VectorXcd::Zero(k_list_last.size()));
        for(int o=0; o<IF1.env_ops[n].size(); o++){
          for(int d_=0; d_<(int)k_list_last.size(); d_++){
            new_ops[n][o](d_) += IF1.env_ops[n][o](k_list_last[d_].first)
                               * IF2.c[n](k_list_last[d_].second);
          }
        }
      }else if(IF1.env_ops.size()<1 || IF1.env_ops[n].size()<1){
        new_ops[n]=std::vector<Eigen::VectorXcd>(IF2.env_ops[n].size(),Eigen::VectorXcd::Zero(k_list_last.size()));
        for(int o=0; o<IF2.env_ops[n].size(); o++){
          for(int d_=0; d_<(int)k_list_last.size(); d_++){
            new_ops[n][o](d_) += IF2.env_ops[n][o](k_list_last[d_].second)
                               * IF1.c[n](k_list_last[d_].first);
          }
        }
      }else{
        if(IF1.env_ops[n].size()!=IF2.env_ops[n].size()){
          std::cerr<<"BlockCombineCompressBack_NF: IF1.env_ops[n].size()!=IF2.env_ops[n].size() ("<<IF1.env_ops[n].size()<<" vs. "<<IF2.env_ops[n].size()<<")!"<<std::endl;
          exit(1);
        }
/*   
        new_ops[n]=std::vector<Eigen::VectorXcd>(1);
        new_ops[n][0]=Eigen::VectorXcd::Zero(k_list_last.size());
        for(int d_=0; d_<(int)k_list_last.size(); d_++){
          new_ops[n][0](d_) = IF1.c[n](k_list_last[d_].first)
                            * IF2.c[n](k_list_last[d_].second);
        }
*/
        new_ops[n]=std::vector<Eigen::VectorXcd>(IF1.env_ops[n].size(),Eigen::VectorXcd::Zero(k_list_last.size()));
        for(int d_=0; d_<(int)k_list_last.size(); d_++){
          new_ops[n][0](d_) += IF1.env_ops[n][0](k_list_last[d_].first)
                             * IF2.env_ops[n][0](k_list_last[d_].second);
        }
        for(int o=1; o<new_ops[n].size(); o++){
          for(int d_=0; d_<(int)k_list_last.size(); d_++){
            new_ops[n][o](d_) += IF1.env_ops[n][o](k_list_last[d_].first)
                               * IF2.env_ops[n][0](k_list_last[d_].second);
            new_ops[n][o](d_) += IF2.env_ops[n][o](k_list_last[d_].second)
                               * IF1.env_ops[n][0](k_list_last[d_].first);
          }
        }
      }
    }
    
    if(n>0){//Compress:
      M.inner_multiply_right(L);
      for(int o=0; o<new_ops[n].size(); o++){
//std::cout<<"Linv.rows()="<<Linv.rows()<<" Linv.cols()="<<Linv.cols();
//std::cout<<" new_ops["<<n<<"].["<<o<<"].rows()="<<new_ops[n][o].rows()<<std::endl;
        new_ops[n][o]=Linv*new_ops[n][o];
      }

      if(noSVD){ 
        L=Eigen::MatrixXcd::Identity(M.dim_d1, M.dim_d1);
        Linv=Eigen::MatrixXcd::Identity(M.dim_d1, M.dim_d1);
        IF1.a[n].swap(M);

      }else{
        Eigen::MatrixXcd A = M.get_Matrix_d1_id2();
        Eigen::JacobiSVD<Eigen::MatrixXcd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
        int maxdim = svd.singularValues().size();
        for(int i=1; i<svd.singularValues().size(); i++){
          if(svd.singularValues()(i)<threshold*svd.singularValues()(0)){
            maxdim=i;
            break;
          }
        }
        double sigma0=svd.singularValues()(0);

        IF1.a[n].set_from_Matrix_d1_id2(  sigma0 *
svd.matrixV().block(0,0,svd.matrixV().rows(), maxdim).adjoint(), M.dim_i);

        L = svd.matrixU().block(0,0,svd.matrixU().rows(), maxdim) 
          * ((1./sigma0)*svd.singularValues().head(maxdim).asDiagonal());
  
        Eigen::VectorXd diag_inv=(1./sigma0)*svd.singularValues().head(maxdim);
        for(int r=0; r<diag_inv.rows(); r++){
          diag_inv(r)=1./diag_inv(r);
        }
        Linv = diag_inv.asDiagonal()*(svd.matrixU().block(0,0,svd.matrixU().rows(), maxdim).adjoint());
      }

    }else{ //last step:
      M.inner_multiply_right(L);
      for(int o=0; o<new_ops[n].size(); o++){
        new_ops[n][o]=Linv*new_ops[n][o];
      }

      IF1.a[n].swap(M);
    }
    k_list_last=k_list;
  }
  
  IF1.dict=newdict;
  IF1.calculate_closures();
//  IF1.env_ops.clear();
  IF1.env_ops=new_ops;
  
}


void sweep_high_to_low(InfluenceFunctional_OD &IF, RankCompressor &compressor, bool print_timesteps){

  double keep_weight=0.;

  for(int n=(int)IF.a.size()-1; n>=1; n--){
    if(print_timesteps){
      if(n==(int)IF.a.size()-1){std::cout<<"Sweep backward: "<<n<<std::flush;}
      else{
        std::stringstream ss_last; ss_last<<n+1;
        for(int i=ss_last.str().length(); i>0; i--)std::cout<<"\b \b";
        std::cout<<n<<std::flush;
      }
      if(n==2)std::cout<<std::endl;
    }
    compressor.sweep_block_high_to_low(n, IF, keep_weight, &IF);
  }
}

void sweep_low_to_high(InfluenceFunctional_OD &IF, RankCompressor &compressor, bool print_timesteps){

  double keep_weight=0.;

  for(int n=1; n<(int)IF.a.size(); n++){
    if(print_timesteps){
      if(n==1){std::cout<<"Sweep forward: "<<n<<std::flush;}
      else{
        std::stringstream ss_last; ss_last<<n-1;
        for(int i=ss_last.str().length(); i>0; i--)std::cout<<"\b";
        std::cout<<n<<std::flush;
      }
      if(n==IF.a.size()-1)std::cout<<std::endl;
    }
    compressor.sweep_block_low_to_high(n-1, IF, keep_weight, &IF);
  }
}


std::vector<Eigen::VectorXd> make_forwardNF(InfluenceFunctional_OD &IF, 
     double threshold, bool print_timesteps, int fix_phase){

  std::vector<Eigen::VectorXd> forwardNF(IF.a.size());
  if(IF.a.size()<1)return forwardNF;

  Eigen::MatrixXcd R = Eigen::MatrixXcd::Identity(IF.a[0].dim_d1, IF.a[0].dim_d1);

  for(int n=0; n<(int)IF.a.size()-1; n++){

    IF.a[n].inner_multiply_left(R);

    Eigen::MatrixXcd A = IF.a[n].get_Matrix_d1i_d2();
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
  
    int maxdim = svd.singularValues().size();
    double maxsv=svd.singularValues()(0);
    for(int i=1; i<svd.singularValues().size(); i++){
      if(svd.singularValues()(i)<threshold*maxsv){
        maxdim=i;
        break;
      }
    }
    forwardNF[n]=svd.singularValues().head(maxdim);
    
    IF.a[n].set_from_Matrix_d1i_d2(maxsv*svd.matrixU().block(0,0,svd.matrixU().rows(), maxdim), IF.a[n].dim_i);

    R = (1./maxsv)*forwardNF[n].asDiagonal()
       * svd.matrixV().block(0,0,svd.matrixV().rows(),maxdim).adjoint();

    if(fix_phase==2){
      Eigen::VectorXcd phases=Eigen::VectorXcd::Ones(maxdim);
      for(int c=0; c<maxdim; c++){
        for(int r=0; r<svd.matrixU().rows(); r++){ 
          std::complex<double> val=svd.matrixU()(r,c);
          if(abs(val)>1e-6){
            phases(c)=val/abs(val);
            break;
          }
        }
      }
      Eigen::MatrixXcd P=phases.asDiagonal();
      IF.a[n].inner_multiply_right(P.adjoint());
      R=P*R;
    }else if(fix_phase==1){
      Eigen::VectorXcd phases(maxdim);
      for(int c=0; c<maxdim; c++){
        std::complex<double> max_val=0.;
        for(int r=0; r<svd.matrixU().rows(); r++){ 
          if(abs(svd.matrixU()(r,c))>abs(max_val))max_val=svd.matrixU()(r,c);
        }
        phases(c)=max_val/abs(max_val);
      }
      Eigen::MatrixXcd P=phases.asDiagonal();
      IF.a[n].inner_multiply_right(P.adjoint());
      R=P*R;
    }

//std::cout<<"R.rows()="<<R.rows()<<" R.cols()="<<R.cols()<<" IF.c[n].rows()="<<IF.c[n].rows()<<std::endl;
//std::cout<<"R.rows()="<<R.rows()<<" R.cols()="<<R.cols()<<" IF.c[n].rows()="<<IF.c[n].rows()<<std::endl;
    IF.c[n]=R*IF.c[n];
    if(IF.env_ops.size()>0)for(int o=0; o<IF.env_ops[n].size(); o++){
      IF.env_ops[n][o]=R*IF.env_ops[n][o];
    }
  }
  IF.back().inner_multiply_left(R);

  forwardNF.back()=Eigen::VectorXd(1); forwardNF.back()<<1.;
  
  return forwardNF;
}


void Block_Combine_and_Sweep_Back(InfluenceFunctional_OD &IF,
                                  InfluenceFunctional_OD &IF2, 
                                  double threshold, bool verbose, 
                                  bool IF_print_timesteps){

  RankCompressor_SVD compressor(threshold);
            
  if(verbose){
    std::cout<<"Before compression: IF.get_max_dim()="<<IF.get_max_dim();
    std::cout<<" IF2.get_max_dim()="<<IF2.get_max_dim()<<std::endl;
  }
  BlockCombine_AlternateOrder(IF, IF2, threshold, IF_print_timesteps);
  if(verbose){
    std::cout<<"After combination: IF.get_max_dim()="<<IF.get_max_dim()<<std::endl;
  }
  sweep_high_to_low(IF, compressor, IF_print_timesteps);
  if(verbose){
    std::cout<<"After back sweep: IF.get_max_dim()="<<IF.get_max_dim()<<std::endl;
  }       
}

void Block_Combine_and_Sweep_NF(  InfluenceFunctional_OD &IF,
                                  InfluenceFunctional_OD &IF2, 
                                  double threshold, bool verbose, 
                                  bool IF_print_timesteps, bool noSVD){

  RankCompressor_SVD compressor(threshold);
            
  if(verbose){
    std::cout<<"Before compression: IF.get_max_dim()="<<IF.get_max_dim();
    std::cout<<" IF2.get_max_dim()="<<IF2.get_max_dim()<<std::endl;
  }

  std::vector<Eigen::VectorXd> NF=make_forwardNF(IF, threshold, IF_print_timesteps);
  if(verbose){
    std::cout<<"forwardNF: IF.get_max_dim()="<<IF.get_max_dim()<<std::endl;
  }
  std::vector<Eigen::VectorXd> NF2=make_forwardNF(IF2, threshold, IF_print_timesteps);
  if(verbose){
    std::cout<<"forwardNF: IF2.get_max_dim()="<<IF2.get_max_dim()<<std::endl;
  }

/*
  int maxdim_intermediate;
  BlockCombineCompress_NF(IF, NF, IF2, NF2, threshold, IF_print_timesteps, &maxdim_intermediate);
  if(verbose){
    std::cout<<"After combination: IF.get_max_dim()="<<maxdim_intermediate<<std::endl;
    std::cout<<"After combination & compression: IF.get_max_dim()="<<IF.get_max_dim()<<std::endl;
  }
  sweep_high_to_low(IF, compressor, IF_print_timesteps);
  if(verbose){
    std::cout<<"After back sweep: IF.get_max_dim()="<<IF.get_max_dim()<<std::endl;
  }  
*/ 
    
  int maxdim_intermediate;
  BlockCombineCompressBack_NF(IF, NF, IF2, NF2, threshold, IF_print_timesteps, &maxdim_intermediate, noSVD);

  if(verbose){
    std::cout<<"After combination: IF.get_max_dim()="<<maxdim_intermediate<<std::endl;
    std::cout<<"After combination & compression: IF.get_max_dim()="<<IF.get_max_dim()<<std::endl;
  }

}


}//namespace
