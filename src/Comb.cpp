#include "Comb.hpp"
#include "TempFileName.hpp"
#include "PassOn.hpp"
#include "CCTEMPO.hpp"
#include "Simulation_CTEMPO.hpp"


namespace ACE{


std::vector<Eigen::VectorXcd> Comb::get_cap_d2(){
  if(size()<1)return std::vector<Eigen::VectorXcd>();

  std::vector<Eigen::VectorXcd> cap(size());
  cap[0]=Eigen::VectorXcd::Ones(get_ro(0,ForwardPreload).f.dim_d2);
  for(int l=0; l<n_tot-1; l++){
    const CombElement & e = get_ro(l,ForwardPreload);
    if(e.closure_d1.rows()!=e.f.dim_d1){
      std::cerr<<"Comb::get_cap_d2: closure_d1["<<l-1<<"].rows()!=mpo["<<l-1<<"].dim_d1 ("<<e.closure_d1.rows()<<" vs. "<<e.f.dim_d1<<")!"<<std::endl;
      throw DummyException();
    }

    cap[l+1]=Eigen::VectorXcd::Zero(e.f.dim_d4);
    int N=sqrt(e.f.dim_d3);
    for(int d1=0; d1<e.f.dim_d1; d1++){
      for(int d2=0; d2<e.f.dim_d2; d2++){
        for(int i=0; i<N; i++){
          for(int d4=0; d4<e.f.dim_d4; d4++){
            cap[l+1](d4)+=e.closure_d1(d1)*cap[l](d2)*e.f(d1, d2, i*N+i, d4);
          }
        }
      }
    }
  }
  return cap;
}
std::vector<Eigen::VectorXcd> Comb::get_cap_d4(){
  if(size()<1)return std::vector<Eigen::VectorXcd>();

  std::vector<Eigen::VectorXcd> cap(size());
  cap.back()=Eigen::VectorXcd::Ones(get_ro(size()-1,BackwardPreload).f.dim_d4);
  for(int l=(int)size()-1; l>0; l--){
    const CombElement & e = get_ro(l,BackwardPreload);
    if(e.closure_d1.rows()!=e.f.dim_d1){
      std::cerr<<"Comb::get_cap_d4: closure_d1[l+1].rows()!=mpo[l+1].dim_d1!"<<std::endl;
      throw DummyException();
    }
    cap[l-1]=Eigen::VectorXcd::Zero(e.f.dim_d2);
    int N=sqrt(e.f.dim_d3);
    for(int d1=0; d1<e.f.dim_d1; d1++){
      for(int d2=0; d2<e.f.dim_d2; d2++){
        for(int i=0; i<N; i++){
          for(int d4=0; d4<e.f.dim_d4; d4++){
            cap[l-1](d2)+=e.closure_d1(d1)*cap[l](d4)*e.f(d1, d2, i*N+i, d4);
          }
        }
      }
    }
  }
  return cap;
}
std::complex<double> Comb::evaluate( const std::vector<int> &list_d3){
  if(list_d3.size()!=size()){
    std::cerr<<"Comb::evaluate: list_d3.size()!=size()!"<<std::endl;
    throw DummyException();
  }
  if(size()<1)return 0;
  if(get_ro(0,ForwardPreload).f.dim_d2!=1){
    std::cerr<<"MPO_op::evaluate: get(0).dim_d2!=1!"<<std::endl;
    throw DummyException();
  }
  Eigen::VectorXcd cap=Eigen::VectorXcd::Ones(1);
  for(int l=0; l<size(); l++){
    const CombElement &e = get_ro(l,ForwardPreload);
    if(e.f.dim_d2!=cap.rows()){
      std::cerr<<"Comb::evaluate: mpo["<<l<<"].dim_d2!=cap.rows()!"<<std::endl;
      throw DummyException();
    }
    if(e.closure_d1.rows()!=e.f.dim_d1){
      std::cerr<<"Comb::evaluate: closure_d1["<<l<<"].rows()!=mpo["<<l<<"].dim_d1 ("<<e.closure_d1.rows()<<" vs. "<<e.f.dim_d1<<")!"<<std::endl;
      throw DummyException();
    }
    if(list_d3[l]>=e.f.dim_d3){
      std::cerr<<"Comb::evaluate: list_d3["<<l<<"]>=mpo["<<l<<"].dim_d3!"<<std::endl;
      throw DummyException();
    }
    Eigen::VectorXcd cap2=Eigen::VectorXcd::Zero(e.f.dim_d4);
    for(int d4=0; d4<e.f.dim_d4; d4++){
      for(int d2=0; d2<e.f.dim_d2; d2++){
        for(int d1=0; d1<e.f.dim_d1; d1++){
          cap2(d4)+=e.closure_d1(d1)*cap(d2)*e.f(d1, d2, list_d3[l], d4);
        }
      }
    }
    cap2.swap(cap);
    if(get_ro(size()-1).f.dim_d4!=1){
      std::cerr<<"Comb::evaluate: get_ro(size()-1).f.dim_d4!=1!"<<std::endl;
      throw DummyException();
    }
  }
  return cap(0);
}

void Comb::join_d3_d1_compress_d2_d4(const MPO & other, const TruncatedSVD &trunc, int verbosity){
  if(size()!=other.size()){
    std::cerr<<"Comb::join_d3_d1_compress_d2_d4: size()!=other.size()!"<<std::endl;
    throw DummyException();
  }
  if(size()<1){return;}

  Eigen::MatrixXcd P;
  {int tmp_dim=get(0,ForwardPreload).f.dim_d2;
   P=Eigen::MatrixXcd::Identity(tmp_dim*other[0].dim_d2,tmp_dim*other[0].dim_d2);
  }
  for(int l=0; l<(int)size()-1; l++){
    CombElement &e=get(l,ForwardPreload);
    
    e.f.join_d3_d1_multiply_d2(other[l], P);
    Eigen::MatrixXcd A = e.f.get_Matrix_d1d2d3_d4();
    TruncatedSVD_Return ret = trunc.compress(A);
//    std::cout<<"["<<l<<"]: singular values: "<<ret.sigma.transpose()<<std::endl;
    double keep=ret.sigma(0);
    if(trunc.keep!=0.)keep=trunc.keep;
#ifdef PRINT_KEEP
std::cout<<"keep: "<<keep<<std::endl;
#endif
    e.f.set_from_Matrix_d1d2d3_d4( ret.U*keep, e.f.dim_d2, e.f.dim_d3);
    P=((ret.sigma/keep).asDiagonal()*ret.Vdagger).transpose();
    if(verbosity>0){std::cout<<"forward sweep of site "<<l<<": "<<A.rows()<<","<<A.cols()<<" -> "<<P.cols()<<std::endl;}
  }
  get(size()-1).f.join_d3_d1_multiply_d2(other.back(), P);
}



void Comb::sweep_d4_d2(const TruncatedSVD &trunc, int verbosity){
  if(n_tot<2)return;
  Eigen::MatrixXcd P;
  for(int l=(int)n_tot-1; l>=0; l--){
    if(verbosity>0){std::cout<<"site "<<l<<" ";}
    get(l,BackwardPreload).sweep_d4_d2(trunc, P, (l==0), verbosity);
  }
}

void Comb::set_teeth(const CCTEMPO_IF_Nodes &nodes){
  if(nodes.size()!=n_tot){
    std::cerr<<"Comb::set_teeth: nodes.size()!=n_tot!"<<std::endl;
    throw DummyException();
  }
  for(int r=0; r<n_tot; r++){
    CombElement & e = get(r, ForwardPreload);
    if(nodes[r].PT && nodes[r].PT->size()>0){
      e.tooth.reset();
    }else if(nodes[r].b.size()<1){ 
      e.tooth.reset(); 
    }else{
      e.tooth = std::make_unique<ProcessTensorBuffer>();

      ProcessTensorElement ref;
      int Nreduced=sqrt(nodes[r].b[0].rows());
      int Nexpand=sqrt(e.f.dim_d3);
      ref.set_trivial(Nreduced);
      ref.accessor.dict.set_default_diag(Nreduced);
      ref.accessor.dict.expand_DiagBB(nodes[r].diagBB);
      ref.M.set_from_Matrix_d1i_d2(Eigen::VectorXcd::Ones(Nreduced*Nreduced), Nreduced*Nreduced);
      ref.env_ops.clear();
      e.tooth->resize(nodes[r].b.size(), ref);
      e.tooth->get(0).accessor.dict.set_default_diag(Nexpand);
      e.tooth->get(0).M.set_from_Matrix_d1i_d2(Eigen::VectorXcd::Zero(Nexpand*Nexpand), Nexpand*Nexpand);
    }
  }
}
/*
void Comb::set_teeth(Parameters & param){
  TimeGrid tgrid(param);
  std::vector<int> dims=CCTEMPO::get_dimensions(param);
  for(int r=0; r<n_tot; r++){
    CombElement & e = get(r, ForwardPreload);
    int N=sqrt(e.f.dim_d3);
    std::string prefix="S"+int_to_string(r)+"_Boson";
    Parameters param2; param2.add_from_prefix(prefix, param);

    if(param2.map.size()<1 || param.get_as_bool("S"+int_to_string(r)+"_use_PT")){
//      e.tooth = std::make_unique<ProcessTensorBuffer>();
//      e.tooth->resize(1);
//      e.tooth->get(0).set_trivial(N); 
//      e.tooth->get(0).accessor.dict.set_default_diag(N); 
//      e.tooth->get(0).env_ops.clear();
      continue;
    }
   
    e.tooth = std::make_unique<ProcessTensorBuffer>();
//std::cout<<"tgrid.n_mem="<<tgrid.n_mem<<std::endl;
    DiagBB diagBB; diagBB.setup(param, prefix);
    Eigen::VectorXcd rho = Eigen::VectorXcd::Zero(N*N); 
    Simulation_CTEMPO::initialize_PTB(*e.tooth, tgrid.n_mem, diagBB, rho);
    e.b=Simulation_CTEMPO::initialize_b(tgrid.n_mem, diagBB, tgrid.dt);
  }

  //PT-MPO:
  for(int r=0; r<n_tot; r++){
    std::string prefix="S"+int_to_string(r);
    Parameters paramS; paramS.add_from_prefix(prefix, param);
    //if Sx_Boson is set but not Sx_use_PT, CTEMPO is used instead of PTs
    if(!paramS.get_as_bool("use_PT")){ 
      paramS.erase_with_prefix("Boson");
    } 
    //Need: copy TimeGrid info:
    std::vector<std::string> keys={"ta","te","dt","t_mem","n_mem","threshold"};
    for(size_t i=0; i<keys.size(); i++){
      if(param.is_specified(keys[i])){
        paramS.add_if_not_specified(keys[i],param.get_as_string(keys[i]));
      }
    }
    std::cout<<"Parameters for PT["<<r<<"]:"<<std::endl; paramS.print();
    CombElement & e = get(r);
    e.PT=std::make_shared<ProcessTensorForwardList>(paramS, dims[r]);
  }

  //Process ...env_same_as_S...
  for(int r=0; r<n_tot; r++){
    std::string key="S"+int_to_string(r)+"_env_sameas_S";
    if(param.is_specified(key)){
      int r2=param.get_as_size_t(key);
      if(r2<0||r2>=size()){
        std::cerr<<key<<": r2="<<r2<<" out of bounds!"<<std::endl;
        throw DummyException();
      }
      std::unique_ptr<ProcessTensorBuffer> tooth_tmp;
      std::vector<Eigen::MatrixXcd> b_tmp;
      { const CombElement & e = get_ro(r2);
        tooth_tmp = std::make_unique<ProcessTensorBuffer>();
        tooth_tmp->copy_content(*e.tooth);
        b_tmp = e.b; }
      { CombElement & e = get(r);
        e.tooth = std::make_unique<ProcessTensorBuffer>();
        e.tooth->copy_content(*tooth_tmp);
        e.b = b_tmp;}
    }
  }
}
*/
void Comb::set_initial_state(const MPO & mpo){
  resize(mpo.size());
  for(int l=0; l<n_tot; l++){
    CombElement & e = get(l, ForwardPreload);
    e.f = mpo.a[l];
    e.closure_d1=Eigen::VectorXcd::Ones(e.f.dim_d1);
  }
}

void Comb::print_dims(std::ostream &os){
  os<<"comb.size()="<<size()<<": ";
  for(int l=0; l<(int)size(); l++){
    get_ro(l,ForwardPreload).f.print_dims(os); os<<"; ";
  }
  os<<std::endl;
}
 
void Comb::setup(const MPO & initial, const CCTEMPO_IF_Nodes &nodes, int buffer_blocksize_, std::string buffer_filename_){
  magicString="COMB0";
  if(buffer_blocksize_>0){
    if(buffer_filename_==""){
      TempFileName tmpname; buffer_filename_=tmpname; tmpname.fname="";
    }
    initialize(buffer_filename_, buffer_blocksize_);
    on_exit=ON_EXIT::DeleteOnDestruction;
  }

  set_initial_state( initial );
  set_teeth( nodes );
}



}//namespace
