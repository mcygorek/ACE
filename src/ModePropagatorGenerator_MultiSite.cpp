#include "ModePropagatorGenerator_MultiSite.hpp"
#include "ModePropagatorGenerator.hpp"
#include "Parameters.hpp"
#include "Operators.hpp"
#include "Operators_Boson.hpp"
#include "Equilibrium.hpp"
#include "ReadTable.hpp"
#include "Reader.hpp"
#include "LiouvilleTools.hpp"
#include "ReadTemperature.hpp"
#include "ReducedLiouvilleBasis_Boson.hpp"
#include "ReducedLiouvilleBasis_Boson_FB.hpp"

namespace ACE{

  //print initial boson number per mode
  void ModePropagatorGenerator_MultiSite::print_initial_n(const std::string &fname){
    if(fname=="")return;
    std::ofstream ofs(fname.c_str());
    for(int k=0; k<get_N_modes(); k++){
      Eigen::MatrixXcd init=get_bath_init(k);
      double n=0; for(int r=0; r<init.rows(); r++)n+=r*init(r,r).real();
      ofs<<get_E(k)<<" "<<n<<" "<<k<<std::endl;
    }
  }

  EnvironmentOperators ModePropagatorGenerator_MultiSite::get_env_ops(int k) const{
    std::vector<Eigen::MatrixXcd> mats;
    mats.push_back(Eigen::MatrixXcd::Identity(M,M));
    mats.push_back(Operators_Boson::n(M));
//    mats.push_back(Operators_Boson::n(M)*Operators_Boson::n(M));
    return EnvironmentOperators(mats);
  }

  void ModePropagatorGenerator_MultiSite::setup(Parameters &param){
    MPG_Discretization_E_g E_g_omega(param, name());
/**  NOTE:
    We discretize J(omega). Each omega has a +k and a -k associated to it.
    => Double the number of modes!
*/
    if(E_g_omega.E.size()>0 && E_g_omega.E[0]<0.){
      std::cerr<<"Please choose a non-negative value for '"<<add_name("omega_min")<<"' or '"<<add_name("E_min")<<"'"<<std::endl;
      exit(1);
    }
    {
      int Nomega = E_g_omega.N;
      E_g.N = Nomega * 2;

      E_g.E.resize(E_g.N);
      E_g.dE.resize(E_g.N);
      E_g.g.resize(E_g.N);
      for(int i=0; i<Nomega; i++){
        E_g.E[2*i]=E_g.E[2*i+1]=E_g_omega.E[i];
        E_g.dE[2*i]=E_g.dE[2*i+1]=E_g_omega.dE[i];
        E_g.g[2*i]=E_g.g[2*i+1]=E_g_omega.g[i]/sqrt(2.);
      }
    }

    set_N_modes(E_g.N);
    setup_skip(param);

    if(!param.is_specified(add_name("sites"))){
      std::cerr<<"Please specify positions of sites (in travel time units) by '"<<add_name("sites")<<"'!"<<std::endl;
      exit(1);
    }


    {
    std::vector<std::string> sites_row=param.get_row(add_name("sites"),0);
    tpos.clear();
    bool do_center=false;
    for(size_t i=0; i<sites_row.size(); i++){
      double d;
      if(!canReadDouble(sites_row[i], d)){
        if(i==0 && sites_row[i]=="center"){
          do_center=true;
          tpos.push_back(0);
        }else if(i==0 && sites_row[i]=="equal"){
          if(sites_row.size()<3){
            std::cerr<<"Usage: "<<add_name("sites")<<" equal dtau Nsites !"<<std::endl;
            exit(1);
          }
          double dtau=readDouble(sites_row[1],"MultiSite_sites: dtau");
          size_t Nsites=readSizeT(sites_row[2],"MultiSite_sites: Nsites");
          if(Nsites<1){
            std::cerr<<"'MultiSite_site equal dtau Nsites' : Nsites must be larger than 0!"<<std::endl; 
          }
          tpos.resize(Nsites, dtau); 
          tpos[0]=0.;
          do_center=true;
          break;
        }else{
          std::cerr<<"Usage: "<<add_name("sites")<<" equal dtau1 dtau2 ... !"<<std::endl;
          exit(1);
        }
      }else{
        tpos.push_back(d);
      }
    }
    if(tpos.size()<2){
      std::cerr<<add_name("sites")<<": please specify at least 2 site positions!"<<std::endl;
      exit(1);
    }
    for(size_t i=1; i<tpos.size(); i++){
      tpos[i]+=tpos[i-1];
    }
    if(do_center){
      double center=(tpos.back()-tpos[0])/2.;
      for(size_t i=0; i<tpos.size(); i++){
        tpos[i]-=center;
      }
    }
    }
      
    std::cout<<"tpos (absolute):";
    for(size_t i=0; i<tpos.size(); i++)std::cout<<" "<<tpos[i];
    std::cout<<std::endl;
    std::cout<<"tpos (relative): 0";
    for(size_t i=1; i<tpos.size(); i++)std::cout<<" "<<tpos[i]-tpos[i-1];
    std::cout<<std::endl;
   


    M=param.get_as_size_t(add_name("M"), 2);

    std::string print_E_g=param.get_as_string(add_name("print_E_g"));
    if(print_E_g!=""){
      std::ofstream ofs(print_E_g.c_str());
      for(int i=0; i<E_g.N; i++){
        ofs<<E_g.get_E(i)<<" "<<E_g.get_g(i)<<std::endl;
      }
    }
    std::string print_omega_g=param.get_as_string(add_name("print_omega_g"));
    if(print_omega_g!=""){
      std::ofstream ofs(print_omega_g.c_str());
      for(int i=0; i<E_g.N; i++){
        ofs<<E_g.get_E(i)/hbar_in_meV_ps<<" "<<E_g.get_g(i)<<std::endl;
      }
    }
 
 
    gparam.add_from_prefix(add_name("Propagator"),param);


    use_initial_coherent=false; use_initial_thermal=false;
    if(param.is_specified(add_name("initial_coherent"))){
      use_initial_coherent=true;
      initial_coherent.real(param.get_as_double(add_name("initial_coherent")));
      initial_coherent.imag(param.get_as_double(add_name("initial_coherent"),0.,0,1));
    }else if(param.is_specified(add_name("temperature")) || 
             param.is_specified("temperature") ||
             param.is_specified("temperature_unitless") ||
             param.is_specified(add_name("temperature_unitless")) ||
             param.is_specified(add_name("EFermi"))){
     
      use_initial_thermal=true;
      temperature=readTemperature(param,name());

      EFermi=param.get_as_double(add_name("EFermi"), 0);
    }

    print_initial_n(param.get_as_string(add_name("print_initial_n")));
    reduced=param.get_as_size_t(add_name("reduced"),0);
    reduced_fb=param.get_as_size_t(add_name("reduced_fb"),0);

    
    Nintermediate=param.get_as_size_t(add_name("Nintermediate"), 0);
    if(param.is_specified(add_name("apply_InteractionPicture"))){
      use_IP=true;
      H_IP=param.get_as_operator(add_name("apply_InteractionPicture"));
    }else{
      use_IP=false;
    }
  }

  ModePropagatorPtr ModePropagatorGenerator_MultiSite::getModePropagator(int k)const{
    if(k<0||k>=get_N_modes()){
      std::cerr<<"ModePropagatorGenerator_MultiSite: k<0||k>=get_N_modes()!"<<std::endl; 
      exit(1);
    }

    int N_sites=get_N_sites();
    Operators ops(N_sites);

    int sysdim=N_sites;
    Eigen::MatrixXcd HB_diag=otimes(
         Eigen::MatrixXcd::Identity(sysdim, sysdim), Operators_Boson::n(M));


    Eigen::MatrixXcd HB_base=Eigen::MatrixXcd::Zero(sysdim*M,sysdim*M);

    for(int i=0; i<N_sites; i++){
      double phase=phase_sign(k)*tpos[i]*get_E(k)/hbar_in_meV_ps;
      HB_base+=otimes( ops.ketbra(i,i), 
             exp(std::complex<double>(0.,phase))*Operators_Boson::adagger(M) 
             + exp(std::complex<double>(0.,-phase))*Operators_Boson::a(M)  );
    }


    Eigen::MatrixXcd HB = get_E(k)*HB_diag 
                         + hbar_in_meV_ps*get_g(k)*HB_base;


    Parameters kparam=gparam;
    ModePropagatorPtr ptr=std::make_shared<ModePropagator>(sysdim, get_bath_init(k));
    ptr->FreePropagator::setup(kparam);
    ptr->Nintermediate=Nintermediate;
    ptr->env_ops=get_env_ops(k);
 
    if(use_IP){
      TimedepMatrixPtr mp=std::make_shared<TimedepMatrix_InteractionPicture>(
                        otimes(H_IP, Eigen::MatrixXcd::Identity(M,M)), HB);
      ptr->add_TimedepMatrix(mp);
    }else{
      ptr->add_Hamiltonian(HB);
    }
 
    if(reduced_fb>0){
      ptr->rBasis=std::make_shared<ReducedLiouvilleBasis_Boson_FB>(get_bath_init(k),reduced_fb);
   
      if(k==0){
        std::cout<<"Using reduced_fb="<<reduced_fb<<" (dim=";
        std::cout<<ptr->rBasis->U.cols()<<"): "<<std::endl;
print_diff_from_ortho(ptr->rBasis->U);
      }
    }else if(reduced>0){
      ptr->rBasis=std::make_shared<ReducedLiouvilleBasis_Boson>(get_bath_init(k),reduced);
   
      if(k==0){
        std::cout<<"Using reduced="<<reduced<<" (dim=";
        std::cout<<ptr->rBasis->U.cols()<<"): "<<std::endl;
print_diff_from_ortho(ptr->rBasis->U);
      }
    }
    return ptr;
  }

  Eigen::MatrixXcd ModePropagatorGenerator_MultiSite::get_bath_init(int k)const{

    if(use_initial_coherent){
      Eigen::VectorXcd v(M);
      std::complex<double> fac=1.;
      double anorm=exp(-abs(initial_coherent)*abs(initial_coherent)/2.);
      double N2=0.;
      for(int n=0; n<v.size(); n++){
        v(n)=anorm*fac;
        N2+=abs(anorm*fac)*abs(anorm*fac);
        fac*=initial_coherent/sqrt(n+1);
      }
      std::cout<<"Initialize with coherent state "<<initial_coherent<<": Loss of norm due to truncation: "<<1.-N2<<std::endl;
      
      //for(int n=0; n<v.size(); n++)v(n)/=sqrt(N2);

      Eigen::MatrixXcd initial(M,M);
      for(int i=0; i<initial.rows(); i++){
        for(int j=0; j<initial.cols(); j++){
          initial(i,j)=std::conj(v(i))*v(j);
        }
      }
      return initial;
    }else if(use_initial_thermal){
      Eigen::MatrixXcd HB_diag=Operators_Boson::n(M);
      HB_diag*=(get_E(k)-EFermi);
      return Boson_Equilibrium(HB_diag, temperature);
    }else{
      return Operators_Boson::vacuum(M);
    }
  }

}//namespace
