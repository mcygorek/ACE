#include "ACE.hpp"
#include "DynamicalMap.hpp"

using namespace ACE;

int main(int args, char **argv){
  Parameters param(args, argv, true);
 
  bool do_calculate_stepwise=param.get_as_bool("do_calculate_stepwise");
  bool do_calculate_TEMPO=param.get_as_bool("do_calculate_TEMPO");
  bool do_calculate=param.get_as_bool("do_calculate", do_calculate_stepwise|do_calculate_TEMPO);
  std::string dynfile=param.get_as_string("dynfile");
 try{
  DynamicalMap dyn;
  if(do_calculate){
    if(do_calculate_TEMPO){
      std::cout<<"do_calculate_TEMPO is set"<<std::endl;
      dyn.calculate_TEMPO(param);
    }else if(do_calculate_stepwise){
      std::cout<<"do_calculate_stepwise is set"<<std::endl;
      dyn.calculate_stepwise(param);
    }else{
      std::cout<<"do_calculate is set"<<std::endl;
      dyn.calculate(param);
    }
    dyn.write(dynfile);
    std::cout<<"DynamicalMap written to file '"<<dynfile<<"'"<<std::endl;

  }else{
    if(dynfile==""){
      std::cerr<<"Set 'do_calculate' or specify 'dynfile'"<<std::endl;
      exit(1);
    }
    dyn.read(dynfile);
    std::cout<<"DynamicalMap read from file '"<<dynfile<<"'"<<std::endl;
  }

 
  std::cout<<"ta: "<<dyn.ta<<" dt: "<<dyn.dt;
  std::cout<<" E.size()="<<(int)dyn.E.size()<<std::endl;


  std::string Richardson_combine=param.get_as_string("Richardson_combine");
  if(Richardson_combine!=""){
    std::cout<<"Richardson-combination with '"<<Richardson_combine<<"'"<<std::endl;
    DynamicalMap other(Richardson_combine);
    dyn.Richardson_combine(other);
  }

  Parameters param_prop=param;
  param_prop.add_to("dt", dyn.dt);
  param_prop.add_to("te", dyn.ta+dyn.E.size()*dyn.dt);

  int coarse_grain=param.get_as_int("coarse_grain",1);
  if(coarse_grain>1){
    std::cout<<"Coarse graining with factor "<<coarse_grain<<std::endl;
    param_prop.override_param("dt", dyn.dt*coarse_grain);
    dyn.dt*=coarse_grain;
    int n_new=( (param.get_as_double("te")-dyn.ta)/dyn.dt );
    param_prop.override_param("te", n_new*dyn.dt+dyn.ta);

    std::vector<Eigen::MatrixXcd> tmp(dyn.E.size()/coarse_grain);
    for(size_t i=0; i<tmp.size(); i++){
      tmp[i]=dyn.E[i*coarse_grain];
    }
    dyn.E.swap(tmp);
  }

  int smoothen_iter=param.get_as_int("smoothen_iter");
  if(smoothen_iter>0){
    std::cout<<"smoothen_iter="<<smoothen_iter<<std::endl;
    if(dyn.E.size()<2){
      std::cerr<<"dyn.E.size()<2"<<std::endl;
      throw DummyException();
    }
    int DL=dyn.get(0).rows();
    for(int iter=0; iter<smoothen_iter; iter++){
      std::vector<Eigen::MatrixXcd> dE2(dyn.E.size(),Eigen::MatrixXcd::Zero(DL, DL));
      for(int j=0; j<dE2.size()-1; j+=2){
        Eigen::MatrixXcd tmp=((dyn.get_dE(j+1)*dyn.get_dE(j)).log()*0.5).exp();
        dE2[j]=tmp;
        dE2[j+1]=tmp;
      }
      for(int j=0; j<dE2.size(); j++){
        if(j==0)dyn.E[0]=dE2[0];
        else dyn.E[j]=dE2[j]*dyn.E[j-1];
      }
      for(int j=1; j<dE2.size()-1; j+=2){
        Eigen::MatrixXcd tmp=((dyn.get_dE(j+1)*dyn.get_dE(j)).log()*0.5).exp();
        dE2[j]=tmp;
        dE2[j+1]=tmp;
      }
      for(int j=0; j<dE2.size(); j++){
        if(j==0)dyn.E[0]=dE2[0];
        else dyn.E[j]=dE2[j]*dyn.E[j-1];
      }
    }
  }

  if(param.is_specified("extrapolate_TT")){
    if(param.get_nr_cols("extrapolate_TT",0)<2){
      std::cerr<<"extrapolate_TT T_TO T_FROM"<<std::endl;
      throw DummyException();
    }
    double te2=param.get_as_double("extrapolate_TT", 0, 0, 0);
    double t_length=param.get_as_double("extrapolate_TT", 0, 0, 1);
    int length=round( (t_length-dyn.ta)/dyn.dt );
    std::cout<<"Transfer tensor extrapolation to time "<<te2<<" length="<<length<<std::endl;
    
    dyn.extrapolate_TT(te2, length);
    param_prop.override_param("te", te2);
  
  }else if(param.is_specified("extrapolate_convolute")){
    if(param.get_nr_cols("extrapolate_convolute",0)<2){
      std::cerr<<"Parameters of 'extrapolate_convolute': PROP_TO FWHM!"<<std::endl;
      throw DummyException();
    }
    double te2 = param.get_as_double("extrapolate_convolute");
    double FWHM = param.get_as_double("extrapolate_convolute",0,0,1);
    std::cout<<"Extrapolating using Gaussian convolution with FWHM "<<FWHM<<" to time "<<te2<<std::endl;
    double sigma=FWHM/sqrt(8.*log(2.));
    int n_new=( (te2-dyn.ta)/dyn.dt );
    param_prop.override_param("te", te2);

    double t_from=param.get_as_double("extrapolate_convolute",0,0,2);
    int n_from=0;
    if(t_from>0.){
      n_from=round( (t_from-dyn.ta)/dyn.dt );
    }
    if(n_from>(int)dyn.E.size()){
      std::cerr<<"extrapolate_convolute: n_from>E.size()!"<<std::endl;
      throw DummyException();
    }
    if(n_from>0 && n_from<(int)dyn.E.size()){
      dyn.E.resize(n_from);
    }

    int DL=dyn.get(0).rows();
    Eigen::MatrixXcd average_dE=Eigen::MatrixXcd::Zero(DL,DL);
    double norm=0.;
    for(int l=0; l<(int)dyn.E.size(); l++){
      double y=(l-((int)dyn.E.size()-1))*dyn.dt/sigma;
      double factor=(2./(sigma*sqrt(2.*M_PI))*exp(-0.5*y*y))*dyn.dt;
      norm+=factor;
      average_dE+=factor*dyn.get_dE(l);
    }
    std::cout<<"norm="<<norm<<std::endl;
    average_dE/=norm;

    for(int n=dyn.E.size(); n<n_new; n++){
      dyn.E.push_back( average_dE*dyn.E.back() );
    }   
    std::cout<<"dyn.E.size()="<<dyn.E.size()<<std::endl;

  }else if(param.is_specified("extrapolate_average")){
    if(param.get_nr_cols("extrapolate_average",0)<1){
      std::cerr<<"Parameters of 'extrapolate_average': PROP_TO [EXTR_FROM] [EXTR_TO]!"<<std::endl;
      throw DummyException();
    }
    double te2=param.get_as_double("extrapolate_average", 0, 0, 0);
    double tstart=param.get_as_double("extrapolate_average", dyn.ta+dyn.dt*dyn.E.size()/2., 0, 1);
    double tend=param.get_as_double("extrapolate_average", dyn.ta+dyn.dt*dyn.E.size(), 0, 2);
    int n_new=( (te2-dyn.ta)/dyn.dt );
    param_prop.override_param("te", te2);
    int n_from=0;
    if(tend>0.){
      n_from=round( (tend-dyn.ta)/dyn.dt );
    }
    if(n_from>(int)dyn.E.size()){
      std::cerr<<"extrapolate_average: n_from>E.size()!"<<std::endl;
      throw DummyException();
    }
    if(n_from>0 && n_from<(int)dyn.E.size()){
      dyn.E.resize(n_from);
    }

    std::cout<<"Extrapolating average eigenvalues ["<<tstart<<":"<<tend<<"] to time "<<te2<<std::endl;
    Eigen::MatrixXcd average_dE=dyn.get_average_dE(tstart, tend);
    for(int n=dyn.E.size(); n<n_new; n++){
      dyn.E.push_back( average_dE*dyn.E.back() );
    }
  }else if(param.is_specified("extrapolate_last")){
    if(param.get_nr_cols("extrapolate_last",0)<2){
      std::cerr<<"extrapolate_last T_TO T_FROM"<<std::endl;
      throw DummyException();
    }
    double te2=param.get_as_double("extrapolate_last");
    int n_new=round( (te2-dyn.ta)/dyn.dt );
    param_prop.override_param("te", te2);
    double t_from=param.get_as_double("extrapolate_last",0,0,1);
    int n_from=0;
    if(t_from>0.){
      n_from=round( (t_from-dyn.ta)/dyn.dt );
    }
    if(n_from>(int)dyn.E.size()){
      std::cerr<<"extrapolate_last: n_from>E.size()!"<<std::endl;
      throw DummyException();
    }
    if(n_from>0 && n_from<(int)dyn.E.size()){
      dyn.E.resize(n_from);
    }
      
    Eigen::MatrixXcd last_dE=dyn.get_dE(dyn.E.size()-1);
    if(param.get_as_bool("make_last_physical")){
      dyn.make_physical_dE(last_dE);
    }
    for(int n=dyn.E.size(); n<n_new; n++){
      dyn.E.push_back( last_dE*dyn.E.back() );
    }
  }else if(param.is_specified("extrapolate_last_ref")){
    if(param.get_nr_cols("extrapolate_last_ref",0)<4){
      std::cout<<"extrapolate_last_ref: T_TO T_FROM EPS TREF"<<std::endl;
    }
    double te2=param.get_as_double("extrapolate_last_ref");
    int n_new=round( (te2-dyn.ta)/dyn.dt );
    param_prop.override_param("te", te2);
    double t_from=param.get_as_double("extrapolate_last_ref",0,0,1);
    int n_from=0;
    if(t_from>0.){
      n_from=round( (t_from-dyn.ta)/dyn.dt );
    }
    if(n_from>(int)dyn.E.size()){
      std::cerr<<"extrapolate_last_ref: n_from>E.size()!"<<std::endl;
      throw DummyException();
    }
    if(n_from>0 && n_from<(int)dyn.E.size()){
      dyn.E.resize(n_from);
    }
    double eps=param.get_as_double("extrapolate_last_ref",0,0,2);
    double tref=param.get_as_double("extrapolate_last_ref",0,0,3);
      
    Eigen::MatrixXcd last_dE=dyn.get_dE_ref(dyn.E.size()-1, eps, tref);
    if(param.get_as_bool("make_last_physical")){
      dyn.make_physical_dE(last_dE);
    }
    for(int n=dyn.E.size(); n<n_new; n++){
      dyn.E.push_back( last_dE*dyn.E.back() );
    }
  }else if(param.is_specified("extrapolate_Pade")){
    if(param.get_nr_cols("extrapolate_Pade",0)<2){
      std::cerr<<"extrapolate_Pade T_TO T_FROM [n_step] [n]"<<std::endl;
      throw DummyException();
    }
    double te2=param.get_as_double("extrapolate_Pade");
    int n_new=round( (te2-dyn.ta)/dyn.dt );
    param_prop.override_param("te", te2);
    double t_from=param.get_as_double("extrapolate_Pade",0,0,1);
    int n_from=0;
    if(t_from>0.){
      n_from=round( (t_from-dyn.ta)/dyn.dt );
    }
    if(n_from>(int)dyn.E.size()){
      std::cerr<<"extrapolate_Pade: n_from>E.size()!"<<std::endl;
      throw DummyException();
    }
    if(n_from>0 && n_from<(int)dyn.E.size()){
      dyn.E.resize(n_from);
    }
 
    int n_step=param.get_as_size_t("extrapolate_Pade",1,0,2);
    int n=param.get_as_size_t("extrapolate_Pade",3,0,3);
    std::vector<int> steps(n);
    for(int i=0; i<n; i++){
      steps[i]=n_from-n_step*(n-i);
    }
    dyn.calculate_Pade(steps);

    if(dyn.Pade.size()<3){
      std::cerr<<"extrapolate_Pade: Pade.size()<3!"<<std::endl;
      throw DummyException();
    }
    int dim=dyn.Pade[0].rows();
    for(int n=dyn.E.size(); n<n_new; n++){
      Eigen::MatrixXcd num=dyn.Pade[0]+n*dyn.Pade[1];
      Eigen::MatrixXcd denom=Eigen::MatrixXcd::Identity(dim,dim);
      denom+=n*dyn.Pade[2];
      Eigen::JacobiSVD<Eigen::MatrixXcd> svd(denom, Eigen::ComputeFullU | Eigen::ComputeFullV);
      Eigen::MatrixXcd thisdE=svd.solve(num);

      dyn.E.push_back( thisdE*dyn.E.back() );
    }
  }else if(param.is_specified("extrapolate_Pade2")){
    if(param.get_nr_cols("extrapolate_Pade2",0)<2){
      std::cerr<<"extrapolate_Pade2 T_TO T_FROM [n_step] [n]"<<std::endl;
      throw DummyException();
    }
    double te2=param.get_as_double("extrapolate_Pade2");
    int n_new=round( (te2-dyn.ta)/dyn.dt );
    param_prop.override_param("te", te2);
    double t_from=param.get_as_double("extrapolate_Pade2",0,0,1);
    int n_from=0;
    if(t_from>0.){
      n_from=round( (t_from-dyn.ta)/dyn.dt );
    }
    if(n_from>(int)dyn.E.size()){
      std::cerr<<"extrapolate_Pade2: n_from>E.size()!"<<std::endl;
      throw DummyException();
    }
    if(n_from>0 && n_from<(int)dyn.E.size()){
      dyn.E.resize(n_from);
    }
 
    int n_step=param.get_as_size_t("extrapolate_Pade2",1,0,2);
    int n=param.get_as_size_t("extrapolate_Pade2",3,0,3);
    std::vector<int> steps(n);
    for(int i=0; i<n; i++){
      steps[i]=n_from-n_step*(n-i);
    }
    dyn.calculate_Pade2(steps);

    if(dyn.Pade.size()<5){
      std::cerr<<"extrapolate_Pade2: Pade.size()<5!"<<std::endl;
      throw DummyException();
    }
    int dim=dyn.Pade[0].rows();
    for(int n=dyn.E.size(); n<n_new; n++){
      Eigen::MatrixXcd num=dyn.Pade[0]+n*dyn.Pade[1]+n*n*dyn.Pade[2];
      Eigen::MatrixXcd denom=Eigen::MatrixXcd::Identity(dim,dim);
      denom.noalias()+=n*dyn.Pade[3]+n*n*dyn.Pade[4];
      Eigen::JacobiSVD<Eigen::MatrixXcd> svd(denom, Eigen::ComputeFullU | Eigen::ComputeFullV);
      Eigen::MatrixXcd thisdE=svd.solve(num);

      dyn.E.push_back( thisdE*dyn.E.back() );
    }
  }



  if(param.get_as_bool("make_physical")){
    dyn.make_physical();
  }

  if(param.get_as_bool("print_E")){
    for(size_t n=0; n<dyn.E.size(); n++){
      std::cout<<"n="<<n<<":"<<std::endl;
      std::cout<<dyn.get(n)<<std::endl;
    }
  }
 
  std::string print_Lindblad_Hall=param.get_as_string("print_Lindblad_Hall");
  if(print_Lindblad_Hall!=""){
    dyn.print_Lindblad_Hall(print_Lindblad_Hall);
  }
  std::string print_Lindblad=param.get_as_string("print_Lindblad");
  if(print_Lindblad!=""){
    dyn.print_Lindblad(print_Lindblad);
  }

  std::string analyze_physicality=param.get_as_string("analyze_physicality");
  if(analyze_physicality!=""){
    dyn.analyze_physicality(analyze_physicality);
  }

  std::string print_normdiff=param.get_as_string("print_normdiff");
  if(print_normdiff!=""){
    dyn.print_normdiff(print_normdiff);
  }
  std::string print_TT_normdiff=param.get_as_string("print_TT_normdiff");
  if(print_TT_normdiff!=""){
    dyn.print_TT_normdiff(print_TT_normdiff);
  }
  std::string print_TT_norm=param.get_as_string("print_TT_norm");
  if(print_TT_norm!=""){
    dyn.print_TT_norm(print_TT_norm);
  }
  std::string print_dE_eigenvalues=param.get_as_string("print_dE_eigenvalues");
  if(print_dE_eigenvalues!=""){
    dyn.print_eigenvalues(print_dE_eigenvalues);
  }
  std::string print_E_eigenvalues=param.get_as_string("print_E_eigenvalues");
  if(print_E_eigenvalues!=""){
    dyn.print_E_eigenvalues(print_E_eigenvalues);
  }
  std::string print_eigenvalues=param.get_as_string("print_eigenvalues");
  if(print_eigenvalues!=""){
    dyn.print_eigenvalues(print_eigenvalues, param.get_as_double("print_eigenvalues_regularize"));
  }
 
  std::string print_eigenvalues_ref=param.get_as_string("print_eigenvalues_ref");
  if(print_eigenvalues_ref!=""){
    double eps=param.get_as_double("print_eigenvalues_ref",0,0,1);
    double tref=param.get_as_double("print_eigenvalues_ref",0,0,2);
    if(eps<=0.||tref<=0.){
      std::cerr<<"print_eigenvalues_ref EPS TREF: arguments must both be larger than zero!"<<std::endl;
      throw DummyException();
    }
    std::cout<<"print_eigenvalues_ref "<<eps<<" "<<tref<<std::endl;
    dyn.print_eigenvalues_ref(print_eigenvalues_ref, eps, tref);
  } 

  std::string print_dE_eigenvalues_ref=param.get_as_string("print_dE_eigenvalues_ref");
  if(print_dE_eigenvalues_ref!=""){
    double eps=param.get_as_double("print_dE_eigenvalues_ref",0,0,1);
    double tref=param.get_as_double("print_dE_eigenvalues_ref",0,0,2);
    if(eps<=0.||tref<=0.){
      std::cerr<<"print_dE_eigenvalues_ref EPS TREF: arguments must both be larger than zero!"<<std::endl;
      throw DummyException();
    }
    std::cout<<"print_dE_eigenvalues_ref "<<eps<<" "<<tref<<std::endl;
    dyn.print_dE_eigenvalues_ref(print_dE_eigenvalues_ref, eps, tref);
  }  

  std::string print_dE_singularvalues=param.get_as_string("print_dE_singularvalues");
  if(print_dE_singularvalues!=""){
    dyn.print_dE_singularvalues(print_dE_singularvalues);
  }  
  std::string print_singularvalues=param.get_as_string("print_singularvalues");
  if(print_singularvalues!=""){
    dyn.print_singularvalues(print_singularvalues);
  } 
 
  std::string print_fullE_eigenvalues=param.get_as_string("print_fullE_eigenvalues");
  if(print_fullE_eigenvalues!=""){
    dyn.print_fullE_eigenvalues(print_fullE_eigenvalues);
  } 

  std::string print_average_eigenvalues=param.get_as_string("print_average_eigenvalues");
  if(print_average_eigenvalues!=""){
    double tstart=param.get_as_double("print_average_eigenvalues",dyn.ta+dyn.E.size()*dyn.dt/2,0,1);
    double tend=param.get_as_double("print_average_eigenvalues",dyn.ta+dyn.E.size()*dyn.dt,0,2);
   
    std::cout<<"Printing average eigenvalues ["<<tstart<<":"<<tend<<"] to file '"<<print_average_eigenvalues<<"'"<<std::endl;
    dyn.print_average_eigenvalues(print_average_eigenvalues, tstart, tend);
  } 

  bool analyze_SVD=param.get_as_bool("analyze_SVD");
  if(analyze_SVD){
    int offset=param.get_as_int("analyze_SVD_offset");
    Eigen::MatrixXcd E_penult=dyn.get(dyn.E.size()-1-offset);
    Eigen::MatrixXcd E_last=dyn.get(dyn.E.size()-1);
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd_penult(E_penult, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd_last(E_last, Eigen::ComputeFullU | Eigen::ComputeFullV);
    int DL=E_last.rows();

    std::cout<<"Analyze SVDs (penultimate vs. last E):"<<std::endl;
    std::cout<<"U overlaps:"<<std::endl;
    for(int i=0; i<DL; i++){
      for(int j=0; j<DL; j++){
        std::cout<<svd_penult.matrixU().col(i).dot( svd_last.matrixU().col(j) )<<" ";
      }
      std::cout<<std::endl;
    }
    std::cout<<"V overlaps:"<<std::endl;
    for(int i=0; i<DL; i++){
      for(int j=0; j<DL; j++){
        std::cout<<svd_penult.matrixV().col(i).dot( svd_last.matrixV().col(j) )<<" ";
      }
      std::cout<<std::endl;
    }
    std::cout<<"sigma (penultimate): "<<svd_penult.singularValues().transpose()<<std::endl;
    std::cout<<"sigma (last): "<<svd_last.singularValues().transpose()<<std::endl;
    std::cout<<"ratios:";
    for(int i=0; i<DL; i++){
      std::cout<<" "<<svd_last.singularValues()(i)/svd_penult.singularValues()(i);
    }std::cout<<std::endl;
    
    std::cout<<"(svd_penult.matrixU().adjoint()*svd_last.matrixU()).trace()="<<(svd_penult.matrixU().adjoint()*svd_last.matrixU()).trace()<<std::endl;
    std::cout<<"(svd_penult.matrixV().adjoint()*svd_last.matrixV()).trace()="<<(svd_penult.matrixV().adjoint()*svd_last.matrixV()).trace()<<std::endl;
    std::cout<<"(svd_penult.matrixU()-svd_last.matrixU()).norm()="<<(svd_penult.matrixU()-svd_last.matrixU()).norm()<<std::endl;
    std::cout<<"(svd_penult.matrixV()-svd_last.matrixV()).norm()="<<(svd_penult.matrixV()-svd_last.matrixV()).norm()<<std::endl;
  }

  if(param.is_specified("outfile")){
//    param_prop.print();
    TimeGrid tgrid(param_prop);
    InitialState initial(param_prop);
    OutputPrinter printer(param_prop);
    std::cout<<"Propagating from initial state:"<<std::endl;
    std::cout<<(Eigen::MatrixXcd)initial<<std::endl;   

    dyn.propagate(param_prop, initial, printer);
//    dyn.propagate_dE(param_prop, initial, printer);
  }
 }catch(std::exception &e){
   throw e;
 }
  return 0;
}
