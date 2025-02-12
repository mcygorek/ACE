#include "ACE.hpp"
#include "Parameters.hpp"
#include "ReadTable.hpp"
#include <Eigen/Dense>
#include "LeastSquares.hpp"

using namespace ACE;

int main(int args, char** argv){

  Parameters param(args, argv);
  std::string infile=param.get_as_string_check("infile");

  int N; int n_tot; double dt; Eigen::VectorXd ref;
  {
    ReadTable tab(infile, 0, 1, 2);
    if(tab.size()<2){
      std::cerr<<"File '"<<infile<<"' has less than 2 usable lines!"<<std::endl;
      exit(1);
    }
    n_tot=tab.size();
    N=tab.size();
    dt=tab[1][0]-tab[0][0];
    ref=Eigen::VectorXd::Zero(2*N);
    for(int i=0; i<N; i++){
      ref(i)=tab[i][1];
      ref(i+N)=tab[i][2];
    }
  }
  
  //Parameters:  a * exp(-(k+i*w)*t) + b * exp(-(k-i*w)*t)
  Eigen::VectorXd guess=Eigen::VectorXd::Zero(4); 
  guess(0)=param.get_as_double_check("gamma"); //k  
  guess(1)=param.get_as_double_check("omega"); //w
  guess(2)=param.get_as_double_check("a"); //a
  guess(3)=param.get_as_double_check("b"); //b

  double scale_exp=param.get_as_double("scale_exp",0);

  auto get_time = [dt] (int n) -> double {
      return n*dt;
    };
  
  auto scale = [scale_exp] (double t) -> double {
      return exp(scale_exp*t);
    };

  Eigen::VectorXd ref_scale=ref;
  for(int n=0; n<N; n++){
    double t=get_time(n);
    ref_scale(n)*=scale(t);
    ref_scale(n+N)*=scale(t);
  }

  auto f = [N,get_time] (Eigen::VectorXd guess) -> Eigen::VectorXd {
      Eigen::VectorXd res=Eigen::VectorXd::Zero(2*N);
      for(int n=0; n<N; n++){
        double t=get_time(n);
        res(n)=(guess(2)+guess(3))*exp(-guess(0)*t)*cos(guess(1)*t);
        res(n+N)=(-guess(2)+guess(3))*exp(-guess(0)*t)*sin(guess(1)*t);
      }
      return res;
    };
  auto f_scale = [N,get_time,f,scale] (Eigen::VectorXd guess) -> Eigen::VectorXd {
      Eigen::VectorXd res=f(guess);
      for(int n=0; n<N; n++){
        double t=get_time(n);
        res(n)*=scale(t);
        res(n+N)*=scale(t);
      }
      return res;
    };
  auto f_with_prefactor = [N,dt] (Eigen::VectorXd guess) -> Eigen::VectorXd {
      Eigen::VectorXd res=Eigen::VectorXd::Zero(2*N);
      double k=guess(0);
      double w=guess(1);
      double g1=sqrt(guess(2))/dt;
      double g2=sqrt(guess(3))/dt;
      std::complex<double> i(0.,1.);
      std::complex<double> c;
      c=g1*g1*(dt/(k+i*w)+(exp(-(k+i*w)*dt)-1.)/(k+i*w)/(k+i*w));
      c+=g2*g2*(dt/(k-i*w)+(exp(-(k-i*w)*dt)-1.)/(k-i*w)/(k-i*w));
      res(0)=c.real();
      res(N)=c.imag();

      for(int n=1; n<N; n++){
        double t=n*dt;
        c=g1*g1*(exp((k+i*w)*dt)+exp(-(k+i*w)*dt)-2.)/(k+i*w)/(k+i*w)*exp(-(k+i*w)*t);
        c+=g2*g2*(exp((k-i*w)*dt)+exp(-(k-i*w)*dt)-2.)/(k-i*w)/(k-i*w)*exp(-(k-i*w)*t);
        res(n)=c.real();
        res(n+N)=c.imag();
      }
      return res;
    };
  auto J = [N,get_time] (Eigen::VectorXd guess) -> Eigen::MatrixXd {
      Eigen::MatrixXd res=Eigen::MatrixXd::Zero(2*N,4);
      for(int n=0; n<N; n++){
        double t=get_time(n);
        res(n,0)=-t*(guess(2)+guess(3))*exp(-guess(0)*t)*cos(guess(1)*t);
        res(n+N,0)=-t*(-guess(2)+guess(3))*exp(-guess(0)*t)*sin(guess(1)*t);
        res(n,1)=-t*(guess(2)+guess(3))*exp(-guess(0)*t)*sin(guess(1)*t);
        res(n+N,1)=t*(-guess(2)+guess(3))*exp(-guess(0)*t)*cos(guess(1)*t);
        res(n,2)=exp(-guess(0)*t)*cos(guess(1)*t);
        res(n,3)=exp(-guess(0)*t)*cos(guess(1)*t);
        res(n+N,2)=-exp(-guess(0)*t)*sin(guess(1)*t);
        res(n+N,3)=exp(-guess(0)*t)*sin(guess(1)*t);
      }
      return res;
    };
  auto J_scale = [N,get_time,J,scale] (Eigen::VectorXd guess) -> Eigen::MatrixXd {
      Eigen::MatrixXd res=J(guess);
      for(int n=0; n<N; n++){
        double t=get_time(n);
        for(int i=0; i<guess.rows(); i++){
          res(n,i)*=scale(t);
          res(n+N,i)*=scale(t);
        }
      }
      return res;
    };


  guess=LeastSquares::GaussNewton(ref_scale, f_scale, J_scale, guess, 1e-8, 100, 1);
  bool with_prefactor=param.get_as_bool("with_prefactor",true);
  Eigen::VectorXd res;
  if(with_prefactor){
    res=f_with_prefactor(guess);
  }else{
    res=f(guess);
  }
  

  std::cout<<"Final: "<<guess.transpose()<<std::endl;

  //-----------------

  //-----------------
  double k=guess(0);
  double w=guess(1);
  double g1=sqrt(guess(2))/dt;
  double g2=sqrt(guess(3))/dt;

  std::string outfile=param.get_as_string("outfile");
  if(outfile!=""){
    std::ofstream ofs(outfile.c_str());
    ofs<<"# fit: "<<guess.transpose()<<std::endl;
    for(int n=0; n<N; n++){
      double t=n*dt;
      ofs<<t<<" "<<res(n)<<" "<<res(n+N)<<std::endl;
    }
  }

  std::string print_diff=param.get_as_string("print_diff","");
  if(print_diff!=""){
    std::ofstream ofs_diff(print_diff.c_str());
    ofs_diff<<"# fit: "<<guess.transpose()<<std::endl;
    for(int n=0; n<N; n++){
      double t=n*dt;
      ofs_diff<<t<<" "<<ref(n)-res(n)<<" "<<ref(n+N)-res(n+N)<<std::endl;
    }
    std::cout<<"printed differences to file '"<<print_diff<<"'."<<std::endl;
  }

  std::string prefix=param.get_as_string("prefix");
  if(prefix!=""){
    int M=param.get_as_size_t("M",10);
    std::string Mstr=int_to_string(M);
    std::string SysOp=param.get_as_string("SysOp", "|1><1|_2");

    for(int i=1; i<=2; i++){
      double sign=1;
      double g=g1;
      if(i==2){ sign=-1; g=g2;}
      std::string prefixmparam=prefix+int_to_string(i)+".mparam";
      std::cout<<prefixmparam<<std::endl;
      {
        std::ofstream ofs(prefixmparam.c_str());
        ofs<<"add_Hamiltonian {hbar*("<<sign*w<<")*(Id_2 otimes n_"<<M<<") }"<<std::endl;
        ofs<<"add_Hamiltonian {hbar*"<<g<<"*({"<<SysOp<<"} otimes (bdagger_"<<M<<"+b_"<<M<<"))}"<<std::endl;
        ofs<<"add_Lindblad {2*"<<k<<"} { Id_2 otimes b_"<<M<<" }"<<std::endl;
      } 

      double threshold=param.get_as_double("threshold",1e-9);

      std::string prefixparam=prefix+int_to_string(i)+".param";
      std::cout<<prefixparam<<std::endl;
      {
        std::ofstream ofs(prefixparam.c_str());
        ofs<<"te "<<dt*n_tot<<std::endl;
        ofs<<"dt "<<dt<<std::endl;
        ofs<<"threshold "<<threshold<<std::endl;
        ofs<<"dict_zero 1e-16"<<std::endl;
        ofs<<"add_single_mode_from_file "<<prefixmparam<<" {|0><0|_"<<M<<"}"<<std::endl;
        ofs<<"initial  {0.5*(Id_2+sigma_x)}"<<std::endl;
        ofs<<"add_Output {|0><1|_2}"<<std::endl;
        ofs<<"write_PT "<<prefix<<i<<".pt"<<std::endl;
        ofs<<"outfile "<<prefix<<i<<".out"<<std::endl;
      }

 
      if(param.get_as_bool("print_generateBCF")){
        std::string prefixBCF=prefix+int_to_string(i)+"_generateBCF.param";
        std::cout<<prefixBCF<<std::endl;
        std::ofstream ofs(prefixBCF.c_str());
        ofs<<"te "<<dt*n_tot<<std::endl;
        ofs<<"dt "<<dt<<std::endl;
        ofs<<"initial {"<<g<<"*(bdagger_"<<M<<"+b_"<<M<<")*|0><0|_"<<M<<"}"<<std::endl;
        ofs<<"add_Hamiltonian {hbar*("<<sign*w<<")*(n_"<<M<<") }"<<std::endl;      
        ofs<<"add_Output {"<<g<<"*(bdagger_"<<M<<"+b_"<<M<<")}"<<std::endl;
        ofs<<"add_Lindblad {2*"<<k<<"} { b_"<<M<<" }"<<std::endl;
        ofs<<"outfile "<<prefix<<i<<"_generateBCF.out"<<std::endl;
      }

    }
  }

  return 0;
}


