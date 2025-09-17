#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>


const double hbar=0.658212; //meVps
const double wn=0.188364997625644; //cm/ps

double J_AR(double w){
  double S = 0.29;
  double s1 = 0.8;
  double s2 = 0.5;
  double w1 = 0.069/hbar;
  double w2 = 0.24/hbar;
  return S/(s1+s2)*(s1/(5040.*2.*pow(w1,4))*pow(w,5)*exp(-sqrt(w/w1)) 
                   +s2/(5040.*2.*pow(w2,4))*pow(w,5)*exp(-sqrt(w/w2)));
}
double J_62(double w, double gamma, const std::vector<std::pair<double, double> > &wS){
  double res=0.;
  for(int i=0; i<wS.size(); i++){
    double wk=wS[i].first;
    double sk=wS[i].second;
    res+=4.*wk*sk*gamma*(wk*wk+gamma*gamma)*w/
(M_PI*(pow(w+wk,2)+pow(gamma,2))*(pow(w-wk,2)+pow(gamma,2)));
  }
  return res;
}
double J(double w, double gamma, const std::vector<std::pair<double, double> > &wS){
  return J_62(w,gamma,wS)+J_AR(w);
}
int main(int args, char **argv){
  double gamma=1;
  std::vector<std::pair<double, double> > wS;


  std::ifstream ifs("62modes.dat");
  std::string line;
  while(std::getline(ifs, line)){
    std::istringstream iss(line);
    double w, S;
    if(!(iss>>w>>S)){break;}
    wS.push_back(std::make_pair(w*wn, S));
  }
  

  double wmax=1000;
  int N = 100000;
  double dw=wmax/N;
  for(int i=0; i<N; i++){
    double w=i*dw;
    std::cout<<w<<" "<<J(w, gamma, wS)<<std::endl;
  }
 
  return 0;
}
