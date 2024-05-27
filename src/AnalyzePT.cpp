#include "AnalyzePT.hpp"
#include "InfluenceFunctional_OD.hpp"
#include "MPS_Matrix.hpp"
#include "PrintPPM.h"

namespace ACE{

  void AnalyzePT::canvas_add_MPS_Matrix(PPM_Canvas &canv, const MPS_Matrix &m, int a, int posx, int posy){
    int w=m.dim_d2;
    int h=m.dim_d1;
    if(posx<0||posy<0||posx+w>canv.get_width()||posy+h>canv.get_height()){
      std::cerr<<"canvas_add_MPS_Matrix: dimensions out of bounds: ";
      std::cerr<<"("<<posx<<".."<<posx+w<<", "<<posy<<".."<<posy+h<<") / ";
      std::cerr<<"("<<canv.get_width()<<", "<<canv.get_height()<<")!"<<std::endl;
      exit(1);
    }
    for(int d1=0; d1<h; d1++){
      for(int d2=0; d2<w; d2++){
        canv(posx+d2,posy+d1)=abs(m(a, d1, d2));
      }
    }
  }
  void AnalyzePT::print_single_pgm(const std::string &file, const InfluenceFunctional_OD &IF, int n, int a){

    if(n<0||n>=(int)IF.a.size()){
      std::cerr<<"print_single_pgm: n<0||n>=IF.a.size()!"<<std::endl;
      exit(1);
    }
    if(a<0||a>=IF.a[n].dim_i){
      std::cerr<<"print_single_pgm: a<0||a>=IF.a[n].dim_a!"<<std::endl;
      exit(1);
    }

    double maxval=0;
    double minval=IF.a[n].max_element_abs();

    PPM_Canvas canv(IF.a[n].dim_d2,IF.a[n].dim_d1);
    canvas_add_MPS_Matrix(canv, IF.a[n], a, 0, 0);
    canv.print_pgm(file, maxval, minval);
  }


  void AnalyzePT::print_pgm(const std::string &file, const InfluenceFunctional_OD &IF, int margins){
    int margin=margins;
    int hmargin=margins;
    
    int nr=IF.a.size();
    if(nr<1){
      std::cerr<<"AnalyzePT::print_ppm: IF empty!"<<std::endl;
      exit(1);
    }
    int dim_a=IF.dict.get_reduced_dim();

    int totalwidth=margin*(nr+1);

    std::vector<int> widths;
    int maxheight=0;
    for(int n=0; n<nr; n++){
      widths.push_back(IF.a[n].dim_d2);
      totalwidth+=widths.back();
      if(IF.a[n].dim_d1>maxheight)maxheight=IF.a[n].dim_d1;
    }

    int totalheight=(dim_a*maxheight)+(dim_a+1)*hmargin;


    double maxval=0;
    double minval=IF.max_element_abs();
    PPM_Canvas canv(totalwidth,totalheight);
    canv.fill(0.2*minval);
    
    for(int a=0; a<dim_a; a++){
      int posy=hmargin*(a+1)+maxheight*a;
      int posx=margin;
      for(int n=0; n<nr; n++){
        canvas_add_MPS_Matrix(canv, IF.a[n], a, posx, posy);
        posx+=IF.a[n].dim_d2+margin;
      }
    }
    canv.print_pgm(file, maxval, minval);
  }

  void AnalyzePT::print_SVD(const std::string &file, const InfluenceFunctional_OD &IF, const std::string &print_PT){

    std::ofstream ofs(file.c_str());

    if(print_PT!=""){ 
      InfluenceFunctional_OD IF2; 
      IF2.set_none(IF.a.size(), IF.get_sys_dim() );
      Eigen::MatrixXcd lastR=Eigen::MatrixXcd::Identity(1,1);
      for(int n=0; n<IF2.a.size()-1; n++){
        MPS_Matrix M = IF.get_a(n);
        M.inner_multiply_left(lastR);
        Eigen::MatrixXcd A = M.get_Matrix_d1i_d2();
        Eigen::JacobiSVD<Eigen::MatrixXcd> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
        ofs<<n<<": "<<svd.singularValues().transpose()<<std::endl;

        Eigen::VectorXcd diag=svd.singularValues();
        if(diag.size()>0) diag /= svd.singularValues()(0);
        lastR = diag.asDiagonal()*svd.matrixV().adjoint();
 
        IF2.a[n].set_from_Matrix_d1i_d2(svd.matrixU()*svd.singularValues()(0), IF.get_a(n).dim_i);
      }
      IF2.a[IF2.a.size()-1] = IF.get_a(IF.a.size()-1); 
      IF2.a[IF2.a.size()-1].inner_multiply_left(lastR);
 
      IF2.dict=IF.dict;
      IF2.calculate_closures();
      IF2.write_binary(print_PT);

    }else{
      Eigen::MatrixXcd lastR=Eigen::MatrixXcd::Identity(1,1);
      for(int n=0; n<IF.a.size(); n++){
        MPS_Matrix M = IF.get_a(n);
        M.inner_multiply_left(lastR);
        Eigen::MatrixXcd A = M.get_Matrix_d1i_d2();
    
        Eigen::JacobiSVD<Eigen::MatrixXcd> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
        ofs<<n<<": "<<svd.singularValues().transpose()<<std::endl;

        Eigen::VectorXcd diag=svd.singularValues();
        if(diag.size()>0) diag /= svd.singularValues()(0);
        lastR = diag.asDiagonal()*svd.matrixV().adjoint();
      }
    }
  }

  void AnalyzePT::print_summary(std::ostream &ofs, const InfluenceFunctional_OD &IF){
    ofs<<"Length: "<<IF.size()<<std::endl;
    ofs<<"Dictionary: "; IF.dict.print_beta(ofs); std::cout<<std::endl;
    ofs<<"Dimensions: "; IF.print_inner_dims(ofs); std::cout<<std::endl;
    ofs<<"Largest element: "<<IF.max_element_abs()<<std::endl;
  }

  void AnalyzePT::print_summary(const std::string &file, const InfluenceFunctional_OD &IF){
    std::ofstream ofs(file.c_str());
    print_summary(ofs, IF);
  }

}//namespace
