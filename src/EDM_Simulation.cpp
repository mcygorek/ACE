#include "EDM_State.hpp"
#include "EDM_Simulation.hpp"
#include "EDM_TwoBodyTerm_FreePropagator.hpp"
#include "EDM_TwoBodyTerm_PT.hpp"
#include "LiouvilleTools.hpp"
#include "ReadExpression.hpp"
#include "ReaderBasics.hpp"
#include "DummyException.hpp"

namespace ACE {
  void EDM_Simulation::extract_Observables(double t, const EDM_State &state){ 
    //check consistency
    for(int o=0; o<oneBodyObservables.size(); o++){
      if(oneBodyObservables[o].second<0 || oneBodyObservables[o].second>=state.Ldim.size()){
        std::cerr<<"extract_Observables: oneBodyObservables[o].second<0 || oneBodyObservables[o].second>=state.Ldim.size()!"<<std::endl;
        throw DummyException();
      }
      if(oneBodyObservables[o].first.rows()!=state.Ldim[oneBodyObservables[o].second]){
        std::cerr<<"extract_Observables: oneBodyObservables[o].first.rows()!=state.Ldim[oneBodyObservables[o].second]!"<<std::endl;
        throw DummyException();
      }
    }

    //check which reduced density matrices are required and calculate them:
    std::vector<bool> which_reduced(state.Ldim.size(), false);
    for(int o=0; o<oneBodyObservables.size(); o++){
      which_reduced[oneBodyObservables[o].second]=true;
    }
    std::vector<Eigen::VectorXcd> reduced(state.Ldim.size(),Eigen::VectorXcd(0));
    for(int r=0; r<state.Ldim.size(); r++){
      if(which_reduced[r]){
        reduced[r]=state.get_reduced(r);
      }
    }

    //evaluate
    if(printer){
      printer.print_time(t);
      for(int o=0; o<oneBodyObservables.size(); o++){
        std::complex<double> c = (oneBodyObservables[o].first.transpose()* 
                               reduced[oneBodyObservables[o].second])(0,0);
        printer.print_value(c);
      }
      printer.print_endl();
    }
  }

  void EDM_Simulation::apply_OneBodyTerms(int n, const TimeGrid &tgrid, EDM_State &state, const EDM_Filter &filter, bool reverse_order){
   if(reverse_order){
    for(int i=(int)oneBodyTerms.size()-1; i>=0; i--){ 
      try{
        oneBodyTerms[i].first->update(tgrid.get_t(n), tgrid.get_dt(n));
        EDM_State out=oneBodyTerms[i].first->apply_filtered(state, oneBodyTerms[i].second, filter);
        state.swap(out);
      }catch (DummyException &e){
        std::cerr<<"called by EDM_Simulation::apply_OneBodyTerms: i="<<i<<" n="<<n<<std::endl; 
        throw e;
      }
    }
   }else{
    for(int i=0; i<oneBodyTerms.size(); i++){ 
      try{
        oneBodyTerms[i].first->update(tgrid.get_t(n), tgrid.get_dt(n));
        EDM_State out=oneBodyTerms[i].first->apply_filtered(state, oneBodyTerms[i].second, filter);
        state.swap(out);
      }catch (DummyException &e){
        std::cerr<<"called by EDM_Simulation::apply_OneBodyTerms: i="<<i<<" n="<<n<<std::endl; 
        throw e;
      }
    }
   }
  }

  void EDM_Simulation::apply_TwoBodyTerms(int n, const TimeGrid &tgrid, EDM_State &state, const EDM_Filter &filter, bool reverse_order){
   if(reverse_order){
    for(int i=(int)twoBodyTerms.size()-1; i>=0; i--){
      try{
        twoBodyTerms[i].first->update(n, tgrid);
        EDM_State out=twoBodyTerms[i].first->apply_filtered(state, twoBodyTerms[i].second, filter);
        state.swap(out);
      }catch (DummyException &e){
        std::cerr<<"called by EDM_Simulation::apply_TwoBodyTerms: i="<<i<<" n="<<n<<std::endl; 
        throw e;
      }
    }
   }else{
    for(int i=0; i<twoBodyTerms.size(); i++){ 
      try{
        twoBodyTerms[i].first->update(n, tgrid);
        EDM_State out=twoBodyTerms[i].first->apply_filtered(state, twoBodyTerms[i].second, filter);
        state.swap(out);
      }catch (DummyException &e){
        std::cerr<<"called by EDM_Simulation::apply_TwoBodyTerms: i="<<i<<" n="<<n<<std::endl; 
        throw e;
      }
    }
   }
  }

  void EDM_Simulation::one_step(int n, const TimeGrid &tgrid, EDM_State &state, const EDM_Filter &filter){
    if(use_symmetric_Trotter && (!propagate_alternate)){
      TimeGrid tgrid2=tgrid; 
      tgrid2.half_dt();
      apply_OneBodyTerms(2*n, tgrid2, state, filter);
      apply_TwoBodyTerms(n, tgrid, state, filter);
      apply_OneBodyTerms(2*n+1, tgrid2, state, filter);
    }else if(n%2==0 || !propagate_alternate){
      apply_TwoBodyTerms(n, tgrid, state, filter);
      apply_OneBodyTerms(n, tgrid, state, filter);
    }else{
      apply_OneBodyTerms(n, tgrid, state, filter, true);
      apply_TwoBodyTerms(n, tgrid, state, filter, true);
    }
  }

  void EDM_Simulation::run(const TimeGrid &tgrid, const EDM_State &initial, const EDM_Filter &filter){
    if(oneBodyObservables.size()<1){
      std::cerr<<"EDM_Simulatin::run: No observable defined!"<<std::endl;
      throw DummyException();
    }

    EDM_State state=initial;
    extract_Observables(tgrid.get_t(0), state);

    for(int n=0; n<tgrid.n_tot; n++){
     try{
      if(print_timesteps){
        double total=1.;
        for(const int &i : state.Ldim.list){ total*=i; }
        std::cout<<"time="<<tgrid.get_t(n)<<": ";
        std::cout<<" coefficients="<<state.coeffs.size();
        std::cout<<" candidates="<<state.candidates.size();
        std::cout<<" of "<<total;
        std::cout<<std::endl;
      }
      one_step(n, tgrid, state, filter);

      print_debug(tgrid.get_t(n+1), state);
      extract_Observables(tgrid.get_t(n+1), state);

     }catch(DummyException &e){
       std::cerr<<"called by EDM_Simulation::run at n="<<n<<std::endl;
       throw e;
     }
    }

/*
    {
     state.reduce(filter);
     std::ofstream ofs("DUMP.dat");
     for(const auto &elem : state.coeffs){
       ofs<<std::abs(elem.second)<<std::endl;
     } 
    }
    {
     std::ofstream ofs("DUMP2.dat");
     for(const auto &elem : state.candidates){
       ofs<<std::abs(elem.second)<<std::endl;
     } 
    }
*/
  }
  

  void EDM_Simulation::print_debug(double t, const EDM_State &state){
    if(debug_printer){
      (*debug_printer.ofs)<<t<<" "<<state.coeffs.size()<<std::endl;
    }
  }

  void EDM_Simulation::clear(){
    oneBodyTerms.clear();
    twoBodyTerms.clear();
    oneBodyObservables.clear();
  }
  void EDM_Simulation::setup(Parameters &param){
    int N_sites=param.get_as_size_t("N_sites",1);
   
    printer.setup(param); 
    debug_printer.setup(param, "", "info_outfile"); 
    print_timesteps=param.get_as_bool("print_timesteps",false);
    use_symmetric_Trotter=param.get_as_bool("use_symmetric_Trotter",false);
    propagate_alternate=param.get_as_bool("propagate_alternate",true);

    //Define observables:
    for(int r=0; r<N_sites; r++){
      std::string key="S"+int_to_string(r)+"_add_Output";
      std::vector<std::vector<std::string> > svv=param.get(key);
      for(int o=0; o<svv.size(); o++){
        Eigen::MatrixXcd op=ReadExpression(svv[o][0]);
        oneBodyObservables.push_back(std::make_pair(
                        H_Matrix_to_L_Vector(op.transpose()),r));
      }
    }

    //OneBodyTerms:
    for(int r=0; r<N_sites; r++){
      std::shared_ptr<EDM_OneBodyTerm> OBT(new EDM_OneBodyTerm_FreePropagator(param, r));
      if(dynamic_cast<EDM_OneBodyTerm_FreePropagator*>(OBT.get())->fprop){
        oneBodyTerms.push_back(std::make_pair(OBT, r));
      }
    }
std::cout<<"oneBodyTerms.size()="<<oneBodyTerms.size()<<std::endl;

    //TwoBodyTerms:
    for(int r0=0; r0<N_sites; r0++){
      for(int r1=0; r1<N_sites; r1++){
        std::shared_ptr<EDM_TwoBodyTerm> TBT(new EDM_TwoBodyTerm_FreePropagator(param, std::make_pair(r0,r1)));
        if(dynamic_cast<EDM_TwoBodyTerm_FreePropagator*>(TBT.get())->fprop){
          twoBodyTerms.push_back(std::make_pair(TBT, std::make_pair(r0,r1)));
        }
      }
    }
std::cout<<"Systems: twoBodyTerms.size()="<<twoBodyTerms.size()<<std::endl;
    //PT-MPOs:
    for(int r0=0; r0<N_sites; r0++){
      for(int r1=0; r1<N_sites; r1++){
        std::shared_ptr<EDM_TwoBodyTerm> TBT(new EDM_TwoBodyTerm_PT(param, std::make_pair(r0,r1)));
        if(dynamic_cast<EDM_TwoBodyTerm_PT*>(TBT.get())->PT){
          twoBodyTerms.push_back(std::make_pair(TBT, std::make_pair(r0,r1)));
        }
      }
    }
std::cout<<"Systems+PT-MPOs: twoBodyTerms.size()="<<twoBodyTerms.size()<<std::endl;

  }
}
