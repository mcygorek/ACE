#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/iostream.h>
#include "Simulation_PT.hpp"
#include "FreePropagator.hpp"
#include "ComplexFunction_Interpolate.hpp"
#include "DynamicalMap.hpp"

namespace py = pybind11;


PYBIND11_MODULE(ACE, m) {
  m.doc() = "pybind11 ACE plugin"; // optional module docstring
  m.attr("hbar") = ACE::hbar_in_meV_ps;

  py::class_<ACE::Parameters>(m, "Parameters")
    .def(py::init<>())
    .def(py::init<const std::string &>())
    .def(py::init([](const py::list &list){ 
      auto param=std::unique_ptr<ACE::Parameters>(new ACE::Parameters()); 
      for(const auto & line: list){
        param->add_from_line(line.cast<std::string>());
      }return param;}))
    .def("add_from_file", &ACE::Parameters::add_from_file)
    .def("get_as_single_string", &ACE::Parameters::get_as_single_string, 
      "get value as a single line", py::arg("key"), py::arg("row")=0)
    .def("get_lines", &ACE::Parameters::get_lines)
    .def("__repr__", [](const ACE::Parameters &param){
      auto list=param.get_lines(); std::string res;
      for(const std::string & line : list){ res+=line+"\n"; } return res;})
    .def("add_from_line", &ACE::Parameters::add_from_line)
    .def("add_from_file", &ACE::Parameters::add_from_file)
    .def("add", [](ACE::Parameters &param, const py::list &list){
      for(const auto & line: list){
        param.add_from_line(line.cast<std::string>());}})
    .def("clear", &ACE::Parameters::clear)
    ;

  py::class_<ACE::FreePropagator>(m, "FreePropagator")
    .def(py::init<>())
    .def(py::init<ACE::Parameters &>())
    .def("get_dim",&ACE::FreePropagator::get_dim)
    .def("set_dim",static_cast<int (ACE::FreePropagator::*)(int,const std::string &)>(&ACE::FreePropagator::set_dim), "set dimension", py::arg("dim"), py::arg("error_comment")="")
    .def("set_Hamiltonian",&ACE::FreePropagator::set_Hamiltonian)
    .def("add_Hamiltonian",&ACE::FreePropagator::add_Hamiltonian)
    .def("add_Pulse",static_cast<void (ACE::FreePropagator::*)(const std::pair<std::vector<double>,std::vector<std::complex<double> > > & shape, const Eigen::MatrixXcd &A)>(&ACE::FreePropagator::add_Pulse))
    .def("add_forward_Pulse",static_cast<void (ACE::FreePropagator::*)(const std::pair<std::vector<double>,std::vector<std::complex<double> > > & shape, const Eigen::MatrixXcd &A)>(&ACE::FreePropagator::add_forward_Pulse))
    .def("add_backward_Pulse",static_cast<void (ACE::FreePropagator::*)(const std::pair<std::vector<double>,std::vector<std::complex<double> > > & shape, const Eigen::MatrixXcd &A)>(&ACE::FreePropagator::add_backward_Pulse))
    .def("add_Lindblad", [](ACE::FreePropagator &prop, double gamma, const Eigen::MatrixXcd & L){prop.add_Lindblad(gamma, L);})
    .def("apply_Operator_left", [](ACE::FreePropagator &prop, double time, const Eigen::MatrixXcd & Op, bool apply_before){prop.add_MultitimeOp(time, Op, Eigen::MatrixXcd::Identity(Op.rows(), Op.cols()), apply_before);}, 
      py::arg("time"), py::arg("Op"), py::arg("apply_before")=false)
    .def("apply_Operator_right", [](ACE::FreePropagator &prop, double time, const Eigen::MatrixXcd & Op, bool apply_before){prop.add_MultitimeOp(time, Eigen::MatrixXcd::Identity(Op.rows(), Op.cols()), Op, apply_before);}, 
      py::arg("time"), py::arg("Op"), py::arg("apply_before")=false)
    .def("get_Htot",&ACE::FreePropagator::get_Htot)
    .def_readwrite("const_H",&ACE::FreePropagator::const_H)
    ;

  py::class_<ACE::ProcessTensorForwardList>(m, "ProcessTensors")
    .def(py::init<>())
    .def(py::init([](ACE::Parameters & param){
      py::scoped_ostream_redirect stream(std::cout,py::module_::import("sys").attr("stdout"));
      return new ACE::ProcessTensorForwardList(param);}))
    ;
//    .def(py::init<ACE::Parameters &>())
   
  py::class_<ACE::InitialState>(m, "InitialState")
    .def(py::init<>())
    .def(py::init<ACE::Parameters &>())
    .def(py::init<const Eigen::MatrixXcd &>())
    .def_readwrite("rho", &ACE::InitialState::rho)
    ;

  py::class_<ACE::TimeGrid>(m, "TimeGrid")
    .def(py::init<>())
    .def(py::init<ACE::Parameters &>())
    .def(py::init([](double ta, double te, double dt){
       ACE::Parameters param;
       param.add_to("ta", ta); 
       param.add_to("te", te); 
       param.add_to("dt", dt); 
      return new ACE::TimeGrid(param);}))
    .def_readwrite("ta", &ACE::TimeGrid::ta)
    .def_readwrite("dt", &ACE::TimeGrid::dt)
    .def_readwrite("n_tot", &ACE::TimeGrid::n_tot)
    .def_readwrite("n_mem", &ACE::TimeGrid::n_mem)
    .def("get_dt", &ACE::TimeGrid::get_dt)
    .def("get_t", &ACE::TimeGrid::get_t)
    .def("get_all", &ACE::TimeGrid::get_all)
    .def("get_closest_n", &ACE::TimeGrid::get_closest_n)
    ;

  py::class_<ACE::OutputPrinter>(m, "OutputPrinter")
    .def(py::init<>())
    .def(py::init<ACE::Parameters &>())
    .def(py::init<const std::string &, const std::vector<Eigen::MatrixXcd>&>())
/*    .def(py::init([](ACE::Parameters param){
       bool do_extract=false;
       if(!param.is_specified("outfile")){
         do_extract=true;
         param.add_to("outfile","/dev/null");
       }
       std::unique_ptr<ACE::OutputPrinter> printer(new ACE::OutputPrinter(param));
       printer->do_extract=do_extract;
       return printer;
     }))*/
     .def_readwrite("do_extract",&ACE::OutputPrinter::do_extract)
     .def("clear",&ACE::OutputPrinter::clear)
     .def("extract",&ACE::OutputPrinter::extract)
    ;

  py::class_<ACE::Simulation_PT>(m, "Simulation")
    .def(py::init<>())
    .def(py::init<ACE::Parameters &>())
//    .def("run", static_cast<Eigen::MatrixXcd (ACE::Simulation_PT::*)(
//      ACE::Propagator &prop, ACE::ProcessTensorForwardList &PT,
//      const Eigen::MatrixXcd & initial, const ACE::TimeGrid &tgrid,
//      ACE::OutputPrinter &printer)>(&ACE::Simulation_PT::run))
    .def("run",&ACE::Simulation_PT::run_)
    ;

  py::class_<ACE::DynamicalMap>(m, "DynamicalMap")
    .def(py::init<>())
    .def(py::init<ACE::Parameters &>())
    .def(py::init<ACE::FreePropagator &, ACE::ProcessTensorForwardList &, ACE::Simulation_PT &, const ACE::TimeGrid &>())
    .def("calculate", static_cast<void (ACE::DynamicalMap::*)(ACE::Parameters &)>(&ACE::DynamicalMap::calculate))
    .def("calculate", static_cast<void (ACE::DynamicalMap::*)(ACE::Propagator &, ACE::ProcessTensorForwardList &, ACE::Simulation_PT &, const ACE::TimeGrid &, int)>(&ACE::DynamicalMap::calculate))
    .def_readwrite("E", &ACE::DynamicalMap::E)
    .def("size", &ACE::DynamicalMap::size)
    .def("get_dE", [](const ACE::DynamicalMap& DM, int i){return DM.get_dE(i);})
//    .def("get_dE", &ACE::DynamicalMap::get_dE, py::arg("i") , py::arg("regularize") = 0., py::arg("os") = (std::ostream*)NULL)
//    .def("get_dE", static_cast<Eigen::MatrixXcd (ACE::DynamicalMap::*)(int)>(&ACE::DynamicalMap::get_dE))
    ;
}
