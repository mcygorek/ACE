/*
 * ACE Python bindings – nanobind version
 *
 * Exposed classes:
 *   Parameters, FreePropagator, ModePropagator,
 *   ProcessTensors (ProcessTensorForwardList),
 *   InitialState, TimeGrid, OutputPrinter,
 *   Simulation  (Simulation_PT),
 *   DynamicalMap
 *
 * GIL notes
 * ─────────
 * Heavy constructors and all run/calculate calls release the GIL so that
 * multiple threads can drive independent simulations simultaneously.
 * Constructors that accept Python-managed shared_ptr<ModePropagator> objects
 * intentionally keep the GIL to ensure safe reference-count manipulation.
 */

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/eigen/dense.h>

#include <memory>

#include "Simulation_PT.hpp"
#include "FreePropagator.hpp"
#include "ModePropagatorGenerator_SingleModes.hpp"
#include "DynamicalMap.hpp"
#include "ReadExpression.hpp"
#include "Constants.hpp"

namespace nb = nanobind;

NB_MODULE(_ACE, m) {
  m.doc() = "nanobind ACE module";
  m.attr("hbar") = ACE::hbar_in_meV_ps;

  m.def("StringToMatrix", [](const std::string& str) {
    nb::gil_scoped_release release;
    return (Eigen::MatrixXcd)ACE::ReadExpression(str);
  });

  // ── Parameters ─────────────────────────────────────────────────────────────
  nb::class_<ACE::Parameters>(m, "Parameters")
    .def(nb::init<>())
    .def(nb::init<const std::string&>(), nb::call_guard<nb::gil_scoped_release>())
    // Construct from a list of parameter lines
    .def("__init__", [](ACE::Parameters* self, const nb::list& lines) {
      new (self) ACE::Parameters();
      for (auto item : lines)
        self->add_from_line(nb::cast<std::string>(item));
    })
    .def("add_from_file",       &ACE::Parameters::add_from_file,
         nb::call_guard<nb::gil_scoped_release>())
    .def("add_from_line",       &ACE::Parameters::add_from_line)
    .def("get_as_single_string",&ACE::Parameters::get_as_single_string,
         nb::arg("key"), nb::arg("row") = 0)
    .def("get_lines",           &ACE::Parameters::get_lines)
    .def("add", [](ACE::Parameters& self, const nb::list& lines) {
      for (auto item : lines)
        self.add_from_line(nb::cast<std::string>(item));
    })
    .def("clear", &ACE::Parameters::clear)
    .def("__repr__", [](const ACE::Parameters& self) {
      std::string res;
      for (const auto& line : self.get_lines()) res += line + "\n";
      return res;
    })
    ;

  // ── FreePropagator ─────────────────────────────────────────────────────────
  nb::class_<ACE::FreePropagator>(m, "FreePropagator")
    .def(nb::init<>(),                  nb::call_guard<nb::gil_scoped_release>())
    .def(nb::init<ACE::Parameters&>(),  nb::call_guard<nb::gil_scoped_release>())
    .def("get_dim",
         &ACE::FreePropagator::get_dim,
         nb::call_guard<nb::gil_scoped_release>())
    .def("set_dim",
         static_cast<int (ACE::FreePropagator::*)(int, const std::string&)>(
           &ACE::FreePropagator::set_dim),
         nb::arg("dim"), nb::arg("error_comment") = "",
         nb::call_guard<nb::gil_scoped_release>())
    .def("set_Hamiltonian", &ACE::FreePropagator::set_Hamiltonian,
         nb::call_guard<nb::gil_scoped_release>())
    .def("add_Hamiltonian",
         &ACE::FreePropagator::add_Hamiltonian,
         nb::call_guard<nb::gil_scoped_release>())
    .def("add_Pulse",
         static_cast<void (ACE::FreePropagator::*)(
           const std::pair<std::vector<double>,
                           std::vector<std::complex<double>>>&,
           const Eigen::MatrixXcd&)>(&ACE::FreePropagator::add_Pulse),
         nb::call_guard<nb::gil_scoped_release>())
    .def("add_Lindblad",
         [](ACE::FreePropagator& self, double gamma, const Eigen::MatrixXcd& L) {
           self.add_Lindblad(gamma, L);
         },
         nb::call_guard<nb::gil_scoped_release>())
    .def("apply_Operator_left",
         [](ACE::FreePropagator& self, double t,
            const Eigen::MatrixXcd& Op, bool before) {
           self.add_MultitimeOp(
             t, Op,
             Eigen::MatrixXcd::Identity(Op.rows(), Op.cols()), before);
         },
         nb::arg("time"), nb::arg("Op"), nb::arg("apply_before") = false,
         nb::call_guard<nb::gil_scoped_release>())
    .def("apply_Operator_right",
         [](ACE::FreePropagator& self, double t,
            const Eigen::MatrixXcd& Op, bool before) {
           self.add_MultitimeOp(
             t,
             Eigen::MatrixXcd::Identity(Op.rows(), Op.cols()), Op, before);
         },
         nb::arg("time"), nb::arg("Op"), nb::arg("apply_before") = false,
         nb::call_guard<nb::gil_scoped_release>())
    .def("get_Htot",
         &ACE::FreePropagator::get_Htot,
         nb::call_guard<nb::gil_scoped_release>())
    .def("update",
         &ACE::FreePropagator::update,
         nb::call_guard<nb::gil_scoped_release>())
    .def_rw("propagate_Taylor",           &ACE::FreePropagator::propagate_Taylor)
    .def_rw("propagate_system_threshold", &ACE::FreePropagator::propagate_system_threshold)
    .def_rw("const_H",                    &ACE::FreePropagator::const_H)
    .def_rw("M",                          &ACE::FreePropagator::M)
    ;

  // ── ModePropagator ─────────────────────────────────────────────────────────
  nb::class_<ACE::ModePropagator, ACE::FreePropagator>(m, "ModePropagator")
    .def(nb::init<>(), nb::call_guard<nb::gil_scoped_release>())
    .def(nb::init<const std::string&, const Eigen::MatrixXcd&>(),
         nb::call_guard<nb::gil_scoped_release>())
    .def(nb::init<ACE::Parameters&,   const Eigen::MatrixXcd&>(),
         nb::call_guard<nb::gil_scoped_release>())
    .def_rw("initial",       &ACE::ModePropagator::bath_init)
    .def("get_N_system",     &ACE::ModePropagator::get_N_system)
    .def("get_N_mode",       &ACE::ModePropagator::get_N_mode)
    .def("get_bath_init",    &ACE::ModePropagator::get_bath_init,
         nb::call_guard<nb::gil_scoped_release>())
    .def("get_initial",      &ACE::ModePropagator::get_initial,
         nb::call_guard<nb::gil_scoped_release>())
    .def("copy", [](const ACE::ModePropagator& self) {
      nb::gil_scoped_release release;
      return ACE::ModePropagator(self);
    })
    ;

  // ── ProcessTensors (ProcessTensorForwardList) ───────────────────────────────
  nb::class_<ACE::ProcessTensorForwardList>(m, "ProcessTensors")
    .def(nb::init<>())
    // Construct from Parameters (heavy I/O – release GIL)
    .def("__init__",
         [](ACE::ProcessTensorForwardList* self, ACE::Parameters& param) {
           new (self) ACE::ProcessTensorForwardList(param);
         },
         nb::call_guard<nb::gil_scoped_release>())
    // Construct from Parameters + list of ModePropagators.
    // Build the generator list while holding the GIL (shared_ptr refcounting),
    // then release it for the heavy PT construction.
    .def("__init__",
         [](ACE::ProcessTensorForwardList* self,
            ACE::Parameters& param,
            std::vector<std::shared_ptr<ACE::ModePropagator>>& mpgs) {
           std::vector<std::shared_ptr<ACE::ModePropagatorGenerator>> conv;
           conv.push_back(
             std::make_shared<ACE::ModePropagatorGenerator_SingleModes>(mpgs));
           //nb::gil_scoped_release release;
           new (self) ACE::ProcessTensorForwardList(param, conv);
         })
    .def("add_PT",
         static_cast<void (ACE::ProcessTensorForwardList::*)(
           const std::string&, int, int)>(
           &ACE::ProcessTensorForwardList::add_PT),
         nb::arg("filename"),
         nb::arg("expand_dim_front") = 0,
         nb::arg("expand_dim_back")  = 0,
         nb::call_guard<nb::gil_scoped_release>())
    ;

  // ── InitialState ───────────────────────────────────────────────────────────
  nb::class_<ACE::InitialState>(m, "InitialState")
    .def(nb::init<>())
    .def(nb::init<ACE::Parameters&>(),         nb::call_guard<nb::gil_scoped_release>())
    .def(nb::init<const Eigen::MatrixXcd&>(),  nb::call_guard<nb::gil_scoped_release>())
    .def_rw("rho", &ACE::InitialState::rho)
    ;

  // ── TimeGrid ───────────────────────────────────────────────────────────────
  nb::class_<ACE::TimeGrid>(m, "TimeGrid")
    .def(nb::init<>())
    .def(nb::init<ACE::Parameters&>(), nb::call_guard<nb::gil_scoped_release>())
    // Convenience constructor: TimeGrid(ta, te, dt)
    .def("__init__", [](ACE::TimeGrid* self,
                        double ta, double te, double dt) {
      ACE::Parameters p;
      p.add_to("ta", ta);
      p.add_to("te", te);
      p.add_to("dt", dt);
      nb::gil_scoped_release release;
      new (self) ACE::TimeGrid(p);
    })
    .def_rw("ta",    &ACE::TimeGrid::ta)
    .def_rw("dt",    &ACE::TimeGrid::dt)
    .def_rw("n_tot", &ACE::TimeGrid::n_tot)
    .def_rw("n_mem", &ACE::TimeGrid::n_mem)
    .def("get_dt",        &ACE::TimeGrid::get_dt)
    .def("get_t",         &ACE::TimeGrid::get_t)
    .def("get_all",       &ACE::TimeGrid::get_all,
         nb::call_guard<nb::gil_scoped_release>())
    .def("get_closest_n", &ACE::TimeGrid::get_closest_n)
    ;

  // ── OutputPrinter ──────────────────────────────────────────────────────────
  nb::class_<ACE::OutputPrinter>(m, "OutputPrinter")
    .def(nb::init<ACE::Parameters&>(),
         nb::call_guard<nb::gil_scoped_release>())
    .def(nb::init<const std::string&, const std::vector<Eigen::MatrixXcd>&>(),
         nb::call_guard<nb::gil_scoped_release>())
    // OutputPrinter(op)  – extract a single observable
    .def("__init__", [](ACE::OutputPrinter* self,
                        const Eigen::MatrixXcd& op) {
      nb::gil_scoped_release release;
      new (self) ACE::OutputPrinter("", std::vector<Eigen::MatrixXcd>{op});
      self->do_extract = true;
    })
    // OutputPrinter([op1, op2, …])  – extract a list of observables
    .def("__init__", [](ACE::OutputPrinter* self,
                        const std::vector<Eigen::MatrixXcd>& ops) {
      nb::gil_scoped_release release;
      new (self) ACE::OutputPrinter("", ops);
      self->do_extract = true;
    })
    // OutputPrinter()  – extract the full density matrix
    .def("__init__", [](ACE::OutputPrinter* self) {
      nb::gil_scoped_release release;
      new (self) ACE::OutputPrinter("", std::vector<Eigen::MatrixXcd>{});
      self->do_extract    = true;
      self->full_densmat  = true;
    })
    .def_rw("do_extract", &ACE::OutputPrinter::do_extract)
    .def("clear",   &ACE::OutputPrinter::clear,
         nb::call_guard<nb::gil_scoped_release>())
    .def("extract", &ACE::OutputPrinter::extract,
         nb::call_guard<nb::gil_scoped_release>())
    ;

  // ── Simulation (Simulation_PT) ─────────────────────────────────────────────
  nb::class_<ACE::Simulation_PT>(m, "Simulation")
    .def(nb::init<>(), nb::call_guard<nb::gil_scoped_release>())
    .def(nb::init<ACE::Parameters&>(), nb::call_guard<nb::gil_scoped_release>())
    // run(prop, PT, InitialState, tgrid, printer)
    .def("run", &ACE::Simulation_PT::run_,
         nb::call_guard<nb::gil_scoped_release>())
    // run(prop, PT, rho_matrix, tgrid, printer)
    .def("run", &ACE::Simulation_PT::run__,
         nb::call_guard<nb::gil_scoped_release>())
    // Simulation(prop, PT, initial, tgrid, printer) – construct-and-run
    .def("__init__",
         [](ACE::Simulation_PT* self,
            ACE::FreePropagator& prop,
            ACE::ProcessTensorForwardList& PT,
            const ACE::InitialState& initial,
            const ACE::TimeGrid& tgrid,
            ACE::OutputPrinter& printer) {
           new (self) ACE::Simulation_PT();
           self->run(prop, PT, initial.rho, tgrid, printer);
         },
         nb::call_guard<nb::gil_scoped_release>())
    .def("__init__",
         [](ACE::Simulation_PT* self,
            ACE::FreePropagator& prop,
            ACE::ProcessTensorForwardList& PT,
            const Eigen::MatrixXcd& rho,
            const ACE::TimeGrid& tgrid,
            ACE::OutputPrinter& printer) {
           new (self) ACE::Simulation_PT();
           self->run(prop, PT, rho, tgrid, printer);
         },
         nb::call_guard<nb::gil_scoped_release>())
    ;

  // ── DynamicalMap ───────────────────────────────────────────────────────────
  nb::class_<ACE::DynamicalMap>(m, "DynamicalMap")
    .def(nb::init<>())
    .def(nb::init<ACE::Parameters&>(), nb::call_guard<nb::gil_scoped_release>())
    .def("__init__",
         [](ACE::DynamicalMap* self,
            ACE::FreePropagator& prop,
            ACE::ProcessTensorForwardList& PT,
            ACE::Simulation_PT& sim,
            const ACE::TimeGrid& tgrid) {
           // Upcast FreePropagator -> Propagator implicitly
           new (self) ACE::DynamicalMap(prop, PT, sim, tgrid);
         },
         nb::call_guard<nb::gil_scoped_release>())
    .def("calculate",
         static_cast<void (ACE::DynamicalMap::*)(ACE::Parameters&)>(
           &ACE::DynamicalMap::calculate),
         nb::call_guard<nb::gil_scoped_release>())
    .def("calculate",
         static_cast<void (ACE::DynamicalMap::*)(
           ACE::Propagator&, ACE::ProcessTensorForwardList&,
           ACE::Simulation_PT&, const ACE::TimeGrid&, int)>(
           &ACE::DynamicalMap::calculate),
         nb::call_guard<nb::gil_scoped_release>())
    .def_rw("E",   &ACE::DynamicalMap::E)
    .def("size",   &ACE::DynamicalMap::size)
    .def("get_dE", [](const ACE::DynamicalMap& self, int i) {
      nb::gil_scoped_release release;
      return self.get_dE(i);
    })
    ;

} // NB_MODULE
