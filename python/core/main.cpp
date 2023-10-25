#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals;

// Everytime add a new h-cpp pair, add the individual module intialization here.
void init_input(py::module &);
void init_modules(py::module &m) {
  init_input(m);
}

PYBIND11_MODULE(_sketching, m) {
  m.doc() = "StripMaker Python Bindings";

  init_modules(m);
}
