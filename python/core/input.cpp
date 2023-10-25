#include "glm/detail/type_vec.hpp"
#include <../src/Util.h>
#include <../stroke_strip_src/Cluster.h>

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>

namespace py = pybind11;
using namespace pybind11::literals;

namespace pybind11 {
namespace detail {
template <>
struct type_caster<glm::dvec2> {
public:
  PYBIND11_TYPE_CASTER(glm::dvec2, _("dvec2"));

  // Python -> C++
  bool load(py::handle src, bool convert) {
    if (!convert && !py::array_t<double>::check_(src))
      return false;

    auto buf =
      py::array_t<double, py::array::c_style | py::array::forcecast>::ensure(
        src);
    if (!buf || buf.ndim() != 1 || buf.size() != 2)
      return false;

    value = glm::dvec2(buf.data()[0], buf.data()[1]);

    return true;
  }

  // C++ -> Python
  static py::handle cast(const glm::dvec2 &src,
                         py::return_value_policy /*policy*/,
                         py::handle /*parent*/) {
    auto a = py::array(py::buffer_info(
      nullptr, /* Pointer to data (nullptr -> ask NumPy to allocate!) */
      sizeof(double), /* Size of one item */
      py::format_descriptor<double>::value, /* Buffer format */
      1, /* How many dimensions? */
      {2}, /* Number of elements for each dimension */
      {sizeof(double)} /* Strides for each dimension */
      ));
    auto buf = a.request();
    double *ptr = (double *)buf.ptr;
    ptr[0] = src.x;
    ptr[1] = src.y;

    return a.release();
  }
};
} // namespace detail
} // namespace pybind11

void init_input(py::module &m) {
  py::class_<Cluster::Stroke>(m, "Stroke")
    .def_readwrite("cluster_ind", &Cluster::Stroke::cluster_ind)
    .def_readwrite("stroke_ind", &Cluster::Stroke::stroke_ind)
    .def_readwrite("points", &Cluster::Stroke::points)
    .def_readwrite("u", &Cluster::Stroke::u);
  py::class_<FittedCurve>(m, "FittedCurve")
    .def_readwrite("centerline", &FittedCurve::centerline)
    .def_readwrite("widths", &FittedCurve::widths);

  py::class_<Cluster>(m, "Cluster")
    .def_readwrite("strokes", &Cluster::strokes)
    .def_readwrite("periodic", &Cluster::periodic)
    .def_readwrite("fit", &Cluster::fit);

  py::class_<Input>(m, "Input")
    .def_readwrite("thickness", &Input::thickness)
    .def_readwrite("orig_thickness", &Input::orig_thickness)
    .def_readwrite("orig_center", &Input::orig_center)
    .def_readwrite("width", &Input::width)
    .def_readwrite("height", &Input::height)
    .def_readwrite("clusters", &Input::clusters)
    .def("__repr__", &Input::repr);

  m.def("read_input", [](const std::string &scap_filename) {
    int width, height;
    double input_thickness;
    Input input;
    bool to_preprocess = false;
    read_input(scap_filename, input, width, height, input_thickness,
               to_preprocess);
    return input;
  });

  m.def(
    "save_fitting",
    [](const std::string &scap_filename, const std::string &output_filename,
       const bool fit_width, const std::string &vis_folder,
       const std::string &cut_filename, const bool fit_single,
       const bool to_spiral, const bool disable_cut_orientation) {
      save_fitting(scap_filename, output_filename, fit_width, vis_folder,
                   cut_filename, fit_single, to_spiral,
                   disable_cut_orientation);
    },
    py::arg("scap_filename"), py::arg("output_filename"),
    py::arg("fit_width") = false, py::arg("vis_folder") = "",
    py::arg("cut_filename") = "", py::arg("fit_single") = false,
    py::arg("to_spiral") = false, py::arg("disable_cut_orientation") = false);

  m.def("example_hash", [](const std::vector<size_t> &stroke_indices1,
                           const std::vector<size_t> &stroke_indices2) {
    return example_hash(stroke_indices1, stroke_indices2);
  });
}
