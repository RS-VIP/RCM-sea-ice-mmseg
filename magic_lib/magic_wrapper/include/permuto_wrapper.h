//
// Created by max on 2021-11-17.
//

#ifndef MAGIC_PORT_PERMUTO_WRAPPER_H
#define MAGIC_PORT_PERMUTO_WRAPPER_H

#include "permutohedral.h"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

py::array_t<float> permuto_gauss_transform(py::array_t<float> sources, py::array_t<float>targets, py::array_t<float> weights, py::array_t<float> bandwidths);

#endif //MAGIC_PORT_PERMUTO_WRAPPER_H
