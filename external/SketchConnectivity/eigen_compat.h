#pragma once

#include "types.h"

#include <Eigen/Core>

namespace sketching {

// (From pybind11.)
using DynStride = Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>;
template <typename MatrixType>
using DynRef = Eigen::Ref<MatrixType, 0, DynStride>;

inline Vec2 from_eigen(const Eigen::Vector2d& v) { return Vec2(v.x(), v.y()); }

inline Eigen::Vector2d to_eigen(Vec2 v) { return Eigen::Vector2d(v.x_, v.y_); }

inline Eigen::Map<Eigen::VectorXd> as_eigen(span<Float> array) {
  return {array.data(), (Eigen::Index)array.size(), 1};
}

inline Eigen::Map<const Eigen::VectorXd> as_eigen(span<const Float> array) {
  return {array.data(), (Eigen::Index)array.size(), 1};
}

} // namespace sketching
