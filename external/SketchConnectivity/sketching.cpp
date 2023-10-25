#include "sketching.h"

#include "detail/util.h"

#include <spdlog/spdlog.h>

namespace sketching {

namespace {

void array_repr(span<const Float> arr, std::stringstream &ss) {
  const auto n = arr.size();
  ss << '[';
  if (n > 0) {
    ss << arr[0];
    for (auto i = decltype(n){1}; i < n; ++i) {
      ss << ", " << arr[i];
    }
  }
  ss << ']';
}

void copy_segment(const Stroke &from, Stroke &to, const double start_index,
                  const double end_index) {
  double from_i_f;
  const auto u = std::modf(start_index, &from_i_f);
  auto from_i = (Index)from_i_f;
  double endf;
  const auto remainder = std::modf(end_index, &endf);

  assert(to.has_time() == from.has_time());
  if (u == 0.0) {
    // Avoid out-of-bounds access if from_i == from.size() - 1.
    to.x(0) = from.x(from_i);
    to.y(0) = from.y(from_i);
    to.width(0) = from.width(from_i);
    if (to.has_time())
      to.time(0) = from.time(from_i);
  } else {
    to.x(0) = (1 - u) * from.x(from_i) + u * from.x(from_i + 1);
    to.y(0) = (1 - u) * from.y(from_i) + u * from.y(from_i + 1);
    to.width(0) = (1 - u) * from.width(from_i) + u * from.width(from_i + 1);
    if (to.has_time())
      to.time(0) = (1 - u) * from.time(from_i) + u * from.time(from_i + 1);
  }

  Index i = 1;
  const auto end = Index(endf - 1);
  for (; from_i <= end; ++i) {
    from_i++;
    to.x(i) = from.x(from_i);
    to.y(i) = from.y(from_i);
    to.width(i) = from.width(from_i);
    if (to.has_time()) {
      to.time(i) = from.time(from_i);
    }
  }

  const auto r = remainder;
  if (r > 0.0) {
    to.x(i) = (1 - r) * from.x(from_i) + r * from.x(from_i + 1);
    to.y(i) = (1 - r) * from.y(from_i) + r * from.y(from_i + 1);
    to.width(i) = (1 - r) * from.width(from_i) + r * from.width(from_i + 1);
    if (to.has_time())
      to.time(i) = (1 - r) * from.time(from_i) + r * from.time(from_i + 1);
    if (to.xy(i).isApprox(to.xy(i - 1))) {
      // Need to avoid creating duplicate vertices when lerping with tiny
      // factors.
      to.resize(to.size() - 1);
    }
  }
}

} // namespace

Stroke::Stroke(Index npoints, bool has_time)
  : memory_(new Float[npoints * (has_time ? 4 : 3)])
  , size_(npoints)
  , capacity_(npoints)
  , x_(memory_.get())
  , y_(&memory_[npoints])
  , width_(&memory_[2 * npoints])
  , time_(has_time ? &memory_[3 * npoints] : nullptr) {}

Stroke::Stroke()
  : memory_(nullptr)
  , size_(0)
  , capacity_(0)
  , x_(nullptr)
  , y_(nullptr)
  , width_(nullptr)
  , time_(nullptr) {}

Stroke Stroke::clone() const {
  auto dst = Stroke(size(), has_time());
  std::copy(x_, x_ + size_, dst.x_);
  std::copy(y_, y_ + size_, dst.y_);
  std::copy(width_, width_ + size_, dst.width_);
  if (has_time())
    std::copy(time_, time_ + size_, dst.time_);
  dst.index = index;
  return dst;
}

void Stroke::compute_arclengths() const {
  const auto n_points = size();
  arclength_.reset(new Float[n_points]);
  if (n_points > 0) {
    auto sofar = Float(0.0);
    arclength_[0] = 0.0;
    bool warned = false;
    for (Index i = 1; i < n_points; ++i) {
      const auto dist =
        std::sqrt(square(x_[i] - x_[i - 1]) + square(y_[i] - y_[i - 1]));
      if (!warned && dist == 0.0) {
        SPDLOG_WARN("Repeated vertex {}", i);
        warned = true;
      }
      sofar += dist;
      arclength_[i] = sofar;
    }
  }
}

void Stroke::invalidate_arclengths() { arclength_.reset(); }

std::pair<Index, Float> Stroke::fractional_index(const Float s) const {
  ensure_arclengths();

  if (s == 0.0) {
    return {0, 0.0};
  } else if (s == length()) {
    return {size_ - 2, 1.0};
  }

  // Essentially copied from VPaint.
  auto i = Index(0);
  auto j = Index(size() - 1);
  auto si = arclength_[i];
  auto sj = arclength_[j];

  while (j - i > 1) {
    // compute an index hopefully close to s
    const auto u = (s - si) / (sj - si);
    auto k = Index(std::floor((1 - u) * i + u * j));
    // make sure i < k < j
    k = std::min(j - 1, std::max(i + 1, k));
    // recurse
    const auto sk = arclength(k);
    if (sk > s) {
      j = k;
      sj = sk;
    } else {
      i = k;
      si = sk;
    }
  }
  assert(j == i + 1);
  assert(arclength(i) - 1e-6 <= s);
  assert(arclength(j) + 1e-6 >= s);
  const auto u = (s - si) / (sj - si);
  assert(u >= -1e-6);
  return {i, u};
}

Float Stroke::avg_sampling() const {
  ensure_arclengths();
  return length() / Float(size() - 1);
}

void Stroke::insert(Index position, Float x, Float y, Float width, Float time) {
  assert(position >= 0 && position <= size_ && "out of bounds");
  if (size_ == capacity_) {
    reserve(std::max(Index(1.5 * capacity_), capacity_ + 2));
  }
  assert(capacity_ > size_);
  for (auto i = size_; i >= position + 1; --i) {
    x_[i] = x_[i - 1];
    y_[i] = y_[i - 1];
    width_[i] = width_[i - 1];
    if (has_time())
      time_[i] = time_[i - 1];
  }
  x_[position] = x;
  y_[position] = y;
  width_[position] = width;
  if (has_time())
    time_[position] = time;
  size_++;
  invalidate_arclengths();
}

void Stroke::push_back(Float x, Float y, Float width, Float time) {
  if (size_ == capacity_) {
    reserve(std::max(Index(1.5 * capacity_), capacity_ + 2));
  }
  assert(capacity_ > size_);
  x_[size_] = x;
  y_[size_] = y;
  width_[size_] = width;
  if (has_time())
    time_[size_] = time;
  size_++;
  invalidate_arclengths();
}

void Stroke::reserve(Index new_capacity) {
  if (new_capacity > capacity_) {
    auto new_memory =
      std::make_unique<Float[]>(new_capacity * (has_time() ? 4 : 3));
    auto *new_x = &new_memory[0];
    std::copy(&x_[0], &x_[size_], new_x);
    auto *new_y = &new_memory[new_capacity];
    std::copy(&y_[0], &y_[size_], new_y);
    auto *new_width = &new_memory[2 * new_capacity];
    std::copy(&width_[0], &width_[size_], new_width);
    auto *new_time = &new_memory[3 * new_capacity];
    if (has_time()) {
      std::copy(&time_[0], &time_[size_], new_time);
    }
    memory_ = std::move(new_memory);
    x_ = new_x;
    y_ = new_y;
    width_ = new_width;
    if (has_time()) {
      time_ = new_time;
    }
    capacity_ = new_capacity;
  }
}

void Stroke::reverse() {
  const auto n2 = Index{size_ / 2};
  for (auto i = 0; i < n2; ++i) {
    std::swap(x_[i], x_[size_ - i - 1]);
  }
  for (auto i = 0; i < n2; ++i) {
    std::swap(y_[i], y_[size_ - i - 1]);
  }
  for (auto i = 0; i < n2; ++i) {
    std::swap(width_[i], width_[size_ - i - 1]);
  }
  if (time_) {
    for (auto i = 0; i < n2; ++i) {
      std::swap(time_[i], time_[size_ - i - 1]);
    }
  }
  invalidate_arclengths();
}

void Stroke::trim(const Index start, const Index stop) {
  assert(start >= 0 && "out of bounds");
  assert(stop <= size_ && "out of bounds");
  assert(start < stop && "invalid range");
  size_ = stop - start;
  if (start == 0) {
    return;
  }
  for (auto i = 0; i < stop - start; ++i) {
    x_[i] = x_[i + start];
  }
  for (auto i = 0; i < stop - start; ++i) {
    y_[i] = y_[i + start];
  }
  for (auto i = 0; i < stop - start; ++i) {
    width_[i] = width_[i + start];
  }
  if (time_) {
    for (auto i = 0; i < stop - start; ++i) {
      time_[i] = time_[i + start];
    }
  }
  invalidate_arclengths();
}

void Stroke::trim(const Float start, const Float stop) {
  assert(start >= 0.0 && "out of bounds");
  assert(start < stop && "invalid range");
  if (start == 0.0 && stop >= length()) {
    return;
  }
  auto [i, u] = fractional_index(start);
  auto dst = 0;
  if (u >= 1.0 - 1e-10) {
    i++;
  } else if (u > 0.0) {
    x_[dst] = lerp(x_, i, u);
    y_[dst] = lerp(y_, i, u);
    width_[dst] = lerp(width_, i, u);
    if (has_time())
      time_[dst] = lerp(time_, i, u);
    i++, dst++;
  }
  for (; i < size_ && arclength_[i] <= stop; ++i, ++dst) {
    x_[dst] = x_[i];
    y_[dst] = y_[i];
    width_[dst] = width_[i];
    if (has_time())
      time_[dst] = time_[i];
  }
  if (i != size_) {
    auto [j, v] = fractional_index(stop);
    if (v > 1e-10) {
      x_[dst] = lerp(x_, j, v);
      y_[dst] = lerp(y_, j, v);
      width_[dst] = lerp(width_, j, v);
      if (has_time())
        time_[dst] = lerp(time_, j, v);
      dst++;
    }
  }
  size_ = dst;
  invalidate_arclengths();
}

Stroke Stroke::fractional_slice(Float start, Float stop) const {
  auto out =
    clone(); // TODO: This is inefficient because it overallocates capacity.
  out.trim(start, stop);
  return out;
}

void Stroke::split(const span<const Float> split_values,
                   std::vector<Stroke> &out_splits) const {
  assert(!(out_splits.data() <= this &&
           this < out_splits.data() + out_splits.size()));
  const auto n_split_values = split_values.size();
#ifndef NDEBUG
  for (size_t i = 0; i < n_split_values; ++i) {
    if (i > 0)
      assert(split_values[i - 1] < split_values[i] &&
             "split values must be sorted and not repeat");
    assert(split_values[i] <= length());
    assert(split_values[i] >= 0);
  }
#endif
  const auto n = size();
  assert(n_split_values >= 2);
  if (n_split_values < 2 || n <= 1 || split_values.front() >= length()) {
    return;
  }

  // Find first vertex.
  auto between = fractional_index(split_values.front());
  if (between.second >= 1.0 - 1e-10) {
    // Avoid creating duplicate vertices through numerical imprecision.
    between.first++;
    between.second = 0.0;
  }
  // `i` will either be at or just before the next split value.
  Index i = between.first;
  if (i >= n - 1) {
    return;
  }
  // Fractional values indicate lerping is needed.
  double start_index = (double)between.first + between.second;

  for (std::size_t split_index = 1; split_index < n_split_values;
       ++split_index) {
    const auto split = split_values[split_index];
    while (i < n - 2 && arclength(i + 1) <= split) {
      i++;
    }
    auto u = (split - arclength(i)) / (arclength(i + 1) - arclength(i));
    const auto end_index = std::min(double(n - 1), i + u);
    auto &new_stroke = out_splits.emplace_back(
      Index(std::ceil(end_index) - std::floor(start_index) + 1), has_time());
    copy_segment(*this, new_stroke, start_index, end_index);
    start_index = end_index;
    if (u >= 1.0 - 1e-10) {
      start_index = ceil(start_index);
    }
  }
}

std::string Stroke::repr() const {
  auto ss = std::stringstream();
  ss.precision(6);
  ss << "Stroke(x=";
  array_repr(x(), ss);
  ss << ", y=";
  array_repr(y(), ss);
  ss << ", width=";
  array_repr(width(), ss);
  if (has_time()) {
    ss << ", time=";
    array_repr(time(), ss);
  }
  ss << ')';
  return ss.str();
}

#if defined(_MSC_VER) && !defined(__clang__)
// In debug mode, MSVC normally doesn't optimize `do { blah; } while(0)` to
// `blah;`.  :(
#pragma optimize("gs", on)
#endif

void copy_run(Stroke &dst, const Index dst_start, const Index dst_end, //
              const Stroke &src, const Index src_start, const Index src_end) {
  assert(std::abs(dst_start - dst_end) == std::abs(src_start - src_end));
  assert(dst.has_time() == src.has_time());

#define LOOP_BODY                                                              \
  do {                                                                         \
    dst.x(i) = src.x(j);                                                       \
    dst.y(i) = src.y(j);                                                       \
    dst.width(i) = src.width(j);                                               \
    if (dst.has_time()) {                                                      \
      dst.time(i) = src.time(j);                                               \
    }                                                                          \
  } while (0)

  if (dst_start <= dst_end && src_start <= src_end) {
    for (Index i = dst_start, j = src_start; i <= dst_end; ++i, ++j) {
      LOOP_BODY;
    }
  } else if (dst_start <= dst_end && src_start > src_end) {
    for (Index i = dst_start, j = src_start; i <= dst_end; ++i, --j) {
      LOOP_BODY;
    }
  } else if (dst_start > dst_end && src_start <= src_end) {
    for (Index i = dst_start, j = src_start; i >= dst_end; --i, ++j) {
      LOOP_BODY;
    }
  } else {
    for (Index i = dst_start, j = src_start; i >= dst_end; --i, --j) {
      LOOP_BODY;
    }
  }

#undef LOOP_BODY
}

#if defined(_MSC_VER) && !defined(__clang__)
#pragma optimize("", on)
#endif

} // namespace sketching
