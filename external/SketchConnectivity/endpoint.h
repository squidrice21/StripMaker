#pragma once

#include <cstddef>
#include <cstdint>
#include <utility>

namespace sketching {

struct Endpoint {
  using IdType = std::int32_t;

  static Endpoint Head(const std::size_t idx) { return Endpoint(IdType(idx)); }

  static Endpoint Tail(const std::size_t idx) { return Endpoint(~IdType(idx)); }

  Endpoint(const std::size_t idx, const bool is_head)
    : m_id(is_head ? IdType(idx) : ~IdType(idx)) {}

  explicit Endpoint(const IdType id)
    : m_id(id) {}

  /**
   * Return the index of the stroke associated with this endpoint.
   */
  [[nodiscard]] std::size_t stroke_idx() const {
    return (m_id >= 0 ? m_id : ~m_id);
  }

  [[nodiscard]] bool is_head() const { return m_id >= 0; }

  [[nodiscard]] bool is_tail() const { return !is_head(); }

  [[nodiscard]] IdType as_int() const { return m_id; }

  bool operator==(const Endpoint &rhs) const { return m_id == rhs.m_id; }
  bool operator!=(const Endpoint &rhs) const { return !operator==(rhs); }
  bool operator<(const Endpoint &rhs) const { return m_id < rhs.m_id; }

  std::pair<int, double> as_pair() const {
    return {(int)stroke_idx(), is_head() ? 0.0 : 1.0};
  }

private:
  IdType m_id;
};

} // namespace sketching

// Make Endpoint hashable.
namespace std {
template <>
struct hash<sketching::Endpoint> {
  size_t operator()(sketching::Endpoint x) const noexcept {
    return std::hash<sketching::Endpoint::IdType>()(x.as_int());
  }
};
} // namespace std
