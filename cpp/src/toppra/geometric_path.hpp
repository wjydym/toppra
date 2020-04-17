#ifndef TOPPRA_GEOMETRIC_PATH_HPP
#define TOPPRA_GEOMETRIC_PATH_HPP

#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <toppra/toppra.hpp>
#include <vector>

namespace toppra {
class GeometricPath {};

/**
 * Piecewise polynomial geometric path.
 */
class PiecewisePolyPath : public GeometricPath {
public:
  PiecewisePolyPath(const Matrices &, const std::vector<value_type> &);
  Vector eval(value_type, int degree = 0);
  Vectors eval(std::vector<value_type>, int degree = 0);

  size_t findSegmentIndex(value_type pos) const;
private:
  Matrices m_coefficients, m_coefficients_1, m_coefficients_2;
  std::vector<value_type> m_breakpoints;
  int m_dof, m_degree;
};


} // namespace toppra

#endif
