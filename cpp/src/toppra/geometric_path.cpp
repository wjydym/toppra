#include "toppra/geometric_path.hpp"
#include "toppra/toppra.hpp"
#include <vector>
namespace toppra {

Vectors GeometricPath::eval(const Vector &positions, int order) const {
  Vectors outputs;
    outputs.resize(positions.size());
    for (size_t i = 0; i < positions.size(); i++) {
      outputs[i] = eval_single(positions(i), order);
    }
    return outputs;
  };


Vector GeometricPath::proposeGridpoints(double maxErrThreshold, int maxIteration, double maxSegLength, int minNbPoints) const {
  std::vector<value_type> gridpoints_vec;
  gridpoints_vec.push_back(pathInterval()[0]);
  gridpoints_vec.push_back(pathInterval()[1]);

  // Add points according to error threshold
  for (auto i=0; i < maxIteration; i++){
    bool add_new_points = false;
    auto current_size = gridpoints_vec.size();
    for (auto j=0; j < current_size - 1; j++){

      value_type p_mid = 0.5 * (gridpoints_vec[j] + gridpoints_vec[j + 1]);
      auto dist = gridpoints_vec[j + 1] - gridpoints_vec[j];

      if (dist > maxSegLength){
        gridpoints_vec.push_back(p_mid);
        add_new_points = true;
        continue;
      }

      auto max_err = (0.5 * eval_single(p_mid, 2) * dist * dist).cwiseAbs().maxCoeff();
      if (max_err > maxErrThreshold){
        gridpoints_vec.push_back(p_mid);
        add_new_points = true;
        continue;
      }
    }

    if (!add_new_points) break;
    sort(gridpoints_vec.begin(), gridpoints_vec.end());
  }

  // Add points according to smallest number of points
  while (gridpoints_vec.size() < minNbPoints){

    auto current_size = gridpoints_vec.size();
    for (auto j=0; j < current_size - 1; j++){
      value_type p_mid = 0.5 * (gridpoints_vec[j] + gridpoints_vec[j + 1]);
      gridpoints_vec.push_back(p_mid);
    }

    sort(gridpoints_vec.begin(), gridpoints_vec.end());
  }

  // Return the Eigen vector
  Vector gridpoints(gridpoints_vec.size());
  for (int i=0; i < gridpoints_vec.size(); i++) gridpoints(i) = gridpoints_vec[i];
  return gridpoints;
}
}