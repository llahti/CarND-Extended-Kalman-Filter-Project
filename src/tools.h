#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

  /**
  * A helper method to calculate Jacobians.
  */
  MatrixXd CalculateJacobian(const VectorXd& x_state);

  /**
   * @brief Polar2Cartesian Converts polar coordinates to cartesian coordinates
   * @param polar vector of rho, theta, rho_dot
   * @return  vector of px, py, vx, vy
   */
  VectorXd Polar2Cartesian(VectorXd &polar);

  /**
   * @brief Cartesian2Polar Converts cartesian coordinate and speed vector (pos x, pos y, vel x, vely) to cartesian coordinate system (distance, angle, rad vel)
   * @param cartesian (px, py, vx, vy
   * @return polar (rho, theta, rho_dot)
   */
  VectorXd Cartesian2Polar(VectorXd &cartesian);
};
#endif /* TOOLS_H_ */
