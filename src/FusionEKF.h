#ifndef FusionEKF_H_
#define FusionEKF_H_

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "kalman_filter.h"
#include "tools.h"
#include "cmath"

class FusionEKF {
public:
  /**
  * Constructor.
  */
  FusionEKF();

  /**
  * Destructor.
  */
  virtual ~FusionEKF();

  /**
   * @brief FirstUpdate This method handles the very first update cycle by initializing current state
   * @param measurement_pack
   */
  void FirstUpdate(const MeasurementPackage &measurement_pack);

  /**
  * Run the whole flow of the Kalman Filter from here.
  */
  void ProcessMeasurement(const MeasurementPackage &measurement_pack);

  /**
  * Kalman Filter update and prediction math lives in here.
  */
  KalmanFilter ekf_;



private:
  // check whether the tracking toolbox was initialized or not (first measurement)
  bool is_initialized_;

  // previous timestamp
  long long previous_timestamp_;

  // tool object used to compute Jacobian and RMSE
  Tools tools;

  // Kalman filter matrices
  Eigen::MatrixXd F_;
  Eigen::MatrixXd R_laser_;
  Eigen::MatrixXd R_radar_;
  Eigen::MatrixXd H_laser_;
  Eigen::MatrixXd Hj_;
  Eigen::MatrixXd I;
  Eigen::MatrixXd P_;
  Eigen::MatrixXd Q_;
  Eigen::MatrixXd x_;

  // Use these to define process noise Q-matrix
  float noise_ax;
  float noise_ay;

  /**
   * @brief This function updates EKF's Q_ with new values
   * @param dt time difference in seconds of this and previous update
   */
  void UpdateQ_(float dt);
};

#endif /* FusionEKF_H_ */
