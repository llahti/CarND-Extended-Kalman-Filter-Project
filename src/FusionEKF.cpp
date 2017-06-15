#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  Hj_ = MatrixXd(3, 4);

  // measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  // measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  // the initial transition matrix F_, Later this will be updated by delta t
  F_ = MatrixXd(4, 4);
  F_ << 1, 0, 1, 0,
        0, 1, 0, 1,
        0, 0, 1, 0,
        0, 0, 0, 1;
  ekf_.F_ = F_;

  // state covariance matrix P
  P_ = MatrixXd(4, 4);
  P_ << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1;
  ekf_.P_ = P_;

  // Process covariance matrix Q. This need to be updated in ProcessMeasurement
  Q_ = MatrixXd(4, 4);
  Q_ << 0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0;
  ekf_.Q_ = Q_;

  // Use these to define noise Q-matrix
  noise_ax = 9;
  noise_ay = 9;

  // Define Identity matrix
  //long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(4, 4);
  ekf_.I_ = I;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}


void FusionEKF::FirstUpdate(const MeasurementPackage &measurement_pack) {

   // Initialize state vector
   VectorXd x_ = VectorXd(4);
   //x_ << 1, 1, 0, 0;

   if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
     /**
     Convert radar from polar to cartesian coordinates and initialize state.
     */
     VectorXd polar = measurement_pack.raw_measurements_;
     x_ =  tools.Polar2Cartesian(polar);
     // Radar do not have enough information to predict vx and vy so initialize them as zero
     x_(2) = 0; x_(3) = 0;
   }
   else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
     /**
     Initialize state.
     */
     //set the state with the initial location and zero velocity, and use default valus for vx and vy
     x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
   }
   previous_timestamp_ = measurement_pack.timestamp_;

   ekf_.x_ = x_;

   // done initializing, no need to predict or update
   is_initialized_ = true;
   return;
}

void FusionEKF::UpdateQ_(float dt)
{
  // Calculate mostly used values of dt for Q matrix
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
  ekf_.Q_ << dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
             0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
             dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
             0, dt_3/2*noise_ay, 0, dt_2*noise_ay;
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
      return FirstUpdate(measurement_pack);
  }

  /*****************************************************************************
   *  Prediction
   *  * Update the state transition matrix F according to the new elapsed time.
   *    - Time is measured in seconds.
   *  * Update the process noise covariance matrix.
   ****************************************************************************/

  //compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	// dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  //Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  //set the process covariance matrix Q
  UpdateQ_(dt);

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   *  * Use the sensor type to perform the update step.
   *  * Update the state and covariance matrices.
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }
}


