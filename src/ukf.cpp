#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  // set status to uninitialized
  is_initialized_ = false;
  // create variable for timestamp
  time_us_ = 0;
  // dimensions state vector
  n_x_ = 5;
  // augmented dimensions state vector
  n_aug_ = 7;
  // lambda
  lambda_ = 3 - n_aug_;
  // number of sigma points
  num_sig_pts = 1 + 2 * n_aug_;
  // avoid additional computation


  // set weights
  weights_ = VectorXd(num_sig_pts);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  float weight = 1 / (2 * (lambda_ + n_aug_));
  for(int j = 1; j < num_sig_pts; j++) {
    weights_(j) = weight;
  }

  // don't initialize until we process our first measurement
  is_initialized_ = false;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if(!is_initialized_) {

    // initialize the filter
    P_ << 1, 0, 0, 0, 0,
            0, 1, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;

    previous_timestamp_ = meas_package.timestamp_;

    if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }

    if(meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      VectorXd pos_cartesian = ConvPolarToCart(meas_package.raw_measurements_);
      float radial_speed = meas_package.raw_measurements_[2];
      x_ << pos_cartesian[0], pos_cartesian[1], radial_speed, 0, 0;
    }

    is_initialized_ = true;
    return;
  }

  float dt = (meas_package.timestamp_ - previous_timestamp_)/1000000.0;
  previous_timestamp_ = meas_package.timestamp_;
  // Predict state
  Prediction(dt);

  // Update state with last measurement
  if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    UpdateLidar(meas_package);
  }

  if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    UpdateRadar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  // augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug << x_, 0, 0;

  // augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) << P_;
  P_aug.bottomRightCorner(2, 2) << (std_a_ * std_a_), 0, 0, (std_yawdd_ * std_yawdd_);

  // sigma point matrix
  double scaler = sqrt(lambda_ + n_aug_);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, num_sig_pts);
  MatrixXd A = P_aug.llt().matrixL();

  // compute augmented sigma points
  MatrixXd first_column  = x_aug;
  MatrixXd second_block(n_aug_-1,n_aug_);
  second_block << x_aug, x_aug, x_aug, x_aug, x_aug, x_aug, x_aug;
  second_block += scaler*A;
  MatrixXd third_block(n_aug_-1,n_aug_);
  third_block << x_aug, x_aug, x_aug, x_aug, x_aug, x_aug, x_aug;
  third_block -= scaler*A;
  Xsig_aug << x_aug, second_block, third_block ;

  // predict sigma points
  for(int i=0; i<num_sig_pts;i++){
    double px = Xsig_aug(0,i);
    double py = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yaw = Xsig_aug(6,i);

    float dt2 = delta_t * delta_t;
    float sin_yaw = sin(yaw);
    float cos_yaw = cos(yaw);

    if(fabs(yawd) > 0.001){
      Xsig_pred_(0,i) = px + (v / yawd ) * (sin(yaw + yawd * delta_t) - sin_yaw) + 0.5 * dt2 * cos_yaw * nu_a;
      Xsig_pred_(1,i) = py + (v / yawd ) * (-cos(yaw + yawd * delta_t) + cos_yaw) + 0.5 * dt2 * sin_yaw * nu_a;
    } else{
      Xsig_pred_(0,i) = px + v * cos_yaw * delta_t + 0.5 * dt2 * cos_yaw * nu_a;
      Xsig_pred_(1,i) = py + v * sin_yaw * delta_t + 0.5 * dt2 * sin_yaw * nu_a;
    }
    Xsig_pred_(2,i) = v + nu_a*delta_t;
    Xsig_pred_(3,i) = yaw + yawd*delta_t + 0.5*delta_t*delta_t*nu_yaw;
    Xsig_pred_(4,i) = yawd + nu_yaw*delta_t;
  }

  // predict state mean
  x_.fill(0.0);
  for(int i=0; i<2*n_aug_+1;i++){
    x_ = x_+ weights_(i)*Xsig_pred_.col(i);
  }
  // predict state covariance
  P_.fill(0.0);
  for(int i=0; i<2*n_aug_+1;i++){
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    while (x_diff(3)>M_PI) x_diff(3) -= 2*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3) += 2*M_PI;
    P_ += weights_(i)*x_diff*x_diff.transpose();
  }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */


void UKF::UpdateLidar(MeasurementPackage meas_package) {

  //set measurement dimension, lidar can measure px, py
  int n_z = 2;

  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, num_sig_pts);
  Zsig.fill(0.0);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);

  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);

  // FIXME an interesting problem happens when
  // these for loops are combined.
  // The RMSE in position jumps up by .2
  // I leave them separate for now due to time but
  // I'm definitely curious what's going on

  // transform sigma points into measurement space
  for(int i = 0; i < num_sig_pts; i++) {
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);

    Zsig(0,i) = p_x;
    Zsig(1,i) = p_y;
  }

  // mean prediction
  for(int i = 0; i < num_sig_pts; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // measurement covariance mtx
  for(int i = 0; i < num_sig_pts; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  S(0,0) += std_laspx_ * std_laspx_;
  S(1,1) += std_laspy_ * std_laspy_;

  // cross correlation mtx
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  for(int i = 0; i < num_sig_pts; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // atan2 handles this for us
    x_diff(3) = atan2(sin(x_diff(3)), cos(x_diff(3))) ;
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // calculate Kalman gain K;
  MatrixXd kalman_gain = Tc * S.inverse();

  // set z to raw measurements
  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_diff = z - z_pred;

  // update state mean and covariance matrix
  x_ = x_ + kalman_gain * z_diff;
  P_ = P_ - kalman_gain *  S * kalman_gain.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  // TODO update Radar and Update Lidar do repeat
  // some small tasks that could be reduced
  // to make the code cleaner. ignoring for time.

  // matrix for sigma points
  MatrixXd Zsig = MatrixXd(n_z, num_sig_pts);
  Zsig.fill(0.0);

  // mean prediction
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);

  // measurement covariance mtx
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);

  // transform sigma points into measurement space
  for(int i = 0; i < num_sig_pts; i++) {
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double psi = Xsig_pred_(3,i);

    Zsig(0,i) = sqrt(p_x * p_x + p_y * p_y);
    Zsig(1,i) = atan2(p_y , p_x);
    Zsig(2,i) = (p_x * cos(psi) * v + p_y * sin(psi) * v) / Zsig(0,i);
  }

  // mean predicted measurement
  for(int i = 0; i < num_sig_pts; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // measurement covariance mtx
  for(int i = 0; i < num_sig_pts; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = atan2(sin(z_diff(1)), cos(z_diff(1))) ;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  S(0,0) += std_radr_ * std_radr_;
  S(1,1) += std_radphi_ * std_radphi_;
  S(2,2) += std_radrd_ * std_radrd_;

  // TODO this could also be put into a separate function
  // ignoring because late submission

  // cross correlation mtx
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  for(int i = 0; i < num_sig_pts; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = atan2(sin(z_diff(1)), cos(z_diff(1))) ;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    x_diff(3) = atan2(sin(x_diff(3)), cos(x_diff(3))) ;
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // kalman gain;
  MatrixXd kalman_gain = Tc * S.inverse();

  // real measured values
  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_diff = z - z_pred;

  z_diff(1) = atan2(sin(z_diff(1)), cos(z_diff(1))) ;

  // update state and covariance mtx
  x_ = x_ + kalman_gain * z_diff;
  P_ = P_ - kalman_gain *  S * kalman_gain.transpose();
}


/**
 * converts a vector from polar coordinates to cartesian coordinates
 */
VectorXd UKF::ConvPolarToCart(const VectorXd& x) {
  float rho     = x[0];
  float phi     = x[1];
  float rho_dot = x[2];
  float px      = rho * cos(phi);
  float py      = rho * sin(phi);
  float vx      = rho_dot * cos(phi);
  float vy      = rho_dot * sin(phi);

  VectorXd cartesian(4);
  cartesian << px, py, vx, vy;

  return cartesian;
}
