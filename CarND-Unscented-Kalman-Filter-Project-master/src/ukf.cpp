#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
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



  /* change ******************************//////////////////////
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.8;

 // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.7;
  

/* change ******************************//////////////////////

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  H_laser = MatrixXd(2, 5);
  R_laser = MatrixXd(2, 2);
  
  H_laser << 1,0,0,0,0,
       0,1,0,0,0;

  
  R_laser << std_laspx_, 0,
       0, std_laspy_;

  n_aug_=7;
  n_x_= 5;
  repeatedEq = 2* n_aug_ +1;

  Xsig_aug = MatrixXd(n_aug_, repeatedEq);
  x_aug = VectorXd(7);
  Xsig_pred_=MatrixXd(n_x_, repeatedEq);
  P_aug = MatrixXd(7, 7);
  lambda_ = 3 - n_aug_;
  
  n_z = 3;
  z_pred = VectorXd(n_z);
  S = MatrixXd(n_z,n_z);
  //create matrix for sigma points in measurement space
  Zsig = MatrixXd(n_z, repeatedEq);

  //create matrix for cross correlation Tc
  Tc = MatrixXd(n_x_, n_z);

  weights_= VectorXd(repeatedEq);
  R_radar = MatrixXd(n_z,n_z);
  R_radar <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
   
  is_initialized_=false;
  
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
  if (!is_initialized_) {
   x_<<  0,0,0,0,0;

   P_<< 0.2,0,0,0,0,
        0,0.2,0,0,0,
        0,0,0.2,0,0,
        0,0,0,0.2,0,
        0,0,0,0,0.2;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
    float rho = meas_package.raw_measurements_[0]; // range
    float phi = meas_package.raw_measurements_[1]; // bearing
    float rho_dot = meas_package.raw_measurements_[2]; // velocity of rho
    // Coordinates convertion from polar to cartesian
    float x = rho * cos(phi); 
    float y = rho * sin(phi);
    x_(0)=x;
    x_(1)=y;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_(0)= meas_package.raw_measurements_(0);
      x_(1)= meas_package.raw_measurements_(1);
    }
    time_us_=meas_package.timestamp_;


    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
float dt= (meas_package.timestamp_ - time_us_)/1000000.0;
time_us_=meas_package.timestamp_;
Prediction(dt);
if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
  UpdateRadar(meas_package);
}
else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
  UpdateLidar(meas_package);
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

  
 
  //1- generate and augmente sigma points
  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }


  //2-predict sigma points
  for (int i = 0; i< repeatedEq; i++)
  {
    //extract values for better readability
    const double p_x = Xsig_aug(0,i);
    const double p_y = Xsig_aug(1,i);
    const double v = Xsig_aug(2,i);
    const double yaw = Xsig_aug(3,i);
    const double yawd = Xsig_aug(4,i);
    const double nu_a = Xsig_aug(5,i);
    const double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }


    //3- the predicted mean and covariance
    // set weights
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<repeatedEq; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }
  x_.fill(0.0);
  //predicted state mean
  for (int i = 0; i < repeatedEq; i++) {  //iterate over sigma points
    x_ =x_+weights_(i) * Xsig_pred_.col(i);
  }
P_.fill(0.0);
  //predicted state covariance matrix
VectorXd x_diff;
  RepeatedLoop(x_diff,P_,Xsig_pred_,x_,x_diff,3);
  

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  

  
  VectorXd z_pred = H_laser * x_;
  VectorXd y = meas_package.raw_measurements_ - z_pred;
  MatrixXd Ht = H_laser.transpose();
  MatrixXd PHt = P_ * Ht;
  MatrixXd S = H_laser * PHt + R_laser;
  MatrixXd Si = S.inverse();
  MatrixXd K = PHt * Si;
  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_laser) * P_;
  
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  //set measurement dimension, radar can measure r, phi, and r_dot
  
  
  //1- predict radar measurement

  //transform sigma points into measurement space
  for (int i = 0; i < repeatedEq; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  
  z_pred.fill(0.0);
  for (int i=0; i < repeatedEq; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //innovation covariance matrix S
  
  S.fill(0.0);
  VectorXd z_diff_temp;
  RepeatedLoop(z_diff_temp,S,Zsig,z_pred,z_diff_temp,1);
  
  //add measurement noise covariance matrix

  S = S + R_radar; 
  //2- Update


  //calculate cross correlation matrix
  Tc.fill(0.0);
  //VectorXd x_diff;
  //RepeatedLoop(x_diff,Tc,Xsig_pred_,x_,z_diff_temp,3);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    NormalizeAngle(z_diff(1));
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    NormalizeAngle(x_diff(3));
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

  //angle normalization
  NormalizeAngle(z_diff(1));
  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();


}
void UKF::NormalizeAngle(double &angle )
{
  while (angle> M_PI) angle-=2.*M_PI;
  while (angle<-M_PI) angle+=2.*M_PI;
}

void UKF::RepeatedLoop(VectorXd &input,MatrixXd &ouput, MatrixXd param1, VectorXd param2,VectorXd &param3, int index)
{

for (int i = 0; i < repeatedEq; i++) {  
    input = param1.col(i) - param2;
    NormalizeAngle(input(index));
    ouput = ouput+weights_(i) * input * param3.transpose() ;
  }
}