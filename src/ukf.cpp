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

  n_x_ = 5;
  // initial state vector
  x_ = VectorXd(n_x_);

  
  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  
  // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 0.2;//1; //change 6 m/s^2 => 3^2 for vehicle -> 2 => 1 => 1  for bike

  // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 0.4;//0.5 * M_PI;//change 0.2
  
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
  
  
  n_aug_ = 7;
  lambda_ = 3-n_aug_; //change
  Xsig_pred_ =  MatrixXd(n_x_, 2 * n_aug_ + 1);
  
  weights_ = VectorXd(2*n_aug_+1);
    
  // set weights
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }
 
  H_laser_ = MatrixXd(2, 5);
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  
  H_laser_ << 1, 0, 0, 0, 0,
              0, 1, 0, 0, 0;
    //measurement covariance matrix - laser
    
 //
   // std_laspy_ = sqrt(0.0225);
  //  std_laspx_ = sqrt(0.0225);
    R_laser_ << std_laspx_*std_laspx_ , 0,
                0, std_laspy_*std_laspy_;
    
    
    ///* Radar measurement noise standard deviation radius in m
   // std_radr_ = 0.5; //sqrt(0.09);
    ///* Radar measurement noise standard deviation angle in rad
   // std_radphi_ = M_PI / 6 ; //sqrt(0.09)
    ///* Radar measurement noise standard deviation radius change in m/s
   // std_radrd_ = M_PI / 6; //sqrt(0.09)
    //measurement noise covariance matrix - radar
    R_radar_ << std_radr_*std_radr_, 0, 0,
                0, std_radphi_*std_radphi_, 0,
                0, 0, std_radrd_*std_radrd_;

    //create sigma point matrix
    Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);


}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  //Initialize
    
    if (!is_initialized_) {
        // first measurement
        
        float bike_speed_avg = 6.9; //m/s = 25 km/h = 6,9 m/s
        x_ <<  0,0,bike_speed_avg/2,0,0;
        
        MatrixXd P_ = MatrixXd(5, 5);
        /*P_ << 1,0,0,0,0, //0.15*0.15 , 1
              0,1,0,0,0, //0.15*0.15
              0,0,10*10,0,0, //difference in v
              0,0,0,M_PI* M_PI/4, 0, //diffrenece in angle pi/2
              0,0,0,0, M_PI * M_PI/4; // difference in yawn*/
        P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;
        
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_laser_) {
            /**
             Convert radar from polar to cartesian coordinates and initialize state.
             */
            x_(0) = meas_package.raw_measurements_(0) * cos(meas_package.raw_measurements_(1));
            x_(1) = meas_package.raw_measurements_(0) * sin(meas_package.raw_measurements_(1));
            double rho_dot = meas_package.raw_measurements_(2);
            double vx = rho_dot * cos(meas_package.raw_measurements_(1));
            double vy = rho_dot * sin(meas_package.raw_measurements_(1));
            x_(2) = sqrt(vx * vx + vy * vy);

        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_radar_) {
            x_(0) = meas_package.raw_measurements_(0);
            x_(1) = meas_package.raw_measurements_(1);

        }
        
        
        time_us_ = meas_package.timestamp_;
        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }
    
    //prediction
    float delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
    Prediction(delta_t);
    
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        UpdateRadar(meas_package);
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
        UpdateLidar(meas_package);
    }
    time_us_ = meas_package.timestamp_;
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
    //create augmented mean state
    VectorXd x_aug = VectorXd(7);
    
    //create augmented state covariance
    MatrixXd P_aug = MatrixXd(7, 7);
    
    
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
        Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
        Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
    }
    
    //predict sigma points
    
    //predict sigma points
    for (int i = 0; i< 2*n_aug_+1; i++)
    {
        //extract values for better readability
        double p_x = Xsig_aug(0,i);
        double p_y = Xsig_aug(1,i);
        double v = Xsig_aug(2,i);
        double yaw = Xsig_aug(3,i);
        double yawd = Xsig_aug(4,i);
        double nu_a = Xsig_aug(5,i);
        double nu_yawdd = Xsig_aug(6,i);
        
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
    
    
    //predict mean and covariance matrix
    //create vector for predicted state
    
    //predicted state mean
    x_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        x_ = x_+ weights_(i) * Xsig_pred_.col(i);
    }
    
    //predicted state covariance matrix
    P_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        
        P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
    }
    
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

    VectorXd y = meas_package.raw_measurements_ - H_laser_ * x_;
    MatrixXd Ht = H_laser_.transpose();
    MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
    MatrixXd Si = S.inverse();
    MatrixXd K =  P_ * Ht * Si;
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    
    P_ = (I - K * H_laser_) * P_;
    
    //lidar limit = 5.99 (df = 2)
    float nis = (y).transpose()* S.inverse()*(y);
    cout << "nis for lidar (limit 5.99)," << nis << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  
    int n_z = 3;
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    //sigma points into measurement space
    for (int i = 0; i< 2*n_aug_+1; i++)
    {
        //extract values for better readability
        double p_x = Xsig_aug(0,i);
        double p_y = Xsig_aug(1,i);
        double v = Xsig_aug(2,i);
        double yaw = Xsig_aug(3,i);
       
        double v1 = cos(yaw)*v;
        double v2 = sin(yaw)*v;
        
        // measurement model
        Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
        //if (p_x == 0.0 && p_y == 0.0)
        //    return;
        Zsig(1,i) = atan2(p_y,p_x);//phi
        //if (sqrt(p_x*p_x + p_y*p_y) < 0.001) {
        //    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);
        //} else {
        Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
    }
    
    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i=0; i < 2*n_aug_+1; i++) {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }
    
    //innovation covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        
        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
        
        S = S + weights_(i) * z_diff * z_diff.transpose();
    }
    
    S = S + R_radar_;
    
    //Update Radar
    
    VectorXd z = meas_package.raw_measurements_;
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    
    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
        
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        //angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
        
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        
        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }
    
    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();
    
    //residual
    VectorXd z_diff = z - z_pred;
    
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    
    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();
    
    //radar limit = 7.815 (df = 3)
    float nis = (z - z_pred).transpose()* S.inverse()*(z - z_pred);
    cout << "nis for radar (limit 7.815)," << nis << endl;
    
}
