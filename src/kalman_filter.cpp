#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd z_predict = H_ * x_;
  VectorXd y = z - z_predict;
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd PH = P_ * H_.transpose();
  MatrixXd K = PH * S.inverse();

  x_ = x_ + (K * y);
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  //recover state parameters
  float pos_x = x_(0);
  float pos_y = x_(1);
  float vel_x = x_(2);
  float vel_y = x_(3);

  // Equations for h_func below
  float h1 = sqrt(pos_x * pos_x + pos_y * pos_y);

  float h2 = atan2(pos_y,pos_x);
  float h3 = (pos_x*vel_x+pos_y*vel_y)/h1;

  //Feed in equations above
  VectorXd H_func(3);
  H_func << h1, h2, h3;

  VectorXd y = z - H_func;

  while (y(1)>M_PI) {
    y(1) -= 2 * M_PI;
  }
  while (y(1)<-M_PI) {
    y(1) += 2 * M_PI;
  }
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd PH = P_ * H_.transpose();
  MatrixXd K = PH * S.inverse();

  x_ = x_ + (K * y);
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  P_ = (I - K * H_) * P_;
}
