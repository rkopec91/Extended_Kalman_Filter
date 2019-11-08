#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  for(unsigned int i=0; i < estimations.size(); ++i){
    VectorXd residuals = estimations[i] - ground_truth[i];
    residuals = residuals.array() * residuals.array();
    rmse += residuals;
  }

  rmse = rmse / estimations.size();

  rmse = sqrt(rmse.array());

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3,4);

  float pos_x = x_state(0);
  float pos_y = x_state(1);
  float vel_x = x_state(2);
  float vel_y = x_state(3);

  float h1 = pos_x * pos_x + pos_y * pos_y;

  float h2 = sqrt(h1);
  float h3 = h1 * h2;

  Hj << pos_x/h2, pos_y/h2, 0, 0,
        -pos_y/h1, pos_x/h1, 0, 0,
        (pos_y*(vel_x*pos_y - vel_y*pos_x))/h3, (pos_x*(vel_y*pos_x - vel_x*pos_y))/h3, pos_x/h2, pos_y/h2;

  return Hj;
}
