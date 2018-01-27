#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check for errors/invalid values in the inputs:
  // estimations and ground_truth must have equal size
  // estimations can't be zero
  if(estimations.size() != ground_truth.size()
     || estimations.size() == 0){
    std::cout << "Error in RMSE data data" << std::endl;
    return rmse;
  }

  //sum squared errors
  for(unsigned int i=0; i < estimations.size(); ++i){

    VectorXd sum_errors = estimations[i] - ground_truth[i];

    //coefficient-wise multiplication
    sum_errors = sum_errors.array()*sum_errors.array();
    rmse += sum_errors;
  }

  //calculate mean and root
  rmse = (rmse/estimations.size()).array().sqrt();

  //return the result
  return rmse;
}
