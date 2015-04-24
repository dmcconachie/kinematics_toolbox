#ifndef KINEMATICS_H
#define KINEMATICS_H

#include <Eigen/Dense>
#include <vector>

namespace kinematics
{

typedef Eigen::Matrix<double, 6, 1> Vector6d;
typedef Eigen::Matrix<double, 6, 6> Matrix6d;
typedef Eigen::Matrix<double, 6, Eigen::Dynamic> Matrix6Xd;

Eigen::Matrix3d skew(const Eigen::Vector3d &w);
Eigen::Vector3d unskew(const Eigen::Matrix3d &w_hat);

Vector6d createTwist(const Eigen::Vector3d &omega, const Eigen::Vector3d &q);
Eigen::Matrix4d twistHat(const Vector6d &xi);
Vector6d twistUnhat(const Eigen::Matrix4d &xi_hat);

Matrix6d adj(const Eigen::Matrix4d &g);
Eigen::Matrix3d expmExact(const Eigen::Matrix3d &w_hat, const double theta);
Eigen::Matrix4d expTwist(const Vector6d &twist, const double theta);

Matrix6Xd bodyJacobian(const std::vector<Vector6d> &twists, const std::vector<double> &theta, const Eigen::Matrix4d &g_theta);


};

#endif
