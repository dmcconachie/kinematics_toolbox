#ifndef KINEMATICS_H
#define KINEMATICS_H

#include <Eigen/Dense>
#include <vector>

namespace kinematics
{

typedef Eigen::Matrix<double, 6, 1> Vector6d;
typedef Eigen::Matrix<double, 6, 6> Matrix6d;
typedef Eigen::Matrix<double, 6, Eigen::Dynamic> Matrix6Xd;

Eigen::Matrix3f rot_z(double angle);

void createTwist(Vector6d *xi, Eigen::Vector3d *omega, Eigen::Vector3d *q);
void skew(Eigen::Matrix3d *w_hat, Eigen::Vector3d *w);
void unskew(Eigen::Vector3d *w, Eigen::Matrix3d *w_hat);
void twist_hat(Eigen::Matrix4d *s_hat, Vector6d *s);
void twist_unhat(Vector6d *s, Eigen::Matrix4d *s_hat);
void adj(Matrix6d *adj_G, Eigen::Matrix4d *G);
void expmExact(Eigen::Matrix3d *expM, Eigen::Matrix3d *w_hat, double theta);
void expTwist(Eigen::Matrix4d *expT, Vector6d *twist, double theta);
void bodyJacobian(std::vector<Vector6d> *twists, std::vector<double> *theta, Eigen::Matrix4d *g_theta);

};

#endif
