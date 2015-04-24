#ifndef KINEMATICS_H
#define KINEMATICS_H

#include <Eigen/Dense>

namespace kinematics
{

typedef Eigen::Matrix<double, 6, 1> Vector6d;
typedef Eigen::Matrix<double, 6, 6> Matrix6d;

Eigen::Matrix3f rot_z(double angle);

void skew(Eigen::Matrix3d *w_hat, Eigen::Vector3d *w);
void unskew(Eigen::Vector3d *w, Eigen::Matrix3d *w_hat);
void twist_hat(Eigen::Matrix4d *s_hat, Vector6d *s);
void twist_unhat(Vector6d *s, Eigen::Matrix4d *s_hat);
void adj(Matrix6d *adj_G, Eigen::Matrix4d *G);
void expm(Eigen::Matrix3d *expM, Eigen::Matrix3d *M);
void expTwist(Eigen::Matrix4d *expT, Vector6d *twist, double theta);

};

#endif
