#ifndef KINEMATICS_H
#define KINEMATICS_H

#include <Eigen/Core>
#include <vector>

namespace kinematics
{

typedef Eigen::Matrix< double, 6, 1 > Vector6d;
typedef Eigen::Matrix< double, 6, 6 > Matrix6d;
typedef Eigen::Matrix< double, 6, Eigen::Dynamic > Matrix6Xd;

typedef struct {
  Vector6d vel;
  bool valid;
} Velocity6d;

////////////////////////////////////////////////////////////////////////////////
// Basic Transforms
////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix4d rotX( const double theta );
Eigen::Matrix4d rotY( const double theta );
Eigen::Matrix4d rotZ( const double theta );
Eigen::Matrix4d translate( const Eigen::Vector3d& p );
Eigen::Matrix4d transX( const double x );
Eigen::Matrix4d transY( const double y );
Eigen::Matrix4d transZ( const double z );

////////////////////////////////////////////////////////////////////////////////
// Skew and Unskew for 3x3 matrices
////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix3d skew( const Eigen::Vector3d& w );
Eigen::Vector3d unskew( const Eigen::Matrix3d& w_hat );

////////////////////////////////////////////////////////////////////////////////
// Twist create, cacluting twists, twist hat and unhat
////////////////////////////////////////////////////////////////////////////////

Vector6d createTwist( const Eigen::Vector3d& omega, const Eigen::Vector3d& q );
std::vector<Vector6d> createTwist( const std::vector< Eigen::Vector3d >& omega,
                                   const std::vector< Eigen::Vector3d >& q );
std::vector<Vector6d> calculateTwists( const Eigen::Matrix4d& g_base,
                                       const std::vector< Eigen::Vector3d >& omega0,
                                       const std::vector< Eigen::Vector3d >& q0 );
Eigen::Matrix4d twistHat( const Vector6d& xi );
Vector6d twistUnhat( const Eigen::Matrix4d& xi_hat );

////////////////////////////////////////////////////////////////////////////////
// Adjoints and twist exponentials
////////////////////////////////////////////////////////////////////////////////

Matrix6d adj( const Eigen::Matrix4d& g );
Eigen::Matrix3d expmExact( const Eigen::Matrix3d& w_hat, const double theta );
Eigen::Matrix4d expTwist( const Vector6d& xi, const double theta );
Eigen::Matrix4d expTwist( const std::vector< Vector6d >& xi,
                          const std::vector< double >& theta );

////////////////////////////////////////////////////////////////////////////////
// Geometric Jacobians
////////////////////////////////////////////////////////////////////////////////

Matrix6Xd spatialJacobian( const std::vector< Vector6d >& xi,
                           const std::vector< double >& theta,
                           const Eigen::Matrix4d& g_theta );
Matrix6Xd bodyJacobian( const std::vector< Vector6d >& xi,
                        const std::vector< double >& theta,
                        const Eigen::Matrix4d& g_zero );

////////////////////////////////////////////////////////////////////////////////
// Other
////////////////////////////////////////////////////////////////////////////////

Vector6d calculateError( const Eigen::Matrix4d& g_current,
                         const Eigen::Matrix4d& g_desired );

};

#endif
