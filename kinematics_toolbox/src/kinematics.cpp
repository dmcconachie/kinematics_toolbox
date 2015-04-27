#include "kinematics_toolbox/kinematics.h"

#include <cmath>
#include <unsupported/Eigen/MatrixFunctions>
#include <iostream>

using namespace kinematics;

////////////////////////////////////////////////////////////////////////////////
// Skew and Unskew for 3x3 matrices
////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix3d kinematics::skew(const Eigen::Vector3d& w)
{
    Eigen::Matrix3d w_hat;

    w_hat <<     0, -w(2),  w(1),
              w(2),     0, -w(0),
             -w(1),  w(0),     0;
    
    return w_hat;
}

Eigen::Vector3d kinematics::unskew(const Eigen::Matrix3d& w_hat)
{
    Eigen::Vector3d w;
    Eigen::Matrix3d w_hat_sym;
    w_hat_sym = (w_hat - w_hat.transpose())/2;
    
    w << w_hat_sym(2,1), w_hat_sym(0,2), w_hat_sym(1,0);

    return w;
}


////////////////////////////////////////////////////////////////////////////////
// Twist create, twist hat and unhat
////////////////////////////////////////////////////////////////////////////////

Vector6d kinematics::createTwist(const Eigen::Vector3d& omega, const Eigen::Vector3d& q)
{
    Vector6d xi;

    if (omega(0) == 0 && omega(1) == 0 && omega(2) == 0)
    {
        xi.segment<3>(0) = q;
        xi.segment<3>(3) = omega;
    }
    else
    {
        xi.segment<3>(0) = -1*omega.cross(q);
        xi.segment<3>(3) = omega;
    }
    
    return xi;
}

std::vector<Vector6d> kinematics::createTwist(const std::vector<Eigen::Vector3d>& omega,
                                              const std::vector<Eigen::Vector3d>& q)
{
    std::vector<Vector6d> xi(omega.size());

    for (unsigned int i = 0; i < omega.size(); i++)
    {
        xi[i] = createTwist(omega[i], q[i]);
    }

    return xi;
}


std::vector<Vector6d> kinematics::calculateTwists(const Eigen::Matrix4d& g_base,
                                                  const std::vector<Eigen::Vector3d>& omega0,
                                                  const std::vector<Eigen::Vector3d>& q0)
{
  std::vector<Vector6d> xi(q0.size());

  for (unsigned int i = 0; i < q0.size(); i++)
  {
    Eigen::Vector4d omega;
    Eigen::Vector4d q;

    omega << omega0[i], 0;
    q     << q0[i], 1;

    omega = g_base*omega;
    q     = g_base*q;

    xi[i] = kinematics::createTwist(omega.segment<3>(0), q.segment<3>(0));
  }

  return xi;
}

Eigen::Matrix4d kinematics::twistHat(const Vector6d& xi)
{
    Eigen::Matrix4d xi_hat = Eigen::Matrix4d::Zero();

    Eigen::Vector3d v = xi.segment<3>(0);
    Eigen::Vector3d w = xi.segment<3>(3);
    Eigen::Matrix3d w_hat = skew(w);
    
    xi_hat.block<3,3>(0,0) = w_hat;
    xi_hat.block<3,1>(0,3) = v;
    
    return xi_hat;
}

Vector6d kinematics::twistUnhat(const Eigen::Matrix4d& xi_hat)
{
    Vector6d xi;

    Eigen::Vector3d v = xi_hat.block<3,1>(0,3);
    Eigen::Matrix3d w_hat = xi_hat.block<3,3>(0,0);
    Eigen::Vector3d w = unskew(w_hat);
    
    xi.segment<3>(0) = v;
    xi.segment<3>(3) = w;
    
    return xi;
}

////////////////////////////////////////////////////////////////////////////////
// Adjoints and twist exponentials
////////////////////////////////////////////////////////////////////////////////

Matrix6d kinematics::adj(const Eigen::Matrix4d& g)
{
    Eigen::Matrix3d R = g.block<3,3>(0,0);
    Eigen::Vector3d p = g.block<3,1>(0,3);
    Eigen::Matrix3d p_hat = skew(p);
    
    Matrix6d adj_g;
    
    adj_g.block<3,3>(0,0) = R;
    adj_g.block<3,3>(0,3) = p_hat*R;
    adj_g.block<3,3>(3,0) = Eigen::Matrix3d::Zero(3,3);
    adj_g.block<3,3>(3,3) = R;
    
    return adj_g;
}

/*Matrix6d kinematics::adjinv(const Eigen::Matrix4d& g)
{
    Eigen::Matrix3d R = g.block<3,3>(0,0);
    Eigen::Vector3d p = g.block<3,1>(0,3);
    Eigen::Matrix3d p_hat = skew(p);
    
    Matrix6d adjinv_g;
    
    adjinv_g.block<3,3>(0,0) = R.transpose();
    adjinv_g.block<3,3>(0,3) = -1*R.trnaspose()*p_hat;
    adjinv_g.block<3,3>(3,0) = Eigen::Matrix3d::Zero(3,3);
    adjinv_g.block<3,3>(3,3) = R.transpose();
    
    return adjinv_g;
}*/

Eigen::Matrix3d kinematics::expmExact(const Eigen::Matrix3d& w_hat, const double theta)
{
    Eigen::Matrix3d eye3 = Eigen::Matrix3d::Identity();
    
    Eigen::Matrix3d expM;
    expM = eye3 + w_hat*std::sin(theta) + w_hat*w_hat*(1-std::cos(theta));
    
    return expM;
}

Eigen::Matrix4d kinematics::expTwist(const Vector6d& xi, const double theta)
{
    Eigen::Vector3d v = xi.segment<3>(0);
    Eigen::Vector3d w = xi.segment<3>(3);
    
    Eigen::Matrix4d expT = Eigen::Matrix4d::Identity();
    
    if (w(0)==0 && w(1)==0 && w(2)==0)
    {
        expT.block<3,1>(0,3) = v*theta;
    }
    else
    {
        Eigen::Matrix3d w_hat = skew(w);
        Eigen::Matrix3d exp_w_hat_theta = expmExact(w_hat, theta);
        Eigen::Matrix3d eye3 = Eigen::Matrix3d::Identity();
        
        expT.block<3,3>(0,0) = exp_w_hat_theta;
        expT.block<3,1>(0,3) = (eye3 - exp_w_hat_theta) * w.cross(v) + w*w.transpose()*v*theta;
    }
    
    return expT;
}

Eigen::Matrix4d kinematics::expTwist(const std::vector<Vector6d>& xi,
                                     const std::vector<double>& theta)
{
    Eigen::Matrix4d g = Eigen::Matrix4d::Identity();

    for (unsigned int i = 0; i < theta.size(); i++)
    {
      g = g * expTwist(xi[i], theta[i]);
    }

    return g;
}
////////////////////////////////////////////////////////////////////////////////
// Geometric Jacobians
////////////////////////////////////////////////////////////////////////////////

Matrix6Xd kinematics::bodyJacobian(const std::vector<Vector6d>& xi,
                                   const std::vector<double>& theta,
                                   const Eigen::Matrix4d& g_theta)
{
    int numTheta = theta.size();
    Matrix6Xd J_s(6,numTheta);
    Matrix6Xd J_b(6,numTheta);
    
    Eigen::Matrix4d g = Eigen::Matrix4d::Identity();

    Vector6d lastTwist;
    Vector6d currentTwist;
    Eigen::Matrix4d expT;
    
    for(int i = 0; i < numTheta; i++)
    {
        if (i == 0)
        {
            currentTwist = xi[i];
            J_s.block<6,1>(0,i) = currentTwist;
        }
        else
        {
            lastTwist = currentTwist;
            currentTwist = xi[i];
            expT = expTwist(lastTwist, theta[i-1]);
            g = g * expT;
            J_s.block<6,1>(0,i) = adj(g)*currentTwist;
        }
    }
    
    J_b = adj(g_theta.inverse()) * J_s;

    return J_b;
}

////////////////////////////////////////////////////////////////////////////////
// Other
////////////////////////////////////////////////////////////////////////////////

Vector6d kinematics::calculateError(const Eigen::Matrix4d& g_current,
                                    const Eigen::Matrix4d& g_desired)
{
    Vector6d xi;

    Eigen::Matrix4d g_diff = g_current.inverse()*g_desired;

//    std::cerr << "g_diff:\n" << g_diff << std::endl;

    xi = twistUnhat(g_diff.log());

    return xi;
}

