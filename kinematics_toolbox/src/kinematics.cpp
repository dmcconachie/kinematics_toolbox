#include "kinematics_toolbox/kinematics.h"

#include <math.h>

using namespace kinematics;

////////////////////////////////////////////////////////////////////////////////
// Skew and Unskew for 3x3 matrices
////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix3d skew(const Eigen::Vector3d w)
{
    Eigen::Matrix3d w_hat;

    w_hat <<     0, -w(2),  w(1),
              w(2),     0, -w(0),
             -w(1),  w(0),     0;
    
    return w_hat;
}

Eigen::Vector3d unskew(const Eigen::Matrix3d &w_hat)
{
    Eigen::Matrix3d w;
    Eigen::Matrix3d w_hat_sym;
    w_hat_sym = ((*w_hat) - (*w_hat).transpose())/2;
    
    w << w_hat_sym(2,1), w_hat_sym(0,2), w_hat_sym(1,0);

    return w;
}


////////////////////////////////////////////////////////////////////////////////
// Twist create, twist hat and unhat
////////////////////////////////////////////////////////////////////////////////

Vector6d createTwist(const Eigen::Vector3d &omega, const Eigen::Vector3d &q)
{
    Vector6d xi;

    if (omega(0) == 0 && omega(1) == 0 && omega(2) == 0)
    {
        xi.segment<3>(0) = q
        xi.segment<3>(3) = omega;
    }
    else
    {
        xi.segment<3>(0) = -1*omega.cross(q);
        xi.segment<3>(3) = omega;
    }
}

Eigen::Matrix4d twistHat(const Vector6d &xi)
{
    Eigen::Matrix4d xi_hat = Eigen::Matrix4d::Zero();

    Eigen::Vector3d v = s.segment<3>(0);
    Eigen::Vector3d w = s.segment<3>(3);
    Eigen::Matrix3d w_hat = skew(w)
    
    xi_hat.block<3,3>(0,0) = w_hat;
    xi_hat.block<3,1>(0,3) = v;
}

Vector6d twistUnhat(const Eigen::Matrix4d &xi_hat)
{
    Vector6d xi;

    Eigen::Vector3d v = xi_hat.block<3,1>(0,3);
    Eigen::Matrix3d w_hat = xi_hat.block<3,3>(0,0);
    Eigen::Vector3d w = unskew(w_hat);
    
    xi.segment<3>(0) = v;
    xi.segment<3>(3) = w;
}

////////////////////////////////////////////////////////////////////////////////
// Adjoints and twist exponentials
////////////////////////////////////////////////////////////////////////////////

Matrix6d adj(const Eigen::Matrix4d &g)
{
    Eigen::Matrix3d R = (*G).block<3,3>(0,0);
    Eigen::Vector3d p = (*G).block<3,1>(0,3);
    Eigen::Matrix3d p_hat;
    skew(&p_hat, &p);
    
    (*adj_G).block<3,3>(0,0) = R;
    (*adj_G).block<3,3>(0,3) = p_hat*R;
    (*adj_G).block<3,3>(3,0) = Eigen::Matrix3d::Zero(3,3);
    (*adj_G).block<3,3>(3,3) = R;
}

Eigen::Matrix3d expmExact(const Eigen::Matrix3d &w_hat, const double theta)
{
    Eigen::Matrix3d eye3 = Eigen::Matrix3d::Identity();
    (*expM) = eye3 + (*w_hat) * sin(theta) + (*w_hat)*(*w_hat)*(1-cos(theta));
}

Eigen::matrix4d expTwist(const Vector6d &twist, const double theta)
{
    Eigen::Vector3d v = (*twist).segment<3>(0);
    Eigen::Vector3d w = (*twist).segment<3>(3);
    
    (*expT) = Eigen::Matrix4d::Zero();
    
    
    if (w(0)==0 && w(1)==0 && w(2)==0)
    {
        (*expT).block<3,1>(0,3) = v;
    }
    else
    {
        Eigen::Matrix3d w_hat;
        skew(&w_hat, &w);
        
        Eigen::Matrix3d exp_w_hat_theta;
        expmExact(&exp_w_hat_theta, &w_hat, theta);
        
        Eigen::Matrix3d eye3 = Eigen::Matrix3d::Identity();
        
        (*expT).block<3,3>(0,0) = exp_w_hat_theta;
        (*expT).block<3,1>(0,3) = (eye3 - exp_w_hat_theta) * w.cross(v) + w*w.transpose()*v*theta;
        (*expT)(3,3) = 1;
    }
}

////////////////////////////////////////////////////////////////////////////////
// Geometric Jacobians
////////////////////////////////////////////////////////////////////////////////

Matrix6Xd bodyJacobian(const std::vector<Vector6d> &twists, const std::vector<double> &theta, const Eigen::Matrix4d &g_theta)
{
    int numTheta = theta->size();
    Matrix6Xd J_s(6,numTheta);
    Matrix6Xd J_b(6,numTheta);
    
    Eigen::Matrix4d g = Eigen::Matrix4d::Identity();

    Vector6d lastTwist;
    Vector6d currentTwist;
    Eigen::Matrix4d expT;
    Matrix6d adjG;
    
    for(int i = 0; i < numTheta; i++)
    {
        if (i == 0)
        {
            currentTwist = twists[i];
            J_s.block<6,1>(0,i) = currentTwist;
        }
        else
        {
            lastTwist = currentTwist;
            currentTwist = twists[i];
            expT = expTwist(lastTwist, theta[i-1]);
            g = g * expT;
            J_s.block<6,1>(0,i) = adj(g).inverse()*currentTwist;
        }
    }
    
    J_b = adj(g_theta).inverse() * J_s;

    return J_b;
}
