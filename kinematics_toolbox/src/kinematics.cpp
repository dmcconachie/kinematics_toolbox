#include "kinematics_toolbox/kinematics.h"

#include <math.h>
#include <vector>

using namespace kinematics;

Eigen::Matrix3f rot_z(double angle)
{
	Eigen::Matrix3f R;

	R << cos(angle), -sin(angle), 0,
		 sin(angle),  cos(angle), 0,
		 0,           0,          1;

	return R;
}

void bodyJacobian(std::vector<Vector6d> *twists, std::vector<double> *theta, Eigen::Matrix4d *g_theta) {
    
    int numTheta = theta->size();
    Matrix6Xd J_s(6,numTheta);
    Matrix6Xd J_b(6,numTheta);
    Eigen::Matrix4d g = Eigen::Matrix4d::Identity();
    Vector6d lastTwist;
    Vector6d currentTwist;
    Eigen::Matrix4d expT;
    Matrix6d adjG;
    
    for(int i = 0; i < numTheta; i++) {
        if (i == 0) {
            currentTwist = (*twists)[i];
            J_s.block<6,1>(0,i) = currentTwist;
        }
        else {
            lastTwist = currentTwist;
            currentTwist = (*twists)[i];
            expTwist(&expT, &lastTwist, (*theta)[i-1]);
            g = g * expT;
            adj(&adjG, &g);
            J_s.block<6,1>(0,i) = adjG*currentTwist;
        }
    }
    
    adj(&adjG, g_theta);
    J_b = adjG.inverse() * J_s;
    
}

void expTwist(Eigen::Matrix4d *expT, Vector6d *twist, double theta) {
    
    Eigen::Vector3d v = (*twist).segment<3>(0);
    Eigen::Vector3d w = (*twist).segment<3>(3);
    
    (*expT) = Eigen::Matrix4d::Zero();
    
    
    if (w(0)==0 && w(1)==0 && w(2)==0) {
        (*expT).block<3,1>(0,3) = v;
    }
    else {
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

void expmExact(Eigen::Matrix3d *expM, Eigen::Matrix3d *w_hat, double theta) {
    
    Eigen::Matrix3d eye3 = Eigen::Matrix3d::Identity();
    (*expM) = eye3 + (*w_hat) * sin(theta) + (*w_hat)*(*w_hat)*(1-cos(theta));
}

void adj(Matrix6d *adj_G, Eigen::Matrix4d *G) {
    
    Eigen::Matrix3d R = (*G).block<3,3>(0,0);
    Eigen::Vector3d p = (*G).block<3,1>(0,3);
    Eigen::Matrix3d p_hat;
    skew(&p_hat, &p);
    
    (*adj_G).block<3,3>(0,0) = R;
    (*adj_G).block<3,3>(0,3) = p_hat*R;
    (*adj_G).block<3,3>(3,0) = Eigen::Matrix3d::Zero(3,3);
    (*adj_G).block<3,3>(3,3) = R;
    
}

void twist_hat(Eigen::Matrix4d *s_hat, Vector6d *s) {
    
    Eigen::Vector3d v = (*s).segment<3>(0);
    Eigen::Vector3d w = (*s).segment<3>(3);
    Eigen::Matrix3d w_hat;
    skew(&w_hat, &w);
    
    (*s_hat) = Eigen::Matrix4d::Zero(4,4);
    (*s_hat).block<3,3>(0,0) = w_hat;
    (*s_hat).block<3,1>(0,3) = v;
    
}

void twist_unhat(Vector6d *s, Eigen::Matrix4d *s_hat) {
    
    Eigen::Vector3d v = (*s_hat).block<3,1>(0,3);
    Eigen::Matrix3d w_hat = (*s_hat).block<3,3>(0,0);
    Eigen::Vector3d w;
    unskew(&w, &w_hat);
    
    (*s).segment<3>(0) = v;
    (*s).segment<3>(3) = w;
}

void skew(Eigen::Matrix3d *w_hat, Eigen::Vector3d *w) {
    
    (*w_hat) <<   0,    -(*w)(2),   (*w)(1),
               (*w)(2),     0,     -(*w)(0),
              -(*w)(1),  (*w)(0),       0;
    
}

void unskew(Eigen::Vector3d *w, Eigen::Matrix3d *w_hat) {
    
    Eigen::Matrix3d w_hat_sym;
    w_hat_sym = ((*w_hat) - (*w_hat).transpose())/2;
    (*w) << w_hat_sym(2,1), w_hat_sym(0,2), w_hat_sym(1,0);
}

void createTwist(Vector6d *xi, Eigen::Vector3d *omega, Eigen::Vector3d *q) {
    
    if ((*omega)(0) == 0 && (*omega)(1) == 0 && (*omega)(2) == 0) {
        (*xi).segment<3>(0) = (*q);
        (*xi).segment<3>(3) = (*omega);
    }
    else {
        (*xi).segment<3>(0) = -1*(*omega).cross(*q);
        (*xi).segment<3>(3) = (*omega);
    }
    
}



