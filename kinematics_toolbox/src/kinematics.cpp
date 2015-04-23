#include "kinematics_toolbox/kinematics.h"

#include <math.h>

using namespace kinematics;

Eigen::Matrix3f rot_z(double angle)
{
	Eigen::Matrix3f R;

	R << cos(angle), -sin(angle), 0,
		 sin(angle),  cos(angle), 0,
		 0,           0,          1;

	return R;
}
