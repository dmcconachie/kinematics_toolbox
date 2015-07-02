#include "kinematics_toolbox/kinematics_ros_types.h"

using namespace kinematics_ros_types;

geometry_msgs::Twist kinematics_ros_types::operator- (
    geometry_msgs::Twist l, const geometry_msgs::Twist& r )
{
  l.linear.x -= r.linear.x;
  l.linear.y -= r.linear.y;
  l.linear.z -= r.linear.z;
  l.angular.x -= r.angular.x;
  l.angular.y -= r.angular.y;
  l.angular.z -= r.angular.z;

  return l;
}

geometry_msgs::Twist kinematics_ros_types::operator+ (
    geometry_msgs::Twist l, const geometry_msgs::Twist& r )
{
  l.linear.x += r.linear.x;
  l.linear.y += r.linear.y;
  l.linear.z += r.linear.z;
  l.angular.x += r.angular.x;
  l.angular.y += r.angular.y;
  l.angular.z += r.angular.z;

  return l;
}

geometry_msgs::Twist kinematics_ros_types::zeroVelocity()
{
  geometry_msgs::Twist vel;
  vel.linear.x = 0;
  vel.linear.y = 0;
  vel.linear.z = 0;
  vel.angular.x = 0;
  vel.angular.y = 0;
  vel.angular.z = 0;
  return vel;
}

/**
 * \brief Determines difference between two SE(2) poses that are in the same
 * plane with the same reference frame.
 *
 * \return A twist based representation of the difference in pose.
 */
geometry_msgs::Twist kinematics_ros_types::calculateError(
    const geometry_msgs::Pose2D& desired, const geometry_msgs::Pose2D& current )
{
  // this is homogeneous(current).inverse * homogeneous(desired)
  geometry_msgs::Twist diff;

  diff.linear.x = ( desired.x - current.x ) * cos( current.theta )
    + ( desired.y - current.y ) * sin( current.theta );

  diff.linear.y = ( desired.y - current.y ) * cos( current.theta )
    - ( desired.y - current.y ) * sin( current.theta );

  diff.linear.z = 0;
  diff.angular.x = 0;
  diff.angular.y = 0;
  diff.angular.z = current.theta - desired.theta;
  return diff;
}
