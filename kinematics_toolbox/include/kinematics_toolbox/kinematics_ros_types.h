#ifndef KINEMATICS_ROS_TYPES_H
#define KINEMATICS_ROS_TYPES_H

#include <geometry_msgs/Pose2D.h>
#include <geometry_msgs/Twist.h>

namespace kinematics_ros_types
{
  geometry_msgs::Twist operator- ( geometry_msgs::Twist l, const geometry_msgs::Twist& r );
  geometry_msgs::Twist operator+ ( geometry_msgs::Twist l, const geometry_msgs::Twist& r );

  geometry_msgs::Twist zeroVelocity();

  /**
   * \brief Determines difference between two SE(2) poses that are in the same
   * plane with the same reference frame.
   *
   * \return A body velocity in the current frame that moves us from \a current
   * to \a desired.
   */
  geometry_msgs::Twist calculateError( const geometry_msgs::Pose2D& desired,
                                       const geometry_msgs::Pose2D& current );

}

#endif
