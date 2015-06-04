#include "kinematics_toolbox/kinematics.h"

using namespace kinematics;

////////////////////////////////////////////////////////////////////////////////
// Basic Transforms
////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix4d kinematics::rotX( const double theta )
{
    Eigen::Matrix4d R;
    R << 1,                  0,                  0, 0,
         0,  std::cos( theta ), -std::sin( theta ), 0,
         0,  std::sin( theta ),  std::cos( theta ), 0,
         0,                  0,                  0, 1;
    return R;
}

Eigen::Matrix4d kinematics::rotY( const double theta )
{
    Eigen::Matrix4d R;
    R <<  std::cos( theta ), 0, std::sin( theta ), 0,
                          0, 1,                 0, 0,
         -std::sin( theta ), 0, std::cos( theta ), 0,
                          0, 0,                 0, 1;
    return R;
}

Eigen::Matrix4d kinematics::rotZ( const double theta )
{
    Eigen::Matrix4d R;
    R << std::cos( theta ), -std::sin( theta ), 0, 0,
         std::sin( theta ),  std::cos( theta ), 0, 0,
                         0,                  0, 1, 0,
                         0,                  0, 0, 1;
    return R;
}

Eigen::Matrix4d kinematics::translate( const Eigen::Vector3d& p )
{
    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    T.block< 3, 1 >( 0, 3 ) = p;
    return T;
}

Eigen::Matrix4d kinematics::transX( const double x )
{
    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    T( 0, 3 ) = x;
    return T;
}

Eigen::Matrix4d kinematics::transY( const double y )
{
    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    T( 1, 3 ) = y;
    return T;
}

Eigen::Matrix4d kinematics::transZ( const double z )
{
    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    T( 2, 3 ) = z;
    return T;
}

////////////////////////////////////////////////////////////////////////////////
// Skew and Unskew for 3x3 matrices
////////////////////////////////////////////////////////////////////////////////

Eigen::Matrix3d kinematics::skew( const Eigen::Vector3d& w )
{
    Eigen::Matrix3d w_hat;

    w_hat <<       0, -w( 2 ),  w( 1 ),
              w( 2 ),       0, -w( 0 ),
             -w( 1 ),  w( 0 ),       0;
    
    return w_hat;
}

Eigen::Vector3d kinematics::unskew( const Eigen::Matrix3d& w_hat )
{
    Eigen::Vector3d w;
    Eigen::Matrix3d w_hat_sym;
    w_hat_sym = ( w_hat - w_hat.transpose() )/2;
    
    w << w_hat_sym( 2, 1 ), w_hat_sym( 0, 2 ), w_hat_sym( 1, 0 );

    return w;
}


////////////////////////////////////////////////////////////////////////////////
// Twist create, cacluting twists, twist hat and unhat
////////////////////////////////////////////////////////////////////////////////

Vector6d kinematics::createTwist( const Eigen::Vector3d& omega, const Eigen::Vector3d& q )
{
    Vector6d xi;

    if ( omega( 0 ) == 0 && omega( 1 ) == 0 && omega( 2 ) == 0 )
    {
        xi.segment< 3 >( 0 ) = q;
        xi.segment< 3 >( 3 ) = omega;
    }
    else
    {
        xi.segment< 3 >( 0 ) = -omega.cross( q );
        xi.segment< 3 >( 3 ) = omega;
    }
    
    return xi;
}

std::vector<Vector6d> kinematics::createTwist( const std::vector< Eigen::Vector3d >& omega,
                                               const std::vector< Eigen::Vector3d >& q )
{
    std::vector<Vector6d> xi( omega.size() );

    for ( unsigned int i = 0; i < omega.size(); i++ )
    {
        xi[i] = createTwist( omega[i], q[i] );
    }

    return xi;
}


std::vector<Vector6d> kinematics::calculateTwists( const Eigen::Matrix4d& g_base,
                                                   const std::vector< Eigen::Vector3d >& omega0,
                                                   const std::vector< Eigen::Vector3d >& q0 )
{
    // TODO: is this the same thing as adj( g_base )*twist?
    std::vector<Vector6d> xi( q0.size() );

    for ( unsigned int i = 0; i < q0.size(); i++ )
    {
        Eigen::Vector4d omega;
        Eigen::Vector4d q;

        omega << omega0[i], 0;
        q     << q0[i], 1;

        omega = g_base*omega;
        q     = g_base*q;

        xi[i] = kinematics::createTwist( omega.segment<3>( 0 ), q.segment<3>( 0 ) );
    }

    return xi;
}

Eigen::Matrix4d kinematics::twistHat( const Vector6d& xi )
{
    Eigen::Matrix4d xi_hat = Eigen::Matrix4d::Zero();

    Eigen::Vector3d v = xi.segment< 3 >( 0 );
    Eigen::Vector3d w = xi.segment< 3 >( 3 );
    Eigen::Matrix3d w_hat = skew( w );
    
    xi_hat.block< 3, 3>( 0, 0 ) = w_hat;
    xi_hat.block< 3, 1>( 0, 3 ) = v;
    
    return xi_hat;
}

Vector6d kinematics::twistUnhat( const Eigen::Matrix4d& xi_hat )
{
    Vector6d xi;

    Eigen::Vector3d v = xi_hat.block< 3, 1>( 0, 3 );
    Eigen::Matrix3d w_hat = xi_hat.block< 3, 3>( 0, 0 );
    Eigen::Vector3d w = unskew( w_hat );
    
    xi.segment< 3 >( 0 ) = v;
    xi.segment< 3 >( 3 ) = w;
    
    return xi;
}

////////////////////////////////////////////////////////////////////////////////
// Adjoints and twist exponentials
////////////////////////////////////////////////////////////////////////////////

Matrix6d kinematics::adj( const Eigen::Matrix4d& g )
{
    Eigen::Matrix3d R = g.block< 3, 3>( 0, 0 );
    Eigen::Vector3d p = g.block< 3, 1>( 0, 3 );
    Eigen::Matrix3d p_hat = skew( p );
    
    Matrix6d adj_g;
    
    adj_g.block< 3, 3>( 0, 0 ) = R;
    adj_g.block< 3, 3>( 0, 3 ) = p_hat*R;
    adj_g.block< 3, 3>( 3, 0 ) = Eigen::Matrix3d::Zero( 3, 3 );
    adj_g.block< 3, 3>( 3, 3 ) = R;
    
    return adj_g;
}

/*Matrix6d kinematics::adjinv( const Eigen::Matrix4d& g )
{
    Eigen::Matrix3d R = g.block< 3, 3>( 0, 0 );
    Eigen::Vector3d p = g.block< 3, 1>( 0, 3 );
    Eigen::Matrix3d p_hat = skew( p );
    
    Matrix6d adjinv_g;
    
    adjinv_g.block< 3, 3 >( 0, 0 ) = R.transpose();
    adjinv_g.block< 3, 3 >( 0, 3 ) = -R.transpose()*p_hat;
    adjinv_g.block< 3, 3 >( 3, 0 ) = Eigen::Matrix3d::Zero( 3, 3 );
    adjinv_g.block< 3, 3 >( 3, 3 ) = R.transpose();
    
    return adjinv_g;
}*/

Eigen::Matrix3d kinematics::expmExact( const Eigen::Matrix3d& w_hat, const double theta )
{
    Eigen::Matrix3d eye3 = Eigen::Matrix3d::Identity();

    // w_hat should be normalized before calling this function
    assert( abs(unskew( w_hat ).norm() - 1) < 1e-10 );
    
    Eigen::Matrix3d expM;
    expM = eye3 
           + w_hat * std::sin( theta )
           + w_hat * w_hat * ( 1 - std::cos( theta ) );
    
    return expM;
}

Eigen::Matrix4d kinematics::expTwist( const Vector6d& xi, double theta )
{
    Eigen::Vector3d v = xi.segment< 3 >( 0 );
    Eigen::Vector3d w = xi.segment< 3 >( 3 );
    
    Eigen::Matrix4d expT = Eigen::Matrix4d::Identity();
    
    if ( w( 0 ) == 0 && w( 1 ) == 0 && w( 2 ) == 0 )
    {
        expT.block< 3, 1 >( 0, 3 ) = v*theta;
    }
    else
    {
        double norm_w = w.norm();

        // rescale to make norm(w) == 1
        theta *= norm_w;
        w /= norm_w;
        v /= norm_w;

        Eigen::Matrix3d w_hat = skew( w );
        Eigen::Matrix3d exp_w_hat_theta = expmExact( w_hat, theta );
        Eigen::Matrix3d eye3 = Eigen::Matrix3d::Identity();
        
        expT.block< 3, 3 >( 0, 0 ) = exp_w_hat_theta;
        expT.block< 3, 1 >( 0, 3 ) = ( eye3 - exp_w_hat_theta ) * w.cross( v ) 
                                     + w * w.transpose() * v * theta;
    }
    
    return expT;
}

Eigen::Matrix4d kinematics::expTwist( const std::vector< Vector6d >& xi,
                                      const std::vector< double >& theta )
{
    Eigen::Matrix4d g = Eigen::Matrix4d::Identity();

    for ( unsigned int i = 0; i < theta.size(); i++ )
    {
      g = g * expTwist( xi[i], theta[i] );
    }

    return g;
}

////////////////////////////////////////////////////////////////////////////////
// Jacobians
////////////////////////////////////////////////////////////////////////////////

Matrix6Xd kinematics::spatialJacobian( const std::vector< Vector6d >& xi,
                                       const std::vector< double >& theta )
{
    int num_theta = theta.size();
    Matrix6Xd J_s( 6, num_theta );
    
    Eigen::Matrix4d g = Eigen::Matrix4d::Identity();
    Eigen::Matrix4d expT;
    
    for( int i = 0; i < num_theta; i++ )
    {
        if ( i == 0 )
        {
            J_s.block< 6, 1 >( 0, i ) = xi[i];
        }
        else
        {
            expT = expTwist( xi[i - 1], theta[i - 1] );
            g = g * expT;
            J_s.block< 6, 1 >( 0, i ) = adj( g ) * xi[i];
        }
    }

    return J_s;
}

Matrix6Xd kinematics::bodyJacobian( const std::vector< Vector6d >& xi,
                                    const std::vector< double >& theta,
                                    const Eigen::Matrix4d& g_zero )
{
    int num_theta = theta.size();
    Matrix6Xd J_b( 6, num_theta );
    
    Eigen::Matrix4d g = g_zero;

    Eigen::Matrix4d expT;
    
    for( int i = num_theta - 1; i >= 0; i-- )
    {
        expT = expTwist( xi[i], theta[i] );
        g = expT * g;
        J_b.block<6,1>( 0,i ) = adj( g.inverse() ) * xi[i];
    }

    return J_b;
}

Eigen::Matrix< double, 1, Eigen::Dynamic > kinematics::calculateJJointLimits(
                                        const std::vector< double >& theta,
                                        const double joint_limit )
{
    unsigned int num_joints = theta.size();
    
    Eigen::Matrix< double, 1, Eigen::Dynamic > J( 1, num_joints );

    // finite difference value
    double eps = 1e-6;

    // find the error before perturbing by eps
    double error = jointLimitError( theta, joint_limit );

    // perturb each joint by eps and see how error changes
    
    for ( unsigned int i = 0; i < num_joints; i++ )
    {
        std::vector< double > theta_tweak = theta;
        theta_tweak[i] += eps;
        double error_tweak = jointLimitError( theta_tweak, joint_limit );
        J(0, i) = ( error_tweak - error ) / eps;
    }

    return J;
}

Eigen::MatrixXd kinematics::dampedPinv6Xd( const kinematics::Matrix6Xd& J,
                                           const std::vector< double >& theta,
                                           const double theta_limit,
                                           const double limit_threshold,
                                           const double manipubility_threshold, 
                                           const double damping_ratio )
{
    unsigned int num_joints = theta.size();

    kinematics::Matrix6d JJtranspose = J * J.transpose();
    kinematics::Matrix6d W_x = kinematics::Matrix6d::Identity();
    Eigen::MatrixXd W_q = Eigen::MatrixXd::Identity( num_joints, num_joints );
    
    // Determine least-squares damping and weights, from:
    // "Robust Inverse Kinematics Using Damped Least Squares
    // with Dynamic Weighting"  NASA 1994 
    for( unsigned int i = 1; i < num_joints; i++ )
    {
        double theta_to_limit = theta_limit - std::abs( theta[i] );
        
        if (  theta_to_limit < limit_threshold )
        {
            W_q( i,i ) = 0.1 + 0.9 * theta_to_limit / limit_threshold;
        } 
    } 
    
    kinematics::Matrix6Xd J_w = W_x * J * W_q;
    
    // find the damping ratio
    // Based on Manipulability, in 'Prior Work' of above paper
    double manipubility = JJtranspose.determinant();
    double damping = 0;
    if ( manipubility < manipubility_threshold )
    {
        damping = damping_ratio * std::pow( 1 - manipubility / manipubility_threshold, 2 );
    }
    
    kinematics::Matrix6d tmp = JJtranspose + damping * kinematics::Matrix6d::Identity();
    Eigen::MatrixXd J_inv = J_w.transpose() * tmp.inverse();
    
    return W_q * J_inv * W_x;
}

Eigen::MatrixXd kinematics::dampedPinvXd( const Eigen::MatrixXd& J,
                                          const std::vector< double >& theta,
                                          const double theta_limit,
                                          const double limit_threshold,
                                          const double manipubility_threshold, 
                                          const double damping_ratio )
{
    unsigned int num_joints = theta.size();
    unsigned int num_velocities = J.rows();

    Eigen::MatrixXd JJtranspose = J * J.transpose();
    Eigen::MatrixXd W_x = Eigen::MatrixXd::Identity( num_velocities, num_velocities );
    Eigen::MatrixXd W_q = Eigen::MatrixXd::Identity( num_joints, num_joints );
    
    // Determine least-squares damping and weights, from:
    // "Robust Inverse Kinematics Using Damped Least Squares
    // with Dynamic Weighting"  NASA 1994 
    for( unsigned int i = 1; i < num_joints; i++ )
    {
        double theta_to_limit = theta_limit - std::abs( theta[i] );
        
        if (  theta_to_limit < limit_threshold )
        {
            W_q( i,i ) = 0.1 + 0.9 * theta_to_limit / limit_threshold;
        } 
    } 
    
    Eigen::MatrixXd J_w = W_x * J * W_q;
    
    // find the damping ratio
    // Based on Manipulability, in 'Prior Work' of above paper
    double manipubility = JJtranspose.determinant();
    double damping = 0;
    if ( manipubility < manipubility_threshold )
    {
        damping = damping_ratio * std::pow( 1 - manipubility / manipubility_threshold, 2 );
    }
    
    Eigen::MatrixXd tmp = JJtranspose + damping * Eigen::MatrixXd::Identity( num_velocities, num_velocities );
    Eigen::MatrixXd J_inv = J_w.transpose() * tmp.inverse();
    
    return W_q * J_inv * W_x;
}

////////////////////////////////////////////////////////////////////////////////
// Other
////////////////////////////////////////////////////////////////////////////////

Vector6d kinematics::calculateError( const Eigen::Matrix4d& g_current,
                                     const Eigen::Matrix4d& g_desired )
{
    Vector6d xi;

    Eigen::Matrix4d g_diff = g_current.inverse() * g_desired;

    xi = twistUnhat( g_diff.log() );

    return xi;
}

double kinematics::jointLimitError( const std::vector< double >& theta,
                                    const double joint_limit )
{
    //TODO: get rid of magic number
    double error_scale = 0.01;
    double error = 0;

    for ( unsigned int i = 0; i < theta.size(); i++)
    {
        error += joint_limit / std::pow( std::abs(theta[i]) - joint_limit, 4);
    }

    return error_scale * error;
}
