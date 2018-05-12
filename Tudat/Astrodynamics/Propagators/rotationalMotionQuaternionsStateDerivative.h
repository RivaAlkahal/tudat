/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_ROTATIONAL_MOTION_QUATERNIONS_STATE_DERIVATIVE_H
#define TUDAT_ROTATIONAL_MOTION_QUATERNIONS_STATE_DERIVATIVE_H

#include <vector>
#include <map>
#include <string>

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include "Tudat/Astrodynamics/Propagators/rotationalMotionStateDerivative.h"

namespace tudat
{

namespace propagators
{

//! Function to obtain the matrix by which a quaternion vector is to be pre-multiplied to obtain this quaternion's time-derivative
/*!
 * Function to obtain the matrix by which a quaternion vector (representing body-fixed to inertial frame rotation) is to be
 * pre-multiplied to obtain this quaternion's time-derivative
 * \param angularVelocityVectorInBodyFixedFrame  Current angular velocity vector of body, expressed in its body-fixed frame
 * \return Matrix by which a quaternion vector (representing body-fixed to inertial frame rotation) is to be
 * pre-multiplied to obtain this quaternion's time-derivative
 */
Eigen::Matrix4d getQuaterionToQuaternionRateMatrix( const Eigen::Vector3d& angularVelocityVectorInBodyFixedFrame );

//! Function to obtain the time derivative of a quaternion (in vector representation) of body-fixed to inertial frame
/*!
 * Function to obtain the time derivative of a quaternion (in vector representation) of body-fixed to inertial frame
 * \param currentQuaternionsToBaseFrame Quaternions (in vector representation) that defined the rotation from body-fixed to inertial
 * frame.
 * \param angularVelocityVectorInBodyFixedFrame Current angular velocity vector of body, expressed in its body-fixed frame
 * \return Time derivative of a quaternion (in vector representation) of body-fixed to inertial frame
 */
Eigen::Vector4d calculateQuaternionsDerivative( const Eigen::Vector4d& currentQuaternionsToBaseFrame,
                                                const Eigen::Vector3d& angularVelocityVectorInBodyFixedFrame );

//! Class for computing the state derivative for rotational dynamics of N bodies.
/*!
 *  Class for computing the state derivative for rotational dynamics of N bodies., using quaternion from body-fixed to inertial
 *  frame (in quaternion format) and angular velocity-vector of body expressed in body-fixed frame as the rotational state of a
 *  single body
 */
template< typename StateScalarType = double, typename TimeType = double >
class RotationalMotionQuaternionsStateDerivative: public RotationalMotionStateDerivative< StateScalarType, TimeType >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param torqueModelsPerBody List of torque models (first map key body undergoing acceleration, second map key body exerting
     * acceleration)
     * \param bodiesToPropagate List of names of bodies for which rotational state is to be propagated
     * \param bodyInertiaTensorFunctions List of functions returning inertia tensors of bodiesToPropagate (in same order)
     * \param bodyInertiaTensorTimeDerivativeFunctions List of functions returning time derivatives of inertia tensors of
     *  bodiesToPropagate (in same order). Default empty, denoting time-invariant inertia tensors.
     */
    RotationalMotionQuaternionsStateDerivative(
            const basic_astrodynamics::TorqueModelMap& torqueModelsPerBody,
            const std::vector< std::string >& bodiesToPropagate,
            std::vector< boost::function< Eigen::Matrix3d( ) > > bodyInertiaTensorFunctions,
            std::vector< boost::function< Eigen::Matrix3d( ) > > bodyInertiaTensorTimeDerivativeFunctions =
            std::vector< boost::function< Eigen::Matrix3d( ) > >( ) ):
        RotationalMotionStateDerivative< StateScalarType, TimeType >(
            torqueModelsPerBody, quaternions, bodiesToPropagate, bodyInertiaTensorFunctions,
            bodyInertiaTensorTimeDerivativeFunctions )
    {

    }

    //! Destructor
    ~RotationalMotionQuaternionsStateDerivative( ){ }

    //! Calculates the state derivative of the rotational motion of the system.
    /*!
     *  Calculates the state derivative of the rotational motion of the system at the given time and rotational state.
     *  \param time Time (seconds since reference epoch) at which the system is to be updated.
     *  \param stateOfSystemToBeIntegrated List of 7 * bodiesToPropagate_.size( ), containing rotation quaternion/angular
     *  velocity of the bodies being propagated. The order of the values is defined by the order of bodies in
     *  bodiesToPropagate_
     *  \param stateDerivative Current state derivative (quaternion rate+angular acceleration) of system of bodies
     *  integrated numerically (returned by reference).
     */
    void calculateSystemStateDerivative(
            const TimeType time,
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& stateOfSystemToBeIntegrated,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > stateDerivative )
    {
        stateDerivative = Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( stateOfSystemToBeIntegrated.rows( ), 1 );
        std::vector< Eigen::Vector3d > torquesActingOnBodies = this->sumTorquesPerBody( );

        for( unsigned int i = 0; i < torquesActingOnBodies.size( ); i++ )
        {
            Eigen::Matrix< StateScalarType, 4, 1 > currentQuaternions = stateOfSystemToBeIntegrated.block( i * 7, 0, 4, 1 );
            Eigen::Matrix< StateScalarType, 3, 1 > currentBodyFixedRotationRate = stateOfSystemToBeIntegrated.block( i * 7 + 4, 0, 3, 1 );

            stateDerivative.block( i * 7, 0, 4, 1 ) = calculateQuaternionsDerivative(
                        currentQuaternions.template cast< double >( ), currentBodyFixedRotationRate.template cast< double >( ) ).
                    template cast< StateScalarType >( );
            stateDerivative.block( i * 7 + 4, 0, 3, 1 ) = evaluateRotationalEquationsOfMotion(
                        this->bodyInertiaTensorFunctions_.at( i )( ), torquesActingOnBodies.at( i ),
                        currentBodyFixedRotationRate.template cast< double >( ),
                        this->bodyInertiaTensorTimeDerivativeFunctions_.at( i )( ) ).template cast< StateScalarType >( );

//            std::cout << "Time: " << time - 236455200 << std::endl;
//            std::cout << "Quat: " << currentQuaternions.transpose( ) << std::endl;
//            std::cout << "Rot: " << currentBodyFixedRotationRate.transpose( ) << std::endl;
//            std::cout << "Torque: " << torquesActingOnBodies.at( i ).transpose( ) << std::endl;
//            std::cout << "Deriv: " << stateDerivative.block( i * 7, 0, 7, 1 ).transpose( ) << std::endl << std::endl;
        }
    }

    //! Function to convert the state in the conventional form to the propagator-specific form.
    /*!
     * Function to convert the state in the conventional form to the propagator-specific form. For this propagator,
     * the two are equivalent, and this function returns the input state.
     * \param outputSolution State in 'conventional form'
     * \param time Current time at which the state is valid (not used in this class).
     * \return State (outputSolution), converted to the 'propagator-specific form' (which is equal to outputSolution).
     */
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > convertFromOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& outputSolution, const TimeType& time )
    {
        return outputSolution;
    }

    //! Function to convert the propagator-specific form of the state to the conventional form.
    /*!
     * Function to convert the propagator-specific form of the state to the conventional form. For the this propagator,
     * the two are equivalent, and this function returns the input state.
     * \param internalSolution State in propagator-specific form (which is equal to outputSolution to conventional form for
     * this propagator)
     * \param time Current time at which the state is valid (not used in this class).
     * \param currentLocalSolution State (internalSolution), converted to the 'conventional form',
     * which is equal to outputSolution for this class (returned by reference).
     */
    void convertToOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& internalSolution, const TimeType& time,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > currentLocalSolution )
    {
        currentLocalSolution = internalSolution;
        currentQuaternionLocalSolution_ = currentLocalSolution;
    }

    //! Function to process the state during propagation.
    /*!
     * Function to process the state during propagation. For quaternions, this function normalizes the quaternion vector
     * in case its magnitude differs from 1.0 by a value larger than the tolerance.
     * \param unprocessedState State computed after propagation.
     * \param startRow Dummy variable added for compatibility issues between Eigen::Matrix and Eigen::Block.
     */
    void postProcessState( Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& unprocessedState,
                           const int startRow )
    {
        // Loop over each body
        const double tolerance = 20.0 * std::numeric_limits< double >::epsilon( );
        Eigen::Matrix< StateScalarType, 4, 1 > quaternionsVector;
        StateScalarType quaternionsMagnitude;
        for( unsigned int i = 0; i < this->bodiesToPropagate_.size( ); i++ )
        {
            // Normalize quaternions
            quaternionsVector = unprocessedState.block( i * 7, 0, 4, 1 );
            quaternionsMagnitude = quaternionsVector.norm( );
            if ( std::fabs( quaternionsMagnitude - 1.0 ) >= tolerance )
            {
                // Normalize
                quaternionsVector /= quaternionsMagnitude;

                // Replace old quaternions with normalized quaternions
                unprocessedState.segment( startRow + i * 7, 4 ) = quaternionsVector;
            }
        }
    }

    //! Function to return whether the state needs to be post-processed.
    /*!
     * Function to return whether the state needs to be post-processed. For quaternions this is true.
     * \return Boolean confirming that the state needs to be post-processed.
     */
    bool isStateToBePostProcessed( )
    {
        return true;
    }

private:

    //! Current full state of the propagated bodies, w.r.t. the central bodies, where the attitude is expressed in quaternions.
    /*!
     *  Current full state of the propagated bodies, w.r.t. the central bodies, where the attitude is expressed in quaternions.
     *  These variables are set when calling the convertToOutputSolution function.
     */
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > currentQuaternionLocalSolution_;

};

}

}

#endif // TUDAT_ROTATIONAL_MOTION_QUATERNIONS_STATE_DERIVATIVE_H
