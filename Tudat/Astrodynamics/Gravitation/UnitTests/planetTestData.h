/*    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120720    K. Kumar          File created.
 *      120725    K. Kumar          Completed Doxygen documentation; added get-functions() for all
 *                                  planet test data.
 *      120801    K. Kumar          Renamed functions to retrieve MATLAB test data.
 *      120815    K. Kumar          Added test data functions for data from (Melman, 2012) and
 *                                  (Ronse, 2012).
 *
 *    References
 *      Dekens, E. Orbit Analysis of a Low Frequency Array for Radio Astronomy, MSc thesis, Delft
 *          University of Technlogy, Delft, The Netherlands, 2012.
 *      MathWorks. gravityzonal, MATLAB 2012b, 2012.
 *      Melman, J. Propagate software, J.C.P.Melman@tudelft.nl, 2012.
 *      Ronse, A. A parametric study of space debris impact footprints, MSc thesis, Delft
 *          University of Technlogy, Delft, The Netherlands, in preparation.
 *
 *    Notes
 *      This file is used by the unitTestGravitationalAcceleration.cpp file. The test data in
 *      included in this file was obtained using the gravityzonal() function in MATLAB
 *      (Mathworks, 2012), and code written as part of MSc thesis work by (Dekens, 2012).
 *
 */

#ifndef TUDAT_PLANET_TEST_DATA_H
#define TUDAT_PLANET_TEST_DATA_H

#include <map>
#include <string>
#include <vector>

#include <Eigen/Core>

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

namespace tudat
{
namespace unit_tests
{

//! Zonal coefficients.
/*!
 * Zonal gravity field coefficients up to J4.
 */
enum ZonalCoefficients { central = 0, j2 = 2, j3 = 3, j4 = 4 };

//! Planets
/*!
 * Planets for which test data is available (Moon is taken as a planet).
 */
enum Planets { mercury, venus, earth, moon, mars, jupiter, saturn, uranus, neptune };

//! Structure containing test data for a given planet.
/*!
 * This structure contains test data for a given planet. The data is used for unit testing of
 * Tudat's gravitational acceleration computation free functions.
 */
struct PlanetTestData
{
public:

    //! Definition of useful typedefs.
    typedef std::vector< Eigen::Vector3d > VectorOfVector3ds;
    typedef std::map< int, double > KeyIntValueDoubleMap;
    typedef std::map< int, std::map< int, std::map< int, Eigen::Vector3d > > > IndexedVector3d;

    //! Default constructor.
    PlanetTestData( const std::string& aPlanetName,
                    const double aGravitationalParameter,
                    const KeyIntValueDoubleMap& someZonalCoefficients,
                    const double anEffectiveRadius,
                    const std::vector< Eigen::Vector3d >& positionsOfBodyExertingAcceleration,
                    const std::vector< Eigen::Vector3d >& positionsOfBodySubjectToAcceleration )
        : planetName( aPlanetName ),
          gravitationalParameter( aGravitationalParameter ),
          zonalCoefficients( someZonalCoefficients ),
          effectiveRadius( anEffectiveRadius ),
          body1Positions( positionsOfBodyExertingAcceleration  ),
          body2Positions( positionsOfBodySubjectToAcceleration  ),
          expectedAcceleration( IndexedVector3d( ) )
    { }

    //! Empty constructor.
    PlanetTestData( )
        : planetName( "" ),
          gravitationalParameter( TUDAT_NAN ),
          zonalCoefficients( KeyIntValueDoubleMap( ) ),
          effectiveRadius( TUDAT_NAN ),
          body1Positions( VectorOfVector3ds( ) ),
          body2Positions( VectorOfVector3ds( ) ),
          expectedAcceleration( IndexedVector3d( ) )
    { }

    //! Planet name.
    /*!
     * The planet's name, given as string, used for verification purposes. Note, that the Moon is
     * considered a "planet" in this context.
     */
    std::string planetName;

    //! Gravitational parameter [m s^3].
    /*!
     * The gravitational parameter of the central body, taken from the documentation for MATLAB's
     * gravityzonal() function (Mathworks, 2012).
     */
    double gravitationalParameter;

    //! Zonal gravity field coefficients.
    /*!
     * Map that indexes the zonal gravity coefficients available for a given central body. These
     * are taken from the documentation for MATLAB's gravityzonal() function (Mathworks, 2012).
     */
    KeyIntValueDoubleMap zonalCoefficients;

    //! Effective (equatorial) radius [m].
    /*!
     * The effective radius, or equatorial radius, is part of the spherical harmonics expansion and
     * is given in the documentation of MATLAB's gravityzonal() function (Mathworks, 2012).
     */
    double effectiveRadius;

    //! Positions of body exerting acceleration.
    /*!
     * Positions of body exerting acceleration (the cases run are for the planets and the Moon).
     * For each body, two positions are taken: one located at the origin and one offset from the
     * origin.
     */
    VectorOfVector3ds body1Positions;

    //! Positions of body subject to acceleration.
    /*!
     * Position of the body subject to acceleration. Currently, two positions are included per
     * planet.
     */
    VectorOfVector3ds body2Positions;

    //! Expected gravitational accelerations [m s^-2].
    /*!
     * Expected gravitational accelerations generated using MATLAB's gravityzonal() function
     * (Mathworks, 2012). They are indexed by the position of the body exerting the acceleration,
     * the position of the body subject to the acceleration, and finally the maximum zonal gravity
     * term included.
     */
    IndexedVector3d expectedAcceleration;

protected:

private:
};

//! Add (raw) expected accelerations to planet data.
void addExpectedAccelerations( PlanetTestData& planetData, const Eigen::MatrixXd& rawData );

//! Get Mercury test data generated using MATLAB (Mathworks, 2012).
/*!
 * Returns Mercury test data, generated with MATLAB's gravityzonal() function (Mathworks, 2012).
 * This is used as benchmark testing of Tudat code.
 * \return Structure containing Mercury test data.
 */
PlanetTestData getMercuryMatlabTestData( );

//! Get Venus test data generated using MATLAB (Mathworks, 2012).
/*!
 * Returns Venus test data, generated with MATLAB's gravityzonal() function (Mathworks, 2012).
 * This is used as benchmark testing of Tudat code.
 * \return Structure containing Venus test data.
 */
PlanetTestData getVenusMatlabTestData( );

//! Get Earth test data generated using MATLAB (Mathworks, 2012).
/*!
 * Returns Earth test data, generated with MATLAB's gravityzonal() function (Mathworks, 2012).
 * This is used as benchmark testing of Tudat code.
 * \return Structure containing Earth test data.
 */
PlanetTestData getEarthMatlabTestData( );

//! Get Moon test data generated using MATLAB (Mathworks, 2012).
/*!
 * Returns Moon test data, generated with MATLAB's gravityzonal() function (Mathworks, 2012).
 * This is used as benchmark testing of Tudat code.
 * \return Structure containing Moon test data.
 */
PlanetTestData getMoonMatlabTestData( );

//! Get Mars test data generated using MATLAB (Mathworks, 2012).
/*!
 * Returns Mars test data, generated with MATLAB's gravityzonal() function (Mathworks, 2012).
 * This is used as benchmark testing of Tudat code.
 * \return Structure containing Mars test data.
 */
PlanetTestData getMarsMatlabTestData( );

//! Get Jupiter test data generated using MATLAB (Mathworks, 2012).
/*!
 * Returns Jupiter test data, generated with MATLAB's gravityzonal() function (Mathworks, 2012).
 * This is used as benchmark testing of Tudat code.
 * \return Structure containing Jupiter test data.
 */
PlanetTestData getJupiterMatlabTestData( );

//! Get Saturn test data generated using MATLAB (Mathworks, 2012).
/*!
 * Returns Saturn test data, generated with MATLAB's gravityzonal() function (Mathworks, 2012).
 * This is used as benchmark testing of Tudat code.
 * \return Structure containing Saturn test data.
 */
PlanetTestData getSaturnMatlabTestData( );

//! Get Uranus test data generated using MATLAB (Mathworks, 2012).
/*!
 * Returns Uranus test data, generated with MATLAB's gravityzonal() function (Mathworks, 2012).
 * This is used as benchmark testing of Tudat code.
 * \return Structure containing Uranus test data.
 */
PlanetTestData getUranusMatlabTestData( );

//! Get Neptune test data generated using MATLAB (Mathworks, 2012).
/*!
 * Returns Neptune test data, generated with MATLAB's gravityzonal() function (Mathworks, 2012).
 * This is used as benchmark testing of Tudat code.
 * \return Structure containing Neptune test data.
 */
PlanetTestData getNeptuneMatlabTestData( );

//! Get planet test data generated using MATLAB (Mathworks, 2012).
/*!
 * Returns vector containing all test data available, generated with MATLAB's gravityzonal()
 * function (Mathworks, 2012). Currently, test data is available for: Mercury, Venus, Earth Moon,
 * Mars, Jupiter, Saturn, Uranus, Neptune.
 * \return Vector of structures containing planet test data.
 */
std::vector< PlanetTestData > getPlanetMatlabTestData( );

//! Get Earth test data generated using code from (Melman, 2012).
/*!
 * Returns Earth test data, generated using code written by (Melman, 2012). This is used as
 * benchmark testing of Tudat code.
 * \return Structure containing Earth test data.
 */
PlanetTestData getEarthMelmanTestData( );

//! Get Earth test data generated using code from (Ronse, 2012).
/*!
 * Returns Earth test data, generated using code written by (Ronse, 2012). This is used as
 * benchmark testing of Tudat code.
 * \return Structure containing Earth test data.
 */
PlanetTestData getEarthRonseTestData( );


} // namespace unit_tests
} // namespace tudat

#endif // TUDAT_PLANET_TEST_DATA_H
