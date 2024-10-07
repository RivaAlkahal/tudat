//
// Created by Riva Alkahal on 04/10/2024.
//
//
// Created by Riva Alkahal on 03/09/2024.
//
//
// Created by Riva Alkahal on 28/08/2024.
//
//
// Created by Riva Alkahal on 27/08/2024.
//
//
// Created by Riva Alkahal on 26/08/2024.
//
//
// Created by Riva Alkahal on 23/08/2024.
//
//
// Created by Riva Alkahal on 19/08/2024.
//
//
// Created by Riva Alkahal on 11/08/2024.
//
//
// Created by Riva Alkahal on 26/07/2024.
//

#include <iostream>
#include <fstream>
#include <limits>
#include "fstream"
#include "iostream"

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/aerodynamics/customAerodynamicCoefficientInterface.h"
#include "tudat/astro/aerodynamics/aerodynamicAcceleration.h"
#include "tudat/astro/reference_frames/aerodynamicAngleCalculator.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/interface/spice/spiceEphemeris.h"
#include "tudat/interface/spice/spiceRotationalEphemeris.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/aerodynamics/marsDtmAtmosphereModel.h"
#include "tudat/io/solarActivityData.h"
#include "tudat/astro/basic_astro/unitConversions.h"

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/estimation.h"
#include "tudat/simulation/estimation_setup.h"

#include "tudat/io/readOdfFile.h"
#include "tudat/io/readTabulatedMediaCorrections.h"
#include "tudat/io/readTabulatedWeatherData.h"
#include "tudat/simulation/estimation_setup/processOdfFile.h"
#include "tudat/simulation/propagation_setup/setNumericallyIntegratedStates.h"

#include <boost/date_time/gregorian/gregorian.hpp>

#include "tudat/astro/ground_stations/transmittingFrequencies.h"


int main( ) {
    using namespace tudat;
    using namespace aerodynamics;
    using namespace simulation_setup;
    using namespace numerical_integrators;
    using namespace simulation_setup;
    using namespace basic_astrodynamics;
    using namespace propagators;
    using namespace estimatable_parameters;
    using namespace observation_models;
    using namespace basic_mathematics;
    using namespace basic_astrodynamics;
    using namespace tudat::spice_interface;
    using namespace tudat::ephemerides;
    using namespace tudat::input_output;
    using namespace tudat::orbit_determination;
    using namespace tudat::interpolators;
    using namespace tudat::orbital_element_conversions;

    spice_interface::loadStandardSpiceKernels( );
    /*spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_map4.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_map5.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_map6.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_map7.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_map8.bsp" );

        spice_interface::loadSpiceKernelInTudat( "/Users/ralkahal/OneDrive - Delft University of Technology/esitmate/sod_assignments/mgs_map4.bsp" );
        spice_interface::loadSpiceKernelInTudat( "/Users/ralkahal/OneDrive - Delft University of Technology/esitmate/sod_assignments/mgs_map5.bsp" );
        spice_interface::loadSpiceKernelInTudat( "/Users/ralkahal/OneDrive - Delft University of Technology/esitmate/sod_assignments/mgs_map6.bsp" );
        spice_interface::loadSpiceKernelInTudat( "/Users/ralkahal/OneDrive - Delft University of Technology/esitmate/sod_assignments/mgs_map7.bsp" );
        spice_interface::loadSpiceKernelInTudat( "/Users/ralkahal/OneDrive - Delft University of Technology/esitmate/sod_assignments/mgs_map8.bsp" );
        */

    std::string saveDirectory = "/Users/ralkahal/OneDrive - Delft University of Technology/new-tudat-tests/";
    std::string fileTag = "nospice";

    // set input options
    double epehemeridesTimeStep = 60.0;
    bool useInterpolatedEphemerides = false;
    double observationsSamplingTime = 60.0;
    double buffer = 20.0 * epehemeridesTimeStep;
    double arcDuration = 3*86400.0;//2.0E4;
    std::string dragEst = "per-halfday";//"per-rev";
    double ndays = 0.5;
    double hoursperdaydrag = 2.0;
    double hoursperday = 10.0;
    int iterationNumber = 6;
    //const double gravitationalParameter = 4.2828378e13;
    //const double planetaryRadius = 3389.5E3;

    double oneWayDopplerNoise = 0.0001;
    double twoWayDopplerNoise = 0.0001;
    double rangeNoise = 0.0001;
    // Select ephemeris time range (based on available data in loaded SPICE ephemeris)
    //Time initialEphemerisTime = Time( 185976000 - 1.0 * 86400.0 ); // 23 November 2005, 0h
    //Time finalEphemerisTime = Time( 186580800 + 1.0 * 86400.0 ); // 30 November 2005, 0h
    // time at 2000-01-01T00:00:00
    Time initialEphemerisTime = Time( 0.0 ) ;
    Time finalEphemerisTime = Time( 86400.0 * 360); // 5 years later
    double totalDuration = finalEphemerisTime - initialEphemerisTime;
    std::cout<<"Total duration: "<<totalDuration<<std::endl;

    std::vector< std::string > bodiesToCreate = {
            "Earth", "Sun", "Mercury", "Venus", "Mars", "Jupiter", "Phobos", "Deimos" };

    std::string baseFrameOrientation = "MARSIAU";
    std::string baseFrameOrigin = "SSB";

    // Specify ephemeris time steps and buffers
    Time ephemerisTimeStepPlanets =  epehemeridesTimeStep;
    std::cout<<"Ephemeris time step planets: "<<epehemeridesTimeStep<<std::endl;
    //Time bufferPlanets = Time( 10.0 * ephemerisTimeStepPlanets );
    Time ephemerisTimeStepSpacecraft = Time( epehemeridesTimeStep );
    //Time bufferSpacecraft = Time( 10.0 * ephemerisTimeStepSpacecraft );

    BodyListSettings bodySettings;
    if ( useInterpolatedEphemerides )
    {
        bodySettings = getDefaultBodySettings(
                bodiesToCreate, initialEphemerisTime - buffer, finalEphemerisTime + buffer,
                baseFrameOrigin, baseFrameOrientation, ephemerisTimeStepPlanets );
    }
    else
    {
        bodySettings = getDefaultBodySettings( bodiesToCreate, baseFrameOrigin, baseFrameOrientation );
    }
    std::cout<<"Body Settings created"<<std::endl;

    bodySettings.at( "Earth" )->groundStationSettings = getDsnStationSettings( );

    std::string filename ="/Users/ralkahal/OneDrive - Delft University of Technology/PhD/Programs/atmodensitydtm/dtm_mars";;
    //std::string filename ="/home/ralkahal/new-tudat-tests/dtm-mars";

    bodySettings.at( "Mars" )->atmosphereSettings = marsDtmAtmosphereSettings( filename, 3378.0E3);
    std::cout<<"atmospheric settings created"<<std::endl;
    // Set spherical harmonics gravity field
    // Create spacecraft
    std::string spacecraftName = "MGS";
    bodySettings.addSettings( spacecraftName );
    if ( useInterpolatedEphemerides )
    {
        bodySettings.at( spacecraftName )->ephemerisSettings =
                std::make_shared< InterpolatedSpiceEphemerisSettings >(
                        initialEphemerisTime - buffer, finalEphemerisTime + buffer,
                        ephemerisTimeStepSpacecraft, baseFrameOrigin, baseFrameOrientation,
                        std::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 ), spacecraftName );
    }
    else
    {
        //bodySettings.at( spacecraftName )->ephemerisSettings =
        //       std::make_shared< DirectSpiceEphemerisSettings >( baseFrameOrigin, baseFrameOrientation );
            bodySettings.at( spacecraftName )->ephemerisSettings = std::make_shared<TabulatedEphemerisSettings>
            (std::map< double, Eigen::Vector6d >(), baseFrameOrigin, baseFrameOrientation);
    }

    std::cout<<"bodySettings of spacecraft created"<<std::endl;
    bodySettings.at( spacecraftName )->constantMass = 700.0;
    bodySettings.at( spacecraftName )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
    std::cout<<"creating gravity field variations"<<std::endl;
    // Set gravity field variations
    std::vector< std::shared_ptr< GravityFieldVariationSettings > > gravityFieldVariations;
    std::cout<<"gravity field variation settings called"<<std::endl;
    // Set solid body tide gravity field variation
    std::vector< std::string > deformingBodies;
    deformingBodies.push_back( "Sun" );

    std::map< int, std::vector< std::complex< double > > > loveNumbers;
    std::vector< std::complex< double > > degreeTwoLoveNumbers_;
    degreeTwoLoveNumbers_.push_back( std::complex< double >( 0.169, 0.0 ) );
    loveNumbers[ 2 ] = degreeTwoLoveNumbers_;

    std::shared_ptr< GravityFieldVariationSettings > singleGravityFieldVariation =
            std::make_shared< BasicSolidBodyGravityFieldVariationSettings >( deformingBodies, loveNumbers );
    gravityFieldVariations.push_back( singleGravityFieldVariation );

    // Set periodic gravity field variation
    std::vector< Eigen::MatrixXd > cosineShAmplitudesCosineTime;
    std::vector< Eigen::MatrixXd > cosineShAmplitudesSineTime;
    std::vector< Eigen::MatrixXd > sineShAmplitudesCosineTime;
    std::vector< Eigen::MatrixXd > sineShAmplitudesSineTime;
    std::vector< double > frequencies;
    std::cout<<"initialization of periodic grav vriaions done"<<std::endl;
    // Values from Mars GMM3_120_SHA, https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/data/shadr/gmm3_120_sha.lbl
    cosineShAmplitudesCosineTime.push_back(
            ( Eigen::MatrixXd( 4, 3 )<<2.39E-9, 0.92E-10, 0.0,
                    1.67E-9, -2.22E-10, 0.0,
                    0.85E-10, 0.0, 0.0,
                    0.38E-9, 0.0, 0.0 ).finished( ) );
    cosineShAmplitudesCosineTime.push_back(
            ( Eigen::MatrixXd( 4, 3 )<<1.23E-9, -0.19E-10, 0.0,
                    0.32E-9, 0.0, 0.0,
                    0.35E-10, 0.0, 0.0,
                    0.15E-9, 0.0, 0.0 ).finished( ) );
    cosineShAmplitudesCosineTime.push_back(
            ( Eigen::MatrixXd( 4, 3 )<<0.53E-9, 0.0, 0.0,
                    0.13E-9, 0.0, 0.0,
                    -0.51E-10, 0.0, 0.0,
                    0.32E-9, 0.0, 0.0).finished( ) );

    cosineShAmplitudesSineTime.push_back(
            ( Eigen::MatrixXd( 4, 3 )<<-0.83E-9, -1.68E-9, 0.0,
                    2.35E-9, 0.48E-10, 0.0,
                    -1.56E-10, 0.0, 0.0,
                    1.30E-9, 0.0, 0.0 ).finished( ) );
    cosineShAmplitudesSineTime.push_back(
            ( Eigen::MatrixXd( 4, 3 )<<0.73E-9, -0.16E-10, 0.0,
                    0.21E-9, 0.0, 0.0,
                    -0.24E-10, 0.0, 0.0,
                    0.42E-9, 0.0, 0.0).finished( ) );
    cosineShAmplitudesSineTime.push_back(
            ( Eigen::MatrixXd( 4, 3 )<<0.46E-9, 0.0, 0.0,
                    0.15E-9, 0.0, 0.0,
                    -0.64E-10, 0.0, 0.0,
                    -0.02E-09, 0.0, 0.0 ).finished( ) );

    sineShAmplitudesCosineTime.push_back(Eigen::MatrixXd::Zero( 4, 3 ));
    sineShAmplitudesCosineTime.push_back(Eigen::MatrixXd::Zero( 4, 3 ));
    sineShAmplitudesCosineTime.push_back(Eigen::MatrixXd::Zero( 4, 3 ));

    sineShAmplitudesSineTime.push_back(Eigen::MatrixXd::Zero( 4, 3 ));
    sineShAmplitudesSineTime.push_back(Eigen::MatrixXd::Zero( 4, 3 ));
    sineShAmplitudesSineTime.push_back(Eigen::MatrixXd::Zero( 4, 3 ));
    frequencies.resize( 3 );
    frequencies = { 2*mathematical_constants::PI/(686.98*86400.0), 4*mathematical_constants::PI/(686.98*86400.0), 6*mathematical_constants::PI/(686.98*86400.0) };
    std::cout<<"assigned values for the amplitudes"<<std::endl;
    std::shared_ptr< GravityFieldVariationSettings > periodicGravityFieldVariations =
            std::make_shared< PeriodicGravityFieldVariationsSettings >(
                    cosineShAmplitudesCosineTime, cosineShAmplitudesSineTime, sineShAmplitudesCosineTime, sineShAmplitudesSineTime,
                    frequencies, 0.0, 2, 0 );

    gravityFieldVariations.push_back( periodicGravityFieldVariations );
    std::cout<<"periodic gravity field variation created"<<std::endl;

    // Set polynomial gravity field variation
    std::map<int, Eigen::MatrixXd> cosineAmplitudes;
    cosineAmplitudes[ 1 ] = Eigen::Matrix< double, 4, 3 >::Zero( );
    /*
    cosineAmplitudes[ 1 ]( 0, 0 ) = 1.91924728977996e-11/(365*24*3600); //3years sam-b

    cosineAmplitudes[1](0,1) =-1.69215684027738e-13/(365*24*3600);
    cosineAmplitudes[1](0,2) = 1.66459256502814e-11/(365*24*3600);
    cosineAmplitudes[1](1,0) = 9.33898280285743e-13/(365*24*3600);
    cosineAmplitudes[1](1,1) = 3.67913654008894e-12/(365*24*3600);
    cosineAmplitudes[1](1,2) = 1.12986679730805e-12/(365*24*3600);
    cosineAmplitudes[1](1,3) = 1.15599079203426e-11/(365*24*3600);
    cosineAmplitudes[1](2,0) = -2.9036065707624e-12/(365*24*3600);
*/
    //scM/10
/*
    cosineAmplitudes[ 1 ]( 0, 0 ) = -1.73871738023640e-10/(365*24*3600);

    cosineAmplitudes[1](0,1) =1.83526161335626e-10/(365*24*3600);
    cosineAmplitudes[1](0,2) = -2.29416809459079e-10/(365*24*3600);
    cosineAmplitudes[1](1,0) = -1.96888593742833e-10/(365*24*3600);
    cosineAmplitudes[1](1,1) = 3.64572599774489e-10/(365*24*3600);
    cosineAmplitudes[1](1,2) = -1.99682320271527e-10/(365*24*3600);
    cosineAmplitudes[1](1,3) = 1.19185562092185e-09/(365*24*3600);
    cosineAmplitudes[1](2,0) = 4.33626663660650e-10/(365*24*3600);

*/

    //scM/0.01

    cosineAmplitudes[ 1 ]( 0, 0 ) += -1.73871738023640e-12/(365*24*3600);

    cosineAmplitudes[1](0,1) +=1.83526161335626e-12/(365*24*3600);
    cosineAmplitudes[1](0,2) += -2.29416809459079e-12/(365*24*3600);
    cosineAmplitudes[1](1,0) += -1.96888593742833e-12/(365*24*3600);
    cosineAmplitudes[1](1,1) += 3.64572599774489e-12/(365*24*3600);
    cosineAmplitudes[1](1,2) += -1.99682320271527e-12/(365*24*3600);
    //cosineAmplitudes[1](1,3) += 1.19185562092185e-11/(365*24*3600);
    cosineAmplitudes[1](2,0) += 4.33626663660650e-12/(365*24*3600);
    std::cout<<"assigned values to cosine amplitde"<<std::endl;
    //samuelb 50years
/*
    cosineAmplitudes[ 1 ]( 0, 0 ) = 2.64611514298544e-10/(365*24*3600);

    cosineAmplitudes[1](0,1) =-7.9324131015657e-12/(365*24*3600);
    cosineAmplitudes[1](0,2) = 3.09051829331777e-10/(365*24*3600);
    cosineAmplitudes[1](1,0) = 1.23467940475808e-11/(365*24*3600);
    cosineAmplitudes[1](1,1) = 5.45964542190114e-11/(365*24*3600);
    cosineAmplitudes[1](1,2) = 1.0859659257389e-11/(365*24*3600);
    cosineAmplitudes[1](1,3) = 1.6671318233878e-10/(365*24*3600);
    cosineAmplitudes[1](2,0) = -2.55988773213135e-11/(365*24*3600);
*/


    //cosineAmplitudes = cosineAmplitudes;
    //stein 50 years
    /*
    cosineAmplitudes[1](0,0) = 6.469163469325e-09/(365*24*3600);
    cosineAmplitudes[1](0,1) = 1.63229205379438e-10/(365*24*3600);
    cosineAmplitudes[1](0,2) = 7.66708257544063e-09/(365*24*3600);
    cosineAmplitudes[1](1,0) = -1.80240472010679e-10/(365*24*3600);
    cosineAmplitudes[1](1,1) = 2.48842335218713e-09/(365*24*3600);
    cosineAmplitudes[1](1,2) =-1.59129711509854e-10/(365*24*3600);
    cosineAmplitudes[1](1,3) = 7.61407887767204e-09/(365*24*3600);
    cosineAmplitudes[1](2,0) = -1.35876027631454e-09/(365*24*3600);
    //cosineAmplitudes[ 2 ] = Eigen::Matrix< double, 3, 5 >::Zero( );
     */
    std::map<int, Eigen::MatrixXd> sineAmplitudes;
    //sineAmplitudes[ 1 ] = Eigen::Matrix< double, 4, 3 >::Zero( );
    /*
    sineAmplitudes[ 1 ]( 0, 1 ) = 2.1923517483097e-12; //3years sam-b
    sineAmplitudes[1](0,2) = 1.93032794374673e-11;
    sineAmplitudes[1](1,1) = -9.09002230132412e-12;
    sineAmplitudes[1](1,2) = 8.10801664003273e-13;
    sineAmplitudes[1](1,3) = -4.02014384239032e-12;
     */
    std::cout<<"creating settings for poly grav"<<std::endl;
    std::shared_ptr< GravityFieldVariationSettings > polynomialGravityFieldVariations =
            std::make_shared< PolynomialGravityFieldVariationsSettings >(
                    cosineAmplitudes, sineAmplitudes, 0.0, 2, 0 );
    gravityFieldVariations.push_back( polynomialGravityFieldVariations );

    std::vector< std::shared_ptr< GravityFieldVariationSettings > > gravityFieldVariationSettings =
            gravityFieldVariations;
    bodySettings.at( "Mars" )->gravityFieldVariationSettings = gravityFieldVariations;

    SystemOfBodies bodies = createSystemOfBodies< long double, Time >( bodySettings );

    std::cout<<"Bodies created"<<std::endl;
    // Create radiation pressure settings
    double referenceAreaRadiation = 4.0;
    double radiationPressureCoefficient = 1.2;
    std::vector< std::string > occultingBodies = { "Mars" };
    std::shared_ptr< RadiationPressureInterfaceSettings > radiationPressureSettings =
            std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                    "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodies.at( spacecraftName )->setRadiationPressureInterface(
            "Sun", createRadiationPressureInterface(
                    radiationPressureSettings, spacecraftName, bodies ) );

    // Create aerodynamic coefficients settings
    Eigen::Vector3d customVector(1.2, 0.0, 0.0);
    std::shared_ptr<AerodynamicCoefficientSettings> aerodynamicCoefficientSettings =
            std::make_shared<ConstantAerodynamicCoefficientSettings>(4.0,  1.2 * Eigen::Vector3d::UnitX( ));
    bodies.at( spacecraftName )->setAerodynamicCoefficientInterface(
            createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, spacecraftName, bodies ) );


    // Set accelerations on Vehicle that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Sun" ].push_back( pointMassGravityAcceleration( ) );
    accelerationsOfVehicle[ "Sun" ].push_back( cannonBallRadiationPressureAcceleration( ) );
//    accelerationsOfVehicle[ "Sun" ].push_back( relativisticAccelerationCorrection(  ) );
    accelerationsOfVehicle[ "Mercury" ].push_back( pointMassGravityAcceleration( ) );
    accelerationsOfVehicle[ "Venus" ].push_back( pointMassGravityAcceleration( ) );
    accelerationsOfVehicle[ "Earth" ].push_back( pointMassGravityAcceleration( ) );
    accelerationsOfVehicle[ "Mars" ].push_back( sphericalHarmonicAcceleration( 50, 50 ) );
    accelerationsOfVehicle[ "Mars" ].push_back( relativisticAccelerationCorrection(  ) );
    accelerationsOfVehicle[ "Mars" ].push_back( aerodynamicAcceleration( ) );
    accelerationsOfVehicle[ "Phobos" ].push_back( pointMassGravityAcceleration( ) );
    accelerationsOfVehicle[ "Deimos" ].push_back( pointMassGravityAcceleration( ) );
    accelerationsOfVehicle[ "Jupiter" ].push_back( pointMassGravityAcceleration( ) );

    accelerationMap[ spacecraftName ] = accelerationsOfVehicle;

    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToIntegrate;
    std::string centralBody = "Mars";
    std::vector< std::string > centralBodies = { centralBody };
    bodiesToIntegrate.push_back( spacecraftName );

    // Create acceleration models
    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToIntegrate, centralBodies );
    std::cout<<"acceleration map created"<<std::endl;

    // Create dependent variables
    std::string headers = "";//"Time [s]     Semi-major axis [m]     Eccentricity [-]    Inclination [rad]   Argument of periapsis [rad]     Longitude of ascending node [rad]   True anomaly [rad]  Aerodynamic coefficients 1 [-]  Aerodynamic coefficients 2 [-]  Aerodynamic coefficients 1 [-]      Drag accelration norm [m/s^2]     Spherical harmonic gravity acceleration norm [m/s^2]      Phobos point mass gravity acceleration norm [m/s^2]     Deimos point mass gravity acceleration norm [m/s^2]     Jupiter point mass gravity acceleration norm [m/s^2]    Sun point mass gravity acceleration norm [m/s^2]    Radiation pressure acceleration norm [m/s^2]    Polynomial gravity variations acceleration x [m/s^2]    Polynomial gravity variations acceleration y [m/s^2]    Polynomial gravity variationsacceleration z [m/s^2]     Periodic gravity variations acceleration x [m/s^2]      Periodic gravity variations acceleration y [m/s^2]      Periodic gravity variations acceleration z [m/s^2]";

    std::vector<std::shared_ptr<SingleDependentVariableSaveSettings> > dependentVariablesToSave;
    dependentVariablesToSave.push_back(
            std::make_shared<SingleDependentVariableSaveSettings>(
                    keplerian_state_dependent_variable, spacecraftName, centralBody));

    dependentVariablesToSave.push_back(std::make_shared< SingleDependentVariableSaveSettings >(
            aerodynamic_force_coefficients_dependent_variable, spacecraftName, centralBody ));

    dependentVariablesToSave.push_back(
            std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    aerodynamic, spacecraftName, centralBody, 1 ) );
    dependentVariablesToSave.push_back(
            std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    spherical_harmonic_gravity, spacecraftName, centralBody, 1 ) );
    dependentVariablesToSave.push_back(
            std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    point_mass_gravity, spacecraftName, "Phobos", 1 ) );
    dependentVariablesToSave.push_back(
            std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    point_mass_gravity, spacecraftName, "Deimos", 1 ) );
    dependentVariablesToSave.push_back(
            std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    point_mass_gravity, spacecraftName, "Jupiter", 1 ) );
    dependentVariablesToSave.push_back(
            std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    point_mass_gravity, spacecraftName, "Sun", 1 ) );
    dependentVariablesToSave.push_back(
            std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    radiation_pressure, spacecraftName, "Sun", 1 ) );

    // Define the required parameters
    //std::vector< std::pair< int, int > > componentIndices = { {2, 0}, {2, 1}, {2, 2} };
    gravitation::BodyDeformationTypes deformationType = gravitation::polynomial_variation;

    // Create an instance of SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings
    auto saveSettings = std::make_shared< propagators::SingleVariationSphericalHarmonicAccelerationSaveSettings >(
            spacecraftName,
            centralBody,
            deformationType
    );

    // Add the instance to the list
    dependentVariablesToSave.push_back(saveSettings);

    gravitation::BodyDeformationTypes deformationTypePer = gravitation::periodic_variation;

    // Create an instance of SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings
    auto saveSettingsPer = std::make_shared< propagators::SingleVariationSphericalHarmonicAccelerationSaveSettings >(
            spacecraftName,
            centralBody,
            deformationTypePer
    );
    dependentVariablesToSave.push_back(saveSettingsPer);
    // Define the required parameters
    /*
    std::vector< std::pair< int, int > > componentIndicesPer = { {2, 0}, {2, 1}, {3, 0}, {4, 0}, {5, 0}};
    gravitation::BodyDeformationTypes deformationTypePer = gravitation::periodic_variation;

    // Create an instance of SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings
    auto saveSettingsPer = std::make_shared< propagators::SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings >(
            spacecraftName,
            centralBody,
            componentIndicesPer,
            deformationTypePer
    );

    // Add the instance to the list
    dependentVariablesToSave.push_back(saveSettingsPer);
    */
    dependentVariablesToSave.push_back(
            std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    point_mass_gravity, spacecraftName, "Earth", 1 ) );

    std::cout<<"dependent variables created"<<std::endl;

    // Create integrator times
    std::vector< double > integrationArcStartTimes;
    std::vector< double > integrationArcEndTimes;
    std::vector< double > integrationArcLimits;

    double integrationStartTime = initialEphemerisTime + 120.0; //1.0E2;
    double integrationEndTime = finalEphemerisTime - 120.0 ; //1.0E2;
    double step_size;
    if (dragEst == "per-rev") {
        step_size = hoursperdaydrag * 3600;
    }
    else
    {
        step_size = ndays * 24 * 3600;
    }
    std::cout<<"step size: "<<step_size<<std::endl;
    //double step_size = ndays * 24 * 3600;
    std::vector< double > initial_times_list_drag;
    // Generate the times for drag coeffs
    for (double time = integrationStartTime; time < integrationEndTime; time += step_size) {
        if (integrationEndTime-time < step_size) {
            break;
        }
        initial_times_list_drag.push_back(time);
        std::cout<<"time drag: "<<time<<std::endl;
    }

    std::cout<<"integration start time: "<<integrationStartTime<<std::endl;
    std::cout<<"integration end time: "<<integrationEndTime<<std::endl;
    double totalDurationIntegration = integrationEndTime - integrationStartTime;

    double arcOverlap = 240.0;//2.0E2;

    double currentStartTime = integrationStartTime, currentEndTime = integrationStartTime + arcDuration;
    //integrationArcLimits.push_back( currentStartTime );
    //integrationArcEndTimes.push_back( currentEndTime );
    //integrationArcStartTimes.push_back( currentStartTime );
    std::cout<<"current start time: "<<currentStartTime<<std::endl;
    std::cout<<"current end time: "<<currentEndTime<<std::endl;
    do
    {
        integrationArcLimits.push_back( currentStartTime );
        integrationArcEndTimes.push_back( currentEndTime );
        integrationArcStartTimes.push_back( currentStartTime );
        currentStartTime = currentEndTime - arcOverlap;
        std::cout<<"current start time: "<<currentStartTime<<std::endl;
        currentEndTime = currentStartTime + arcDuration;
        std::cout<<"current end time: "<<currentEndTime<<std::endl;
    }while( currentEndTime <= integrationEndTime );

    integrationArcLimits.push_back( currentStartTime + arcOverlap );
    std::cout<<"arc times created"<<std::endl;
    // Select observation times
    // Compute observed observations. NOTE: don't move this to after the creation of the OrbitDeterminationManager!
    std::vector< Time > observationTimes;
    std::vector< Eigen::Matrix< long double, Eigen::Dynamic, 1 > > observations;
    for ( Time t = integrationStartTime + 600; t < integrationEndTime - 600; t += observationsSamplingTime )
    {
        try
        {
            observations.push_back( bodies.getBody( spacecraftName )->getStateInBaseFrameFromEphemeris< long double, Time >( t ).segment( 0, 3 ) );
            observationTimes.push_back( t );
        }
        catch( ... )
        { }
    }
    // Retrieve state history from SPICE
/*
    std::map< long double, Eigen::Matrix < long double, Eigen::Dynamic, 1 > > spiceStateHistory;
    for ( Time t : observationTimes )
    {
        spiceStateHistory[ t.getSeconds< long double >() ] =
                bodies.getBody( spacecraftName )->getStateInBaseFrameFromEphemeris< long double, Time >( t ) -
                bodies.getBody( centralBody )->getStateInBaseFrameFromEphemeris< long double, Time >( t );
    }
    writeDataMapToTextFile( spiceStateHistory, "stateHistorySpice_" + fileTag + ".txt", saveDirectory,
                            "", 18, 18 );
*/
    // Define integrator settings.
    //std::shared_ptr< IntegratorSettings< double > > integratorSettings =
    //        rungeKutta4Settings< double >( 60.0 );
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
                    ( rungeKutta4, integrationArcStartTimes.at(0), 60.0 );

    std::cout<<"Integration settings created"<<std::endl;

    // start global propagation
    Eigen::Matrix< double, 6, 1 > spacecraftInitialState;
    spacecraftInitialState << -1058004.86077880859,116366.372375488281, -3595420.33442687988, -3113.27424413069639, 950.899628500708786, 938.480820752637783;
           // [-1058004.86077880859,116366.372375488281	 -3595420.33442687988	 -3113.27424413069639	 950.899628500708786	 938.480820752637783];
  /*  Eigen::Matrix< double, 6, 1 > spacecraftInitialState =
            bodies.getBody( spacecraftName )->getStateInBaseFrameFromEphemeris< double, Time >( integrationArcStartTimes[0] ) -
            bodies.getBody( centralBody )->getStateInBaseFrameFromEphemeris< double, Time >( integrationArcStartTimes[0] );
*/
    std::cout<<"spacecraft initial  state created"<<std::endl;

    // Create termination settings
    std::shared_ptr< PropagationTerminationSettings > terminationSettings = propagationTimeTerminationSettings(
            integrationEndTime );
    std::cout<<"termination setting created"<<std::endl;
    // Create propagation settings
    std::shared_ptr< TranslationalStatePropagatorSettings< double> > propagatorSettings = translationalStatePropagatorSettings< double>( centralBodies, accelerationModelMap, bodiesToIntegrate,
                                                                                                                                                          spacecraftInitialState, integrationArcStartTimes[0], integratorSettings, terminationSettings, cowell, dependentVariablesToSave);
    std::cout<<"propagation settings created"<<std::endl;

    SingleArcDynamicsSimulator< > dynamicsSimulator(
            bodies, propagatorSettings );
    std::cout<<"dynamic simulator"<<std::endl;

    std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableResult = dynamicsSimulator.getDependentVariableHistory( );

    writeDataMapToTextFile( integrationResult, "stateHistoryPropagation_" + fileTag + ".txt", saveDirectory,
                            "", 18, 18 );
    writeDataMapToTextFile( dependentVariableResult, "dependentVariablesPropagation_" + fileTag + ".txt", saveDirectory,
                            headers, 18, 18 );

    // create multi arc propagation settings
    int numberOfIntegrationArcs = integrationArcStartTimes.size( );
    for( unsigned int i = 0; i < numberOfIntegrationArcs; i++ ) {
        std::cout << "integration arc start times"<<integrationArcStartTimes.at(i) << std::endl;
        std::cout << "integration arc end times"<<integrationArcEndTimes.at(i) << std::endl;
    }
    std::cout<<"number of integration arcs: "<<numberOfIntegrationArcs<<std::endl;
    //std::vector< Eigen::VectorXd > systemInitialStates(numberOfIntegrationArcs, Eigen::VectorXd(6));
    std::map< double, Eigen::VectorXd > stateVectorsAtStartTimes;

    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > arcPropagationSettingsList;
    for( unsigned int i = 0; i < numberOfIntegrationArcs; i++ )
    {
        std::cout<<"iteration started"<<std::endl;
        std::cout<<i<<std::endl;
        std::cout<<bodiesToIntegrate[ 0 ]<<std::endl;
        std::cout<<integrationArcStartTimes.at(i)<<std::endl;
        // get system initial states from the global propagation
        // Check if the start time exists in the integration results
        if (integrationResult.find( integrationArcStartTimes.at(i)) != integrationResult.end())
        {
            // Store the state vector at the start time
            stateVectorsAtStartTimes[ integrationArcStartTimes.at(i)] = integrationResult[ integrationArcStartTimes.at(i)];
        }
        else
        {
            std::cerr << "Start time " <<  integrationArcStartTimes.at(i) << " not found in integration results." << std::endl;
        }
        //systemInitialStates[ i ]  = spice_interface::getBodyCartesianStateAtEpoch(
        //        bodiesToIntegrate[ 0 ], "Mars", "MARSIAU", "NONE", integrationArcStartTimes.at(i));
        std::cout<<"system initial states created"<<std::endl;
        arcPropagationSettingsList.push_back(
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                        ( centralBodies, accelerationModelMap, bodiesToIntegrate,
                          stateVectorsAtStartTimes[ integrationArcStartTimes.at(i)], integrationArcEndTimes.at( i ), cowell, dependentVariablesToSave, TUDAT_NAN ) );
    }

    std::cout<<"single arc propagation done"<<std::endl;
    std::shared_ptr< MultiArcPropagatorSettings< double > > multiArcPropagatorSettings =
            validateDeprecatedMultiArcSettings< double, double >(
                    integratorSettings, std::make_shared< MultiArcPropagatorSettings< double > >( arcPropagationSettingsList ),
                    integrationArcStartTimes, false, true );



}

