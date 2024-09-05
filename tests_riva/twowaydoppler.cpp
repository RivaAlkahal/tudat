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
    spice_interface::loadSpiceKernelInTudat( "/Users/ralkahal/estimate/sod_assignments/mgs_map4.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/Users/ralkahal/estimate/sod_assignments/mgs_map5.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/Users/ralkahal/estimate/sod_assignments/mgs_map6.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/Users/ralkahal/estimate/sod_assignments/mgs_map7.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/Users/ralkahal/estimate/sod_assignments/mgs_map8.bsp" );

    std::string saveDirectory = "/Users/ralkahal/OneDrive - Delft University of Technology/PhD/Programs/testingestimation/";
    std::string fileTag = "polygrav";

    // set input options
    double epehemeridesTimeStep = 60.0;
    bool useInterpolatedEphemerides = true;
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
    Time finalEphemerisTime = Time( 86400.0 * 120.0 ); // 5 years later
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
    bodySettings.at( "Mars" )->atmosphereSettings = marsDtmAtmosphereSettings( filename, 3378.0E3);

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
        bodySettings.at( spacecraftName )->ephemerisSettings =
                std::make_shared< DirectSpiceEphemerisSettings >( baseFrameOrigin, baseFrameOrientation );
    }
    bodySettings.at( spacecraftName )->constantMass = 700.0;
    bodySettings.at( spacecraftName )->ephemerisSettings->resetMakeMultiArcEphemeris( true );

    // Set gravity field variations
    std::vector< std::shared_ptr< GravityFieldVariationSettings > > gravityFieldVariations;

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
            ( Eigen::MatrixXd( 4, 3 )<<0.53E-9, 0.13E-9, -0.51E-10,
                    0.32E-9).finished( ) );

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
            ( Eigen::MatrixXd( 4, 1 )<<0.46E-9, 0.15E-9, -0.64E-10,
                    -0.02E-09 ).finished( ) );

    sineShAmplitudesCosineTime.push_back(Eigen::MatrixXd::Zero( 4, 1 ));
    sineShAmplitudesCosineTime.push_back(Eigen::MatrixXd::Zero( 4, 1 ));
    sineShAmplitudesCosineTime.push_back(Eigen::MatrixXd::Zero( 4, 1 ));

    sineShAmplitudesSineTime.push_back(Eigen::MatrixXd::Zero( 4, 1 ));
    sineShAmplitudesSineTime.push_back(Eigen::MatrixXd::Zero( 4, 1 ));
    sineShAmplitudesSineTime.push_back(Eigen::MatrixXd::Zero( 4, 1 ));
    frequencies.resize( 3 );
    frequencies = { 2*mathematical_constants::PI/(686.98*86400.0), 4*mathematical_constants::PI/(686.98*86400.0), 6*mathematical_constants::PI/(686.98*86400.0) };

    std::shared_ptr< GravityFieldVariationSettings > periodicGravityFieldVariations =
            std::make_shared< PeriodicGravityFieldVariationsSettings >(
                    cosineShAmplitudesCosineTime, cosineShAmplitudesSineTime, sineShAmplitudesCosineTime, sineShAmplitudesSineTime,
                    frequencies, 0.0, 2, 0 );

    gravityFieldVariations.push_back( periodicGravityFieldVariations );
    std::cout<<"periodic gravity field variation created"<<std::endl;

    // Set polynomial gravity field variation
    std::map<int, Eigen::MatrixXd> cosineAmplitudes;
    cosineAmplitudes[ 1 ] = Eigen::Matrix< double, 4, 3 >::Zero( );
    //cosineAmplitudes[1](0,0) += 6.469163469325e-09;
    //cosineAmplitudes[1](0,1) += 1.63229205379438e-10;
    //cosineAmplitudes[1](0,2) += 7.66708257544063e-09;
    //cosineAmplitudes[ 2 ] = Eigen::Matrix< double, 3, 5 >::Zero( );
    std::map<int, Eigen::MatrixXd> sineAmplitudes;
    //sineAmplitudes[ 1 ] = Eigen::Matrix< double, 4, 3 >::Zero( );
    std::shared_ptr< GravityFieldVariationSettings > polynomialGravityFieldVariations =
            std::make_shared< PolynomialGravityFieldVariationsSettings >(
                    cosineAmplitudes, sineAmplitudes, 0.0, 2, 0 );
    gravityFieldVariations.push_back( polynomialGravityFieldVariations );

    std::vector< std::shared_ptr< GravityFieldVariationSettings > > gravityFieldVariationSettings =
            gravityFieldVariations;
    bodySettings.at( "Mars" )->gravityFieldVariationSettings = gravityFieldVariations;

    SystemOfBodies bodies = createSystemOfBodies< long double, Time >( bodySettings );

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
    std::vector< std::pair< int, int > > componentIndices = { {2, 0}, {2, 1}, {2, 2} };
    gravitation::BodyDeformationTypes deformationType = gravitation::polynomial_variation;

    // Create an instance of SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings
    auto saveSettings = std::make_shared< propagators::SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings >(
            spacecraftName,
            centralBody,
            componentIndices,
            deformationType
    );

    // Add the instance to the list
    dependentVariablesToSave.push_back(saveSettings);

    // Define the required parameters
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
    std::map< long double, Eigen::Matrix < long double, Eigen::Dynamic, 1 > > spiceStateHistory;
    for ( Time t : observationTimes )
    {
        spiceStateHistory[ t.getSeconds< long double >() ] =
                bodies.getBody( spacecraftName )->getStateInBaseFrameFromEphemeris< long double, Time >( t ) -
                bodies.getBody( centralBody )->getStateInBaseFrameFromEphemeris< long double, Time >( t );
    }
    writeDataMapToTextFile( spiceStateHistory, "stateHistorySpice_" + fileTag + ".txt", saveDirectory,
                            "", 18, 18 );

    // Define integrator settings.
    //std::shared_ptr< IntegratorSettings< double > > integratorSettings =
    //        rungeKutta4Settings< double >( 60.0 );
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
                    ( rungeKutta4, integrationArcStartTimes.at(0), 60.0 );

    std::cout<<"Integration settings created"<<std::endl;

    // start global propagation
    Eigen::Matrix< double, 6, 1 > spacecraftInitialState =
            bodies.getBody( spacecraftName )->getStateInBaseFrameFromEphemeris< double, Time >( integrationArcStartTimes[0] ) -
            bodies.getBody( centralBody )->getStateInBaseFrameFromEphemeris< double, Time >( integrationArcStartTimes[0] );

    // Create termination settings
    std::shared_ptr< PropagationTerminationSettings > terminationSettings = propagationTimeTerminationSettings(
            integrationEndTime );

    // Create propagation settings
    std::shared_ptr< TranslationalStatePropagatorSettings< double, double> > propagatorSettings = translationalStatePropagatorSettings< double, double >( centralBodies, accelerationModelMap, bodiesToIntegrate,
                                                                                                                                                          spacecraftInitialState, integrationArcStartTimes[0], integratorSettings, terminationSettings, cowell, dependentVariablesToSave);

    SingleArcDynamicsSimulator< > dynamicsSimulator(
            bodies, propagatorSettings );

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


    //std::shared_ptr< PropagationPrintSettings > multiArcPrintSettings =
    //        std::make_shared< PropagationPrintSettings >( );
    //multiArcPrintSettings->reset(
    //        true, true, TUDAT_NAN, 0, true, true, true, true, true, true );

    // multiArcPropagatorSettings->getOutputSettings( )->resetAndApplyConsistentSingleArcPrintSettings(
    //        multiArcPrintSettings );

    // MultiArcDynamicsSimulator< > dynamicsSimulator(
    //         bodies, multiArcPropagatorSettings );

    //std::map< long double, Eigen::Matrix < long double, Eigen::Dynamic, 1 > > propagatedStateHistory;
    //for ( Time t : observationTimes )
    //{
    //    propagatedStateHistory[ t.getSeconds< long double >() ] =
    //           bodies.getBody( spacecraftName )->getStateInBaseFrameFromEphemeris< long double, Time >( t ) -
    //           bodies.getBody( centralBody )->getStateInBaseFrameFromEphemeris< long double, Time >( t );
    //}
    //writeDataMapToTextFile( propagatedStateHistory, "stateHistoryPropagation_" + fileTag + ".txt", saveDirectory,
    //                        "", 18, 18 );


    //initial state?
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
            getInitialMultiArcParameterSettings< double, double  >( multiArcPropagatorSettings, bodies, integrationArcStartTimes );
    //polynomial gravity field
    std::map<int, std::vector<std::pair<int, int> > > cosineBlockIndicesPerPower;
    cosineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 2, 0 ) );
    std::map<int, std::vector<std::pair<int, int> > > sineBlockIndicesPerPower;
    //sineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 2, 1 ) );
    parameterNames.push_back( std::make_shared< PolynomialGravityFieldVariationEstimatableParameterSettings >(
            "Mars", cosineBlockIndicesPerPower, sineBlockIndicesPerPower ) );


    std::map<int, std::vector<std::pair<int, int> > > cosineBlockIndicesPerPeriod;
    //periodic gravity field
    cosineBlockIndicesPerPeriod[ 0 ].push_back( std::make_pair( 2, 0) );
    cosineBlockIndicesPerPeriod[ 0 ].push_back( std::make_pair( 2, 1 ) );
    cosineBlockIndicesPerPeriod[ 0 ].push_back( std::make_pair( 3, 0 ) );
    cosineBlockIndicesPerPeriod[ 0 ].push_back( std::make_pair( 4, 0 ) );
    cosineBlockIndicesPerPeriod[ 0 ].push_back( std::make_pair( 5, 0 ) );
    std::map<int, std::vector<std::pair<int, int> > > sineBlockIndicesPerPeriod;
    parameterNames.push_back( std::make_shared< PeriodicGravityFieldVariationEstimatableParameterSettings >(
            centralBody, cosineBlockIndicesPerPeriod, sineBlockIndicesPerPeriod ) );


    parameterNames.push_back(std::make_shared< ArcWiseDragCoefficientEstimatableParameterSettings >(spacecraftName, initial_times_list_drag ));
    //parameterNames.push_back(std::make_shared< EstimatableParameterSettings >( spacecraftName, constant_drag_coefficient ));

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( parameterNames, bodies );

    std::cout<<"parameters to estimate created"<<std::endl;

    // observations times
    int days = static_cast<int>(totalDuration / physical_constants::JULIAN_DAY );
    std::cout<<"days: "<<days<<std::endl;
    std::cout<<integrationArcStartTimes.size( )<<std::endl;
    std::vector<double> initial_times_list = { integrationArcStartTimes[0] + buffer};
    for (int day = 1; day < days; ++day) {
        initial_times_list.push_back(integrationArcStartTimes[0] + buffer + day * 24 * 3600);
    }
    std::vector<double> final_times_list = { integrationArcStartTimes[0] + buffer + hoursperday * 3600 };
    for (int day = 1; day < days; ++day) {
        final_times_list.push_back(integrationArcStartTimes[0] + day * 24 * 3600 + hoursperday * 3600);
    }
    std::vector<double> observationTimesList;
    for (int i = 0; i < days; ++i) {
        double start = initial_times_list[i];
        std::cout<<"start time: "<<start<<std::endl;
        double end = final_times_list[i];
        std::cout<<"end time: "<<end<<std::endl;
        for (double time = start; time < end; time += 10) {
            observationTimesList.push_back(time);
        }
    }

    std::cout<< "Observation times created" << std::endl;

    // Create link ends

    // Create list of link ends where the ground station is the transmitter and the spacecraft is the receiver
    std::vector< LinkEnds > stationTransmitterLinkEnds;
    std::vector< LinkEnds > downlinkLinkEnds_;
    std::vector< LinkEnds > uplinkLinkEnds_;
    //linkEnds[ receiver ] = spacecraftName;
    std::vector< std::string > GroundStations = {  "DSS-26" , "DSS-42", "DSS-61"};
    for ( std::string groundStation : GroundStations ) {
        //    linkEnds[ transmitter ] = LinkEndId( "Earth", groundStation );
        //    linkEnds[ receiver ] = spacecraftName;
        //
        //}
        // Define link ends for observations.
        LinkEnds linkEnds;
        linkEnds[transmitter] = LinkEndId("Earth", groundStation);
        linkEnds[reflector1] = spacecraftName;
        linkEnds[receiver] = LinkEndId("Earth", groundStation);
        stationTransmitterLinkEnds.push_back( linkEnds );

        LinkEnds uplinkLinkEnds;
        uplinkLinkEnds[transmitter] = LinkEndId("Earth", groundStation);
        uplinkLinkEnds[receiver] = spacecraftName;
        uplinkLinkEnds_.push_back( uplinkLinkEnds );

        LinkEnds downlinkLinkEnds;
        downlinkLinkEnds[receiver] = LinkEndId("Earth", groundStation);
        downlinkLinkEnds[transmitter] = spacecraftName;

        downlinkLinkEnds_.push_back( downlinkLinkEnds );
    }
    // Define (arbitrary) link ends for each observable
    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
    linkEndsPerObservable[ one_way_range ].push_back( downlinkLinkEnds_[ 0 ] );
    linkEndsPerObservable[ one_way_range ].push_back( uplinkLinkEnds_[ 0 ] );
    linkEndsPerObservable[ one_way_range ].push_back( downlinkLinkEnds_[ 1 ] );

    linkEndsPerObservable[ one_way_doppler ].push_back( downlinkLinkEnds_[ 1 ] );
    linkEndsPerObservable[ one_way_doppler ].push_back( uplinkLinkEnds_[ 2 ] );

    linkEndsPerObservable[ two_way_doppler ].push_back( stationTransmitterLinkEnds[ 0 ] );
    linkEndsPerObservable[ two_way_doppler ].push_back( stationTransmitterLinkEnds[ 1 ] );

    std::cout<<"link ends created"<<std::endl;

    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;

        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            observationSettingsList.push_back(
                    std::make_shared< ObservationModelSettings >(
                            currentObservable, currentLinkEndsList.at( i ), std::shared_ptr< LightTimeCorrectionSettings >( ) ) );
        }
    }

    std::cout<<"observation settings created"<<std::endl;

    // Create orbit determination object.
    OrbitDeterminationManager< double, double > orbitDeterminationManager =
            OrbitDeterminationManager< double, double >(
                    bodies, parametersToEstimate,
                    observationSettingsList, multiArcPropagatorSettings );
    std::cout<<"orbit determination manager created"<<std::endl;


    // Create observation viability settings and calculators
    std::vector< std::shared_ptr< ObservationViabilitySettings > > observationViabilitySettings;
    for ( std::string groundStation : GroundStations )
    {
        observationViabilitySettings.push_back( std::make_shared< ObservationViabilitySettings >(
                minimum_elevation_angle, std::make_pair( "Earth", groundStation ), "",
                unit_conversions::convertDegreesToRadians( 15.0 ) ) );
        observationViabilitySettings.push_back( std::make_shared<ObservationViabilitySettings>(
                body_occultation,
                std::make_pair(groundStation, spacecraftName),
                "Mars",
                TUDAT_NAN
        ));
    }
    std::cout<<"Observation viability settings created"<<std::endl;

    std::vector< std::shared_ptr< ObservationViabilityCalculator > > viabilityCalculators;
    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;

        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            std::vector< std::shared_ptr< ObservationViabilityCalculator > > calculators =
                    createObservationViabilityCalculators(bodies, currentLinkEndsList.at(i), currentObservable, observationViabilitySettings);
            viabilityCalculators.insert(viabilityCalculators.end(), calculators.begin(), calculators.end());
        }
    }
    std::cout<<"Observation viability calculator created"<<std::endl;

    // Create noise
    // Create noise functions per observable
/*
    std::map< ObservableType, std::function< double( const double ) > > noiseFunctions;
    noiseFunctions[ one_way_doppler ] =
            std::bind( &utilities::evaluateFunctionWithoutInputArgumentDependency< double, const double >,
                       createBoostContinuousRandomVariableGeneratorFunction(
                               tudat::statistics::normal_boost_distribution, { 0.0, oneWayDopplerNoise }, 0.0 ), std::placeholders::_1 );

    noiseFunctions[ one_way_range ] =
            std::bind( &utilities::evaluateFunctionWithoutInputArgumentDependency< double, const double >,
                       createBoostContinuousRandomVariableGeneratorFunction(
                               tudat::statistics::normal_boost_distribution, { 0.0, rangeNoise }, 0.0 ), std::placeholders::_1 );

    noiseFunctions[ two_way_doppler ] =
            std::bind( &utilities::evaluateFunctionWithoutInputArgumentDependency< double, const double >,
                       createBoostContinuousRandomVariableGeneratorFunction(
                               tudat::statistics::normal_boost_distribution, { 0.0, twoWayDopplerNoise }, 0.0 ), std::placeholders::_1 );
*/
    // Define the noise functions map with the required type
    std::map< ObservableType, std::function< Eigen::VectorXd( const double ) > > noiseFunctions;

// Create the noise functions that return Eigen::VectorXd
    noiseFunctions[ one_way_doppler ] =
            [=](const double input) -> Eigen::VectorXd {
                // Call the original function that returns a double
                double noiseValue = utilities::evaluateFunctionWithoutInputArgumentDependency< double, const double >(
                        createBoostContinuousRandomVariableGeneratorFunction(
                                tudat::statistics::normal_boost_distribution, { 0.0, oneWayDopplerNoise }, 0.0
                        ), input
                );
                // Convert the double to Eigen::VectorXd
                Eigen::VectorXd result(1);
                result(0) = noiseValue;
                return result;
            };

    noiseFunctions[ one_way_range ] =
            [=](const double input) -> Eigen::VectorXd {
                // Call the original function that returns a double
                double noiseValue = utilities::evaluateFunctionWithoutInputArgumentDependency< double, const double >(
                        createBoostContinuousRandomVariableGeneratorFunction(
                                tudat::statistics::normal_boost_distribution, { 0.0, rangeNoise }, 0.0
                        ), input
                );
                // Convert the double to Eigen::VectorXd
                Eigen::VectorXd result(1);
                result(0) = noiseValue;
                return result;
            };

    noiseFunctions[ two_way_doppler ] =
            [=](const double input) -> Eigen::VectorXd {
                // Call the original function that returns a double
                double noiseValue = utilities::evaluateFunctionWithoutInputArgumentDependency< double, const double >(
                        createBoostContinuousRandomVariableGeneratorFunction(
                                tudat::statistics::normal_boost_distribution, { 0.0, twoWayDopplerNoise }, 0.0
                        ), input
                );
                // Convert the double to Eigen::VectorXd
                Eigen::VectorXd result(1);
                result(0) = noiseValue;
                return result;
            };

    std::cout<<"noise functions created"<<std::endl;
    // Wrapper function to convert double to Eigen::VectorXd
    /*std::function<Eigen::VectorXd(const double)> noiseFunctionWrapper =
            [noiseFunctions](const double time) -> Eigen::VectorXd
            {
                double noiseValue = noiseFunctions.at(one_way_doppler)(time);
                Eigen::VectorXd noiseVector(1);
                noiseVector(0) = noiseValue;
                return noiseVector;
            };
    */
    //std::cout<<"noise functions converted to Eigen::VectorXd"<<std::endl;

    Eigen::Matrix< double, Eigen::Dynamic, 1 > initialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< double >( );

    // Create observation model settings
    //std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationModelSettingsList;
    //observationModelSettingsList.push_back( positionObservableSettings( linkEnds ) );

    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput;
    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;
        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
        //std::function< double( const double ) > noiseFunction = noiseFunctions[currentObservable];
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            measurementSimulationInput.push_back(
                    std::make_shared< TabulatedObservationSimulationSettings< > >(
                            currentObservable, currentLinkEndsList[ i ], observationTimesList, receiver, observationViabilitySettings, noiseFunctions[currentObservable]) );
        }
    }

    // Simulate observations.
    std::shared_ptr< ObservationCollection< > > observationsAndTimes = simulateObservations< double, double >(
            measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );


    //maybe add the viability settings here
    //check out the measurment simulation settings


        //std::shared_ptr< ObservationCollection< double, double > > observationsAndTimes = simulateObservationsWithNoise< double, double >(
    //        measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), noiseFunctions, viabilityCalculators );
    //PodInputDataType observationsAndTimes = simulateObservationsWithNoise< double, double >(
    //        measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), noiseFunctions,
    //        viabilityCalculators );
    std::cout<<"observations and times created"<<std::endl;

    Eigen::Matrix< double, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;
    int numberOfParameters = initialParameterEstimate.rows( );
    std::cout<<"number of parameters: "<<numberOfParameters<<std::endl;
    // Perturb initial states
    for( unsigned int i = 0; i < integrationArcStartTimes.size( ); i++ )
    {
        initialParameterEstimate[ 0 + 6 * i ] += 100.0E0;
        initialParameterEstimate[ 1 + 6 * i ] += 100.0E0;
        initialParameterEstimate[ 2 + 6 * i ] += 100.0E0;
        initialParameterEstimate[ 3 + 6 * i ] += 0.005;
        initialParameterEstimate[ 4 + 6 * i ] += 0.005;
        initialParameterEstimate[ 5 + 6 * i ] += 0.005;
    }
    std::cout<<"initial parameters perturbed"<<std::endl;
    //Perturb remaining parameters
    /*for( unsigned int i = integrationArcStartTimes.size( );
         i < static_cast< unsigned int >( initialParameterEstimate.rows( ) ); i++ )
    {
        initialParameterEstimate[ i ] *= ( 1.0 + 1.0E-6 );
    }
*/
    parametersToEstimate->resetParameterValues( initialParameterEstimate );
    printEstimatableParameterEntries( parametersToEstimate );

    // set a priori to the drag coefficients
    const int DIAGONALS = 251;
    double aprioriuncertainty =1.0/(0.1*0.1);

    // Create a 2D vector (matrix) filled with zeros
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(numberOfParameters, numberOfParameters);
    // Fill the matrix with the values of the diagonal
    for (int i = DIAGONALS; i < numberOfParameters; ++i) {
        matrix(i,i) = aprioriuncertainty;
    }
    //print matrix
    std::cout<<matrix<<std::endl;

    // Define estimation input
    std::shared_ptr< EstimationInput< double, double  > > estimationInput =
            std::make_shared< EstimationInput< double, double > >(
                    observationsAndTimes,
                    matrix,
                    std::make_shared< EstimationConvergenceChecker >( iterationNumber ) );
    // Call the function with reintegrateVariationalEquations set to true
    estimationInput->defineEstimationSettings(
            true,  // reintegrateEquationsOnFirstIteration
            true,  // reintegrateVariationalEquations
            true,  // saveDesignMatrix
            true,  // printOutput
            true,  // saveResidualsAndParametersFromEachIteration
            true, // saveStateHistoryForEachIteration
            1.0E8, // limitConditionNumberForWarning
            true   // conditionNumberWarningEachIteration
    );
    std::map< observation_models::ObservableType, double > weightPerObservable;
    weightPerObservable[ one_way_doppler ] = std::pow(oneWayDopplerNoise, -2);
    weightPerObservable[ one_way_range ] = std::pow(rangeNoise, -2);
    weightPerObservable[ two_way_doppler ] = std::pow(twoWayDopplerNoise, -2);

    estimationInput->setConstantPerObservableWeightsMatrix( weightPerObservable );
    std::cout<<"estimation input created"<<std::endl;



    std::shared_ptr< CovarianceAnalysisInput< double, double > > covarianceInput =
            std::make_shared< CovarianceAnalysisInput< double, double > >(
                    observationsAndTimes, matrix );
    std::cout<<"covariance input created"<<std::endl;

    std::shared_ptr< CovarianceAnalysisOutput< double, double > > covarianceOutput = orbitDeterminationManager.computeCovariance(
            covarianceInput );
    std::cout<<"covariance output created"<<std::endl;

    // Perform estimation
    std::shared_ptr< EstimationOutput< double, double > > estimationOutput = orbitDeterminationManager.estimateParameters(
            estimationInput );
    std::cout<<"estimation performed"<<std::endl;

    Eigen::Matrix< double, Eigen::Dynamic, 1 > finalParameters = parametersToEstimate->template getFullParameterValues< double >( );


    // retrieve simulate observations
    std::ofstream outputfile1(saveDirectory + "concatenatedlinkedIdnames" + fileTag + ".txt");
    outputfile1 << std::setprecision(17);
    std::map<tudat::observation_models::ObservableType, std::map<int, std::vector<std::shared_ptr<SingleObservationSet < double, double>>>>> sortedObservationSets = observationsAndTimes->getSortedObservationSets();
    std::vector< double > concatenatedObservationTimes = observationsAndTimes->getConcatenatedTimeVector();
    Eigen::Matrix< double, Eigen::Dynamic, 1 > observationVector = observationsAndTimes->getObservationVector();
    std::vector< LinkEnds >  concatenatedLinkEndIdNames = observationsAndTimes->getConcatenatedLinkEndIdNames();
    std::cout << "Observation Times size:" << concatenatedObservationTimes.size() << std::endl;
    std::cout << "Observation Vector size:" << concatenatedObservationTimes.size() << std::endl;
    std::cout<<"concatenated link end ids: " << std::endl;
    for (const auto &linkEnds : concatenatedLinkEndIdNames) {
        for (const auto &linkEnd : linkEnds) {
            outputfile1 << "Link End ID: " << linkEnd.first << std::endl;
            outputfile1 << "Body: " << linkEnd.second.bodyName_ << ", Point: " << linkEnd.second.stationName_ << std::endl;

        }
    }
    outputfile1.close();
    // Open file to save observations and times
    std::ofstream outputFile(saveDirectory + "observations_and_times_" + fileTag + ".txt");
    outputFile << std::setprecision(17);
    outputFile << "station_id,observable_type,time,observation\n";
    std::cout<<"output file created, starting loop over the sorted observations"<<std::endl;

    for (const auto &observableType: sortedObservationSets) {
        std::cout << "Observable type: " << observableType.first << std::endl;
        for (const auto &stationId: observableType.second) {
            std::cout << "Station ID: " << stationId.first << std::endl;
            for (const auto &obsSetPtr: stationId.second) {
                auto time = obsSetPtr->getObservationTimes();
                auto observation = obsSetPtr->getObservationsVector();

                std::cout << "Observation time: " << time.size() << std::endl;
                std::cout << "Observation: " << observation.size() << std::endl;
                // Save times and observations to file
                for (size_t i = 0; i < time.size(); ++i)
                {
                    outputFile << stationId.first << "," << observableType.first << "," << time[i] << "," << observation[i] << std::endl;
                }
            }
        }
    }
    outputFile.close();

    // save true errors
    Eigen::Matrix< double, Eigen::Dynamic, 1 > TrueError = ( estimationOutput->parameterEstimate_ - truthParameters ).transpose( );
    std::cout<< "normalized design matrix:"<< std::endl;
    // Define the output file
    std::ofstream fe(saveDirectory + "matrix.csv");
    // Write the matrix to the file
    fe << estimationOutput->normalizedDesignMatrix_;
    // Close the file
    fe.close();

    //std::cout<< estimationOutput->normalizedDesignMatrix_ << std::endl;
    //TrueError = ( estimationOutput->parameterEstimate_ - truthParameters ).transpose( );
    std::cout<<"True error: "<<( estimationOutput->parameterEstimate_ - truthParameters ).transpose( )<<std::endl;
    // save formal errors
    Eigen::Matrix<double, Eigen::Dynamic, 1> FormalError = estimationOutput->getFormalErrorVector( );
    std::cout<<"Formal error: "<<estimationOutput->getFormalErrorVector( ).transpose( )<<std::endl;
    std::cout<<"Error ratio: "<<( ( 1.0E-3 * estimationOutput->getFormalErrorVector( ).segment( 0, numberOfParameters ) ).cwiseQuotient(
            estimationOutput->parameterEstimate_ - truthParameters ) ).transpose( )<<std::endl;

    std::cout<<"True and formal errors saved"<<std::endl;

    // retrieve residuals
    Eigen::MatrixXd residualHistory = estimationOutput->getResidualHistoryMatrix( );
    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > parameterHistory = estimationOutput->getParameterHistoryMatrix( );

    Eigen::MatrixXd residualsWithTime;
    residualsWithTime.resize( residualHistory.rows( ), residualHistory.cols( ) + 1 );
    residualsWithTime.rightCols( residualHistory.cols( ) ) = residualHistory;

    for ( unsigned int i = 0; i < observationVector.size( ); ++i )
    {
        residualsWithTime( i, 0 ) = static_cast< Time >( concatenatedObservationTimes.at( i )
        ).getSeconds< long double >();
    }

    std::ofstream file(saveDirectory + "residuals_" + fileTag + ".txt");
    file << std::setprecision( 17 ) << residualsWithTime;
    file.close();

    std::cout<<"residuals saved"<<std::endl;

    // Retrieve estimated state history
    // Retrieve the best iteration simulation results
    //std::shared_ptr< propagators::SimulationResults< long double, double > > postFitSimulationResults =
    //        estimationOutput->getBestIterationSimulationResults( );
    // parametersToEstimate->resetParameterValues(estimationOutput->parameterHistory_.at(estimationOutput->bestIteration_));
    //finalParameters = parametersToEstimate->template getFullParameterValues< double >( );

    std::shared_ptr<tudat::propagators::SimulationResults<double, double>> bestIterationOutput = estimationOutput->simulationResultsPerIteration_.back();//estimationOutput->bestIteration_ );
    //std::cout << "Type of bestIterationOutput: " << typeid(*bestIterationOutput).name() << std::endl;

    auto multiArcResults = std::dynamic_pointer_cast<tudat::propagators::MultiArcSimulationResults<tudat::propagators::SingleArcVariationalSimulationResults, double, double>>(bestIterationOutput);
    auto singleArcResults = multiArcResults->getSingleArcResults();

    // std::vector< std::shared_ptr< propagators::MultiArcSimulationResults<SingleArcSimulationResults, double, Time >>> bestIterationSimulationResultsPerArc = bestIterationOutput->getSingleArcResults();

    std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > concatenatedStateHistoryPostFitDynamic;
    std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > concatenatedDependentVariablesPostFitDynamic;
    std::map< double, Eigen::Matrix< double, 6, 1 > > concatenatedStateHistoryPostFit;
    std::cout<<singleArcResults.size()<<std::endl;
    for ( unsigned int arcIndex = 0; arcIndex < singleArcResults.size(); ++arcIndex )
    {
        // Retrieve the best iteration simulation results for the current arc
        //std::shared_ptr< propagators::SimulationResults< long double, Time > > bestIterationSimulationResults =
        //        singleArcResults[arcIndex];
        //singleArcResults.push_back( bestIterationSimulationResults );

        // Retrieve the state history for the current arc
        std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > arcStateHistoryPostFitDynamic = singleArcResults[arcIndex]->getDynamicsResults( )->getEquationsOfMotionNumericalSolution();
        // Retrieve the dependent variables history for the current arc
        std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > arcDependentVariablesPostFitDynamic = singleArcResults[arcIndex]->getDynamicsResults( )->getDependentVariableHistory( );
        // Retrieve the state history for the current arc
        //  std::map< Time, Eigen::Matrix< long double, Eigen::Dynamic, 1 > > arcStateHistoryPostFitDynamic =
        //        std::dynamic_pointer_cast< SingleArcVariationalSimulationResults< long double, Time > >(
        //               bestIterationSimulationResults )->getDynamicsResults( )->getEquationsOfMotionNumericalSolution();

        // Append the state history to the concatenated map
        for ( auto it = arcStateHistoryPostFitDynamic.begin(); it != arcStateHistoryPostFitDynamic.end(); ++it )
        {
            concatenatedStateHistoryPostFitDynamic[ it->first ] = it->second;
            //std::cout<<"concatenated state history post fit dynamic: "<<concatenatedStateHistoryPostFitDynamic[ it->first ]<<std::endl;

        }
        // Append the dependent variables history to the concatenated map
        for ( auto it = arcDependentVariablesPostFitDynamic.begin(); it != arcDependentVariablesPostFitDynamic.end(); ++it )
        {
            concatenatedDependentVariablesPostFitDynamic[ it->first ] = it->second;
        }
    }

    // Convert the dynamic state history to fixed size
    for ( auto it = concatenatedStateHistoryPostFitDynamic.begin(); it != concatenatedStateHistoryPostFitDynamic.end(); ++it )
    {
        concatenatedStateHistoryPostFit[ it->first ] = it->second;
        //std::cout<<"concatenated state history post fit: "<<concatenatedStateHistoryPostFit[ it->first ]<<std::endl;
    }
    writeDataMapToTextFile( concatenatedStateHistoryPostFit, "stateHistoryPropagatedPostFit_" + fileTag + ".txt", saveDirectory,
                            "", 18, 18 );

    // Write dependent variables to file
    writeDataMapToTextFile( concatenatedDependentVariablesPostFitDynamic, "dependentVariablesPropagatedPostFit_" + fileTag + ".txt", saveDirectory,
                            headers, 18, 18 );

    /*

    std::shared_ptr< interpolators::OneDimensionalInterpolator< Time, Eigen::Matrix< long double, 6, 1 > > > postFitStateInterpolator =
            propagators::createStateInterpolator< Time, long double >( concatenatedStateHistoryPostFit );
    std::map< long double, Eigen::Matrix < long double, Eigen::Dynamic, 1 > > propagatedStateHistoryPostFitToWrite;
    for ( Time t : observationTimes )
    {
        propagatedStateHistoryPostFitToWrite[ t.getSeconds< long double >() ] = postFitStateInterpolator->interpolate( t );
    }
    writeDataMapToTextFile( propagatedStateHistoryPostFitToWrite, "stateHistoryPropagatedPostFit_" + fileTag + ".txt", saveDirectory,
                            "", 18, 18 );
*/
    //parametersToEstimate->resetParameterValues( estimationOutput->parameterHistory_.at( estimationOutput->bestIteration_ ) );
/*
    // Retrieve post-fit state history
    std::shared_ptr< propagators::SimulationResults< double, double > > postFitSimulationResults =
            estimationOutput->getBestIterationSimulationResults( );
    std::map < double, Eigen::Matrix < double, Eigen::Dynamic, 1 > > propagatedStateHistoryPostFitDynamic =
            std::dynamic_pointer_cast< SingleArcVariationalSimulationResults< double, double > >(
                    postFitSimulationResults )->getDynamicsResults( )->getEquationsOfMotionNumericalSolution( );
    std::map < double, Eigen::Matrix < double, 6, 1 > > propagatedStateHistoryPostFit;
    for ( auto it = propagatedStateHistoryPostFitDynamic.begin( ); it != propagatedStateHistoryPostFitDynamic.end( ); ++it )
    {
        propagatedStateHistoryPostFit[ it->first ] = it->second;
    }
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 6, 1 > > > postFitStateInterpolator =
            propagators::createStateInterpolator< double, double >( propagatedStateHistoryPostFit );
    std::map< double, Eigen::Matrix < double, Eigen::Dynamic, 1 > > propagatedStateHistoryPostFitToWrite;
    for ( Time t : observationTimes )
    {
        propagatedStateHistoryPostFitToWrite[ t.getSeconds< double >() ] = postFitStateInterpolator->interpolate( t );
    }
    writeDataMapToTextFile( propagatedStateHistoryPostFitToWrite, "stateHistoryPropagatedPostFit_" + fileTag + ".txt", saveDirectory,
                            "", 18, 18 );

*/
    Eigen::MatrixXd covarianceOutputEstimation = estimationOutput->getCorrelationMatrix();

    // Save covariance matrix
    std::ofstream fileCov(saveDirectory + "covariance_" + fileTag + ".txt");
    fileCov << std::setprecision( 17 ) << covarianceOutputEstimation;
    fileCov.close();



}

