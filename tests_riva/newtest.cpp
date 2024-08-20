//
// Created by Riva Alkahal on 06/08/2024.
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

    spice_interface::loadStandardSpiceKernels();
    spice_interface::loadSpiceKernelInTudat("/Users/ralkahal/estimate/sod_assignments/mgs_map4.bsp");
    spice_interface::loadSpiceKernelInTudat("/Users/ralkahal/estimate/sod_assignments/mgs_map5.bsp");
    spice_interface::loadSpiceKernelInTudat("/Users/ralkahal/estimate/sod_assignments/mgs_map6.bsp");
    spice_interface::loadSpiceKernelInTudat("/Users/ralkahal/estimate/sod_assignments/mgs_map7.bsp");
    spice_interface::loadSpiceKernelInTudat("/Users/ralkahal/estimate/sod_assignments/mgs_map8.bsp");

    // Set output directory
    std::string saveDirectory = "/Users/ralkahal/OneDrive - Delft University of Technology/PhD/Programs/testingestimation/";
    std::string fileTag = "polygrav";

    // Set simulation time settings
    Time initialEphemerisTime = Time(1.0E7);//Time( 0.0 );
    Time finalEphemerisTime = Time(3.0E7);//Time( 86400.0 * 8.0 ); // 5 years later
    double totalDuration = finalEphemerisTime - initialEphemerisTime;
    double ephemeridesTimeStep = 3600;//60.0;
    bool useInterpolatedEphemerides = true;
    double maximumTimeStep = 3600;//60.0;
    double hoursperday = 10.0;
    Time buffer = Time(10.0 * maximumTimeStep);

    // Set arcs, integrator, and observation times
    std::vector< double > integrationArcStartTimes;
    std::vector< double > integrationArcEndTimes;
    std::vector< double > integrationArcLimits;

    double integrationStartTime = initialEphemerisTime + 1.0E4;
    double integrationEndTime = finalEphemerisTime - 1.0E4;
    double arcDuration = 1.6E7; //86400.0*4.0;
    double arcOverlap = 2.0E4;

    // arc times
    double currentStartTime = integrationStartTime, currentEndTime = integrationStartTime + arcDuration;
    do
    {
        integrationArcLimits.push_back( currentStartTime );
        integrationArcEndTimes.push_back( currentEndTime );
        integrationArcStartTimes.push_back( currentStartTime );
        currentStartTime = currentEndTime - arcOverlap;
        currentEndTime = currentStartTime + arcDuration;
    }
    while( currentEndTime < integrationEndTime );
    integrationArcLimits.push_back( currentStartTime + arcOverlap );

    // observations times
    int days = static_cast<int>(totalDuration / physical_constants::JULIAN_DAY );
    std::vector<double> initial_times_list = { integrationArcStartTimes[0] };
    for (int day = 1; day < days; ++day) {
        initial_times_list.push_back(integrationArcStartTimes[0] + day * 24 * 3600);
    }
    std::vector<double> final_times_list = { integrationArcStartTimes[0] + hoursperday * 3600 };
    for (int day = 1; day < days; ++day) {
        final_times_list.push_back(integrationArcStartTimes[0] + day * 24 * 3600 + hoursperday * 3600);
    }
    std::vector<double> observationTimesList;
    for (int i = 0; i < days; ++i) {
        double start = initial_times_list[i];
        double end = final_times_list[i];
        for (double time = start; time < end; time += 3600.0) {
            observationTimesList.push_back(time);
        }
    }
    std::cout<< "Observation times created" << std::endl;
    // print observation times
    std::cout << "Observation Times:" << std::endl;
    for (const auto& time : observationTimesList) {
        std::cout << time << std::endl;
    }
    // Set up the environment settings
    // Create bodies
    std::vector< std::string > bodiesToCreate = {
            "Earth", "Sun", "Mercury", "Venus", "Mars", "Jupiter", "Phobos", "Deimos" };
    std::string baseFrameOrientation = "MARSIAU";
    std::string baseFrameOrigin = "SSB";

    // Create body settings
    BodyListSettings bodySettings = getDefaultBodySettings(
            bodiesToCreate, initialEphemerisTime - buffer, finalEphemerisTime + buffer,
            baseFrameOrigin, baseFrameOrientation, ephemeridesTimeStep );
    //BodyListSettings bodySettings =
    //        getDefaultBodySettings( bodiesToCreate, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    // at Earth
    bodySettings.at( "Earth" )->groundStationSettings = getDsnStationSettings( );
    // at Mars
    std::string filename ="/Users/ralkahal/OneDrive - Delft University of Technology/PhD/Programs/atmodensitydtm/dtm_mars";;
    bodySettings.at( "Mars" )->atmosphereSettings = marsDtmAtmosphereSettings( filename, 3378.0E3);
    std::cout << "Mars atmosphere settings created" << std::endl;
    // Create spacecraft
    std::string spacecraftName = "MGS";
    bodySettings.addSettings( spacecraftName );
    bodySettings.at( spacecraftName )->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
            initialEphemerisTime - buffer, finalEphemerisTime + buffer,
            maximumTimeStep, baseFrameOrigin, baseFrameOrientation,
            std::make_shared< interpolators::LagrangeInterpolatorSettings >( 6 ), spacecraftName );

    bodySettings.at( spacecraftName )->constantMass = 700.0;
    std::cout << "Spacecraft settings created" << std::endl;
    // Create radiation pressure settings
    double referenceAreaRadiation = 20.0;
    double radiationPressureCoefficient = 1.0;
    std::vector< std::string > occultingBodies = { "Mars" };
    bodySettings.at( spacecraftName )->radiationPressureSettings[ "Sun" ] =
            cannonBallRadiationPressureSettings( "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );
    std::cout<< "Radiation pressure settings created" << std::endl;
    // Create aerodynamic settings
    std::shared_ptr<AerodynamicCoefficientSettings> aerodynamicCoefficientSettings =
            std::make_shared<ConstantAerodynamicCoefficientSettings>(
                    2.0, 4.0, Eigen::Vector3d::Zero( ), Eigen::Vector3d::UnitX( ), Eigen::Vector3d::Zero( ),
                    negative_aerodynamic_frame_coefficients, negative_aerodynamic_frame_coefficients );
    bodySettings.at( spacecraftName )->aerodynamicCoefficientSettings = aerodynamicCoefficientSettings;

    // Create gravity field variations
    std::vector< std::shared_ptr< GravityFieldVariationSettings > > gravityFieldVariations;
    std::map<int, Eigen::MatrixXd> cosineAmplitudes;
    cosineAmplitudes[ 1 ] = Eigen::Matrix< double, 3, 5 >::Zero( );
    std::map<int, Eigen::MatrixXd> sineAmplitudes;
    sineAmplitudes[ 1 ] = Eigen::Matrix< double, 3, 5 >::Zero( );
    std::shared_ptr< GravityFieldVariationSettings > polynomialGravityFieldVariations =
            std::make_shared< PolynomialGravityFieldVariationsSettings >(
                    cosineAmplitudes, sineAmplitudes, -1.0E4, 2, 0 );
    gravityFieldVariations.push_back( polynomialGravityFieldVariations );
    std::vector< std::shared_ptr< GravityFieldVariationSettings > > gravityFieldVariationSettings =
            gravityFieldVariations;
    bodySettings.at( "Mars" )->gravityFieldVariationSettings = gravityFieldVariations;

    // Create body map
    SystemOfBodies bodies = createSystemOfBodies< long double, Time >( bodySettings );
    std::cout<< "Bodies created" << std::endl;
    // Create acceleration map
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Sun" ].push_back( pointMassGravityAcceleration( ) );
    accelerationsOfVehicle[ "Sun" ].push_back( cannonBallRadiationPressureAcceleration( ) );
    accelerationsOfVehicle[ "Sun" ].push_back( relativisticAccelerationCorrection(  ) );
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
    bodiesToIntegrate.push_back( spacecraftName );
    unsigned int numberOfNumericalBodies = bodiesToIntegrate.size( );
    std::vector< std::string > bodiesToEstimate;
    bodiesToEstimate.push_back( spacecraftName );
    std::string centralBody = "Mars";
    std::vector< std::string > centralBodies = { centralBody };

    // Create acceleration models
    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToIntegrate, centralBodies );
    std::cout<< "Acceleration models created" << std::endl;
    // Propagation settings
    // Define integrator settings
    std::shared_ptr< IntegratorSettings< double > > integratorSettings =
            rungeKutta4Settings< double >( 3600.0 );

    // Define propagator settings
    // Get initial state
    double marsGravitationalParameter =  bodies.at( "Mars" )->getGravityFieldModel( )->getGravitationalParameter( );
    //std::vector< Eigen::VectorXd > systemInitialStates;
    std::vector< Eigen::Vector6d > initialKeplerElements;
    std::cout<<integrationArcStartTimes[0]<<std::endl;
    // Retrieve state history from SPICE
    std::map< long double, Eigen::Matrix < long double, Eigen::Dynamic, 1 > > spiceStateHistory;
    for ( Time t : observationTimesList )
    {
        spiceStateHistory[ t.getSeconds< long double >() ] =
                bodies.getBody( spacecraftName )->getStateInBaseFrameFromEphemeris< long double, Time >( t ) -
                bodies.getBody( centralBody )->getStateInBaseFrameFromEphemeris< long double, Time >( t );
    }
    writeDataMapToTextFile( spiceStateHistory, "stateHistorySpice_" + fileTag + ".txt", saveDirectory,
                            "", 18, 18 );
    // Create propagator settings
   // Eigen::Matrix< long double, 6, 1 > systemInitialState = bodies.getBody( spacecraftName )->getStateInBaseFrameFromEphemeris< long double, double >( integrationArcStartTimes[0] ) -
    //        bodies.getBody( centralBody )->getStateInBaseFrameFromEphemeris< long double, Time >( integrationArcStartTimes[0] );
    //propagatorSettingsList.push_back(
    //> systemInitialState = getBodyCartesianStateAtEpoch(spacecraftName, "Mars", "MARSIAU", "NONE", integrationArcStartTimes[0] );
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double, double > > > propagatorSettingsList;
    for( unsigned int i = 0; i < integrationArcStartTimes.size( ); i++ ) {
        std::cout<< i << std::endl;
        Eigen::Matrix< double, 6, 1 > systemInitialStates = bodies.getBody( spacecraftName )->getStateInBaseFrameFromEphemeris< double, double >( integrationArcStartTimes[i] ) -
                bodies.getBody( centralBody )->getStateInBaseFrameFromEphemeris< double, double >( integrationArcStartTimes[i] );
        //systemInitialStates[i] = spice_interface::getBodyCartesianStateAtEpoch(
        //        spacecraftName, "Mars", "MARSIAU", "NONE", integrationArcStartTimes[i] );
        std::cout<< "Initial state created" << std::endl;

        propagatorSettingsList.push_back(
                std::make_shared< TranslationalStatePropagatorSettings< double, double > >
                        ( centralBodies, accelerationModelMap, bodiesToIntegrate,
                          systemInitialStates,
                          integrationArcStartTimes[i],
                          integratorSettings,
                          propagationTimeTerminationSettings( integrationArcEndTimes[i])));
           }

    std::cout<< "Propagator settings created" << std::endl;

    std::shared_ptr< MultiArcPropagatorSettings< double, double > > propagatorSettings =
           std::make_shared< MultiArcPropagatorSettings< double, double > >( propagatorSettingsList, 0 );
    /*std::shared_ptr< MultiArcPropagatorSettings< double > > multiArcPropagatorSettings =
            validateDeprecatedMultiArcSettings< double, double >(
                    integratorSettings, std::make_shared< MultiArcPropagatorSettings< double > >( propagatorSettingsList ),
                    integrationArcStartTimes, true, true );
                    */
    std::cout<< "Multi arc propagator settings created" << std::endl;
    std::shared_ptr< PropagationPrintSettings > multiArcPrintSettings =
            std::make_shared< PropagationPrintSettings >( );
    multiArcPrintSettings->reset(
            true, true, TUDAT_NAN, 0, true, true, true, true, true, true );

    propagatorSettings->getOutputSettings( )->resetAndApplyConsistentSingleArcPrintSettings(
            multiArcPrintSettings );
    std::cout<< "Multi arc print settings created" << std::endl;
    //multiArcPropagatorSettings->getOutputSettings( )->resetAndApplyConsistentSingleArcPrintSettings(
    //        multiArcPrintSettings );

    MultiArcDynamicsSimulator< > dynamicsSimulator(
            bodies, propagatorSettings);

    // Retrieve state history
    std::map< long double, Eigen::Matrix < long double, Eigen::Dynamic, 1 > > propagatedStateHistory;
    for ( Time t : observationTimesList )
    {
        propagatedStateHistory[ t.getSeconds< long double >() ] =
                bodies.getBody( spacecraftName )->getStateInBaseFrameFromEphemeris< long double, Time >( t ) -
                bodies.getBody( centralBody )->getStateInBaseFrameFromEphemeris< long double, Time >( t );
    }
    writeDataMapToTextFile( propagatedStateHistory, "stateHistoryPropagation_" + fileTag + ".txt", saveDirectory,
                            "", 18, 18 );


    }