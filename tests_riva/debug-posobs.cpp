//
// Created by Riva Alkahal on 08/10/2024.
//
//
// Created by Riva Alkahal on 07/10/2024.
//
//
// Created by Riva Alkahal on 05/10/2024.
//
//
// Created by Riva Alkahal on 31/08/2024.
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
    spice_interface::loadSpiceKernelInTudat( "/Users/ralkahal/OneDrive - Delft University of Technology/esitmate/sod_assignments/mgs_map4.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/Users/ralkahal/OneDrive - Delft University of Technology/esitmate/sod_assignments/mgs_map5.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/Users/ralkahal/OneDrive - Delft University of Technology/esitmate/sod_assignments/mgs_map6.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/Users/ralkahal/OneDrive - Delft University of Technology/esitmate/sod_assignments/mgs_map7.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/Users/ralkahal/OneDrive - Delft University of Technology/esitmate/sod_assignments/mgs_map8.bsp" );


    spice_interface::loadSpiceKernelInTudat( "/Users/ralkahal/OneDrive - Delft University of Technology/esitmate/sod_assignments/mgs_map3_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/Users/ralkahal/OneDrive - Delft University of Technology/esitmate/sod_assignments/mgs_map4_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/Users/ralkahal/OneDrive - Delft University of Technology/esitmate/sod_assignments/mgs_map5_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/Users/ralkahal/OneDrive - Delft University of Technology/esitmate/sod_assignments/mgs_map6_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/Users/ralkahal/OneDrive - Delft University of Technology/esitmate/sod_assignments/mgs_map7_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/Users/ralkahal/OneDrive - Delft University of Technology/esitmate/sod_assignments/mgs_map8_ipng_mgs95j.bsp" );

    std::string saveDirectory = "/Users/ralkahal/OneDrive - Delft University of Technology/new-tudat-tests/";
    std::string fileTag = "observspice-RK78";

    // set input options
    double epehemeridesTimeStep = 60.0;
    bool useInterpolatedEphemerides = true;
    double observationsSamplingTime = 60.0;
    double buffer = 20.0 * epehemeridesTimeStep;
    double arcDuration = 3.0 * 86400.0;//2.0E4;
    std::string dragEst = "per-halfday";//"per-rev";
    double ndays = 1.0;
    double hoursperdaydrag = 1.0;
    double hoursperday = 10.0;
    int iterationNumber = 5;
    //const double gravitationalParameter = 4.2828378e13;
    //const double planetaryRadius = 3389.5E3;

    double oneWayDopplerNoise = 0.0001;
    double twoWayDopplerNoise = 0.0001;
    double rangeNoise = 0.0001;
    Time initialEphemerisTime = Time(0.0);
    Time finalEphemerisTime = Time(86400.0 * 60.0);
    double totalDuration = finalEphemerisTime - initialEphemerisTime;
    std::cout << "Total duration: " << totalDuration << std::endl;

    std::vector<std::string> bodiesToCreate = {
            "Earth", "Sun", "Mercury", "Venus", "Mars", "Jupiter", "Phobos", "Deimos"};

    std::string baseFrameOrientation = "MARSIAU";
    std::string baseFrameOrigin = "SSB";

    // Specify ephemeris time steps and buffers
    Time ephemerisTimeStepPlanets = epehemeridesTimeStep;
    std::cout << "Ephemeris time step planets: " << epehemeridesTimeStep << std::endl;
    Time ephemerisTimeStepSpacecraft = Time(epehemeridesTimeStep);

    BodyListSettings bodySettings;
    if (useInterpolatedEphemerides) {
        bodySettings = getDefaultBodySettings(
                bodiesToCreate, initialEphemerisTime - buffer, finalEphemerisTime + buffer,
                baseFrameOrigin, baseFrameOrientation, ephemerisTimeStepPlanets);
    } else {
        bodySettings = getDefaultBodySettings(bodiesToCreate, baseFrameOrigin, baseFrameOrientation);
    }
    std::cout << "Body Settings created" << std::endl;

    bodySettings.at("Earth")->groundStationSettings = getDsnStationSettings();

    std::string filename = "/Users/ralkahal/OneDrive - Delft University of Technology/PhD/Programs/atmodensitydtm/dtm_mars";;
    bodySettings.at("Mars")->atmosphereSettings = marsDtmAtmosphereSettings(filename, 3378.0E3);

    // Set spherical harmonics gravity field
    // Create spacecraft
    std::string spacecraftName = "MGS";
    bodySettings.addSettings(spacecraftName);
    if (useInterpolatedEphemerides) {
        bodySettings.at(spacecraftName)->ephemerisSettings =
                std::make_shared<InterpolatedSpiceEphemerisSettings>(
                        initialEphemerisTime - buffer, finalEphemerisTime + buffer,
                        ephemerisTimeStepSpacecraft, baseFrameOrigin, baseFrameOrientation,
                        std::make_shared<interpolators::LagrangeInterpolatorSettings>(8), spacecraftName);
    } else {
        bodySettings.at(spacecraftName)->ephemerisSettings =
                std::make_shared<DirectSpiceEphemerisSettings>(baseFrameOrigin, baseFrameOrientation);
    }
    bodySettings.at(spacecraftName)->constantMass = 700.0;
    //bodySettings.at(spacecraftName)->ephemerisSettings->resetMakeMultiArcEphemeris(false);

    // Set gravity field variations
    std::vector<std::shared_ptr<GravityFieldVariationSettings> > gravityFieldVariations;

    // Set solid body tide gravity field variation
    std::vector<std::string> deformingBodies;
    deformingBodies.push_back("Sun");

    std::map<int, std::vector<std::complex<double> > > loveNumbers;
    std::vector<std::complex<double> > degreeTwoLoveNumbers_;
    degreeTwoLoveNumbers_.push_back(std::complex<double>(0.169, 0.0));
    loveNumbers[2] = degreeTwoLoveNumbers_;

    std::shared_ptr<GravityFieldVariationSettings> singleGravityFieldVariation =
            std::make_shared<BasicSolidBodyGravityFieldVariationSettings>(deformingBodies, loveNumbers);
    gravityFieldVariations.push_back(singleGravityFieldVariation);

    // Set periodic gravity field variation
    /*std::vector<Eigen::MatrixXd> cosineShAmplitudesCosineTime;
    std::vector<Eigen::MatrixXd> cosineShAmplitudesSineTime;
    std::vector<Eigen::MatrixXd> sineShAmplitudesCosineTime;
    std::vector<Eigen::MatrixXd> sineShAmplitudesSineTime;
    std::vector<double> frequencies;

    // Values from Mars GMM3_120_SHA, https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/data/shadr/gmm3_120_sha.lbl
    cosineShAmplitudesCosineTime.push_back(
            (Eigen::MatrixXd(4, 1) << 2.39E-9, 1.67E-9, 0.85E-10,
                    0.38E-9).finished());
    cosineShAmplitudesCosineTime.push_back(
            (Eigen::MatrixXd(4, 1) << 1.23E-9, 0.32E-9, 0.35E-10,
                    0.15E-9).finished());
    cosineShAmplitudesCosineTime.push_back(
            (Eigen::MatrixXd(4, 1) << 0.53E-9, 0.13E-9, -0.51E-10,
                    0.32E-9).finished());

    cosineShAmplitudesSineTime.push_back(
            (Eigen::MatrixXd(4, 1) << -0.83E-9, 2.35E-9, -1.56E-10,
                    1.30E-9).finished());
    cosineShAmplitudesSineTime.push_back(
            (Eigen::MatrixXd(4, 1) << 0.73E-9, 0.21E-9, -0.24E-10,
                    0.42E-9).finished());
    cosineShAmplitudesSineTime.push_back(
            (Eigen::MatrixXd(4, 1) << 0.46E-9, 0.15E-9, -0.64E-10,
                    -0.02E-09).finished());

    sineShAmplitudesCosineTime.push_back(Eigen::MatrixXd::Zero(4, 1));
    sineShAmplitudesCosineTime.push_back(Eigen::MatrixXd::Zero(4, 1));
    sineShAmplitudesCosineTime.push_back(Eigen::MatrixXd::Zero(4, 1));

    sineShAmplitudesSineTime.push_back(Eigen::MatrixXd::Zero(4, 1));
    sineShAmplitudesSineTime.push_back(Eigen::MatrixXd::Zero(4, 1));
    sineShAmplitudesSineTime.push_back(Eigen::MatrixXd::Zero(4, 1));
    frequencies.resize(3);
    frequencies = {2 * mathematical_constants::PI / (686.98 * 86400.0),
                   4 * mathematical_constants::PI / (686.98 * 86400.0),
                   9 * mathematical_constants::PI / (686.98 * 86400.0)};

    std::shared_ptr<GravityFieldVariationSettings> periodicGravityFieldVariations =
            std::make_shared<PeriodicGravityFieldVariationsSettings>(
                    cosineShAmplitudesCosineTime, cosineShAmplitudesSineTime, sineShAmplitudesCosineTime,
                    sineShAmplitudesSineTime,
                    frequencies, 0.0, 2, 0);

    gravityFieldVariations.push_back(periodicGravityFieldVariations);
    std::cout << "periodic gravity field variation created" << std::endl;
*/
    // Set polynomial gravity field variation
    /*std::map<int, Eigen::MatrixXd> cosineAmplitudes;
    cosineAmplitudes[1] = Eigen::Matrix<double, 4, 3>::Zero();
    //cosineAmplitudes[ 2 ] = Eigen::Matrix< double, 3, 5 >::Zero( );
    std::map<int, Eigen::MatrixXd> sineAmplitudes;
    sineAmplitudes[1] = Eigen::Matrix<double, 4, 3>::Zero();
    std::shared_ptr<GravityFieldVariationSettings> polynomialGravityFieldVariations =
            std::make_shared<PolynomialGravityFieldVariationsSettings>(
                    cosineAmplitudes, sineAmplitudes, 0.0, 2, 0);
    gravityFieldVariations.push_back(polynomialGravityFieldVariations);
*/
    std::vector<std::shared_ptr<GravityFieldVariationSettings> > gravityFieldVariationSettings =
            gravityFieldVariations;
    bodySettings.at("Mars")->gravityFieldVariationSettings = gravityFieldVariations;

    SystemOfBodies bodies = createSystemOfBodies<long double, Time>(bodySettings);

    // Create radiation pressure settings
    double referenceAreaRadiation = 4.0;
    double radiationPressureCoefficient = 1.2;
    std::vector<std::string> occultingBodies = {"Mars"};
    std::shared_ptr<RadiationPressureInterfaceSettings> radiationPressureSettings =
            std::make_shared<CannonBallRadiationPressureInterfaceSettings>(
                    "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies);

    // Create and set radiation pressure settings
    bodies.at(spacecraftName)->setRadiationPressureInterface(
            "Sun", createRadiationPressureInterface(
                    radiationPressureSettings, spacecraftName, bodies));

    // Create aerodynamic coefficients settings
    Eigen::Vector3d customVector(1.2, 0.0, 0.0);
    std::shared_ptr<AerodynamicCoefficientSettings> aerodynamicCoefficientSettings =
            std::make_shared<ConstantAerodynamicCoefficientSettings>(4.0, 1.2 * Eigen::Vector3d::UnitX());
    bodies.at(spacecraftName)->setAerodynamicCoefficientInterface(
            createAerodynamicCoefficientInterface(aerodynamicCoefficientSettings, spacecraftName, bodies));


    // Set accelerations on Vehicle that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map<std::string, std::vector<std::shared_ptr<AccelerationSettings> > > accelerationsOfVehicle;
    accelerationsOfVehicle["Sun"].push_back(pointMassGravityAcceleration());
    accelerationsOfVehicle["Sun"].push_back(cannonBallRadiationPressureAcceleration());
    accelerationsOfVehicle["Mercury"].push_back(pointMassGravityAcceleration());
    accelerationsOfVehicle["Venus"].push_back(pointMassGravityAcceleration());
    accelerationsOfVehicle["Earth"].push_back(pointMassGravityAcceleration());
    accelerationsOfVehicle["Mars"].push_back(sphericalHarmonicAcceleration(50, 50));
    //accelerationsOfVehicle["Mars"].push_back(relativisticAccelerationCorrection());
    accelerationsOfVehicle["Mars"].push_back(aerodynamicAcceleration());
    accelerationsOfVehicle["Phobos"].push_back(pointMassGravityAcceleration());
    accelerationsOfVehicle["Deimos"].push_back(pointMassGravityAcceleration());
    accelerationsOfVehicle["Jupiter"].push_back(pointMassGravityAcceleration());

    accelerationMap[spacecraftName] = accelerationsOfVehicle;

    // Set bodies for which initial state is to be estimated and integrated.
    std::vector<std::string> bodiesToIntegrate;
    std::string centralBody = "Mars";
    std::vector<std::string> centralBodies = {centralBody};
    bodiesToIntegrate.push_back(spacecraftName);

    // Create acceleration models
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(bodies, accelerationMap, bodiesToIntegrate,
                                                                       centralBodies);
    std::cout << "acceleration map created" << std::endl;

    // Create dependent variables
    std::vector<std::shared_ptr<SingleDependentVariableSaveSettings> > dependentVariablesToSave;
    dependentVariablesToSave.push_back(
            std::make_shared<SingleDependentVariableSaveSettings>(
                    keplerian_state_dependent_variable, spacecraftName, centralBody));
        dependentVariablesToSave.push_back(std::make_shared<SingleDependentVariableSaveSettings>(
            aerodynamic_force_coefficients_dependent_variable, spacecraftName, centralBody));

    std::cout << "dependent variables created" << std::endl;

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
    for (double time = integrationStartTime + 10000 ; time < integrationEndTime + 10000; time += step_size) {
        //if (integrationEndTime-time < step_size) {
        //    break;
        //}
        initial_times_list_drag.push_back(time);
        std::cout<<"time drag: "<<time<<std::endl;
    }

    std::cout<<"integration start time: "<<integrationStartTime<<std::endl;
    std::cout<<"integration end time: "<<integrationEndTime<<std::endl;
    double totalDurationIntegration = integrationEndTime - integrationStartTime;

    double arcOverlap = 240.0;//2.0E2;

    double currentStartTime = integrationStartTime, currentEndTime = integrationStartTime + arcDuration;

    std::cout<<"current start time: "<<currentStartTime<<std::endl;
    std::cout<<"current end time: "<<currentEndTime<<std::endl;
    do
    {
        integrationArcLimits.push_back( currentStartTime );
        integrationArcEndTimes.push_back( currentEndTime );
        integrationArcStartTimes.push_back( currentStartTime );
        currentStartTime = currentEndTime - arcOverlap;
        currentEndTime = currentStartTime + arcDuration;
    }while( currentEndTime <= integrationEndTime );

    integrationArcLimits.push_back( currentStartTime + arcOverlap );
    std::cout<<"arc times created"<<std::endl;
    // Select observation times
    // Compute observed observations. NOTE: don't move this to after the creation of the OrbitDeterminationManager!
    std::vector< double > observationTimes;
    std::vector< Eigen::Matrix< long double, Eigen::Dynamic, 1 > > observations;
    for ( Time t = integrationStartTime + 600; t < integrationEndTime - 600; t += observationsSamplingTime )
    {
        try
        {
            observations.push_back( bodies.getBody( spacecraftName )->getStateInBaseFrameFromEphemeris< long double, double >( t ).segment( 0, 3 ) );
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
                bodies.getBody( spacecraftName )->getStateInBaseFrameFromEphemeris< long double, double >( t ) -
                bodies.getBody( centralBody )->getStateInBaseFrameFromEphemeris< long double, double >( t );
    }
    writeDataMapToTextFile( spiceStateHistory, "stateHistorySpice_" + fileTag + ".txt", saveDirectory,
                            "", 18, 18 );

    // Define integrator settings
    //std::shared_ptr< IntegratorSettings< > > integratorSettings =
    //        std::make_shared< IntegratorSettings< > >
    //                ( rungeKutta4, integrationStartTime + 600, 10.0 );
    std::shared_ptr<IntegratorSettings<> >integratorSettings =
            std::make_shared<RungeKuttaFixedStepSizeSettings<> >( 30, CoefficientSets::rungeKutta87DormandPrince );

    std::cout<<"Integration settings created"<<std::endl;
    int numberOfIntegrationArcs = integrationArcStartTimes.size( );

    // start global propagation
    Eigen::Matrix<long double, 6, 1 > spacecraftInitialState =
            bodies.getBody( spacecraftName )->getStateInBaseFrameFromEphemeris< long double, Time >( integrationStartTime + 600 ) -
            bodies.getBody( centralBody )->getStateInBaseFrameFromEphemeris< long double, Time >( integrationStartTime + 600 );

    // Create termination settings
    std::shared_ptr< PropagationTerminationSettings > terminationSettings = propagationTimeTerminationSettings(
            integrationEndTime - 600 );

    // Create propagation settings
    std::shared_ptr< TranslationalStatePropagatorSettings< long double, double> > propagatorSettings =
        translationalStatePropagatorSettings< long double, double >( centralBodies, accelerationModelMap, bodiesToIntegrate,
            spacecraftInitialState, integrationStartTime + 600, integratorSettings, terminationSettings, cowell, dependentVariablesToSave);


    //SingleArcDynamicsSimulator< > dynamicsSimulator(
    //        bodies, propagatorSettings );

    //std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    //std::map< double, Eigen::VectorXd > dependentVariableResult = dynamicsSimulator.getDependentVariableHistory( );

    //writeDataMapToTextFile( integrationResult, "stateHistoryPropagation_" + fileTag + ".txt", saveDirectory,
    //                        "", 18, 18 );
    //writeDataMapToTextFile( dependentVariableResult, "dependentVariablesPropagation_" + fileTag + ".txt", saveDirectory,
    //                        "", 18, 18 );

    // Select parameters to estimate
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< long double > >(
            spacecraftName, spacecraftInitialState, centralBody, baseFrameOrientation ) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( spacecraftName, radiation_pressure_coefficient ) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( spacecraftName, constant_drag_coefficient ) );

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > > parametersToEstimate =
            createParametersToEstimate< long double, double >( parameterNames, bodies );

    std::cout<<"parameters to estimate created"<<std::endl;

    // Create link ends
    LinkEnds linkEnds;
    linkEnds[ observed_body ] = spacecraftName;

    // Create observation model settings
    std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationModelSettingsList;
    observationModelSettingsList.push_back( positionObservableSettings( linkEnds ) );

    std::cout<<"link ends created"<<std::endl;

    // Create orbit determination object.
    OrbitDeterminationManager< long double, double > orbitDeterminationManager =
            OrbitDeterminationManager< long double, double >(
                    bodies, parametersToEstimate,
                    observationModelSettingsList, propagatorSettings );
    std::cout<<"orbit determination manager created"<<std::endl;


    // Retrieve state history
    std::map< long double, Eigen::Matrix < long double, Eigen::Dynamic, 1 > > propagatedStateHistory;
    for ( Time t : observationTimes )
    {
        propagatedStateHistory[ t.getSeconds< long double >() ] =
                bodies.getBody( spacecraftName )->getStateInBaseFrameFromEphemeris< long double, Time >( t ) -
                    bodies.getBody( centralBody )->getStateInBaseFrameFromEphemeris< long double, Time >( t );
    }

    writeDataMapToTextFile( propagatedStateHistory, "stateHistoryPropagatedPreFit_" + fileTag + ".txt", saveDirectory,
                            "", 18, 18 );
    writeDataMapToTextFile( spiceStateHistory, "stateHistorySpice_" + fileTag + ".txt", saveDirectory,
                            "", 18, 18 );



    std::vector< std::shared_ptr< SingleObservationSet< long double, double > > > observationSetList;
    observationSetList.push_back(
            std::make_shared< SingleObservationSet<long double, double > >(
                    position_observable, linkEnds, observations, observationTimes, observed_body ) );
    std::shared_ptr< ObservationCollection< long double, double > > observedObservationCollection =
            std::make_shared< ObservationCollection< long double, double > >( observationSetList );

    std::cout<<"observations and times created"<<std::endl;

  // Define estimation input
    std::shared_ptr< EstimationInput< long double, double  > > estimationInput =
            std::make_shared< EstimationInput< long double, double > >(
                    observedObservationCollection,
                    Eigen::MatrixXd::Zero( 0, 0 ),
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

    // Perform estimation
    std::shared_ptr< EstimationOutput< long double, double > > estimationOutput = orbitDeterminationManager.estimateParameters(
                estimationInput );

    // Retrieve post-fit state history
    std::shared_ptr< propagators::SimulationResults< long double, double > > postFitSimulationResults =
            estimationOutput->getBestIterationSimulationResults( );
    std::map < double, Eigen::Matrix < long double, Eigen::Dynamic, 1 > > propagatedStateHistoryPostFitDynamic =
            std::dynamic_pointer_cast< SingleArcVariationalSimulationResults< long double, double > >(
                    postFitSimulationResults )->getDynamicsResults( )->getEquationsOfMotionNumericalSolution( );
    std::map < double, Eigen::Matrix < long double, 6, 1 > > propagatedStateHistoryPostFit;
    for ( auto it = propagatedStateHistoryPostFitDynamic.begin( ); it != propagatedStateHistoryPostFitDynamic.end( ); ++it )
    {
        propagatedStateHistoryPostFit[ it->first ] = it->second;
    }
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< long double, 6, 1 > > > postFitStateInterpolator =
            propagators::createStateInterpolator< double, long double >( propagatedStateHistoryPostFit );
    std::map< long double, Eigen::Matrix < long double, Eigen::Dynamic, 1 > > propagatedStateHistoryPostFitToWrite;
    for ( Time t : observationTimes )
    {
        propagatedStateHistoryPostFitToWrite[ t.getSeconds< long double >() ] = postFitStateInterpolator->interpolate( t );
    }
    //    for ( auto it = propagatedStateHistoryPostFit.begin( ); it != propagatedStateHistoryPostFit.end( ); ++it )
    //    {
    //        propagatedStateHistoryPostFitToWrite[ it->first.getSeconds< long double >() ] = it->second;
    //    }
    writeDataMapToTextFile( propagatedStateHistoryPostFitToWrite, "stateHistoryPropagatedPostFit_" + fileTag + ".txt", saveDirectory,
                            "", 18, 18 );

    // Retrieve residuals and set them in matrix
    Eigen::MatrixXd residualHistory = estimationOutput->getResidualHistoryMatrix( );
    Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic > parameterHistory = estimationOutput->getParameterHistoryMatrix( );

    Eigen::MatrixXd residualsWithTime;
    residualsWithTime.resize( residualHistory.rows( ), residualHistory.cols( ) + 1 );
    residualsWithTime.rightCols( residualHistory.cols( ) ) = residualHistory;

    for ( unsigned int i = 0; i < observedObservationCollection->getObservationVector( ).size( ); ++i )
    {
        residualsWithTime( i, 0 ) = static_cast< Time >( observedObservationCollection->getConcatenatedTimeVector( ).at( i )
                ).getSeconds< long double >();
    }

    std::ofstream file(saveDirectory + "residuals_" + fileTag + ".txt");
    file << std::setprecision( 17 ) << residualsWithTime;
    file.close();

    std::ofstream file3(saveDirectory + "parameters_" + fileTag + ".txt");
    file3 << std::setprecision( 21 ) << parameterHistory;
    file3.close();

    // Retrieve covariance matrix
    Eigen::MatrixXd normalizedCovarianceMatrix = estimationOutput->getNormalizedCovarianceMatrix( );
    Eigen::MatrixXd unnormalizedCovarianceMatrix = estimationOutput->getUnnormalizedCovarianceMatrix( );
    std::ofstream file4(saveDirectory + "covariance_" + fileTag + ".txt");
    file4 << std::setprecision( 17 ) << "Normalized covariance matrix: " << std::endl << normalizedCovarianceMatrix;
    file4 << std::endl << std::endl;
    file4 << "Unnormalized covariance matrix: " << std::endl << unnormalizedCovarianceMatrix;
    file4.close( );

    std::ofstream file2(saveDirectory + "observationsStartAndSize_" + fileTag + ".txt");
    std::map< ObservableType, std::map< int, std::vector< std::pair< int, int > > > > observationSetStartAndSize =
            observedObservationCollection->getObservationSetStartAndSizePerLinkEndIndex();
    for ( auto it = observationSetStartAndSize.begin(); it != observationSetStartAndSize.end(); ++it )
    {
        ObservableType observable = it->first;
        for ( auto it2 = it->second.begin(); it2 != it->second.end(); ++it2 )
        {
            int linkEnd = it2->first;
            for ( unsigned int i = 0; i < it2->second.size(); ++i )
            {
                file2 << std::setprecision( 15 ) << observable << " " << linkEnd << " " << it2->second.at( i ).first <<
                    " " << it2->second.at( i ).second << std::endl;
            }
        }
    }
    file2.close();



}