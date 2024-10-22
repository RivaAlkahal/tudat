// last run on 2024-10-17
// Created by Riva Alkahal on 09/10/2024.
//
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
    std::string fileTag = "observspice-RK78-multiarc-drag-noap-per-arc-emprconst";

    // set input options
    double epehemeridesTimeStep = 60.0;
    bool useInterpolatedEphemerides = true;
    double observationsSamplingTime = 60.0;
    double buffer = 20.0 * epehemeridesTimeStep;
    double arcDuration = 3.0 * 86400.0;//2.0E4;
    std::string dragEst = "per-arc";//"per-rev";
    double ndays = 3.0;
    double hoursperdaydrag = 2.0;
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
    bodySettings.at(spacecraftName)->ephemerisSettings->resetMakeMultiArcEphemeris(true);

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
    double referenceAreaRadiation = 10.0;
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
            std::make_shared<ConstantAerodynamicCoefficientSettings>(10.0, 1.2 * Eigen::Vector3d::UnitX());
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
    accelerationsOfVehicle["Mars"].push_back(relativisticAccelerationCorrection());
    accelerationsOfVehicle["Mars"].push_back(aerodynamicAcceleration());
    accelerationsOfVehicle["Phobos"].push_back(pointMassGravityAcceleration());
    accelerationsOfVehicle["Deimos"].push_back(pointMassGravityAcceleration());
    accelerationsOfVehicle["Jupiter"].push_back(pointMassGravityAcceleration());

    //double empiricalAccelerationNorm = 1.0E-6;
    accelerationsOfVehicle[ "Mars" ].push_back( std::make_shared< EmpiricalAccelerationSettings >(
                                                             Eigen::Vector3d::Zero( ),
                                                             Eigen::Vector3d::Zero( ),
                                                             Eigen::Vector3d::Zero( )) );

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
    dependentVariablesToSave.push_back(std::make_shared<SingleDependentVariableSaveSettings>(
            local_density_dependent_variable, spacecraftName, centralBody));
    dependentVariablesToSave.push_back(std::make_shared<SingleAccelerationDependentVariableSaveSettings>(empirical_acceleration, spacecraftName, centralBody));
    dependentVariablesToSave.push_back(std::make_shared<SingleDependentVariableSaveSettings>(
            radiation_pressure_coefficient_dependent_variable, spacecraftName, centralBody));

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
    for (double time = integrationStartTime ; time < integrationEndTime; time += step_size) {
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
    for ( Time t = integrationStartTime + 600; t < integrationArcEndTimes.back(); t += observationsSamplingTime )
    {
        try
        {
            observations.push_back( bodies.getBody( spacecraftName )->getStateInBaseFrameFromEphemeris< long double, double >( t ).segment( 0, 3 ) );
            observationTimes.push_back( t );
        }
        catch( ... )
        { }
    }
    std::cout<< observationTimes.back()<<std::endl;
    std::cout<<"observations and times created"<<std::endl;
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
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
                    ( rungeKutta4, integrationStartTime + 600, 30.0 );
   // std::shared_ptr<IntegratorSettings<> >integratorSettings =
    //        std::make_shared<RungeKuttaFixedStepSizeSettings<> >( 30, CoefficientSets::rungeKutta87DormandPrince );

    std::cout<<"Integration settings created"<<std::endl;

    int numberOfIntegrationArcs = integrationArcStartTimes.size( );

    std::vector< Eigen::Matrix < long double, Eigen::Dynamic, 1 > >  systemInitialStates(numberOfIntegrationArcs);

    std::vector< std::shared_ptr< SingleArcPropagatorSettings< long double > > > arcPropagationSettingsList;
    for( unsigned int i = 0; i < numberOfIntegrationArcs; i++ )
        {
            std::cout<<"iteration started"<<std::endl;
            std::cout<<i<<std::endl;
            std::cout<<bodiesToIntegrate[ 0 ]<<std::endl;
            std::cout<<integrationArcStartTimes.at(i)<<std::endl;
            std::cout<<integrationArcEndTimes.at(i)<<std::endl;
            systemInitialStates[ i ]  =bodies.getBody( spacecraftName )->getStateInBaseFrameFromEphemeris< long double, double >( integrationArcStartTimes.at(i) ) -
            bodies.getBody( centralBody )->getStateInBaseFrameFromEphemeris< long double, double >( integrationArcStartTimes.at(i) );
            std::cout<<"system initial states created"<<std::endl;
            arcPropagationSettingsList.push_back(std::make_shared< TranslationalStatePropagatorSettings< long double > >
                                ( centralBodies, accelerationModelMap, bodiesToIntegrate,
                                  systemInitialStates.at( i ), integrationArcEndTimes.at( i ), cowell, dependentVariablesToSave, TUDAT_NAN ) );
        }

    std::shared_ptr< MultiArcPropagatorSettings< long double > > multiArcPropagatorSettings =
            validateDeprecatedMultiArcSettings< long double, double >(
            integratorSettings, std::make_shared< MultiArcPropagatorSettings< long double > >( arcPropagationSettingsList ),
            integrationArcStartTimes, false, true );


    // Select parameters to estimate
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
                getInitialMultiArcParameterSettings< long double, double  >( multiArcPropagatorSettings, bodies, integrationArcStartTimes );
    //parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( spacecraftName, radiation_pressure_coefficient ) );
    parameterNames.push_back(std::make_shared< ArcWiseDragCoefficientEstimatableParameterSettings >(spacecraftName, initial_times_list_drag ));
    parameterNames.push_back(std::make_shared< ArcWiseRadiationPressureCoefficientEstimatableParameterSettings >(spacecraftName, integrationArcStartTimes ));

    std::map< EmpiricalAccelerationComponents, std::vector< EmpiricalAccelerationFunctionalShapes > > empiricalAccelerationComponents;
    empiricalAccelerationComponents[ across_track_empirical_acceleration_component ].push_back( constant_empirical );
    //empiricalAccelerationComponents[ across_track_empirical_acceleration_component ].push_back( sine_empirical );
    empiricalAccelerationComponents[ along_track_empirical_acceleration_component ].push_back( constant_empirical );
    //empiricalAccelerationComponents[ along_track_empirical_acceleration_component ].push_back( sine_empirical );
    empiricalAccelerationComponents[ radial_empirical_acceleration_component ].push_back( constant_empirical );
    //empiricalAccelerationComponents[ radial_empirical_acceleration_component ].push_back( sine_empirical );

    parameterNames.push_back( std::make_shared< ArcWiseEmpiricalAccelerationEstimatableParameterSettings >(
            spacecraftName,"Mars", empiricalAccelerationComponents, integrationArcStartTimes ) );

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > > parametersToEstimate =
            createParametersToEstimate< long double, double >( parameterNames, bodies, multiArcPropagatorSettings );
    printEstimatableParameterEntries( parametersToEstimate );


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
                    observationModelSettingsList, multiArcPropagatorSettings );
    std::cout<<"orbit determination manager created"<<std::endl;


    // Retrieve state history
    //std::map< long double, Eigen::Matrix < long double, Eigen::Dynamic, 1 > > propagatedStateHistory;
    /*for ( Time t : observationTimes )
    {
        propagatedStateHistory[ t.getSeconds< long double >() ] =
                bodies.getBody( spacecraftName )->getStateInBaseFrameFromEphemeris< long double, Time >( t ) -
                    bodies.getBody( centralBody )->getStateInBaseFrameFromEphemeris< long double, Time >( t );
    }

    writeDataMapToTextFile( propagatedStateHistory, "stateHistoryPropagatedPreFit_" + fileTag + ".txt", saveDirectory,
                            "", 18, 18 );
    writeDataMapToTextFile( spiceStateHistory, "stateHistorySpice_" + fileTag + ".txt", saveDirectory,
                            "", 18, 18 );
*/


    std::vector< std::shared_ptr< SingleObservationSet< long double, double > > > observationSetList;
    observationSetList.push_back(
            std::make_shared< SingleObservationSet<long double, double > >(
                    position_observable, linkEnds, observations, observationTimes, observed_body ) );
    std::shared_ptr< ObservationCollection< long double, double > > observedObservationCollection =
            std::make_shared< ObservationCollection< long double, double > >( observationSetList );

    std::cout<<"observations and times created"<<std::endl;

    // set a priori to the drag coefficients
    const int DIAGONALS = numberOfIntegrationArcs*6;
    double aprioriuncertainty =1.0/(0.001*0.001);

    Eigen::Matrix< double, Eigen::Dynamic, 1 > initialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< double >( );
    int numberOfParameters = initialParameterEstimate.rows( );
    // Create a 2D vector (matrix) filled with zeros
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(numberOfParameters, numberOfParameters);
    // Fill the matrix with the values of the diagonal
    for (int i = DIAGONALS; i < numberOfParameters - 4*numberOfIntegrationArcs ; ++i) {
        matrix(i,i) = aprioriuncertainty;
    }
  // Define estimation input
    std::shared_ptr< EstimationInput< long double, double  > > estimationInput =
            std::make_shared< EstimationInput< long double, double > >(
                    observedObservationCollection,
                    Eigen::MatrixXd::Zero(0,0),//matrix,
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


    std::shared_ptr< CovarianceAnalysisInput< long double, double > > covarianceInput =
        std::make_shared< CovarianceAnalysisInput< long double, double > >(
                observedObservationCollection);

    covarianceInput->getConsiderCovariance();
    covarianceInput->defineCovarianceSettings(true, false);

    std::cout<<"covariance input created"<<std::endl;

    std::shared_ptr< CovarianceAnalysisOutput< long double, double > > covarianceOutput = orbitDeterminationManager.computeCovariance(
            covarianceInput );

    Eigen::MatrixXd propcorrelationMatrix = covarianceOutput->getCorrelationMatrix( );

    std::ofstream file11(saveDirectory + "propcorrelations_" + fileTag + ".txt");
    file11 << std::setprecision( 21 ) << propcorrelationMatrix ;
    file11.close();

    std::cout<<"covariance output created"<<std::endl;

    // Perform estimation
    std::shared_ptr< EstimationOutput< long double, double > > estimationOutput = orbitDeterminationManager.estimateParameters(
                estimationInput );

    int itr_number = 0;
    for (const auto& iterationOutput : estimationOutput->getSimulationResults()){ //simulationResultsPerIteration_) {
        auto multiArcResults = std::dynamic_pointer_cast<tudat::propagators::MultiArcSimulationResults<tudat::propagators::SingleArcVariationalSimulationResults, long double, double>>(iterationOutput);
        auto singleArcResults = multiArcResults->getSingleArcResults();
        std::map< double, Eigen::Matrix< long double, Eigen::Dynamic, 1 > > concatenatedStateHistoryPostFitDynamic;
        std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > concatenatedDependentVariablesPostFitDynamic;
        std::map< double, Eigen::Matrix< long double, 6, 1 > > concatenatedStateHistoryPostFit;
        //std::cout<<singleArcResults.size()<<std::endl;
        for ( unsigned int arcIndex = 0; arcIndex < singleArcResults.size(); ++arcIndex )
        {
            // Retrieve the state history for the current arc
            std::map< double, Eigen::Matrix< long double, Eigen::Dynamic, 1 > > arcStateHistoryPostFitDynamic = singleArcResults[arcIndex]->getDynamicsResults( )->getEquationsOfMotionNumericalSolution();

            // Retrieve the dependent variables history for the current arc
            std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > arcDependentVariablesPostFitDynamic = singleArcResults[arcIndex]->getDynamicsResults( )->getDependentVariableHistory( );
            // Append the state history to the concatenated map
            for ( auto it = arcStateHistoryPostFitDynamic.begin(); it != arcStateHistoryPostFitDynamic.end(); ++it )
            {
                concatenatedStateHistoryPostFitDynamic[ it->first ] = it->second;
                //std::cout<<concatenatedStateHistoryPostFitDynamic[ it->first ]<<std::endl;
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
        writeDataMapToTextFile( concatenatedStateHistoryPostFit, "stateHistoryPropagatedPostFit_" + fileTag + "itr" + std::to_string(itr_number) + ".txt", saveDirectory,
                                "", 18, 18 );

        // Write dependent variables to file
        writeDataMapToTextFile( concatenatedDependentVariablesPostFitDynamic, "dependentVariablesPropagatedPostFit_" + fileTag + "itr" + std::to_string(itr_number) + ".txt", saveDirectory,
                                "", 18, 18 );

        itr_number++;
    }
    // Retrieve residuals and set them in matrix
    Eigen::MatrixXd residualHistory = estimationOutput->getResidualHistoryMatrix( );
    Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic > parameterHistory = estimationOutput->getParameterHistoryMatrix( );

    Eigen::MatrixXd correlationMatrix = estimationOutput->getCorrelationMatrix( );

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

    std::ofstream file10(saveDirectory + "correlations_" + fileTag + ".txt");
    file10 << std::setprecision( 21 ) << correlationMatrix ;
    file10.close();


}