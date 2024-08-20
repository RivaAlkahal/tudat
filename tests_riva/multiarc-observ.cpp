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
    double buffer = 10.0 * epehemeridesTimeStep;
    //const double gravitationalParameter = 4.2828378e13;
    //const double planetaryRadius = 3389.5E3;

    // Select ephemeris time range (based on available data in loaded SPICE ephemeris)
    //Time initialEphemerisTime = Time( 185976000 - 1.0 * 86400.0 ); // 23 November 2005, 0h
    //Time finalEphemerisTime = Time( 186580800 + 1.0 * 86400.0 ); // 30 November 2005, 0h
    // time at 2000-01-01T00:00:00
    Time initialEphemerisTime = Time( 0.0 ) ;
    Time finalEphemerisTime = Time( 86400.0 * 6.0 ); // 5 years later
    double totalDuration = finalEphemerisTime - initialEphemerisTime;
    std::cout<<"Total duration: "<<totalDuration<<std::endl;
    double hoursperday = 10.0;

    // Define bodies to use.
    //std::vector< std::string > bodiesToCreate = {
    //        "Earth", "Sun", "Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Phobos", "Deimos",
    //        "Io", "Ganymede", "Callisto", "Europa", "Titan" };
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

    // Create radiation pressure settings
    double referenceAreaRadiation = 4.0;
    double radiationPressureCoefficient = 1.2;
    std::vector< std::string > occultingBodies = { "Mars" };
    bodySettings.at( spacecraftName )->radiationPressureSettings[ "Sun" ] =
            cannonBallRadiationPressureSettings( "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    std::shared_ptr<AerodynamicCoefficientSettings> aerodynamicCoefficientSettings =
            std::make_shared<ConstantAerodynamicCoefficientSettings>(
                    2.0, 4.0, Eigen::Vector3d::Zero( ), Eigen::Vector3d::UnitX( ), Eigen::Vector3d::Zero( ),
                    negative_aerodynamic_frame_coefficients, negative_aerodynamic_frame_coefficients );
    bodySettings.at( spacecraftName )->aerodynamicCoefficientSettings = aerodynamicCoefficientSettings;
    bodySettings.at( spacecraftName )->ephemerisSettings->resetMakeMultiArcEphemeris( true );

    std::vector< std::shared_ptr< GravityFieldVariationSettings > > gravityFieldVariations;

    std::map<int, Eigen::MatrixXd> cosineAmplitudes;
    cosineAmplitudes[ 1 ] = Eigen::Matrix< double, 3, 5 >::Zero( );
    //cosineAmplitudes[ 2 ] = Eigen::Matrix< double, 3, 5 >::Zero( );
    std::map<int, Eigen::MatrixXd> sineAmplitudes;
    sineAmplitudes[ 1 ] = Eigen::Matrix< double, 3, 5 >::Zero( );
    std::shared_ptr< GravityFieldVariationSettings > polynomialGravityFieldVariations =
            std::make_shared< PolynomialGravityFieldVariationsSettings >(
                    cosineAmplitudes, sineAmplitudes, -1.0E4, 2, 0 );
    gravityFieldVariations.push_back( polynomialGravityFieldVariations );

    std::vector< std::shared_ptr< GravityFieldVariationSettings > > gravityFieldVariationSettings =
            gravityFieldVariations;
    bodySettings.at( "Mars" )->gravityFieldVariationSettings = gravityFieldVariations;

    SystemOfBodies bodies = createSystemOfBodies< long double, Time >( bodySettings );

    //Time initialPropagationTime = initialEphemerisTime + 600.0;
    //Time finalPropagationTime = finalEphemerisTime - 600.0;



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
    // Create integrator settings
    std::vector< double > integrationArcStartTimes;
    std::vector< double > integrationArcEndTimes;
    std::vector< double > integrationArcLimits;

    double integrationStartTime = initialEphemerisTime + 120.0; //1.0E2;
    double integrationEndTime = finalEphemerisTime - 120.0 ; //1.0E2;
    double totalDurationIntegration = integrationEndTime - integrationStartTime;
    double arcDuration = 86400.0;//2.0E4;
    double arcOverlap = 60.0;//2.0E2;

    double currentStartTime = integrationStartTime, currentEndTime = integrationStartTime + arcDuration;

    do
    {
        integrationArcLimits.push_back( currentStartTime );
        integrationArcEndTimes.push_back( currentEndTime );
        std::cout<<"integration arc start times: "<<currentStartTime<<std::endl;
        std::cout<<"integration arc end times: "<<currentEndTime<<std::endl;
        integrationArcStartTimes.push_back( currentStartTime );
        currentStartTime = currentEndTime - arcOverlap;
        std::cout<<"current start time: "<<currentStartTime<<std::endl;
        currentEndTime = currentStartTime + arcDuration;
        std::cout<<"current end time: "<<currentEndTime<<std::endl;
    }
    while( currentEndTime < integrationEndTime );
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

    // Create integrator settings
    /* std::pair< double, double > integrationMinMaxStep = std::make_pair( 1e-16, 1e16 );

     double initialStep = 5.0;
     if ( integrationMinMaxStep.first == integrationMinMaxStep.second )
     {
         initialStep = integrationMinMaxStep.first;
     }
     std::shared_ptr< IntegratorSettings< Time > > integratorSettings = rungeKuttaVariableStepSettingsScalarTolerances< Time >(
             initialStep, rungeKutta87DormandPrince, integrationMinMaxStep.first,
             integrationMinMaxStep.second, 1e-10, 1e-10 );*/

    /*std::shared_ptr< IntegratorSettings< Time > > integratorSettings =
            rungeKutta4Settings< Time >( 60.0 );*/
    // Set initial state from ephemerides
    /* Eigen::Matrix< long double, 6, 1 > spacecraftInitialState =
             bodies.getBody( spacecraftName )->getStateInBaseFrameFromEphemeris< long double, Time >( initialPropagationTime ) -
             bodies.getBody( centralBody )->getStateInBaseFrameFromEphemeris< long double, Time >( initialPropagationTime );
 */


    int numberOfIntegrationArcs = integrationArcStartTimes.size( );
    for( unsigned int i = 0; i < numberOfIntegrationArcs; i++ ) {
        std::cout << "integration arc start times"<<integrationArcStartTimes.at(i) << std::endl;
        std::cout << "integration arc end times"<<integrationArcEndTimes.at(i) << std::endl;
    }
    std::cout<<"number of integration arcs: "<<numberOfIntegrationArcs<<std::endl;
    std::vector< Eigen::VectorXd > systemInitialStates(numberOfIntegrationArcs, Eigen::VectorXd(6));

    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > arcPropagationSettingsList;
    for( unsigned int i = 0; i < numberOfIntegrationArcs; i++ )
    {
        std::cout<<"iteration started"<<std::endl;
        std::cout<<i<<std::endl;
        std::cout<<bodiesToIntegrate[ 0 ]<<std::endl;
        std::cout<<integrationArcStartTimes.at(i)<<std::endl;
        systemInitialStates[ i ]  = spice_interface::getBodyCartesianStateAtEpoch(
                bodiesToIntegrate[ 0 ], "Mars", "MARSIAU", "NONE", integrationArcStartTimes.at(i));
        std::cout<<"system initial states created"<<std::endl;
        arcPropagationSettingsList.push_back(
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                        ( centralBodies, accelerationModelMap, bodiesToIntegrate,
                          systemInitialStates.at( i ), integrationArcEndTimes.at( i ) ) );
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

    MultiArcDynamicsSimulator< > dynamicsSimulator(
            bodies, multiArcPropagatorSettings );

    std::map< long double, Eigen::Matrix < long double, Eigen::Dynamic, 1 > > propagatedStateHistory;
    for ( Time t : observationTimes )
    {
        propagatedStateHistory[ t.getSeconds< long double >() ] =
                bodies.getBody( spacecraftName )->getStateInBaseFrameFromEphemeris< long double, Time >( t ) -
                bodies.getBody( centralBody )->getStateInBaseFrameFromEphemeris< long double, Time >( t );
    }
    writeDataMapToTextFile( propagatedStateHistory, "stateHistoryPropagation_" + fileTag + ".txt", saveDirectory,
                            "", 18, 18 );


    //initial state?
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
            getInitialMultiArcParameterSettings< double, double  >( multiArcPropagatorSettings, bodies, integrationArcStartTimes );
    //polynomial gravity field
    std::map<int, std::vector<std::pair<int, int> > > cosineBlockIndicesPerPower;
    cosineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 2, 0 ) );
    std::map<int, std::vector<std::pair<int, int> > > sineBlockIndicesPerPower;
    sineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 3, 2 ) );
    parameterNames.push_back( std::make_shared< PolynomialGravityFieldVariationEstimatableParameterSettings >(
            "Mars", cosineBlockIndicesPerPower, sineBlockIndicesPerPower ) );

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
        final_times_list.push_back(integrationArcStartTimes[0] + buffer + day * 24 * 3600 + hoursperday * 3600);
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
    LinkEnds linkEnds;
    // Create list of link ends where the ground station is the transmitter and the spacecraft is the receiver
    std::vector< LinkEnds > stationTransmitterLinkEnds;
    //linkEnds[ receiver ] = spacecraftName;
    std::vector< std::string > GroundStations = {  "DSS-26" , "DSS-42", "DSS-61"};
    for ( std::string groundStation : GroundStations )
    {
        linkEnds[ transmitter ] = LinkEndId( "Earth", groundStation );
        linkEnds[ receiver ] = spacecraftName;
        stationTransmitterLinkEnds.push_back( linkEnds );
    }
    // Define (arbitrarily) link ends to be used for the chosen observable type (e.g., one-way range)
    //std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
    ObservableType selectedObservable = one_way_range;  // Define the observable type (replace with desired type)
    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
    // Populate the map with all the transmitter link ends
    for( unsigned int i = 0; i < stationTransmitterLinkEnds.size( ); i++ )
    {
        //linkEndsPerObservable[selectedObservable].push_back( stationTransmitterLinkEnds[ i ] );
        observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
                one_way_range, stationTransmitterLinkEnds.at( i ) ) );
    }

    std::cout<<"link ends created"<<std::endl;

    // Create observation settings
    //std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
    // observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
    //        one_way_range, linkEnds ) );
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

    //std::vector< std::shared_ptr< ObservationViabilityCalculator > > viabilityCalculators;
   // for (unsigned int i = 0; i<stationTransmitterLinkEnds.size(); i++) {

    //    viabilityCalculators.push_back(
   //             std::vector< std::shared_ptr< ObservationViabilityCalculator > >(createObservationViabilityCalculators( bodies, stationTransmitterLinkEnds[i], one_way_range, observationViabilitySettings )));

   // }
    std::vector< std::shared_ptr< ObservationViabilityCalculator > > viabilityCalculators;
    for (unsigned int i = 0; i < stationTransmitterLinkEnds.size(); i++) {
        std::vector< std::shared_ptr< ObservationViabilityCalculator > > calculators =
                createObservationViabilityCalculators(bodies, stationTransmitterLinkEnds[i], one_way_range, observationViabilitySettings);
        viabilityCalculators.insert(viabilityCalculators.end(), calculators.begin(), calculators.end());
    }
    // PerObservableObservationViabilityCalculatorList viabilityCalculators = createObservationViabilityCalculators(
    //         bodies, stationTransmitterLinkEnds, observationViabilitySettings );

    std::cout<<"Observation viability calculator created"<<std::endl;

    // Create noise
    // Create noise functions per observable
    double rangeNoise = 0.1;
    std::map< ObservableType, std::function< double( const double ) > > noiseFunctions;
    noiseFunctions[ one_way_range ] =
            std::bind( &utilities::evaluateFunctionWithoutInputArgumentDependency< double, const double >,
                       createBoostContinuousRandomVariableGeneratorFunction(
                               tudat::statistics::normal_boost_distribution, { 0.0, rangeNoise }, 0.0 ), std::placeholders::_1 );

    std::cout<<"noise functions created"<<std::endl;
    // Wrapper function to convert double to Eigen::VectorXd
    std::function<Eigen::VectorXd(const double)> noiseFunctionWrapper =
            [noiseFunctions](const double time) -> Eigen::VectorXd
            {
                double noiseValue = noiseFunctions.at(one_way_range)(time);
                Eigen::VectorXd noiseVector(1);
                noiseVector(0) = noiseValue;
                return noiseVector;
            };
    std::cout<<"noise functions converted to Eigen::VectorXd"<<std::endl;

    Eigen::Matrix< double, Eigen::Dynamic, 1 > initialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< double >( );

    // Create observation model settings
    //std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationModelSettingsList;
    //observationModelSettingsList.push_back( positionObservableSettings( linkEnds ) );


// Observation simulation settings
    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput;
    for (unsigned int r = 0; r < stationTransmitterLinkEnds.size(); r++) {
        measurementSimulationInput.push_back(
                std::make_shared< TabulatedObservationSimulationSettings< double > >(
                        one_way_range, stationTransmitterLinkEnds[r], observationTimesList, receiver, observationViabilitySettings, noiseFunctionWrapper));
    }
    //maybe add the viability settings here
    //check out the measurment simulation settings


    std::shared_ptr< ObservationCollection< double, double > > observationsAndTimes = simulateObservations< double, double >(
            measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies  );
    //std::shared_ptr< ObservationCollection< double, double > > observationsAndTimes = simulateObservationsWithNoise< double, double >(
    //        measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), noiseFunctions, viabilityCalculators );
    //PodInputDataType observationsAndTimes = simulateObservationsWithNoise< double, double >(
    //        measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), noiseFunctions,
    //        viabilityCalculators );
    std::cout<<"observations and times created"<<std::endl;

    Eigen::Matrix< double, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;
    int numberOfParameters = initialParameterEstimate.rows( );
    // Perturb initial states
    for( unsigned int i = 0; i < integrationArcStartTimes.size( ); i++ )
    {
        initialParameterEstimate[ 0 + 6 * i ] += 1.0E0;
        initialParameterEstimate[ 1 + 6 * i ] += 1.0E0;
        initialParameterEstimate[ 2 + 6 * i ] += 1.0E0;
        initialParameterEstimate[ 3 + 6 * i ] += 1.0E-5;
        initialParameterEstimate[ 4 + 6 * i ] += 1.0E-5;
        initialParameterEstimate[ 5 + 6 * i ] += 1.0E-5;
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

    // Define estimation input
    std::shared_ptr< EstimationInput< double, double  > > estimationInput =
            std::make_shared< EstimationInput< double, double > >(
                    observationsAndTimes,
                    Eigen::MatrixXd::Zero( 0, 0 ),
                    std::make_shared< EstimationConvergenceChecker >( 5 ) );
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
    std::cout<<"estimation input created"<<std::endl;

    //std::shared_ptr< EstimationInput< double, double > > estimationInput =
    //        std::make_shared< EstimationInput< double, double > >(
    //                observationsAndTimes );
    std::shared_ptr< CovarianceAnalysisInput< double, double > > covarianceInput =
            std::make_shared< CovarianceAnalysisInput< double, double > >(
                    observationsAndTimes );
    std::cout<<"covariance input created"<<std::endl;

    // Perform estimation
    std::shared_ptr< EstimationOutput< double, double > > estimationOutput = orbitDeterminationManager.estimateParameters(
            estimationInput );
    std::cout<<"estimation performed"<<std::endl;

    Eigen::Matrix< double, Eigen::Dynamic, 1 > finalParameters = parametersToEstimate->template getFullParameterValues< double >( );

    // save true errors
    Eigen::Matrix< double, Eigen::Dynamic, 1 > TrueError = ( estimationOutput->parameterEstimate_ - truthParameters ).transpose( );
    //TrueError = ( estimationOutput->parameterEstimate_ - truthParameters ).transpose( );
    std::cout<<"True error: "<<( estimationOutput->parameterEstimate_ - truthParameters ).transpose( )<<std::endl;
    // save formal errors
    Eigen::Matrix<double, Eigen::Dynamic, 1> FormalError = estimationOutput->getFormalErrorVector( );
    std::cout<<"Formal error: "<<estimationOutput->getFormalErrorVector( ).transpose( )<<std::endl;
    std::cout<<"Error ratio: "<<( ( 1.0E-3 * estimationOutput->getFormalErrorVector( ).segment( 0, numberOfParameters ) ).cwiseQuotient(
            estimationOutput->parameterEstimate_ - truthParameters ) ).transpose( )<<std::endl;

    // Retrieve estimated state history
    // Retrieve the best iteration simulation results
    parametersToEstimate->resetParameterValues(estimationOutput->parameterHistory_.at(estimationOutput->bestIteration_));

    std::vector< std::shared_ptr< propagators::SimulationResults< long double, Time > > > bestIterationSimulationResultsPerArc;
    std::map< Time, Eigen::Matrix< long double, Eigen::Dynamic, 1 > > concatenatedStateHistoryPostFitDynamic;
    std::map< Time, Eigen::Matrix< long double, 6, 1 > > concatenatedStateHistoryPostFit;

    for ( unsigned int arcIndex = 0; arcIndex < bestIterationSimulationResultsPerArc.size(); ++arcIndex )
    {
        // Retrieve the best iteration simulation results for the current arc
        std::shared_ptr< propagators::SimulationResults< long double, Time > > bestIterationSimulationResults =
                bestIterationSimulationResultsPerArc[arcIndex];
        bestIterationSimulationResultsPerArc.push_back( bestIterationSimulationResults );

        // Retrieve the state history for the current arc
        std::map< Time, Eigen::Matrix< long double, Eigen::Dynamic, 1 > > arcStateHistoryPostFitDynamic =
                std::dynamic_pointer_cast< SingleArcVariationalSimulationResults< long double, Time > >(
                        bestIterationSimulationResults )->getDynamicsResults( )->getEquationsOfMotionNumericalSolution();

        // Append the state history to the concatenated map
        for ( auto it = arcStateHistoryPostFitDynamic.begin(); it != arcStateHistoryPostFitDynamic.end(); ++it )
        {
            concatenatedStateHistoryPostFitDynamic[ it->first ] = it->second;
        }
    }

// Convert the dynamic state history to fixed size
    for ( auto it = concatenatedStateHistoryPostFitDynamic.begin(); it != concatenatedStateHistoryPostFitDynamic.end(); ++it )
    {
        concatenatedStateHistoryPostFit[ it->first ] = it->second;
    }

    //parametersToEstimate->resetParameterValues( estimationOutput->parameterHistory_.at( estimationOutput->bestIteration_ ) );

    std::shared_ptr< CovarianceAnalysisOutput< double, double > > covarianceOutput = orbitDeterminationManager.computeCovariance(
            covarianceInput );

    //compareEstimationAndCovarianceResults( estimationOutput, covarianceOutput );


    //return ( finalParameters - truthParameters ).template cast< double >( );




    // Create link ends
    /*LinkEnds linkEnds;
    linkEnds[ observed_body ] = spacecraftName;

    // Create observation model settings
    std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationModelSettingsList;
    observationModelSettingsList.push_back( positionObservableSettings( linkEnds ) );

    // Select parameters to estimate
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< long double > >(
            spacecraftName, spacecraftInitialState, centralBody, baseFrameOrientation ) );
    //parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( spacecraftName, radiation_pressure_coefficient ) );
    //parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( spacecraftName, constant_drag_coefficient ) );

    std::map<int, std::vector<std::pair<int, int> > > cosineBlockIndicesPerPower;
    cosineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 2, 0 ) );
    std::map<int, std::vector<std::pair<int, int> > > sineBlockIndicesPerPower;
    sineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 3, 2 ) );
    parameterNames.push_back( std::make_shared< PolynomialGravityFieldVariationEstimatableParameterSettings >(
            "Mars", cosineBlockIndicesPerPower, sineBlockIndicesPerPower ) );

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > > parametersToEstimate =
            createParametersToEstimate< long double, Time >( parameterNames, bodies );

    // Create orbit determination object.
    OrbitDeterminationManager< long double, Time > orbitDeterminationManager =
            OrbitDeterminationManager< long double, Time >(
                    bodies, parametersToEstimate,
                    observationModelSettingsList, propagatorSettings, true );

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

    std::vector< std::shared_ptr< SingleObservationSet< long double, Time > > > observationSetList;
    observationSetList.push_back(
            std::make_shared< SingleObservationSet< long double, Time > >(
                    position_observable, linkEnds, observations, observationTimes, observed_body ) );
    std::shared_ptr< ObservationCollection< long double, Time > > observedObservationCollection =
            std::make_shared< ObservationCollection< long double, Time > >( observationSetList );

    // Define estimation input
    std::shared_ptr< EstimationInput< long double, Time  > > estimationInput =
            std::make_shared< EstimationInput< long double, Time > >(
                    observedObservationCollection,
                    Eigen::MatrixXd::Zero( 0, 0 ),
                    std::make_shared< EstimationConvergenceChecker >( 5 ) );
    estimationInput->saveStateHistoryForEachIteration_ = true;

    // Perform estimation
    std::shared_ptr< EstimationOutput< long double, Time > > estimationOutput = orbitDeterminationManager.estimateParameters(
            estimationInput );

    // Retrieve post-fit state history
    std::shared_ptr< propagators::SimulationResults< long double, Time > > postFitSimulationResults =
            estimationOutput->getBestIterationSimulationResults( );
    std::map < Time, Eigen::Matrix < long double, Eigen::Dynamic, 1 > > propagatedStateHistoryPostFitDynamic =
            std::dynamic_pointer_cast< SingleArcVariationalSimulationResults< long double, Time > >(
                    postFitSimulationResults )->getDynamicsResults( )->getEquationsOfMotionNumericalSolution( );
    std::map < Time, Eigen::Matrix < long double, 6, 1 > > propagatedStateHistoryPostFit;
    for ( auto it = propagatedStateHistoryPostFitDynamic.begin( ); it != propagatedStateHistoryPostFitDynamic.end( ); ++it )
    {
        propagatedStateHistoryPostFit[ it->first ] = it->second;
    }
    std::shared_ptr< interpolators::OneDimensionalInterpolator< Time, Eigen::Matrix< long double, 6, 1 > > > postFitStateInterpolator =
            propagators::createStateInterpolator< Time, long double >( propagatedStateHistoryPostFit );
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
    */

}

