//
// Created by ralkahal on 3-6-25.
//
//
// Created by Riva Alkahal on 07/11/2024.
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

// Function to compute the cross product of position and velocity and store it in a map
std::map<double, Eigen::VectorXd> computeCrossProduct(const std::map<double, Eigen::VectorXd>& stateHistory) {
        std::map<double, Eigen::VectorXd> crossProductMap;

        for (const auto& [time, stateVector] : stateHistory) {
                // Ensure the state vector has exactly 6 elements (3 for position, 3 for velocity)
                if (stateVector.size() == 6) {
                        // Extract position and velocity vectors
                        Eigen::Vector3d position = stateVector.head(3);
                        Eigen::Vector3d velocity = stateVector.tail(3);

                        // Compute the cross product
                        Eigen::Vector3d crossProduct = position.cross(velocity);

                        // Store the result in the map with the time as the key
                        crossProductMap[time] = crossProduct;
                } else {
                        std::cerr << "Warning: State vector at time " << time << " does not have exactly 6 elements." << std::endl;
                        crossProductMap[time] = Eigen::Vector3d::Zero();  // Placeholder if state vector is not size 6
                }
        }

        return crossProductMap;
}
// Example function to compute dot products for a map of vectors
std::map<double, double> computeDotProductMap(const std::map<double, Eigen::VectorXd>& vectorMap1, const std::map<double, Eigen::VectorXd>& vectorMap2) {
        std::map<double, double> dotProductMap;

        for (const auto& [time, vector1] : vectorMap1) {
                // Ensure the time key exists in both maps
                if (vectorMap2.find(time) != vectorMap2.end()) {
                        const Eigen::VectorXd& vector2 = vectorMap2.at(time);
                        // Compute the dot product
                        dotProductMap[time] = vector1.dot(vector2);

                } else {
                        std::cerr << "Warning: Time key " << time << " not found in both maps." << std::endl;
                }
        }

        return dotProductMap;
}
std::map<double, double> computeNorms(const std::map<double, Eigen::VectorXd>& relativePos) {
        std::map<double, double> norms;

        for (const auto& [key, vector] : relativePos) {
                // Calculate the norm of the vector (last 3 columns)
                double norm = vector.norm();
                norms[key] = norm;
        }

        return norms;
}
void arcLengthRuns( double hoursperday, double initialTime, double finalTime, int arcLength, int iterationNumber, double perturbPos, double perturbVel, int number_of_arcs, int itotalDuration, int startTime, int intperturbPos, int intperturbVel , int ihoursperday, bool performEst)
{
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
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_map1.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_map2.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_map3.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_map4.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_map5.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_map6.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_map7.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_map8.bsp" );

    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_map4_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_map5_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_map6_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_map7_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_map8_ipng_mgs95j.bsp" );

    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext1.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext2.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext3.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext4.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext5.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext6.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext7.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext8.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext9.bsp" );

    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext5_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext6_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext7_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext8_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext9_ipng_mgs95j.bsp" );


//    std::string saveDirectory = "/Users/ralkahal/OneDrive - Delft University of Technology/new-tudat-tests/propagateDynamics/";
    std::string saveDirectory = "/home/ralkahal/nnew-tudat-tests/arcLengths_paperTests/";


    //std::string fileTag = "arcLengths_160darc_startat320_0per_10hobs_6itr";
    std::string fileTag = "RK78_60S_arcLengths_" +  std::to_string(arcLength) + std::to_string(itotalDuration) +  "darc_startat" + std::to_string(startTime)
                          + "_" + std::to_string(intperturbPos) + "perpos_" + std::to_string(intperturbVel) + "_pervel_" + std::to_string(ihoursperday) + "hobs_" + std::to_string(iterationNumber) + "itr" ;

    // set input options
    double epehemeridesTimeStep = 60.0;
    bool useInterpolatedEphemerides = true;
    double observationsSamplingTime = 60.0;
    double buffer = 30.0 * epehemeridesTimeStep;
    double arcDuration = arcLength*86400.0;//2.0E4;
    std::string dragEst = "per-halfday";//"per-rev";
    double ndays = 0.5;
    double hoursperdaydrag = 2.0;
    //double hoursperday = 10.0;
    //int iterationNumber = 6;
    //const double gravitationalParameter = 4.2828378e13;
    //const double planetaryRadius = 3389.5E3;

    double oneWayDopplerNoise = 0.0001;
    double twoWayDopplerNoise = 0.0001;
    double rangeNoise = 1.0;
    // Select ephemeris time range (based on available data in loaded SPICE ephemeris)
    //Time initialEphemerisTime = Time( 185976000 - 1.0 * 86400.0 ); // 23 November 2005, 0h
    //Time finalEphemerisTime = Time( 186580800 + 1.0 * 86400.0 ); // 30 November 2005, 0h
    // time at 2000-01-01T00:00:00
    Time initialEphemerisTime = Time(initialTime); //Time( 86400.0 * 350 ) ;
    Time finalEphemerisTime = Time(finalTime);//Time( 86400.0 * 350 + 86400.0 * 60.0 ); // 5 years later
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

    std::string filename = "/home/ralkahal/new-tudat-tests/dtm-mars";
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
    bodySettings.at( spacecraftName )->constantMass = 1030.5;
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
    std::vector<Eigen::MatrixXd> cosineShAmplitudesCosineTime;
    std::vector<Eigen::MatrixXd> cosineShAmplitudesSineTime;
    std::vector<Eigen::MatrixXd> sineShAmplitudesCosineTime;
    std::vector<Eigen::MatrixXd> sineShAmplitudesSineTime;
    std::vector<double> frequencies;
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
    cosineAmplitudes[ 1 ] = Eigen::Matrix< double, 4, 5 >::Zero( );
        //nVec Root 800 km depth
        cosineAmplitudes[ 1 ]( 0, 0 ) += -7.00583559071078e-13/(365*24*3600);
        cosineAmplitudes[1](0,1) +=-5.44909982178362e-14/(365*24*3600);
        cosineAmplitudes[1](0,2) += -8.31134553095663e-13/(365*24*3600);
        cosineAmplitudes[1](1,0) += 1.15640729781333e-13/(365*24*3600);
        cosineAmplitudes[1](1,1) += -2.61531330180095e-13/(365*24*3600);
        cosineAmplitudes[1](1,2) += 1.02212332919433e-13/(365*24*3600);
        cosineAmplitudes[1](1,3) += -8.00674595992919e-13/(365*24*3600);
        cosineAmplitudes[1](2,0) += 4.32941757959663e-13/(365*24*3600);
        cosineAmplitudes[1](2,1) += 5.98812345051549e-14/(365*24*3600);
        cosineAmplitudes[1](2,2) += 4.4199871466477e-13/(365*24*3600);
        cosineAmplitudes[1](2,3) += 1.25159028954031e-13/(365*24*3600);
        cosineAmplitudes[1](2,4) += -5.3726758654485e-14/(365*24*3600);

    std::map<int, Eigen::MatrixXd> sineAmplitudes;
        std::cout<<"creating settings for poly grav"<<std::endl;
    std::shared_ptr< GravityFieldVariationSettings > polynomialGravityFieldVariations =
            std::make_shared< PolynomialGravityFieldVariationsSettings >(
                    cosineAmplitudes, sineAmplitudes, 0.0, 2, 0 );

    gravityFieldVariations.push_back(polynomialGravityFieldVariations);

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
    accelerationsOfVehicle["Mars"].push_back(sphericalHarmonicAcceleration(95, 95));
    accelerationsOfVehicle["Mars"].push_back(relativisticAccelerationCorrection());
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

        dependentVariablesToSave.push_back(std::make_shared< SingleDependentVariableSaveSettings >(
                relative_position_dependent_variable, "Earth", "Mars" ));
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
    for (double time = integrationStartTime + 10000 ; time < integrationEndTime + 10000 ; time += step_size) {
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
        // observations times
        int days = static_cast<int>(totalDuration / physical_constants::JULIAN_DAY );
        std::cout<<"days: "<<days<<std::endl;
        std::cout<<integrationArcStartTimes.size( )<<std::endl;
        std::vector<double> initial_times_list = { integrationArcStartTimes[0] + buffer};
        std::vector<double> filtered_initial_times_list = { integrationArcStartTimes[0] + buffer};
        if (hoursperday > 20.0) { days = days-1; }

        for (int day = 1; day < days  ; ++day) {
                initial_times_list.push_back(integrationArcStartTimes[0] + buffer + day * 24 * 3600);
        }

        std::vector<double> final_times_list = { initial_times_list[0] + hoursperday * 3600 };
        for (int day = 1; day < days ; ++day) {
                double finalTime = initial_times_list[day] + hoursperday * 3600;
                if (finalTime > integrationArcEndTimes.back()) {
                        break;
                }
                final_times_list.push_back(finalTime);
                filtered_initial_times_list.push_back(initial_times_list[day]);
        }
        std::cout<< " final times list size: " <<final_times_list.size()<<std::endl;
        std::vector<double> observationTimesList;
        for (int i = 0; i < final_times_list.size(); ++i) {
                double start = filtered_initial_times_list[i];
                std::cout<<"start time: "<<start<<std::endl;
                double end = final_times_list[i];
                std::cout<<"end time: "<<end<<std::endl;
                for (double time = start; time < end ; time += 10) {
                        observationTimesList.push_back(time);
                }
        }
        std::cout<< "Observation times created" << std::endl;

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

        // define integrator settings
    std::shared_ptr<IntegratorSettings<> >integratorSettings =
       std::make_shared<RungeKuttaFixedStepSizeSettings<> >( 60, CoefficientSets::rungeKutta87DormandPrince );
        std::cout<<"Integration settings created"<<std::endl;
    int numberOfIntegrationArcs = integrationArcStartTimes.size( );

        std::cout<<"number of integration arcs: "<<numberOfIntegrationArcs<<std::endl;
        std::vector< Eigen::VectorXd > systemInitialStates(numberOfIntegrationArcs, Eigen::VectorXd(6));

        // create multi arc propagation settings
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
                                  systemInitialStates.at(i), integrationArcEndTimes.at( i ), cowell, dependentVariablesToSave, TUDAT_NAN ) );
        }

        std::cout<<"single arc propagation done"<<std::endl;
        std::shared_ptr< MultiArcPropagatorSettings< double > > multiArcPropagatorSettings =
                validateDeprecatedMultiArcSettings< double, double >(
                        integratorSettings, std::make_shared< MultiArcPropagatorSettings< double > >( arcPropagationSettingsList ),
                        integrationArcStartTimes, false, true );

        MultiArcDynamicsSimulator< > dynamicsSimulator(
            bodies, multiArcPropagatorSettings );

        std::map< double, Eigen::Matrix< double,Eigen::Dynamic,1> > integrationResult;
        std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > dependentVariableResult;
        for ( unsigned int arcIndex = 0; arcIndex < numberOfIntegrationArcs; ++arcIndex ) {
                auto singleArcResult = dynamicsSimulator.getMultiArcPropagationResults()->getSingleArcResults( ).at(
                                    arcIndex )->getEquationsOfMotionNumericalSolution( );
                integrationResult.insert(singleArcResult.begin(), singleArcResult.end());
                auto singleArcDepVars = dynamicsSimulator.getMultiArcPropagationResults()->getSingleArcResults( ).at(
                                    arcIndex )->getDependentVariableHistory( );
                dependentVariableResult.insert(singleArcDepVars.begin(), singleArcDepVars.end());
                // std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
                // std::map< double, Eigen::VectorXd > dependentVariableResult = dynamicsSimulator.getDependentVariableHistory( );
        }
        writeDataMapToTextFile( integrationResult, "stateHistoryPropagation_" + fileTag + ".txt", saveDirectory,
                                "", 18, 18 );
        writeDataMapToTextFile( dependentVariableResult, "dependentVariablesPropagation_" + fileTag + ".txt", saveDirectory,
                                "", 18, 18 );



        // save the last column of dependentVariableResult
        std::map<double, Eigen::VectorXd> relativePosEarthtoMars;
        for (const auto& [key, vector] : dependentVariableResult) {
                // Check if the vector has at least 3 elements
                if (vector.size() >= 3) {
                        // Extract the last three elements and store them in the new map
                        relativePosEarthtoMars[key] = vector.tail(3);
                } else {
                        std::cerr << "Warning: Vector at key " << key << " has less than 3 elements." << std::endl;
                        // Optional: you can choose to add an empty or zero vector if needed
                        relativePosEarthtoMars[key] = Eigen::VectorXd::Zero(3);  // Fills with zeros as placeholder
                }
        }

        // Compute the norms of the relative positions
        std::map<double, double> normsEarthtoMars = computeNorms(relativePosEarthtoMars);
        // Compute the cross product of the position and velocity vectors
        std::map<double, Eigen::VectorXd> relativePosMarstoSat = computeCrossProduct(integrationResult);
        // Compute the norm of the cross product
        std::map<double, double> normsMarstoSat = computeNorms(relativePosMarstoSat);
        // Compute the dot product of the norms
        std::map<double, double> dotProducts = computeDotProductMap(relativePosEarthtoMars, relativePosMarstoSat);
        // Compute the angle between the two vectors
        std::map<double, double> angles;
        for (const auto& [key, dotProduct] : dotProducts) {
                // Calculate the angle in radians
                double angle = std::acos(dotProduct / (normsEarthtoMars[key] * normsMarstoSat[key]));
                // Convert to degrees and store in the map
                angles[key] = angle * 180 / mathematical_constants::PI;
        }

        // save angle to file
        writeDataMapToTextFile( angles, "BetaAngles_" + fileTag + ".txt", saveDirectory,
                                "", 18, 18 );


        std::map< double, Eigen::VectorXd > stateVectorsAtStartTimes;

        // Create parameters to estimate
        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
          getInitialMultiArcParameterSettings< double, double  >( multiArcPropagatorSettings, bodies, integrationArcStartTimes );
        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
                    createParametersToEstimate< double, double >( parameterNames, bodies );

        std::cout<<"parameters to estimate created"<<std::endl;


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
        linkEndsPerObservable[ two_way_doppler ].push_back( stationTransmitterLinkEnds[ 0 ] );
        linkEndsPerObservable[ two_way_doppler ].push_back( stationTransmitterLinkEnds[ 1 ] );
        linkEndsPerObservable[ two_way_doppler ].push_back( stationTransmitterLinkEnds[ 2 ] );
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
        // Define the noise functions map with the required type
    std::map< ObservableType, std::function< Eigen::VectorXd( const double ) > > noiseFunctions;

// Create the noise functions that return Eigen::VectorXd
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
        Eigen::Matrix< double, Eigen::Dynamic, 1 > initialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< double >( );

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
        std::cout<<"observations and times created"<<std::endl;

        Eigen::Matrix< double, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;
        int numberOfParameters = initialParameterEstimate.rows( );
        std::cout<<"number of parameters: "<<numberOfParameters<<std::endl;
        printEstimatableParameterEntries( parametersToEstimate );
        // set a priori to the drag coefficients
    //const int DIAGONALS = numberOfIntegrationArcs*6;
    //double aprioriuncertainty =1.0/(0.1*0.1);

    // Create a 2D vector (matrix) filled with zeros
    //Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(numberOfParameters, numberOfParameters);
    // Fill the matrix with the values of the diagonal
    //for (int i = DIAGONALS; i < numberOfParameters; ++i) {
    //    matrix(i,i) = aprioriuncertainty;
    //}
    //print matrix
    //std::cout<<matrix<<std::endl;
// retrieve simulate observations
        std::ofstream outputfile1(saveDirectory + "concatenatedlinkedId-stationNames-Observations" + fileTag + ".txt");
        outputfile1 << std::setprecision(17);
        std::map<tudat::observation_models::ObservableType, std::map<int, std::vector<std::shared_ptr<SingleObservationSet < double, double>>>>> sortedObservationSets = observationsAndTimes->getSortedObservationSets();
        std::vector< double > concatenatedObservationTimes = observationsAndTimes->getConcatenatedTimeVector();
        Eigen::Matrix< double, Eigen::Dynamic, 1 > observationVector = observationsAndTimes->getObservationVector();
        std::vector< LinkEnds >  concatenatedLinkEndIdNames = observationsAndTimes->getConcatenatedLinkEndIdNames();
        std::cout << "Observation Times size:" << concatenatedObservationTimes.size() << std::endl;
        std::cout << "Observation Vector size:" << observationVector.size() << std::endl;
        std::cout<<"concatenated link end ids: " << std::endl;
        if (concatenatedObservationTimes.size() != observationVector.size() || concatenatedObservationTimes.size() != concatenatedLinkEndIdNames.size()) {
                throw std::runtime_error("Error: Observation times, observations vector and link end IDs have different sizes.");
        }

        for (size_t i = 0; i < concatenatedObservationTimes.size(); ++i) {
                double time = concatenatedObservationTimes[i];
                double observation = observationVector[i];
                LinkEnds station = concatenatedLinkEndIdNames[i];

                std::string stationName;
                int linkEndID = -1;
                for (const auto &linkEnd : station) {
                        if (linkEnd.first == 5) {
                                linkEndID = static_cast<int>(linkEnd.first);
                                stationName = linkEnd.second.stationName_;
                                break;
                        }

                }
                if (linkEndID == -1) {
                        for (const auto &linkEnd : station) {
                                if (linkEnd.first == 0) {
                                        linkEndID = static_cast<int>(linkEnd.first);
                                        stationName = linkEnd.second.stationName_;
                                        break;
                                }
                        }
                }
                if (linkEndID == -1) {
                        for (const auto &linkEnd : station) {
                                linkEndID = static_cast<int>(linkEnd.first);
                                stationName = linkEnd.second.stationName_;
                                break;

                        }
                }
                outputfile1 << stationName << "," << linkEndID << "," << time << "," << observation << std::endl;
        }
        //std::map<tudat::observation_models::ObservableType, std::map<int, std::vector<std::shared_ptr<SingleObservationSet < double, double>>>>> sortedObservationSets = observationsAndTimes->getSortedObservationSets();

        std::ofstream outputFile(saveDirectory + "observations_and_times_" + fileTag + ".txt");
        outputFile << std::setprecision(17);
        outputFile << "linkendtype,observable_type,time,observation\n";
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
    // Define estimation input
    std::shared_ptr< EstimationInput< double, double  > > estimationInput =
            std::make_shared< EstimationInput< double, double > >(
                    observationsAndTimes,
                    Eigen::MatrixXd::Zero(0,0),
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
    weightPerObservable[ two_way_doppler ] = std::pow(twoWayDopplerNoise, -2);

    estimationInput->setConstantPerObservableWeightsMatrix( weightPerObservable );
    std::cout<<"estimation input created"<<std::endl;



    std::shared_ptr< CovarianceAnalysisInput< double, double > > covarianceInput =
            std::make_shared< CovarianceAnalysisInput< double, double > >(
                    observationsAndTimes );
    std::cout<<"covariance input created"<<std::endl;
        covarianceInput->setConstantPerObservableWeightsMatrix( weightPerObservable );

    std::shared_ptr< CovarianceAnalysisOutput< double, double > > covarianceOutput = orbitDeterminationManager.computeCovariance(
            covarianceInput );
    std::cout<<"covariance output created"<<std::endl;

        Eigen::MatrixXd correlationMatrix = covarianceOutput->getCorrelationMatrix( );
        std::ofstream file10(saveDirectory + "correlations_" + fileTag + ".txt");
        file10 << std::setprecision( 21 ) << correlationMatrix ;
        file10.close();

        Eigen::MatrixXd covarianceMatrix = covarianceOutput->getUnnormalizedCovarianceMatrix( );
        std::ofstream file11(saveDirectory + "covanCovariances_" + fileTag + ".txt");
        file11 << std::setprecision( 21 ) << covarianceMatrix ;
        file11.close();


    // Perform estimation
        if (performEst) {
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
                std::ofstream fe(saveDirectory + "matrix_"+ fileTag + ".csv");
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

                std::ofstream errorFile(saveDirectory+"TrueError_"+fileTag+".txt");
                errorFile<<std::setprecision(17);
                errorFile<<TrueError<<std::endl;
                errorFile.close();
                std::ofstream formalErrorFile(saveDirectory+"FormalError_"+fileTag+".txt");
                formalErrorFile<<std::setprecision(17);
                formalErrorFile<<FormalError<<std::endl;
                formalErrorFile.close();

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
                                        "", 18, 18 );
        }
}
int main() {
        int iterationNumber = 6;
        //std::vector<int> arcLengths = {1, 3, 5, 10, 30};
        std::vector<int> arcLengths = {3, 5};
	//std::vector<int> number_of_arcs = {60, 20, 12, 6, 2};
        std::vector<int> number_of_arcs = {5, 3};
	std::vector<double>  hoursperday = {10.0,  24.0};
        std::vector<int> ihoursperday = {10, 24};
        //std::vector<double>  initialTimes = {0.0,  86400.0 * 360};
        std::vector<double>  initialTimes = {-240.0*86400.0, -180*86400.0, 0.0, 180.0*86400.0 , 86400.0 * 360};
	std::vector<int> startTime = {-240, -180, 0, 180 ,360};
        //std::vector<double>  perturbPos = {100.0, 1000.0};
        std::vector<double>  perturbPos = {100.0};
	std::vector<int> intperturbPos = {100};
        //std::vector<double>  perturbVel = {0.001, 0.01};
        std::vector<double>  perturbVel = {0.001};
	std::vector<int> intperturbVel = {0001};
        bool performEst = false;
        for (int i = 0;i<arcLengths.size(); i++) {
                for (int hours = 0; hours < hoursperday.size(); hours++) {
                        for (int initialTime = 0; initialTime < initialTimes.size(); initialTime++) {
                                double finalTime = initialTimes[initialTime] + 86400.0*15.0;
                                for (int intperturb = 0; intperturb < perturbPos.size(); intperturb++) {
                                        arcLengthRuns( hoursperday[hours],  initialTimes[initialTime],  finalTime,
                                                arcLengths[i],  iterationNumber, perturbPos[intperturb], perturbVel[intperturb],
                                                number_of_arcs[i],  15,  startTime[initialTime],  intperturbPos[intperturb],
                                                intperturbVel[intperturb],ihoursperday[hours],performEst );
                                }
                        }
                }
        }
}