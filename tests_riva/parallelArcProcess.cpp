//
// Created by ralkahal on 06-08-25.
//

#include <iostream>
#include <thread>
#include <chrono>
#include <vector>
#include <sstream>
#include <fstream>
#include <Eigen/Dense>
#include <cstdlib>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <boost/test/unit_test.hpp>
#include <stdexcept>
#include "nlohmann/json.hpp"

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

struct ArcConfig
{
    int arcIndex;
    double arcStart;
    double arcEnd;
    std::vector<double> obsStartTimes;
    std::vector<double> obsEndTimes;
};
ArcConfig loadArcConfig(const std::string& configFilePath) {
    std::ifstream configFile(configFilePath);
    if (!configFile.is_open()) {
        throw std::runtime_error("Could not open config file: " + configFilePath);
    }

    nlohmann::json configJson;
    configFile >> configJson;

    ArcConfig arcConfig;
    arcConfig.arcIndex = configJson["arc_index"];
    arcConfig.arcStart = configJson["arc_start"];
    arcConfig.arcEnd = configJson["arc_end"];
    arcConfig.obsStartTimes = configJson["obs_start_times"].get<std::vector<double>>();
    arcConfig.obsEndTimes = configJson["obs_end_times"].get<std::vector<double>>();

    return arcConfig;
}

inline void saveParamCounts(const std::string& outDir,
                            int numberOfLocalParameters,
                            int numberOfGlobalParameters) // pass 0 if none / not needed
{
        std::ofstream f(outDir + "/param_counts.txt");
        if (!f) throw std::runtime_error("Cannot open param_counts.txt for writing");
        // Format: one integer per line
        f << numberOfLocalParameters << "\n"
          << numberOfGlobalParameters << "\n";
}
struct GravityCoefficient {
        int degree = 0;    // Degree (n)
        int order = 0;     // Order (m)
        double Cnm = 0.0;  // Cosine coefficient
        double Snm = 0.0;  // Sine coefficient
        double CnmErr = 0.0; // Error in Cnm
        double SnmErr = 0.0; // Error in Snm
};

void loadGravityFieldFile(const std::string& filename,
                          std::vector<GravityCoefficient>& coefficients) {
        std::ifstream file(filename);
        if (!file.is_open()) {
                throw std::runtime_error("Unable to open file: " + filename);
        }
        // Skip the header line
        std::string header;
        if (!std::getline(file, header)) {
                throw std::runtime_error("File is empty or header is missing.");
        }

        // Read coefficients
        std::string line;
        while (std::getline(file, line)) {
                std::replace(line.begin(), line.end(), ',', ' ');

                std::istringstream lineStream(line);
                GravityCoefficient coeff;

                lineStream >> coeff.degree >> coeff.order
                           >> coeff.Cnm >> coeff.Snm
                           >> coeff.CnmErr >> coeff.SnmErr;

                if (lineStream.fail()) {
                        throw std::runtime_error("File format is incorrect in the coefficients section.");
                }

                coefficients.push_back(coeff);
        }

        file.close();
}


void extractErrorsWithinRange(
    const std::vector<GravityCoefficient>& coefficients,
    int minDegree, int maxDegree,
    std::vector<double>& cnmErrors,
    std::vector<double>& snmErrors
) {
        for (const auto& coeff : coefficients) {
                if (coeff.degree >= minDegree && coeff.degree <= maxDegree && coeff.order <= maxDegree) {
                        cnmErrors.push_back(coeff.CnmErr);
                        if (coeff.SnmErr != 0.0) {
                                snmErrors.push_back(coeff.SnmErr);
                        }
                }
        }
}


void ensureDirectoryExists(const std::string& path)
{
    struct stat info;
    if(stat(path.c_str(), &info)!=0)
    {
	if (mkdir(path.c_str(), 0755)!=0)
	{
	    throw std::runtime_error("Failed to create output directory: " + path);
	}
    }
}
void saveMatrixBinary(const std::string& filename, const Eigen::MatrixXd& matrix)
{
    std::ofstream out(filename, std::ios::binary);
    if(!out)
    {
	throw std::runtime_error("Cannot open file: " + filename);
    }

    Eigen::Index rows = matrix.rows(),cols = matrix.cols();
    out.write(reinterpret_cast<const char*>(&rows), sizeof(Eigen::Index));
    out.write(reinterpret_cast<const char*>(&cols), sizeof(Eigen::Index));
    out.write(reinterpret_cast<const char*>(matrix.data()),rows*cols*sizeof(double));
    out.close();
}
void saveVectorBinary(const std::string& filename, const Eigen::VectorXd& vector)
{
        std::ofstream out(filename, std::ios::binary);
        if (!out)
                throw std::runtime_error("Cannot open file: " + filename);

        Eigen::Index size = vector.size();
        out.write(reinterpret_cast<const char*>(&size), sizeof(Eigen::Index));
        out.write(reinterpret_cast<const char*>(vector.data()), size * sizeof(double));
}

void runCovarianceAnalysisForArc(const ArcConfig& config,  std::string saveDirectory, Eigen::MatrixXd& covarianceMatrix, Eigen::MatrixXd& unnormalizedDesignMatrix,  Eigen::MatrixXd& normalizedDesignMatrix,Eigen::VectorXd& normalizationFactor, Eigen::VectorXd& weightMatrixDiagonal, Eigen::MatrixXd& P0_matrix)
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

    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext1.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext2.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext3.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext4.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext5.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext6.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext7.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext8.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext9.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext10.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext11.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext12.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext13.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext14.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext15.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext16.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext17.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext18.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext19.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext20.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext21.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext22.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext23.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext24.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext25.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext26.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_50year_nominal.bsp");

    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_map4_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_map5_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_map6_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_map7_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_map8_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext1_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext2_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext3_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext4_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext5_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext6_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext7_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext8_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext9_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext10_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext11_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext12_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext13_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext14_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext15_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext16_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext17_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext18_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext19_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext20_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext21_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext22_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext23_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext24_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext25_ipng_mgs95j.bsp" );
    spice_interface::loadSpiceKernelInTudat( "/home/ralkahal/new-tudat-tests/mgs_ext26_ipng_mgs95j.bsp" );
    int arcIndex = config.arcIndex;
    double arcStart = config.arcStart;
    double arcEnd = config.arcEnd;
    std::vector<double> obsStartTimes = config.obsStartTimes;
    std::vector<double> obsEndTimes = config.obsEndTimes;
    double arcLength = arcEnd - arcStart;
    double arcDuration = arcLength * 86400.0; // Convert days to seconds
    double twoWayDopplerNoise = 0.0001;
    double initialEphemerisTime = arcStart-120.0;
    double finalEphemerisTime = arcEnd+120.0;
    double epehemeridesTimeStep = 60.0;
    double ephemerisTimeStepPlanets =  epehemeridesTimeStep;
    double ephemerisTimeStepSpacecraft = epehemeridesTimeStep ;
    double buffer = 30.0 * epehemeridesTimeStep;
        double step_size;
        string dragEst = "per-rev";
        double hoursperdaydrag = 2.0;
        int ndays = 1;
        if (dragEst == "per-rev") {
                step_size = hoursperdaydrag * 3600;
        }
        else
        {
                step_size = ndays * 24 * 3600;
        }
        std::vector<double> initial_times_list_emp;
        for (double time =arcStart +buffer; time <= arcEnd-buffer; time += step_size) {
                if (arcEnd-time < step_size) {
                        break;
                }
                initial_times_list_emp.push_back(time);
        }

        std::string fileTag = "accumul_InverseAprALLGlobalPars_drag_5" + std::to_string(arcIndex) +  "darc_startat" + std::to_string(arcStart) + "_" + std::to_string(arcEnd);

    // bodies to create
    std::vector< std::string > bodiesToCreate = {
        "Earth", "Sun", "Mercury", "Venus", "Mars", "Jupiter", "Phobos", "Deimos" };

    std::string baseFrameOrientation = "MARSIAU";
    std::string baseFrameOrigin = "SSB";

    BodyListSettings bodySettings;
    bodySettings = getDefaultBodySettings(
                    bodiesToCreate, initialEphemerisTime - buffer, finalEphemerisTime + buffer,
                    baseFrameOrigin, baseFrameOrientation, ephemerisTimeStepPlanets );
    std::cout<<"Body Settings created"<<std::endl;
    bodySettings.at( "Earth" )->groundStationSettings = getDsnStationSettings( );
    std::string filename = "/home/ralkahal/new-tudat-tests/dtm-mars";
    bodySettings.at( "Mars" )->atmosphereSettings = marsDtmAtmosphereSettings( filename, 3378.0E3);

    std::string spacecraftName = "MGS";
    bodySettings.addSettings( spacecraftName );
    bodySettings.at( spacecraftName )->ephemerisSettings =
                std::make_shared< InterpolatedSpiceEphemerisSettings >(
                        initialEphemerisTime - buffer, finalEphemerisTime + buffer,
                        ephemerisTimeStepSpacecraft, baseFrameOrigin, baseFrameOrientation,
                        std::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 ), spacecraftName );
    bodySettings.at( spacecraftName )->constantMass = 1030.5;
    // Set gravity field variations
    std::vector< std::shared_ptr< GravityFieldVariationSettings > > gravityFieldVariations;

    // Set solid body tide gravity field variation
    std::vector< std::string > deformingBodies;
    deformingBodies.push_back( "Sun" );
    deformingBodies.push_back( "Phobos" );
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
    sineAmplitudes[ 1 ] = Eigen::Matrix< double, 4, 5 >::Zero( );
    sineAmplitudes[1](0,1) +=1.25916272835878e-13/(365*24*3600);
    sineAmplitudes[1](0,2) += -8.84999663538266e-13/(365*24*3600);
    sineAmplitudes[1](1,1) += 6.0445217093141e-13/(365*24*3600);
    sineAmplitudes[1](1,2) += 1.08851817438134e-13/(365*24*3600);
    sineAmplitudes[1](1,3) += 2.88244592826493e-13/(365*24*3600);
    sineAmplitudes[1](2,1) += -1.38365049897319e-13/(365*24*3600);
    sineAmplitudes[1](2,2) += 4.70640224945299e-13/(365*24*3600);
    sineAmplitudes[1](2,3) += -4.50562269498905e-14/(365*24*3600);
    sineAmplitudes[1](2,4) += 8.5339608889136e-13/(365*24*3600);

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
    double referenceAreaRadiation = 15.0;
    double radiationPressureCoefficient = 2.1;
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
                std::make_shared<ConstantAerodynamicCoefficientSettings>(15.0, 2.1 * Eigen::Vector3d::UnitX());
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
        dependentVariablesToSave.push_back(std::make_shared<SingleDependentVariableSaveSettings>(
                local_density_dependent_variable,spacecraftName, centralBody));


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

        //read in observation times
        std::vector< double > observationTimesList;
        for (int t = 0; t<obsEndTimes.size(); ++t)
        {
                double start = obsStartTimes[t];
                double end = obsEndTimes[t];
                for (double time = start; time < end; time += 10.0) // 10 seconds interval
                {
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
        std::shared_ptr<IntegratorSettings<> >integratorSettings =
               std::make_shared<RungeKuttaFixedStepSizeSettings<> >( 30, CoefficientSets::rungeKutta87DormandPrince );
        std::cout<<"Integration settings created"<<std::endl;

        const std::string filenameGrav = "/home/ralkahal/nnew-tudat-tests/jgmro_120d_sha.tab";
        std::vector<GravityCoefficient> coefficients;
        // Filter and extract errors up to degree and order 8
        int maxDegree = 18;
        int minDegree = 2;
        std::vector<double> cnmErrors;
        std::vector<double> snmErrors;

        try {
                // Load the gravity field file
                loadGravityFieldFile(filenameGrav, coefficients);

                extractErrorsWithinRange(coefficients, minDegree, maxDegree, cnmErrors, snmErrors);
                std::cout<<"size of Cnm errors: " << cnmErrors.size() << std::endl;
                std::cout<<"size of Snm errors: " << snmErrors.size() << std::endl;
        } catch (const std::exception& e) {
                std::cerr << "Error: " << e.what() << std::endl;
                exit(1);
        }

        // Create link ends

        // Create list of link ends where the ground station is the transmitter and the spacecraft is the receiver
        std::vector< LinkEnds > stationTransmitterLinkEnds;
        std::vector< LinkEnds > downlinkLinkEnds_;
        std::vector< LinkEnds > uplinkLinkEnds_;
        std::vector< std::string > GroundStations = {  "DSS-26" , "DSS-42", "DSS-61"};
        for ( std::string groundStation : GroundStations ) {
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

        // Create observation viability settings and calculators
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

        Eigen::VectorXd systemInitialStates;
        systemInitialStates = spice_interface::getBodyCartesianStateAtEpoch(
                         bodiesToIntegrate[ 0 ], "Mars", "MARSIAU", "NONE", arcStart);
        std::cout<<"system initial states created"<<std::endl;
        std::shared_ptr< PropagationTerminationSettings > terminationSettings = propagationTimeTerminationSettings(
                        arcEnd);
        std::shared_ptr< TranslationalStatePropagatorSettings< double, double> > propagatorSettings = translationalStatePropagatorSettings< double, double >( centralBodies, accelerationModelMap, bodiesToIntegrate,
                                                                                                                                                      systemInitialStates, arcStart, integratorSettings, terminationSettings, cowell, dependentVariablesToSave);
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                                bodies, propagatorSettings );
        std::map< double, Eigen::VectorXd > stateHistory = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
        writeDataMapToTextFile( stateHistory, "stateHistoryPropagation_arc_" + fileTag + ".txt", saveDirectory,
                        "", 18, 18 );
        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
                            getInitialStateParameterSettings< double, double  >( propagatorSettings, bodies);
        parameterNames.push_back(std::make_shared< EstimatableParameterSettings >(spacecraftName,constant_drag_coefficient));

        parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                         2, 0, 18, 18, "Mars", spherical_harmonics_cosine_coefficient_block ) );
        parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                              2, 1, 18, 18, "Mars", spherical_harmonics_sine_coefficient_block ) );

        std::map<int, std::vector<std::pair<int, int> > > cosineBlockIndicesPerPeriod;
                //periodic gravity field
        cosineBlockIndicesPerPeriod[ 0 ].push_back( std::make_pair( 2, 0) );
        cosineBlockIndicesPerPeriod[ 0 ].push_back( std::make_pair( 3, 0 ) );
        cosineBlockIndicesPerPeriod[ 0 ].push_back( std::make_pair( 4, 0 ) );
        cosineBlockIndicesPerPeriod[ 0 ].push_back( std::make_pair( 5, 0 ) );

        cosineBlockIndicesPerPeriod[ 1 ].push_back( std::make_pair( 2, 0) );
        cosineBlockIndicesPerPeriod[ 1 ].push_back( std::make_pair( 3, 0 ) );
        cosineBlockIndicesPerPeriod[ 1 ].push_back( std::make_pair( 4, 0 ) );
        cosineBlockIndicesPerPeriod[ 1 ].push_back( std::make_pair( 5, 0 ) );

        cosineBlockIndicesPerPeriod[ 2 ].push_back( std::make_pair( 2, 0) );
        cosineBlockIndicesPerPeriod[ 2 ].push_back( std::make_pair( 3, 0 ) );
        cosineBlockIndicesPerPeriod[ 2 ].push_back( std::make_pair( 4, 0 ) );
        cosineBlockIndicesPerPeriod[ 2 ].push_back( std::make_pair( 5, 0 ) );
        std::cout<<" cosine block indices per period size"<<cosineBlockIndicesPerPeriod[0].size()<<std::endl;

        std::map<int, std::vector<std::pair<int, int> > > sineBlockIndicesPerPeriod;
        parameterNames.push_back( std::make_shared< PeriodicGravityFieldVariationEstimatableParameterSettings >(
                                centralBody, cosineBlockIndicesPerPeriod, sineBlockIndicesPerPeriod ) );

        std::map<int, std::vector<std::pair<int, int> > > cosineBlockIndicesPerPower;

        cosineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 2, 0 ) );
        cosineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 2, 1 ) );
        cosineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 2, 2 ) );
        cosineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 3, 0 ) );
        cosineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 3, 1 ) );
        cosineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 3, 2 ) );
        cosineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 3, 3 ) );
        cosineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 4, 0 ) );
        cosineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 4, 1 ) );
        cosineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 4, 2 ) );
        cosineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 4, 3 ) );
        cosineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 4, 4 ) );
	
        std::map<int, std::vector<std::pair<int, int> > > sineBlockIndicesPerPower;
        sineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 2, 1 ) );
        sineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 2, 2 ) );
        sineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 3, 1 ) );
        sineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 3, 2 ) );
        sineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 3, 3 ) );
        sineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 4, 1 ) );
        sineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 4, 2 ) );
        sineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 4, 3 ) );
        sineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 4, 4 ) );

        parameterNames.push_back( std::make_shared< PolynomialGravityFieldVariationEstimatableParameterSettings >(
                "Mars", cosineBlockIndicesPerPower, sineBlockIndicesPerPower ) );


        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
                    createParametersToEstimate< double, double >( parameterNames, bodies, propagatorSettings );

        std::cout<<"parameters to estimate created"<<std::endl;

        // Create orbit determination object.
        OrbitDeterminationManager< double, double > orbitDeterminationManager =
                OrbitDeterminationManager< double, double >(
                        bodies, parametersToEstimate,
                        observationSettingsList, propagatorSettings,true );
        std::cout<<"orbit determination manager created"<<std::endl;
        Eigen::Matrix< double, Eigen::Dynamic, 1 > initialParameterEstimate =
                    parametersToEstimate->template getFullParameterValues< double >( );

        std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput;
        for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
                linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
        {
                ObservableType currentObservable = linkEndIterator->first;
                std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
                for( unsigned int currLinkEnd = 0; currLinkEnd < currentLinkEndsList.size( ); currLinkEnd++ )
                {
                        measurementSimulationInput.push_back(
                                std::make_shared< TabulatedObservationSimulationSettings< > >(
                                                currentObservable, currentLinkEndsList[ currLinkEnd ], observationTimesList, receiver, observationViabilitySettings, noiseFunctions[currentObservable]) );
                }
        }

        // Simulate observations.
        std::shared_ptr< ObservationCollection< > > observationsAndTimes = simulateObservations< double, double >(
        measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );
        std::cout<<"observations and times created"<<std::endl;

        Eigen::Matrix< double, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;
        int numberOfParameters = initialParameterEstimate.rows( );

        printEstimatableParameterEntries( parametersToEstimate );
        int lengthOfTimeListEmp = initial_times_list_emp.size();
        int numberOfLocalParameters = 7;

        const int DIAGONALS = numberOfLocalParameters;
        Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(numberOfParameters, numberOfParameters);
        double aprioriuncertainty =1.0/(10*10);
        matrix(6,6) = aprioriuncertainty;
        
        for (int j = 0; j < cnmErrors.size(); ++j) {
                matrix(j+DIAGONALS,j+DIAGONALS) = 1.0/(cnmErrors[j]*cnmErrors[j]);
        }
        std::cout<<"matrix filled with Cnm errors size: "<< cnmErrors.size()+DIAGONALS <<std::endl;
        for (int j= 0; j < snmErrors.size(); ++j) {
                matrix(j+DIAGONALS+cnmErrors.size(),j+DIAGONALS+cnmErrors.size()) = 1.0/(snmErrors[j]*snmErrors[j]);
        }

        matrix(DIAGONALS+cnmErrors.size()+snmErrors.size(), DIAGONALS+cnmErrors.size()+snmErrors.size()) = 1.0/(0.016E-09*0.016E-09);
        matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+1, DIAGONALS+cnmErrors.size()+snmErrors.size()+1) = 1.0/(0.016E-09*0.016E-09);
        matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+2, DIAGONALS+cnmErrors.size()+snmErrors.size()+2) = 1.0/(0.011E-09*0.011E-09);
        matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+3, DIAGONALS+cnmErrors.size()+snmErrors.size()+3) = 1.0/(0.011E-09*0.011E-09);
        matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+4, DIAGONALS+cnmErrors.size()+snmErrors.size()+4) = 1.0/(0.101E-10*0.101E-10);
        matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+5, DIAGONALS+cnmErrors.size()+snmErrors.size()+5) = 1.0/(0.101E-10*0.101E-10);
        matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+6, DIAGONALS+cnmErrors.size()+snmErrors.size()+6) = 1.0/(0.010E-09*0.010E-09);
        matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+7, DIAGONALS+cnmErrors.size()+snmErrors.size()+7) = 1.0/(0.010E-09*0.010E-09);

        matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+8, DIAGONALS+cnmErrors.size()+snmErrors.size()+8) = 1.0/(0.016E-09*0.016E-09);
        matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+9, DIAGONALS+cnmErrors.size()+snmErrors.size()+9) = 1.0/(0.016E-09*0.016E-09);
        matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+10, DIAGONALS+cnmErrors.size()+snmErrors.size()+10) = 1.0/(0.011E-09*0.011E-09);
        matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+11, DIAGONALS+cnmErrors.size()+snmErrors.size()+11) = 1.0/(0.011E-09*0.011E-09);
        matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+12, DIAGONALS+cnmErrors.size()+snmErrors.size()+12) = 1.0/(0.101E-10*0.101E-10);
        matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+13, DIAGONALS+cnmErrors.size()+snmErrors.size()+13) = 1.0/(0.101E-10*0.101E-10);
        matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+14, DIAGONALS+cnmErrors.size()+snmErrors.size()+14) = 1.0/(0.010E-09*0.010E-09);
        matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+15, DIAGONALS+cnmErrors.size()+snmErrors.size()+15) = 1.0/(0.010E-09*0.010E-09);

        matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+16, DIAGONALS+cnmErrors.size()+snmErrors.size()+16) = 1.0/(0.016E-09*0.016E-09);
        matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+17, DIAGONALS+cnmErrors.size()+snmErrors.size()+17) = 1.0/(0.016E-09*0.016E-09);
        matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+18, DIAGONALS+cnmErrors.size()+snmErrors.size()+18) = 1.0/(0.011E-09*0.011E-09);
        matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+19, DIAGONALS+cnmErrors.size()+snmErrors.size()+19) = 1.0/(0.011E-09*0.011E-09);
        matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+20, DIAGONALS+cnmErrors.size()+snmErrors.size()+20) = 1.0/(0.101E-10*0.101E-10);
        matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+21, DIAGONALS+cnmErrors.size()+snmErrors.size()+21) = 1.0/(0.101E-10*0.101E-10);
        matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+22, DIAGONALS+cnmErrors.size()+snmErrors.size()+22) = 1.0/(0.010E-09*0.010E-09);
        matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+23, DIAGONALS+cnmErrors.size()+snmErrors.size()+23) = 1.0/(0.010E-09*0.010E-09);

        matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+24, DIAGONALS+cnmErrors.size()+snmErrors.size()+24) = 0.0;
        matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+25, DIAGONALS+cnmErrors.size()+snmErrors.size()+25) = 0.0;
        matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+26, DIAGONALS+cnmErrors.size()+snmErrors.size()+26) = 0.0;
        matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+27, DIAGONALS+cnmErrors.size()+snmErrors.size()+27) = 0.0;
        matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+28, DIAGONALS+cnmErrors.size()+snmErrors.size()+28) = 0.0;
        
	P0_matrix = matrix;
        saveParamCounts(saveDirectory, numberOfLocalParameters, 0);
	std::shared_ptr< CovarianceAnalysisInput< double, double > > covarianceInput =
                                std::make_shared< CovarianceAnalysisInput< double, double > >(
                                observationsAndTimes,matrix );
        std::cout<<"covariance input created"<<std::endl;
        std::map< observation_models::ObservableType, double > weightPerObservable;
        weightPerObservable[ two_way_doppler ] = std::pow(twoWayDopplerNoise, -2);

        covarianceInput->setConstantPerObservableWeightsMatrix( weightPerObservable );


        std::shared_ptr< CovarianceAnalysisOutput< double, double > > covarianceOutput = orbitDeterminationManager.computeCovariance(
                    covarianceInput );
        std::cout<<"covariance output created"<<std::endl;
        unnormalizedDesignMatrix = covarianceOutput->getUnnormalizedDesignMatrix( );
        normalizedDesignMatrix = covarianceOutput->getNormalizedDesignMatrix( );
        weightMatrixDiagonal = covarianceOutput->weightsMatrixDiagonal_;

	normalizationFactor = covarianceOutput->designMatrixTransformationDiagonal_;
        covarianceMatrix = covarianceOutput->getUnnormalizedCovarianceMatrix( );;

}

int main(int argc, char* argv[]){

    if(argc < 2){
	std::cerr<< "Usage: " << argv[0] << "--arc-index <index> [--output-dir <dir>]" << std::endl;
	return EXIT_FAILURE;
    }

    std::string configFilePath;
    std::string outputDir = ".";

    for (int i = 1; i<argc; ++i)
    {
	    std::string arg = argv[i];
	    if (arg == "--config" && i + 1 <argc)
	    {
	        configFilePath = argv[++i];
	    }
	    else if (arg == "--output-dir" && i + 1 <argc)
	    {
	        outputDir = argv[++i];
	    }
    }
    try
    {
        ensureDirectoryExists(outputDir);
    }
    catch(const std::exception& ex)
    {
        std::cerr << "Error: " << ex.what() << std::endl;
        return EXIT_FAILURE;
    }
    if (configFilePath.empty())
    {
        std::cerr << "Error: --config must be provided\n";
        return EXIT_FAILURE;
    }
    ArcConfig config;
    try
    {
        config = loadArcConfig(configFilePath);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Failed to load config: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    std::string mkdirCmd = "mkdir -p " + outputDir;
    std::system(mkdirCmd.c_str());
    int arcIndex = config.arcIndex;
    std::cout<< "Running arc index: "<< arcIndex <<std::endl;

    Eigen::MatrixXd covarianceMatrix, unnormalizedDesignMatrix, normalizedDesignMatrix, P0_matrix;
    Eigen::VectorXd weightMatrixDiagonal, normalizationFactor;
    runCovarianceAnalysisForArc(config, outputDir, covarianceMatrix, unnormalizedDesignMatrix, normalizedDesignMatrix,normalizationFactor, weightMatrixDiagonal, P0_matrix);


    // save results
    try
    {
	saveMatrixBinary(outputDir + "/covariance_matrix.bin", covarianceMatrix);
	saveMatrixBinary(outputDir + "/unnormalized_Design_Matrix.bin", unnormalizedDesignMatrix);
        saveMatrixBinary(outputDir + "/normalized_Design_Matrix.bin", normalizedDesignMatrix);
        saveMatrixBinary(outputDir + "/P0_matrix.bin", P0_matrix);
        saveVectorBinary(outputDir + "/weight_matrix_diagonal.bin", weightMatrixDiagonal);
	saveVectorBinary(outputDir + "/normalization_factor.bin", normalizationFactor);

    }
    catch (const std::exception& ex)
    {
	std::cerr << "Error saving matrices: " << ex.what() << std::endl;
	return EXIT_FAILURE;
    }

    std::cout<< "Arc  "<< arcIndex << " completed successfully!" <<std::endl;
    return EXIT_SUCCESS;
}
