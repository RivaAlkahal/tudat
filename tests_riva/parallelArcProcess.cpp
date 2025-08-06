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

struct GravityCoefficient {
        int degree = 0;    // Degree (n)
        int order = 0;     // Order (m)
        double Cnm = 0.0;  // Cosine coefficient
        double Snm = 0.0;  // Sine coefficient
        double CnmErr = 0.0; // Error in Cnm
        double SnmErr = 0.0; // Error in Snm
};

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
void runCovarianceAnalysisForArc(int arcIndex, Eigen::MatrixXd& covarianceMatrix, Eigen::MatrixXd& designMatrix)
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

        //std::string fileTag = "DEBUGnoapr234001poly_aprperiodicand18StaticandDragperarc_" +  std::to_string(arcLength) + std::to_string(itotalDuration) +  "darc_startat" + std::to_string(startTime)
     //                     + "_" + std::to_string(finalTime) + "_" + std::to_string(intperturbPos) + "perpos_" + std::to_string(intperturbVel) + "_pervel_" + std::to_string(ihoursperday) + "hobs_" + std::to_string(iterationNumber) + "itr_spiceprop" ;
    std::string fileTag = "accumul_InverseAprALLGlobalPars_drag_" +  std::to_string(arcLength) + std::to_string(itotalDuration) +  "darc_startat" + std::to_string(startTime)
                                  + "_" + std::to_string(finalTime);


    covarianceMatrix = Eigen::MatrixXd::Identity(6,6) * arcIndex;
    designMatrix = Eigen::MatrixXd::Random(6,10);

}

int main(int argc, char* argv[]){
    // 1) Read dimension (or default to 5000)
    
    if(argc < 2){
	std::cerr<< "Usage: " << argv[0] << "--arc-index <index> [--output-dir <dir>]" << std::endl;
	return EXIT_FAILURE;
    }

    int arcIndex = -1;
    std::string outputDir = ".";

    for (int i = 1; i<argc; ++i)
    {
	std::string arg = argv[i];
	if (arg == "--arc-index" && i + 1 <argc)
	{
	    arcIndex = std::atoi(argv[++i]);
	}
	else if (arg == "--output-dir" && i + 1 <argc)
	{
	    outputDir = argv[++i];
	}
    }

    if (arcIndex < 0)
    {
	std::cerr<< "Error: --arc-index must be provided and >=0" << std::endl;
	return EXIT_FAILURE;
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
    std::cout<< "Running arc index: "<< arcIndex <<std::endl;

    Eigen::MatrixXd covarianceMatrix, designMatrix;
    runCovarianceAnalysisForArc(arcIndex, covarianceMatrix, designMatrix);


    // save results
    try
    {
	saveMatrixBinary(outputDir + "/covariance_matrix.bin", covarianceMatrix);
	saveMatrixBinary(outputDir + "/design_matrix.bin", designMatrix);
    }
    catch (const std::exception& ex)
    {
	std::cerr << "Error saving matrices: " << ex.what() << std::endl;
	return EXIT_FAILURE;
    }

    std::cout<< "Arc  "<< arcIndex << " completed successfully!" <<std::endl;
    return EXIT_SUCCESS;
}
