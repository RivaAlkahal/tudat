//
// Created by ralkahal on 23-4-25.
//
//
// Created by ralkahal on 9-4-25.
//
//
// Created by ralkahal on 17-2-25.
//
//
// Created by ralkahal on 22-1-25.
//
//
// Created by ralkahal on 23-12-24.
//
//
// Created by Riva Alkahal on 20/12/2024.
//
//
// Created by Riva Alkahal on 10/12/2024.
//
//
// Created by Riva Alkahal on 22/11/2024.
//
//
// Created by Riva Alkahal on 07/11/2024.
//

#include <iostream>
#include <fstream>
#include <limits>
#include "fstream"
#include "iostream"
#include <Eigen/Dense>

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



std::string saveDirectory = "/home/ralkahal/nnew-tudat-tests/covAn/twowaydoppler_apriori/accumsResultsHigherOrder/";

struct GravityCoefficient {
        int degree = 0;    // Degree (n)
        int order = 0;     // Order (m)
        double Cnm = 0.0;  // Cosine coefficient
        double Snm = 0.0;  // Sine coefficient
        double CnmErr = 0.0; // Error in Cnm
        double SnmErr = 0.0; // Error in Snm
};

void computeCovarianceMatrix(int startTime, int nIndex, int numberOfGlobalParameters, std::vector<Eigen::MatrixXd> normalizedDesignMatrices,std::vector<Eigen::MatrixXd> unnormalizedDesignMatrices,std::vector<Eigen::VectorXd> normalizationFactors, std::vector<Eigen::VectorXd> weightDiagonals,std::vector<Eigen::MatrixXd> P0_matrices) {

	int itotalDuration = nIndex * 5;
	double totalDur = nIndex * 5.0;
	double finalTime = totalDur * 86400.0;
	std::string fileTag = "accumul_InverseAprALLGlobalPars_drag_5" + std::to_string(itotalDuration) +  "darc_startat" + std::to_string(startTime) + "_" + std::to_string(finalTime);
	//fileTag = nIndex + "nArcs" + fileTag;
	std::cout<<"covariance analysis for" << nIndex <<" arcs done, starting summing them up..."<<std::endl;
        std::cout<< "number of global parameters: "<< numberOfGlobalParameters<<std::endl;
        int stateVectorSize =6;
        int nLocal=1;
        int numberOfLocalParameters = stateVectorSize + nLocal;
	int numCols = numberOfGlobalParameters;
        int nArcs = normalizedDesignMatrices.size();
        int arcMatrixSize = normalizedDesignMatrices[0].rows();
        int total_size = nArcs*(numberOfLocalParameters) + numberOfGlobalParameters;
        Eigen::MatrixXd P_global = Eigen::MatrixXd::Zero(total_size,total_size);

        // Fill the matrix with the values of the diagonal
        int DIAGONALS = numberOfLocalParameters;

        // merge normalization factors for all matrices, except for the global parameters
        Eigen::VectorXd normalizationFactorsMerged = Eigen::VectorXd::Zero(total_size);
        //Eigen::VectorXd localNormalizationFactors= Eigen::VectorXd::Zero(nArcs * nLocal);;
	size_t stateOffst= 0;
	size_t localOffst = nArcs*6;
        for (int i = 0; i < nArcs; i++){
            normalizationFactorsMerged.segment(stateOffst + i*(6),6) = normalizationFactors[i].segment(0,6);
            std::cout<<"first normalization factor merged"<<std::endl;
            normalizationFactorsMerged.segment(localOffst+i*nLocal,nLocal) = normalizationFactors[i].segment(6,nLocal);
            std::cout << P0_matrices[i].diagonal().segment(6,nLocal).array() << std::endl;
                std::cout<< normalizationFactors[i].segment(6,nLocal).array() << std::endl;
            P0_matrices[i].diagonal().segment(6,nLocal).array() = P0_matrices[i].diagonal().segment(6,nLocal).array()/normalizationFactors[i].segment(6,nLocal).array().square();
            std::cout << P0_matrices[i].diagonal().segment(6,nLocal).array() << std::endl;
            std::cout<<"second normalization factor merged"<<std::endl;
        }
        //normalizationFactorsMerged.segment(nArcs*(numberOfLocalParameters-1),nLocal*nArcs) = localNormalizationFactors;
        std::cout<<"Normalization factors merged!"<<std::endl;
	//int totalRows = unnormalizedDesignMatrices.size() * unnormalizedDesignMatrices[0].rows();
	int totalRows = std::accumulate(unnormalizedDesignMatrices.begin(), unnormalizedDesignMatrices.end(),0,[](int sum, const Eigen::MatrixXd& mat) {return sum + mat.rows();});
	std::cout<<"totalRows"<< totalRows<<std::endl;
	Eigen::MatrixXd stackedGlobalParsDesignMatrix(totalRows,numCols);

	// Extract and stack columns
	int position = 0;
	for (const auto &mat: unnormalizedDesignMatrices){
	    for (size_t j = 0; j<numberOfGlobalParameters; j++){
		//std::cout<<"j"<<j<<std::endl;
		//std::cout<<numberOfLocalParameters+j<<std::endl;
		//std::cout<<mat.rows()<<std::endl;
		stackedGlobalParsDesignMatrix.block(position,j,mat.rows(),1) = mat.col(numberOfLocalParameters+j);

	    }
	    position += mat.rows();
	}
	std::cout<<"Stacked unnormalized design matrix of global parameters!"<<std::endl;
	std::cout<<"Now performing normalization..."<<std::endl;

	// Normalize the columns of the stacked global parameters of the design matrices
	Eigen::VectorXd normalizationFactorsGP(stackedGlobalParsDesignMatrix.cols());

	for (int k = 0; k < stackedGlobalParsDesignMatrix.cols(); k++){
	    Eigen::Index maxIndex;
	    double maxAbsVal = stackedGlobalParsDesignMatrix.col(k).cwiseAbs().maxCoeff(&maxIndex);
	    double maxVal = stackedGlobalParsDesignMatrix(maxIndex,k);
	    if (maxVal != 0.0){
		stackedGlobalParsDesignMatrix.col(k) /= maxVal;
		normalizationFactorsGP[k] = maxVal;
	    }
	}
        normalizationFactorsMerged.segment(nArcs*numberOfLocalParameters,numberOfGlobalParameters) = normalizationFactorsGP;
        std::cout<<"Normalization factors size: "<<normalizationFactorsMerged.size()<<std::endl;

        //for (int arc = 0; arc < nArcs; arc++){
        //    P0_matrices[arc].segment(numberOfLocalParameters,numberOfLocalParameters+numberOfGlobalParameters).array() /= normalizationFactorsGP.array();
        //}
        std::cout<< normalizationFactorsGP.array() << std::endl;
        std::cout<< normalizationFactorsGP.size() << std::endl;
        std::cout<< P0_matrices[0].diagonal().segment(numberOfLocalParameters,numberOfGlobalParameters).size()  << std::endl;
	std::cout<<"Done with normalization!"<<std::endl;
	std::cout<<"Now bringing back the values to the normalized design matrices..."<<std::endl;
	// Reconstructing normalized matrix
	position = 0;
	int R = 0;
	for (auto &mat: normalizedDesignMatrices){
	    for (size_t j = 0; j<numberOfGlobalParameters; j++){
	            mat.col(numberOfLocalParameters+j) = stackedGlobalParsDesignMatrix.block(position,j,mat.rows(),1);
	    }
	    position += mat.rows();
	    //std::ofstream feR(saveDirectory + "updateNormalizedDesignMatrix" + std::to_string(R)  + fileTag + ".txt");
            // Write the matrix to the file
            //feR << std::setprecision(32) <<mat;
            // Close the file
            //feR.close();
	    P0_matrices[R].diagonal().segment(numberOfLocalParameters,numberOfGlobalParameters).array() = P0_matrices[R].diagonal().segment(numberOfLocalParameters,numberOfGlobalParameters).array()/normalizationFactorsGP.array().square();
	    //std::ofstream feR(saveDirectory + "normalizedAprioriCovInv" + std::to_string(R)  + fileTag + ".txt");
            // Write the matrix to the file
            //feR << std::setprecision(32) <<P0_matrices[R];
            // Close the file
            //feR.close();
	    R +=1;
	}


	std::cout<<"Updated normalized design matrices!"<<std::endl;
	std::cout<<"Now computing the covariance matrices..."<<std::endl;
	std::vector<Eigen::MatrixXd> resultNormalizedInvCovMatrices;
	for (size_t i = 0; i < normalizedDesignMatrices.size(); i++){
	    Eigen::DiagonalMatrix<double,Eigen::Dynamic> W(weightDiagonals[i]);
	    Eigen::MatrixXd P_im = (normalizedDesignMatrices[i].transpose() * W * normalizedDesignMatrices[i]);
	    Eigen::MatrixXd resultNormalizedInvCovMatrix = P_im + P0_matrices[i];
	    resultNormalizedInvCovMatrices.push_back(resultNormalizedInvCovMatrix);
	}
	std::cout<<"Normalized inverse covariance matrices computed!"<<std::endl;
	std::cout<<"Now summing up the values for the global parameters..."<<std::endl;
	//Eigen::MatrixXd sumCovarianceMatrix = Eigen::MatrixXd::Zero(numberOfGlobalParameters, numberOfGlobalParameters);
        //Eigen::MatrixXd inverseUnnCovarianceMatrixSummed = Eigen::MatrixXd::Zero(numberOfGlobalParameters, numberOfGlobalParameters);
        Eigen::MatrixXd inverseNormalizedCovarianceMatrixSummedGP = Eigen::MatrixXd::Zero(numberOfGlobalParameters, numberOfGlobalParameters);
	for (int i = 0; i < resultNormalizedInvCovMatrices.size(); i++){
                //inverseUnnCovarianceMatrixSummed += inverseUnnormalizedCovarianceMatrices[i].block(7,7,numberOfGlobalParameters,numberOfGlobalParameters);
                //sumCovarianceMatrix += covarianceMatrices[i].block(7,7,numberOfGlobalParameters,numberOfGlobalParameters);
	        if (i == 0){
	        inverseNormalizedCovarianceMatrixSummedGP += resultNormalizedInvCovMatrices[i].block(numberOfLocalParameters,numberOfLocalParameters,numberOfGlobalParameters,numberOfGlobalParameters);
                } else {
                        inverseNormalizedCovarianceMatrixSummedGP += resultNormalizedInvCovMatrices[i].block(numberOfLocalParameters,numberOfLocalParameters,numberOfGlobalParameters,numberOfGlobalParameters) - P0_matrices[i].block(numberOfLocalParameters,numberOfLocalParameters,numberOfGlobalParameters,numberOfGlobalParameters);
	        }
	        }
	std::cout<<"Global parameters values summed up!"<<std::endl;
        std::cout<<"Now assembling the global normalized inverse covariance matrix..."<<std::endl;

        // fill in the global matrix
        int globalOffset = nArcs*(stateVectorSize+nLocal);

        for (int i = 0; i<nArcs;i++) {
                int localOffset = nArcs*stateVectorSize+i*nLocal;
                //arc covariance matrix
                P_global.block(i*stateVectorSize,i*stateVectorSize,stateVectorSize,stateVectorSize) = resultNormalizedInvCovMatrices[i].block(0,0,stateVectorSize,stateVectorSize);
                //arc-local covariance matrix
                P_global.block(i*stateVectorSize,localOffset,stateVectorSize,nLocal) = resultNormalizedInvCovMatrices[i].block(0,stateVectorSize,stateVectorSize,nLocal);
                //local-arc covariance matrix
                P_global.block(localOffset,i*stateVectorSize,nLocal,stateVectorSize) = resultNormalizedInvCovMatrices[i].block(0,stateVectorSize,stateVectorSize,nLocal).transpose();
                //local-local covariance matrix
                P_global.block(localOffset,localOffset,nLocal,nLocal) = resultNormalizedInvCovMatrices[i].block(stateVectorSize,stateVectorSize,nLocal,nLocal);

                //arc-global covariance matrix
                P_global.block(i*stateVectorSize,globalOffset,stateVectorSize,numberOfGlobalParameters) = resultNormalizedInvCovMatrices[i].block(0,stateVectorSize+nLocal,stateVectorSize,numberOfGlobalParameters);
                //global-arc covariance matrix
                P_global.block(globalOffset,i*stateVectorSize,numberOfGlobalParameters,stateVectorSize) = resultNormalizedInvCovMatrices[i].block(0,stateVectorSize+nLocal,stateVectorSize,numberOfGlobalParameters).transpose();
                //local-global covariance matrix
                P_global.block(localOffset,globalOffset,nLocal,numberOfGlobalParameters) = resultNormalizedInvCovMatrices[i].block(stateVectorSize,stateVectorSize+nLocal,nLocal,numberOfGlobalParameters);
                //global-local covariance matrix
                P_global.block(globalOffset,localOffset,numberOfGlobalParameters,nLocal) = resultNormalizedInvCovMatrices[i].block(stateVectorSize,stateVectorSize+nLocal,nLocal,numberOfGlobalParameters).transpose();
        }
        // global-global covariance matrix
        P_global.block(globalOffset,globalOffset,numberOfGlobalParameters,numberOfGlobalParameters) = inverseNormalizedCovarianceMatrixSummedGP;
        std::cout<<"Global normalized inverse covariance matrix assembled!"<<std::endl;
        std::cout<<"Now updating the normalization factors"<<std::endl;
        Eigen::VectorXd normalizationFactorsGP2 = Eigen::VectorXd::Zero(total_size);
        normalizationFactorsGP2.head(stackedGlobalParsDesignMatrix.cols()) = normalizationFactorsGP;
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(P_global, Eigen::ComputeThinU | Eigen::ComputeThinV);
        double tol = 1e-13;
	double cond = svd.singularValues()(0)/svd.singularValues().tail(1)(0);
	std::cout << "condition number: " << cond <<std::endl;

        std::cout<<"Now inverting the inverse covariance matrices..."<<std::endl;
        Eigen::MatrixXd NormalizedCovarianceMatrixSummed = P_global.completeOrthogonalDecomposition().pseudoInverse();
        /*
        Eigen::MatrixXd NormalizedCovarianceMatrixSummed;
        Eigen::VectorXd S_inv = svd.singularValues();
        for (int i = 0; i < S_inv.size(); ++i) {
            if (S_inv(i) > tol) {
                S_inv(i) = 1.0 / S_inv(i);
            } else {
                S_inv(i) = 0.0;
            }
        }
        NormalizedCovarianceMatrixSummed = svd.matrixV() * S_inv.asDiagonal() * svd.matrixU().transpose();
*/
        //std::ofstream file17 (saveDirectory + "sumNormalizedCovarianceMatrix" + fileTag + ".txt");
        //file17 << std::setprecision(32) << NormalizedCovarianceMatrixSummed;
        //file17.close( );
	std::cout<<"Normalized covariance matrix inverted!"<<std::endl;
	std::cout<<"Now unnormalizing it back again..."<<std::endl;
        //Eigen::MatrixXd normalizDiagonals = normalizationFactorsMerged.asDiagonal();
        Eigen::MatrixXd resultedUnnormalizedCovarianceMatrix = Eigen::MatrixXd::Zero(total_size,total_size);
        for( int i = 0; i < normalizationFactorsMerged.rows( ); i++ )
        {
                for( int j = 0; j < normalizationFactorsMerged.rows( ); j++ )
                {
                        resultedUnnormalizedCovarianceMatrix(i,j) = NormalizedCovarianceMatrixSummed( i, j ) /  (normalizationFactorsMerged( i ) *  normalizationFactorsMerged( j ));

                }
        }
	//Eigen::MatrixXd resultedUnnormalizedCovarianceMatrix = normalizDiagonals*NormalizedCovarianceMatrixSummed * normalizDiagonals;
	std::cout<<"Unnormalization done!"<<std::endl;
	std::cout<<"Now saving the matrices..."<<std::endl;

	std::ofstream file13 (saveDirectory + "resultedUnnormalizedCovarianceMatrix" + fileTag + ".txt");
        file13 << std::setprecision(32) << resultedUnnormalizedCovarianceMatrix ;
        file13.close( );
        //std::ofstream file14 (saveDirectory + "normalizationFactors" + fileTag + ".txt");
        //file14 << std::setprecision(32) << normalizationFactorsMerged ;
        //file14.close( );
        //std::ofstream file15 (saveDirectory + "sumCovarianceMatrix" + fileTag + ".txt");
        //file15 << std::setprecision(32) << sumCovarianceMatrix ;
        //file15.close( );
	std::ofstream file16 (saveDirectory + "suminverseNormalizedCovarianceMatrix" + fileTag + ".txt");
        file16 << std::setprecision(32) << P_global;
        file16.close( );

	return;
}

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
void arcLengthRuns( double hoursperday, double initialTime, double finalTime, int arcLength, int iterationNumber, double perturbPos, double perturbVel, int number_of_arcs, int itotalDuration, int startTime, int intperturbPos, int intperturbVel, int ihoursperday, int batchSize, bool performEst, bool filterArcTimes )
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

    // set input options
    double epehemeridesTimeStep = 60.0;
    bool useInterpolatedEphemerides = true;
    double observationsSamplingTime = 60.0;
    double buffer = 30.0 * epehemeridesTimeStep;
    double arcDuration = arcLength*86400.0;//2.0E4;
    std::string dragEst = "per-arc";//"per-rev";
    double ndays = 1.0;
    double hoursperdaydrag = 2.0;

    //double hoursperday = 10.0;
    //int iterationNumber = 6;
    //const double gravitationalParameter = 4.2828378e13;
    //const double planetaryRadius = 3389.5E3;

    double oneWayDopplerNoise = 0.0001;
    double twoWayDopplerNoise = 0.0001;
    double rangeNoise = 1;
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
    //bodySettings.at( spacecraftName )->ephemerisSettings->resetMakeMultiArcEphemeris( true );

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
    cosineAmplitudes[ 1 ] = Eigen::Matrix< double, 3, 5 >::Zero( );
    //scM/10
/*
    cosineAmplitudes[ 1 ]( 0, 0 ) += -1.73871738023640e-10/(365*24*3600);
    cosineAmplitudes[1](0,1) +=1.83526161335626e-10/(365*24*3600);
    cosineAmplitudes[1](0,2) += -2.29416809459079e-10/(365*24*3600);
    cosineAmplitudes[1](1,0) += -1.96888593742833e-10/(365*24*3600);
    cosineAmplitudes[1](1,1) += 3.64572599774489e-10/(365*24*3600);
    cosineAmplitudes[1](1,2) += -1.99682320271527e-10/(365*24*3600);
    cosineAmplitudes[1](2,0) += 4.33626663660650e-10/(365*24*3600);

        cosineAmplitudes[1](2,1) += 5.187166893375834e-11/(365*24*3600);
        cosineAmplitudes[1](2,2) += 4.704102072862383e-10/(365*24*3600);
        cosineAmplitudes[1](2,3) += 1.325699794642290e-10/(365*24*3600);
        cosineAmplitudes[1](2,4) += 1.175466344751454e-10/(365*24*3600);
        */

    //scM/0.01
     cosineAmplitudes[ 1 ]( 0, 0 ) += -1.73871738023640e-12/(365*24*3600);
     cosineAmplitudes[1](0,1) +=1.83526161335626e-12/(365*24*3600);
     cosineAmplitudes[1](0,2) += -2.29416809459079e-12/(365*24*3600);
     cosineAmplitudes[1](1,0) += -1.96888593742833e-12/(365*24*3600);
     cosineAmplitudes[1](1,1) += 3.64572599774489e-12/(365*24*3600);
     cosineAmplitudes[1](1,2) += -1.99682320271527e-12/(365*24*3600);
     cosineAmplitudes[1](1,3) += 1.19185562092185e-11/(365*24*3600);
     cosineAmplitudes[1](2,0) += 4.33626663660650e-12/(365*24*3600);
     cosineAmplitudes[1](2,1) += 5.187166893375834e-13/(365*24*3600);
     cosineAmplitudes[1](2,2) += 4.704102072862383e-12/(365*24*3600);
     cosineAmplitudes[1](2,3) += 1.325699794642290e-12/(365*24*3600);
     cosineAmplitudes[1](2,4) += 1.175466344751454e-12/(365*24*3600);

        /*cosineAmplitudes[ 1 ]( 0, 0 ) += -1.73871738023640e-12/(365*24*3600);
        cosineAmplitudes[1](0,1) +=1.83526161335626e-12/(365*24*3600);
        cosineAmplitudes[1](0,2) += -2.29416809459079e-12/(365*24*3600);
        cosineAmplitudes[1](1,0) += -1.96888593742833e-12/(365*24*3600);
        cosineAmplitudes[1](1,1) += 3.64572599774489e-12/(365*24*3600);
        cosineAmplitudes[1](1,2) += -1.99682320271527e-12/(365*24*3600);
        //cosineAmplitudes[1](1,3) += 1.19185562092185e-11/(365*24*3600);
        cosineAmplitudes[1](2,0) += 4.33626663660650e-12/(365*24*3600);
*/
    std::map<int, Eigen::MatrixXd> sineAmplitudes;
    //sineAmplitudes[ 1 ] = Eigen::Matrix< double, 3, 5 >::Zero( );
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
    // Create integrator times
    std::vector< double > integrationArcStartTimes;
    std::vector< double > integrationArcEndTimes;
    std::vector< double > integrationArcLimits;

    std::vector< double > filteredArcStartTimes;
    std::vector< double > filteredArcEndTimes;

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

    double accumulated_time = 0.0;
    int count = 0;
    double one_month_duration = 30.0*86400.0;


        std::vector<std::vector<double>> observationTimesPerArc;
        std::vector<double> allObservationTimes;
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
            /*
            // create observation times per arc
            double arcDur = currentEndTime - currentStartTime;
            int days = static_cast<int>(arcDur / physical_constants::JULIAN_DAY );
            std::vector<double> initial_times_list = { currentStartTime + buffer};
            std::vector<double> final_times_list = { initial_times_list[0] + hoursperday * 3600 };
            for (int day = 1; day < days  ; ++day) {
                    initial_times_list.push_back(currentStartTime+ buffer + day * 24 * 3600);
                    double finalTime = initial_times_list[day] + hoursperday * 3600;
                    if (finalTime > integrationArcEndTimes.back()) {
                            break;
                    }
                    final_times_list.push_back(finalTime);
            }
            std::vector<double> observationTimesList;
            for (int i = 0; i < final_times_list.size(); ++i) {
                    double start = initial_times_list[i];
                    double end = final_times_list[i];
                    for (double time = start; time < end ; time += 10) {
                            observationTimesList.push_back(time);
                            allObservationTimes.push_back(time);
                    }
            }
            observationTimesPerArc.push_back(observationTimesList);
            */
	//Filtered Arc Times
	if (count<2){
	    filteredArcStartTimes.push_back(currentStartTime);
	    filteredArcEndTimes.push_back(currentEndTime);
	}
	accumulated_time += (currentEndTime-currentStartTime);
	count++;
	if (accumulated_time>= one_month_duration ){
	    count = 0;
	    accumulated_time = 0.0;
	}
	//Store integration arc times
	currentStartTime = currentEndTime - arcOverlap;
        currentEndTime = currentStartTime + arcDuration;
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
       std::make_shared<RungeKuttaFixedStepSizeSettings<> >( 30, CoefficientSets::rungeKutta87DormandPrince );
        std::cout<<"Integration settings created"<<std::endl;
	int numberOfIntegrationArcs;
	if (filterArcTimes){
	    numberOfIntegrationArcs = filteredArcStartTimes.size( );
	}else{
	    numberOfIntegrationArcs = integrationArcStartTimes.size( );
	}

        std::cout<<"number of integration arcs: "<<numberOfIntegrationArcs<<std::endl;
        std::vector< Eigen::VectorXd > systemInitialStates(numberOfIntegrationArcs, Eigen::VectorXd(6));
	std::vector<double> arcTimesToEstimate;
	std::vector<double> arcEndTimesToEstimate;
	if (filterArcTimes){
	    arcTimesToEstimate = filteredArcStartTimes;
	    arcEndTimesToEstimate = filteredArcEndTimes;
	}else{
	    arcTimesToEstimate = integrationArcStartTimes;
	    arcEndTimesToEstimate = integrationArcEndTimes;
	}


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

        std::vector< Eigen::MatrixXd > covarianceMatrices;
        //std::vector<Eigen::MatrixXd> inverseUnnormalizedCovarianceMatrices;
        std::vector<Eigen::VectorXd> weightDiagonals;
        std::vector<Eigen::MatrixXd> unnormalizedDesignMatrices;
	std::vector<Eigen::MatrixXd> normalizedDesignMatrices;
        int numberOfGlobalParameters;
	std::vector<Eigen::MatrixXd> inverseNormalizedCovarianceMatrices;
        std::vector<Eigen::VectorXd> normalizationFactors;
        std::vector<Eigen::MatrixXd> P0_matrices;
        // create multi arc propagation settings
        std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > arcPropagationSettingsList;
        for( unsigned int i = 0; i < numberOfIntegrationArcs; i++ )
        {
                std::cout<<"iteration started"<<std::endl;
                std::cout<<i<<std::endl;
                std::cout<<bodiesToIntegrate[ 0 ]<<std::endl;
                std::vector<double> filteredObservationTimes;
		std::ofstream feo(saveDirectory + "observationTimes_arc_" + std::to_string(i) +  fileTag + ".txt");
                for (const auto& time : observationTimesList) {
                    if (time >= arcTimesToEstimate.at(i) + buffer && time <= arcEndTimesToEstimate.at(i)-buffer) {
                        filteredObservationTimes.push_back(time);
			feo <<  std::setprecision(32) <<time <<"\n";
                    }
                }
                feo.close();
		observationTimesPerArc.push_back(filteredObservationTimes);
                //std::cout<<integrationArcStartTimes.at(i)<<std::endl;
                systemInitialStates[ i ]  = spice_interface::getBodyCartesianStateAtEpoch(
                        bodiesToIntegrate[ 0 ], "Mars", "MARSIAU", "NONE", arcTimesToEstimate.at(i));
                std::cout<<systemInitialStates[i]<<std::endl;
                std::cout<<"system initial states created"<<std::endl;
                // Create termination settings
                std::shared_ptr< PropagationTerminationSettings > terminationSettings = propagationTimeTerminationSettings(
                        arcEndTimesToEstimate.at(i) );

                std::shared_ptr< TranslationalStatePropagatorSettings< double, double> > propagatorSettings = translationalStatePropagatorSettings< double, double >( centralBodies, accelerationModelMap, bodiesToIntegrate,
                                                                                                                                                      systemInitialStates[i], arcTimesToEstimate.at(i), integratorSettings, terminationSettings, cowell, dependentVariablesToSave);

                SingleArcDynamicsSimulator< > dynamicsSimulator(
                        bodies, propagatorSettings );

                std::map< double, Eigen::VectorXd > stateHistory = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
                writeDataMapToTextFile( stateHistory, "stateHistoryPropagation_arc_" + std::to_string(i) + fileTag + ".txt", saveDirectory,
                                "", 18, 18 );

                std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
    	        getInitialStateParameterSettings< double, double  >( propagatorSettings, bodies);
                parameterNames.push_back(std::make_shared< EstimatableParameterSettings >(spacecraftName,constant_drag_coefficient));// initial_times_list_drag ));
                //parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                //                                         2, 0, 4, 0, "Mars", spherical_harmonics_cosine_coefficient_block ) );

                parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                         2, 0, 18, 18, "Mars", spherical_harmonics_cosine_coefficient_block ) );
                parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                              2, 1, 18, 18, "Mars", spherical_harmonics_sine_coefficient_block ) );

                std::map<int, std::vector<std::pair<int, int> > > cosineBlockIndicesPerPeriod;
                //periodic gravity field
                cosineBlockIndicesPerPeriod[ 0 ].push_back( std::make_pair( 2, 0) );
                //cosineBlockIndicesPerPeriod[ 0 ].push_back( std::make_pair( 2, 1 ) );
                cosineBlockIndicesPerPeriod[ 0 ].push_back( std::make_pair( 3, 0 ) );
                cosineBlockIndicesPerPeriod[ 0 ].push_back( std::make_pair( 4, 0 ) );
                cosineBlockIndicesPerPeriod[ 0 ].push_back( std::make_pair( 5, 0 ) );

	        cosineBlockIndicesPerPeriod[ 1 ].push_back( std::make_pair( 2, 0) );
                //cosineBlockIndicesPerPeriod[ 0 ].push_back( std::make_pair( 2, 1 ) );
                cosineBlockIndicesPerPeriod[ 1 ].push_back( std::make_pair( 3, 0 ) );
                cosineBlockIndicesPerPeriod[ 1 ].push_back( std::make_pair( 4, 0 ) );
                cosineBlockIndicesPerPeriod[ 1 ].push_back( std::make_pair( 5, 0 ) );

	        cosineBlockIndicesPerPeriod[ 2 ].push_back( std::make_pair( 2, 0) );
                //cosineBlockIndicesPerPeriod[ 0 ].push_back( std::make_pair( 2, 1 ) );
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
                 /*
                //std::map<int, std::vector<std::pair<int, int> > > cosineBlockIndicesPerPower;
                cosineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 2, 0 ) );
                cosineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 2, 1 ) );
                cosineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 2, 2 ) );
                cosineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 3, 0 ) );
                cosineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 4, 0 ) );
                //cosineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 4, 0 ) );
*/
                std::map<int, std::vector<std::pair<int, int> > > sineBlockIndicesPerPower;

                //sineBlockIndicesPerPower[ 1 ].push_back( std::make_pair( 2, 1 ) );
                parameterNames.push_back( std::make_shared< PolynomialGravityFieldVariationEstimatableParameterSettings >(
                        "Mars", cosineBlockIndicesPerPower, sineBlockIndicesPerPower ) );


                std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
                            createParametersToEstimate< double, double >( parameterNames, bodies );

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
                        //std::function< double( const double ) > noiseFunction = noiseFunctions[currentObservable];
                        for( unsigned int currLinkEnd = 0; currLinkEnd < currentLinkEndsList.size( ); currLinkEnd++ )
                        {
                                measurementSimulationInput.push_back(
                                        std::make_shared< TabulatedObservationSimulationSettings< > >(
                                                currentObservable, currentLinkEndsList[ currLinkEnd ], observationTimesPerArc[i], receiver, observationViabilitySettings, noiseFunctions[currentObservable]) );
                        }
                }

                // Simulate observations.
                std::shared_ptr< ObservationCollection< > > observationsAndTimes = simulateObservations< double, double >(
                        measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );
                std::cout<<"observations and times created"<<std::endl;

                Eigen::Matrix< double, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;
                int numberOfParameters = initialParameterEstimate.rows( );
                numberOfGlobalParameters = numberOfParameters - 7;
                std::cout<<"number of parameters: "<<numberOfParameters<<std::endl;
                printEstimatableParameterEntries( parametersToEstimate );

                // Create a 2D vector (matrix) filled with zeros
	        // since it is one arc together with the one drag coefficient, the first diagonal matix is 7x7
	        const int DIAGONALS = 7;

                Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(numberOfParameters, numberOfParameters);

                double aprioriuncertainty =1.0/(10*10);
                matrix(DIAGONALS-1,DIAGONALS-1) = aprioriuncertainty;


                // Print the extracted errors
                // Fill the matrix with the values of the diagonal
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

                matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+24, DIAGONALS+cnmErrors.size()+snmErrors.size()+24) = 0.0;///(0.1E-19*0.1E-19);
                matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+25, DIAGONALS+cnmErrors.size()+snmErrors.size()+25) = 0.0;///(0.1E-19*0.1E-19);
                matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+26, DIAGONALS+cnmErrors.size()+snmErrors.size()+26) = 0.0;///(0.1E-19*0.1E-19);
                matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+27, DIAGONALS+cnmErrors.size()+snmErrors.size()+27) = 0.0;///(0.1E-19*0.1E-19);
	        matrix(DIAGONALS+cnmErrors.size()+snmErrors.size()+28, DIAGONALS+cnmErrors.size()+snmErrors.size()+28) = 0.0;///(0.1E-19*0.1E-19);
                P0_matrices.push_back(matrix);
                //std::ofstream fe(saveDirectory + "unnormalizedAprioriCovInv_arc_" + std::to_string(i) +  fileTag + ".txt");
                //fe <<  std::setprecision(32) <<matrix;
                //fe.close();

		std::shared_ptr< CovarianceAnalysisInput< double, double > > covarianceInput =
                        std::make_shared< CovarianceAnalysisInput< double, double > >(
                        observationsAndTimes,matrix );
                std::cout<<"covariance input created"<<std::endl;
                std::map< observation_models::ObservableType, double > weightPerObservable;
                //weightPerObservable[ one_way_doppler ] = std::pow(oneWayDopplerNoise, -2);
                // weightPerObservable[ one_way_range ] = std::pow(rangeNoise, -2);
                weightPerObservable[ two_way_doppler ] = std::pow(twoWayDopplerNoise, -2);

                covarianceInput->setConstantPerObservableWeightsMatrix( weightPerObservable );


                std::shared_ptr< CovarianceAnalysisOutput< double, double > > covarianceOutput = orbitDeterminationManager.computeCovariance(
                            covarianceInput );
                std::cout<<"covariance output created"<<std::endl;

		Eigen::MatrixXd unnormalizedDesignMatrix = covarianceOutput->getUnnormalizedDesignMatrix( );
                //std::ofstream fe(saveDirectory + "unnormalizedDesignMatrix_arc_" + std::to_string(i) +  fileTag + ".txt");
                //fe <<  std::setprecision(32) <<unnormalizedDesignMatrix;
                //fe.close();
		unnormalizedDesignMatrices.push_back(unnormalizedDesignMatrix);

		Eigen::MatrixXd normalizedDesignMatrix = covarianceOutput->getNormalizedDesignMatrix( );
                //std::ofstream fe2(saveDirectory + "normalizedDesignMatrix" + std::to_string(i)  + fileTag + ".txt");
                // Write the matrix to the file
                //fe2 << std::setprecision(32) <<normalizedDesignMatrix;
                // Close the file
                //fe2.close();
                normalizedDesignMatrices.push_back(normalizedDesignMatrix);

                Eigen::VectorXd weightMatrixDiagonal = covarianceOutput->weightsMatrixDiagonal_;
                weightDiagonals.push_back(weightMatrixDiagonal);
                //std::ofstream fe3(saveDirectory + "weightsMatrixDiagonal_arc_" + std::to_string(i) +fileTag + ".txt");
                // Write the matrix to the file
                //fe3 <<  std::setprecision(32) << weightMatrixDiagonal;
                // Close the file
                //fe3.close();
                // Eigen::MatrixXd inverseUnnormalizedCovarianceMatrix = covarianceOutput->getUnnormalizedInverseCovarianceMatrix( );
                //inverseUnnormalizedCovarianceMatrices.push_back(inverseUnnormalizedCovarianceMatrix);
                //std::ofstream fein(saveDirectory + "inverseUnormalizedCovarianceMatrix_arc_" + std::to_string(i) + fileTag + ".txt");
                //fein <<  std::setprecision(32) <<inverseUnnormalizedCovarianceMatrix;
                //fein.close();

                Eigen::MatrixXd inverseNormalizedCovarianceMatrix = covarianceOutput->getNormalizedInverseCovarianceMatrix( );
                inverseNormalizedCovarianceMatrices.push_back(inverseNormalizedCovarianceMatrix);
                //std::ofstream fean(saveDirectory + "inverseNormalizedCovarianceMatrix_arc_" + std::to_string(i) + fileTag + ".txt");
                //fean <<  std::setprecision(32) <<inverseNormalizedCovarianceMatrix;
                //fean.close();

                ////std::ofstream fe4(saveDirectory + "normalizedCovarianceMatrix" + fileTag + ".txt");
                // Write the matrix to the file
                ////fe4 << covarianceOutput->normalizedCovarianceMatrix_;
                // Close the file
                ////fe4.close();
                Eigen::VectorXd normalizationFactor = covarianceOutput->designMatrixTransformationDiagonal_;
                normalizationFactors.push_back(normalizationFactor);
                //std::ofstream fe5(saveDirectory + "designMatrixTransformationDiagonal_arc_" + std::to_string(i)  + fileTag + ".txt");
                //fe5 <<  std::setprecision(32) <<normalizationFactor;
                //fe5.close();

                Eigen::MatrixXd correlationMatrix = covarianceOutput->getCorrelationMatrix( );
                //std::ofstream file10 (saveDirectory + "correlationMatrix_arc_" + std::to_string(i) + fileTag + ".txt");
                //file10 << std::setprecision(32) << correlationMatrix ;
                //file10.close( );

                Eigen::MatrixXd covarianceMatrix = covarianceOutput->getUnnormalizedCovarianceMatrix( );
                covarianceMatrices.push_back(covarianceMatrix);
                //std::ofstream file11 (saveDirectory + "conCovarianceMatrix_arc_" + std::to_string(i) + fileTag + ".txt");
                //file11 << std::setprecision(32) << covarianceMatrix ;
                //file11.close( );

                Eigen::Matrix<double, Eigen::Dynamic, 1> FormalError = covarianceOutput->getFormalErrorVector( );
                std::cout<<"formal error: "<<covarianceOutput->getFormalErrorVector( ).transpose( )<<std::endl;
                //std::ofstream file12 (saveDirectory + "FormalError_arc_" + std::to_string(i) + fileTag + ".txt");
                //file12 << std::setprecision(32) << FormalError ;
                //file12.close( );


		if ((i + 1) % batchSize ==0 || i == numberOfIntegrationArcs -1){
		    computeCovarianceMatrix(startTime,i+1,numberOfGlobalParameters,normalizedDesignMatrices,unnormalizedDesignMatrices,normalizationFactors,weightDiagonals,P0_matrices);
		}
        }


        //std::cout<<"single arc propagation done"<<std::endl;

/*
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

*/


// Create the noise functions that return Eigen::VectorXd

    // noiseFunctions[ one_way_range ] =
    //         [=](const double input) -> Eigen::VectorXd {
    //             // Call the original function that returns a double
    //             double noiseValue = utilities::evaluateFunctionWithoutInputArgumentDependency< double, const double >(
    //                     createBoostContinuousRandomVariableGeneratorFunction(
    //                             tudat::statistics::normal_boost_distribution, { 0.0, rangeNoise }, 0.0
    //                     ), input
    //             );
    //             // Convert the double to Eigen::VectorXd
    //             Eigen::VectorXd result(1);
    //             result(0) = noiseValue;
    //             return result;
    //         };






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
/*
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
*/



}
int main() {
        int iterationNumber = 5;
	int batchSize = 210;
        std::vector<int> arcLengths= {5};
        std::vector<int> number_of_arcs= {210};//72,137,210};//{210,283,356};
        std::vector<double> hoursperday = {10.0};
        std::vector<int> ihoursperday = {10};
        //std::vector<double> initialTimes = {-240.0*86400.0, -180.0*86400.0, -60.0*86400.0,0.0, 60*86400.0, 180.0*86400.0, 240.0*86400.0};
        std::vector<double> initialTimes = {0.0*86400.0};
        //std::vector<int> startTime = {-240, -180, -60, 0, 60, 180, 240};
        std::vector<int> startTime = {0};
        std::vector<double> perturbPos = {100.0};
        std::vector<int> intperturbPos = {100};
        std::vector<double> perturbVel = {0.001};
        std::vector<int> intperturbVel = {0001};
        std::vector<int> totalDuration = {1050};//360};//,685,1050}; //{1050,1415,1780};
        std::vector<double> finalTimes= {86400.0*1050.0};//,86400.0*360.0};//,86400.0*685.0,86400.0*1050.0}; //{86400.0*1050.0,86400.0*1415.0,86400.0*1780.0};
        bool performEst = false;
	bool filterArcTimes = false;
        for (int i = 0; i<arcLengths.size();i++) {
                for (int hours = 0; hours<hoursperday.size();hours++) {
                        for (int initialTime = 0; initialTime<finalTimes.size(); initialTime++) {
                                double finalTime = initialTimes[initialTime] + finalTimes.at(initialTime);
                                for (int intperturb = 0; intperturb<perturbPos.size(); intperturb++) {
                                        arcLengthRuns( hoursperday[hours],  initialTimes[0],  finalTimes.at(initialTime), arcLengths[i],  iterationNumber, perturbPos[intperturb], perturbVel[intperturb], number_of_arcs[i],  totalDuration.at(initialTime),  startTime[0],  intperturbPos[intperturb],  intperturbVel[intperturb], ihoursperday[hours],batchSize, performEst,filterArcTimes);

                                }

                        }
                }
        }

}
