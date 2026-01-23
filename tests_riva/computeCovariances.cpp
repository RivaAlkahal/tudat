//
// Created by ralkahal on 7-8-25.
//
#include <iostream>
#include <fstream>
#include <limits>
#include "fstream"
#include "iostream"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/Core>
#ifdef _OPENMP
# include <omp.h>
#endif

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
template <class Sparse>
	double minDiagonal(const Sparse& A) {
	double m = std::numeric_limits<double>::infinity();
	for (int k = 0; k < A.outerSize(); ++k) {
		for (typename Sparse::InnerIterator it(A, k); it; ++it) {
			if (it.row() == it.col()) m = std::min(m, (double)it.value());
		}
	}
	return (std::isinf(m) ? 0.0 : m);
}
Eigen::MatrixXd loadMatrixBinary(const std::string& filename)
{
    std::ifstream in(filename,std::ios::binary);
    if(!in) throw std::runtime_error("Cannot open file:" + filename);
    Eigen::Index rows, cols;
    in.read(reinterpret_cast<char*>(&rows), sizeof(Eigen::Index));
    in.read(reinterpret_cast<char*>(&cols), sizeof(Eigen::Index));
    Eigen::MatrixXd matrix(rows,cols);
    in.read(reinterpret_cast<char*>(matrix.data()),rows*cols*sizeof(double));
    if (!in) throw std::runtime_error("Cannot read data from :" + filename);
    return matrix;
}


Eigen::VectorXd loadVectorBinary(const std::string& filename)
{
    std::ifstream in(filename,std::ios::binary);
    if(!in) throw std::runtime_error("Cannot open file:" + filename);
    Eigen::Index size;
    in.read(reinterpret_cast<char*>(&size), sizeof(Eigen::Index));
    Eigen::VectorXd vector(size);
    in.read(reinterpret_cast<char*>(vector.data()),size*sizeof(double));
    if (!in) throw std::runtime_error("Cannot read data from :" + filename);
    return vector;
}

void saveMatrixBinary(const std::string& filename, const Eigen::MatrixXd& M)
{
    std::ofstream out(filename,std::ios::binary);
    if(!out) throw std::runtime_error("Cannot open file:" + filename);

    Eigen::Index rows= M.rows();
    Eigen::Index cols = M.cols();
    out.write(reinterpret_cast<const char*>(&rows), sizeof(Eigen::Index));
    out.write(reinterpret_cast<const char*>(&cols), sizeof(Eigen::Index));
//    Eigen::MatrixXd matrix(rows,cols);
    out.write(reinterpret_cast<const char*>(M.data()),rows*cols*sizeof(double));
    //if (!out) throw std::runtime_error("Cannot read data from :" + filename);

}
bool fileExists(const std::string& path){
    std::ifstream f(path);
    return static_cast<bool>(f);
}

inline int loadArcResultsRange(const std::string& baseDir,
				int numArcs,
				std::vector<Eigen::MatrixXd>* normalizedDesignMatrices,
				std::vector<Eigen::MatrixXd>* unnormalizedDesignMatrices,
				std::vector<Eigen::VectorXd>* normalizationFactors,
				std::vector<Eigen::VectorXd>* weightDiagonals,
				std::vector<Eigen::MatrixXd>* P0_matrices,
				const std::string& normalizedDesignName = "normalized_Design_Matrix.bin",
				const std::string& unnormalizedDesignName = "unnormalized_Design_Matrix.bin",
				const std::string& normalizationFactorsName = "normalization_factor.bin",
				const std::string& weightDiagonalsName = "weight_matrix_diagonal.bin",
				const std::string& P0Name = "P0_matrix.bin"){

    int used = 0;

    for (int i = 0; i<numArcs; ++i){
	const std::string arcDir = baseDir + "/output_arc_" + std::to_string(i);
	bool ok = false;
	const std::string path = arcDir + "/" + normalizedDesignName;
	if (fileExists(path))
	{
	    try{normalizedDesignMatrices->push_back(loadMatrixBinary(path)); ok =true;}
	    catch (const std::exception& e) { std::cerr <<"Skipped " << path << " : " << e.what() << "\n";}
	}
	else { std::cerr <<"Missing normalized design matrix " << path << "\n";}

	const std::string path1 = arcDir + "/" + unnormalizedDesignName;
	if (fileExists(path1))
	{
	    try{unnormalizedDesignMatrices->push_back(loadMatrixBinary(path1)); ok =true;}
	    catch (const std::exception& e) { std::cerr <<"Skipped " << path1 << " : " << e.what() << "\n";}
	}
	else { std::cerr <<"Missing unnormalized design matrix " << path1 << "\n";}

	const std::string path2 = arcDir + "/" + normalizationFactorsName;
	if (fileExists(path2))
	{
	    try{normalizationFactors->push_back(loadVectorBinary(path2)); ok =true;}
	    catch (const std::exception& e) { std::cerr <<"Skipped " << path2 << " : " << e.what() << "\n";}
	}
	else { std::cerr <<"Missing normalizationFactors " << path2 << "\n";}

	const std::string path3 = arcDir + "/" + weightDiagonalsName;
	if (fileExists(path))
	{
	    try{weightDiagonals->push_back(loadVectorBinary(path3)); ok =true;}
	    catch (const std::exception& e) { std::cerr <<"Skipped " << path3 << " : " << e.what() << "\n";}
	}
	else { std::cerr <<"Missing weightDiagonals " << path3 << "\n";}

	const std::string path4 = arcDir + "/" + P0Name;
	if (fileExists(path4))
	{
	    try{P0_matrices->push_back(loadMatrixBinary(path4)); ok =true;}
	    catch (const std::exception& e) { std::cerr <<"Skipped " << path4 << " : " << e.what() << "\n";}
	}
	else { std::cerr <<"Missing P0  matrix " << path4 << "\n";}

	if (ok) ++used;
    }
    return used;
}

struct combinedResults {
    Eigen::MatrixXd resultedUnnormalizedCovarianceMatrix;
    Eigen::MatrixXd P_global;
};

combinedResults computeCovarianceMatrix( int nIndex, std::vector<Eigen::MatrixXd> normalizedDesignMatrices,std::vector<Eigen::MatrixXd> unnormalizedDesignMatrices,std::vector<Eigen::VectorXd> normalizationFactors, std::vector<Eigen::VectorXd> weightDiagonals,std::vector<Eigen::MatrixXd> P0_matrices) 
{

        Eigen::initParallel();
        Eigen::setNbThreads(64);

	int itotalDuration = nIndex * 5;
	double totalDur = nIndex * 5.0;
	double finalTime = totalDur * 86400.0;
	//std::string fileTag = "accumul_InverseAprALLGlobalPars_drag_5" + std::to_string(itotalDuration);
	//fileTag = nIndex + "nArcs" + fileTag;
	std::cout<<"covariance analysis for" << nIndex <<" arcs done, starting summing them up..."<<std::endl;

	int numberOfParameters = P0_matrices[0].rows();
	std::cout<<P0_matrices.size()<<std::endl;
	std::cout<<P0_matrices[0].size()<<std::endl;
	std::cout<<normalizedDesignMatrices[0].size()<<std::endl;
	std::cout<<unnormalizedDesignMatrices[0].size()<<std::endl;
	std::cout<<normalizationFactors[0].size()<<std::endl;
        int stateVectorSize =6;
        int nLocal=1;
        int numberOfLocalParameters = stateVectorSize + nLocal;
	int numberOfGlobalParameters = numberOfParameters - numberOfLocalParameters;
	std::cout<< "number of global parameters: "<< numberOfGlobalParameters<<std::endl;
	int numCols = numberOfGlobalParameters;
        int nArcs = normalizedDesignMatrices.size();
	if (nArcs != nIndex) throw std::runtime_error("number of arcs do not match!");
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
	// Eigen::JacobiSVD<Eigen::MatrixXd> svd(P_global, Eigen::ComputeThinU | Eigen::ComputeThinV);
 //        double tol = 1e-13;
	// double cond = svd.singularValues()(0)/svd.singularValues().tail(1)(0);
	// std::cout << "condition number: " << cond <<std::endl;

        std::cout<<"Now inverting the inverse covariance matrices..."<<std::endl;

	std::cerr << "min diagonal: " << minDiagonal(P_global) << " (should be > 0 for SPD)\n";
        //IT USED TO START HERE!!
	/*using SpMat = Eigen::SparseMatrix<double>;
        using Solver = Eigen::SparseLU<SpMat>;

        SpMat A = P_global.sparseView();
        Solver solver;
        solver.analyzePattern(A);
        solver.factorize(A);
        if (solver.info() != Eigen::Success) throw std::runtime_error("Error during factorization of the covariance matrix.");
        const int n = A.rows();
        Eigen::MatrixXd invA(n, n);
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n, n);

        int BlockSize = 8; // Adjust this value based on your system's memory capacity

        //#pragma omp parallel for
        /*for (int col = 0; col<n; ++col)
        {
            Eigen::VectorXd e = I.col(col);
            Eigen::VectorXd x = solver.solve(e);
            if (solver.info() != Eigen::Success) {
                std::cerr << "Error during solving for column " << col << " of the inverse covariance matrix." << std::endl;
                return;
            }
            invA.col(col) = x;
        } //INSERT COMMENT CLOSING HERE
        // Using a block-wise approach to handle large matrices
        for (int col = 0; col < n; col += BlockSize) {
                int bs = std::min(BlockSize, n - col);
                // Process a block of columns
                invA.middleCols(col,bs) = solver.solve(I.middleCols(col, bs));
            /*int blockEnd = std::min(col + BlockSize, n);
            Eigen::MatrixXd block = I.block(0, col, n, blockEnd - col);
            Eigen::MatrixXd x = solver.solve(block);
            if (solver.info() != Eigen::Success) {
                std::cerr << "Error during solving for block starting at column " << col << " of the inverse covariance matrix." << std::endl;
                return;
            }
            invA.block(0, col, n, blockEnd - col) = x; //Insert closing the comment here
        }
        Eigen::MatrixXd NormalizedCovarianceMatrixSummed = invA;
*/
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

	//std::ofstream file13 (saveDirectory + "resultedUnnormalizedCovarianceMatrix" + fileTag + ".txt");
        //file13 << std::setprecision(32) << resultedUnnormalizedCovarianceMatrix ;
        //file13.close( );
        //std::ofstream file14 (saveDirectory + "normalizationFactors" + fileTag + ".txt");
        //file14 << std::setprecision(32) << normalizationFactorsMerged ;
        //file14.close( );
        //std::ofstream file15 (saveDirectory + "sumCovarianceMatrix" + fileTag + ".txt");
        //file15 << std::setprecision(32) << sumCovarianceMatrix ;
        //file15.close( );
	//std::ofstream file16 (saveDirectory + "suminverseNormalizedCovarianceMatrix" + fileTag + ".txt");
        //file16 << std::setprecision(32) << P_global;
        //file16.close( );
	resultedUnnormalizedCovarianceMatrix= 0.5 * (resultedUnnormalizedCovarianceMatrix + resultedUnnormalizedCovarianceMatrix.transpose());
	return combinedResults{resultedUnnormalizedCovarianceMatrix, P_global};
}

static void printUsage(const char* exe)
{
    std::cerr
	<< "Usage:\n"
	<< "  " << exe << "	 [--baseDir ouput --num-arcs N] \n";
}

int main(int argc, char* argv[])
{
    //std::string outputFile;
    std::string baseDir;
    int numArcs = -1;
    for (int i=1; i<argc;++i){
	std::string a = argv[i];
	auto next = [&](std::string& dst){ if (++i>=argc) throw std::runtime_error("Missing value after "+a); dst = argv[i]; };

	if (a== "--baseDir") next(baseDir);
	else if (a== "--num-arcs") { std::string t; next(t); numArcs = std::atoi(t.c_str());}
	else { std::cerr << "Unknown arg: " << a << "\n"; printUsage(argv[0]); return 1; }
    }
    if (baseDir.empty() || numArcs <= 0)
    { std:: cerr << "Provide --base-dir AND --num-arcs\n"; printUsage(argv[0]); return 1;}
    //if (outputFile.empty()) { std::cerr << "Missing --ouput\n"; printUsage(argv[0]); return 1;}

    std::vector<Eigen::MatrixXd> normalizedDesignMatrices;
    std::vector<Eigen::MatrixXd> unnormalizedDesignMatrices;
    std::vector<Eigen::VectorXd> normalizationFactors;
    std::vector<Eigen::VectorXd> weightDiagonals;
    std::vector<Eigen::MatrixXd> P0_matrices;
    int itotalDuration = numArcs * 5;
    std::string fileTag = "accumul_InverseAprALLGlobalPars_drag_5" + std::to_string(itotalDuration);
    int used = loadArcResultsRange(baseDir,numArcs, &normalizedDesignMatrices,&unnormalizedDesignMatrices,&normalizationFactors,&weightDiagonals,&P0_matrices);
    if (used ==0) throw std::runtime_error("No matrices are loaded!");

    Eigen::MatrixXd resultedUnnormalizedCovarianceMatrix, P_global;

    combinedResults res = computeCovarianceMatrix( numArcs, normalizedDesignMatrices, unnormalizedDesignMatrices, normalizationFactors, weightDiagonals,P0_matrices);

    saveMatrixBinary(baseDir + "/resultedUnnormalizedCovarianceMatrix_COD" + fileTag + ".bin",res.resultedUnnormalizedCovarianceMatrix);
    //std::ofstream file1("resultedUnnormalizedCovarianceMatrix_"+ fileTag + ".txt");
    //file1<< std::setprecision(32) << res.resultedUnnormalizedCovarianceMatrix;
    //file1.close();
    saveMatrixBinary(baseDir + "/suminverseNormalizedCovarianceMatrixCOD" + fileTag + ".bin",res.P_global);


}