// Created by ralkahal on 2-9-25.
//
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

struct EigenThreadsGuard {
	int prev_{1};
	explicit EigenThreadsGuard(int set_to = 1) {
		prev_ = Eigen::nbThreads();   // save current
		Eigen::setNbThreads(set_to);  // set temporary value
	}
	~EigenThreadsGuard() {
		Eigen::setNbThreads(prev_);   // restore
	}
};

#ifdef _OPENMP
#include <omp.h>
struct OmpThreadsGuard {
	int prev_{1};
	explicit OmpThreadsGuard(int set_to = -1) {
		prev_ = omp_get_max_threads();
		if (set_to > 0) omp_set_num_threads(set_to);
	}
	~OmpThreadsGuard() {
		
	}
};
#endif

#include <Eigen/Sparse>
#include <numeric>
#include <vector>
#include <stdexcept>

// Toggle verbose logs
#ifndef COV_VERBOSE
#define COV_VERBOSE 0
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

// --- Debug switch ---
#ifndef COV_DEBUG
#define COV_DEBUG 1   // set to 0 to disable all debug prints/checks
#endif

#if COV_DEBUG
#include <limits>
#include <iostream>

template <class M>
bool allFinite(const M& A) { return A.array().isFinite().all(); }

template <class V>
bool allNonNegative(const V& v) { return (v.array() >= 0.0).all(); }

template <class Sparse>
double symmetryResidual(const Sparse& A) {
	// ||A - A^T|| / ||A||
	auto AT = Sparse(A.transpose());
	Sparse D = A - AT;
	double num = D.norm();
	double den = A.norm();
	return (den > 0.0) ? num / den : num;
}

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

template <class Dense>
double maxAbsDiff(const Dense& A, const Dense& B) {
	return (A - B).cwiseAbs().maxCoeff();
}
#endif

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
inline bool loadParamCounts(const std::string& arcDir,
							int& numberOfLocalParameters,
							int& numberOfGlobalParameters)
{
	std::ifstream f(arcDir + "/param_counts.txt");
	if (!f) return false;
	if (!(f >> numberOfLocalParameters)) return false;
	if (!(f >> numberOfGlobalParameters)) numberOfGlobalParameters = 0; // tolerate 1-line files
	return true;
}




inline int loadArcResultsRange(const std::string& baseDir,
				int start,
				int step,
				int take,
				int skip,
				bool hasStep,
				bool hasWindow,
				int numArcs,
				std::vector<double>* numberOfLocalParametersAll,
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
    auto process_arc=[&](int i){

	const std::string arcDir = baseDir + "/output_arc_" + std::to_string(i);
    	int numberOfLocalParameters;
    	int nGlobal = 0;
	if (!loadParamCounts(arcDir, numberOfLocalParameters, nGlobal)) {
		// decide: skip arc, or fallback to inference, or throw
		throw std::runtime_error("Missing/invalid param_counts.txt in " + baseDir);
	}
	numberOfLocalParametersAll->push_back(numberOfLocalParameters);
	bool ok = false;
	const std::string path = arcDir + "/" + normalizedDesignName;
	if (fileExists(path))
	{
	    try{normalizedDesignMatrices->push_back(loadMatrixBinary(path));}
	    catch (const std::exception& e) { std::cerr <<"Skipped " << path << " : " << e.what() << "\n";}
	}
	else { std::cerr <<"Missing normalized design matrix " << path << "\n";}

	const std::string path1 = arcDir + "/" + unnormalizedDesignName;
	if (fileExists(path1))
	{
	    try{unnormalizedDesignMatrices->push_back(loadMatrixBinary(path1));}
	    catch (const std::exception& e) { std::cerr <<"Skipped " << path1 << " : " << e.what() << "\n";}
	}
	else { std::cerr <<"Missing unnormalized design matrix " << path1 << "\n";}

	const std::string path2 = arcDir + "/" + normalizationFactorsName;
	if (fileExists(path2))
	{
	    try{normalizationFactors->push_back(loadVectorBinary(path2));}
	    catch (const std::exception& e) { std::cerr <<"Skipped " << path2 << " : " << e.what() << "\n";}
	}
	else { std::cerr <<"Missing normalizationFactors " << path2 << "\n";}

	const std::string path3 = arcDir + "/" + weightDiagonalsName;
	if (fileExists(path))
	{
	    try{weightDiagonals->push_back(loadVectorBinary(path3));}
	    catch (const std::exception& e) { std::cerr <<"Skipped " << path3 << " : " << e.what() << "\n";}
	}
	else { std::cerr <<"Missing weightDiagonals " << path3 << "\n";}

	const std::string path4 = arcDir + "/" + P0Name;
	if (fileExists(path4))
	{
	    try{P0_matrices->push_back(loadMatrixBinary(path4));}
	    catch (const std::exception& e) { std::cerr <<"Skipped " << path4 << " : " << e.what() << "\n";}
	}
	else { std::cerr <<"Missing P0  matrix " << path4 << "\n";}

    };

    if (hasStep){
	for (int i = start; i< numArcs; i+=step){
	    process_arc(i);
	}
    } else if (hasWindow){
	for (int base =start; base <numArcs; base += (take+skip)){
	    for (int j=0;j<take && base + j <numArcs; ++j){
	    process_arc(base+j);
	    }
	}
    } else {
	throw std::runtime_error("invalid selection mode, choose either step, or take and skip");
    }

    return 1;
}

struct combinedResults {
    Eigen::MatrixXd resultedUnnormalizedCovarianceMatrix;
    Eigen::SparseMatrix<double> P_global_sparse;
};


// Helper: append dense block into triplet list at (r0,c0)
static inline void addBlockTriplets(int r0, int c0,
	const Eigen::Ref<const Eigen::MatrixXd>& B,
	std::vector<Eigen::Triplet<double>>& trips)
{
	const int R = B.rows(), C = B.cols();
	for (int r = 0; r < R; ++r) {
		const double* brow = B.data() + r;         // column-major
		for (int c = 0; c < C; ++c) {
			double v = *(brow + c*R);
			if (v != 0.0) trips.emplace_back(r0 + r, c0 + c, v);
		}
	}
}


combinedResults computeCovarianceMatrix(std::vector<double> numberOfLocalParametersAll, int nIndex, std::vector<Eigen::MatrixXd> normalizedDesignMatrices,std::vector<Eigen::MatrixXd> unnormalizedDesignMatrices,std::vector<Eigen::VectorXd> normalizationFactors, std::vector<Eigen::VectorXd> weightDiagonals,std::vector<Eigen::MatrixXd> P0_matrices)
{

	const int nArcs = static_cast<int>(normalizedDesignMatrices.size());
	std::cout<<"nArcs " << nArcs << ", nIndex "<< nIndex << std::endl;

	int itotalDuration = nIndex * 5;
	double totalDur = nIndex * 5.0;
	double finalTime = totalDur * 86400.0;
	std::cout<<"covariance analysis for" << nArcs <<" samples out of " << nIndex << " arcs done, starting summing them up..."<<std::endl;

	int numberOfParameters = static_cast<int>(P0_matrices[0].rows());
	int stateVectorSize =6;
	
	int numberOfGlobalParameters = numberOfParameters - numberOfLocalParametersAll[0];
	int numCols = numberOfGlobalParameters;

	int totalNumberOfLocalParameters = std::accumulate(numberOfLocalParametersAll.begin(), numberOfLocalParametersAll.end(), 0);
	// total parameter size across all arcs (locals) + globals
	int total_size = totalNumberOfLocalParameters + numberOfGlobalParameters;

	//Eigen::MatrixXd P_global = Eigen::MatrixXd::Zero(total_size,total_size);

	// --- Merge normalization factors for all params (locals for each arc + globals) ---
        // merge normalization factors for all matrices, except for the global parameters
    Eigen::VectorXd normalizationFactorsMerged = Eigen::VectorXd::Zero(total_size);
	size_t stateOffst= 0;
	size_t localOffst = static_cast<size_t>(nArcs) * stateVectorSize;
	for (int i = 0; i < nArcs; i++){
		int nLocal =numberOfLocalParametersAll[i] - stateVectorSize ;
		// copy state (6)
		normalizationFactorsMerged.segment(static_cast<Eigen::Index>(stateOffst) + i*stateVectorSize,stateVectorSize) = normalizationFactors[i].segment(0,stateVectorSize);
		// copy local extras (nLocal)
		normalizationFactorsMerged.segment(static_cast<Eigen::Index>(localOffst)+i*nLocal,nLocal) = normalizationFactors[i].segment(stateVectorSize,nLocal);
		// Update P0 local-diagonal with normalization
		P0_matrices[i].diagonal().segment(stateVectorSize,nLocal).array() = P0_matrices[i].diagonal().segment(stateVectorSize,nLocal).array()/normalizationFactors[i].segment(stateVectorSize,nLocal).array().square();
	}
    // --- Stack global columns from unnormalized H and compute their normalization ---
    std::cout<<"Normalization factors merged!"<<std::endl;

	int totalRows = std::accumulate(unnormalizedDesignMatrices.begin(), unnormalizedDesignMatrices.end(),0,[](int sum, const Eigen::MatrixXd& mat) {return sum + mat.rows();});
	Eigen::MatrixXd stackedGlobalParsDesignMatrix(totalRows,numCols);

	// Extract and stack columns
	{
		int position = 0;
		int Arcs = 0;
		for (const auto &mat: unnormalizedDesignMatrices){
			int numberOfLocalParameters = numberOfLocalParametersAll[Arcs];
			for (size_t j = 0; j<numberOfGlobalParameters; j++){
				stackedGlobalParsDesignMatrix.block(position,j,mat.rows(),1) = mat.col(numberOfLocalParameters+j);
			}
			position += mat.rows();
			Arcs +=1;
		}
	}
	std::cout<<"Stacked unnormalized design matrix of global parameters!"<<std::endl;
	std::cout<<"Now performing normalization..."<<std::endl;

	// Normalize the columns of the stacked global parameters of the design matrices
	Eigen::VectorXd normalizationFactorsGP(stackedGlobalParsDesignMatrix.cols());
	normalizationFactorsGP.setOnes();
	for (int k = 0; k < stackedGlobalParsDesignMatrix.cols(); k++){
	    Eigen::Index maxIndex;
	    double maxAbsVal = stackedGlobalParsDesignMatrix.col(k).cwiseAbs().maxCoeff(&maxIndex);
	    double maxVal = stackedGlobalParsDesignMatrix(maxIndex,k);
	    if (maxVal != 0.0){
			stackedGlobalParsDesignMatrix.col(k) /= maxVal;
			normalizationFactorsGP[k] = maxVal;
	    } else {
		    normalizationFactorsGP[k] = 1.0; // If max is zero, keep normalization factor as 1.0
	    }
	}
	normalizationFactorsMerged.segment(totalNumberOfLocalParameters,numberOfGlobalParameters) = normalizationFactorsGP;
	std::cout<<"Done with normalization!"<<std::endl;
	std::cout<<"Now bringing back the values to the normalized design matrices..."<<std::endl;
	// Reconstructing normalized matrix
	{
		int position = 0;
		int R = 0;
		for (auto &mat: normalizedDesignMatrices){
			for (size_t j = 0; j<numberOfGlobalParameters; j++){
				mat.col(numberOfLocalParametersAll[R]+j) = stackedGlobalParsDesignMatrix.block(position,j,mat.rows(),1);
			}
			position += mat.rows();
			P0_matrices[R].diagonal().segment(numberOfLocalParametersAll[R],numberOfGlobalParameters).array() = P0_matrices[R].diagonal().segment(numberOfLocalParametersAll[R],numberOfGlobalParameters).array()/normalizationFactorsGP.array().square();
			++R;
		}
	}

	#if COV_DEBUG
	for (int i = 0; i < nArcs; ++i) {
		if (!allFinite(normalizedDesignMatrices[i])) std::cerr << "NaN/Inf in normalizedDesignMatrices["<<i<<"]\n";
		if (!allFinite(unnormalizedDesignMatrices[i])) std::cerr << "NaN/Inf in unnormalizedDesignMatrices["<<i<<"]\n";
		if (!allFinite(P0_matrices[i])) std::cerr << "NaN/Inf in P0_matrices["<<i<<"]\n";
		if (!allFinite(weightDiagonals[i])) std::cerr << "NaN/Inf in weightDiagonals["<<i<<"]\n";
		if (!allNonNegative(weightDiagonals[i])) std::cerr << "Negative weight(s) in weightDiagonals["<<i<<"]\n";
		if (!allFinite(normalizationFactors[i])) std::cerr << "NaN/Inf in normalizationFactors["<<i<<"]\n";
		// sanity: local priors diagonal positive?
		auto dloc = P0_matrices[i].diagonal().segment(6, /*nLocal=*/1).array();
		if ( (dloc <= 0.0).any() ) std::cerr << "Non-positive local prior diag in P0_matrices["<<i<<"]\n";
	}
	#endif

	std::cout<<"Updated normalized design matrices!"<<std::endl;
	std::cout<<"Now computing the covariance matrices..."<<std::endl;
	// --- Per-arc (Hᵀ W H + P0) and accumulate global-global sum (exclude repeated priors) ---
	std::vector<Eigen::MatrixXd> resultNormalizedInvCovMatrices; resultNormalizedInvCovMatrices.resize(nArcs);
	EigenThreadsGuard eigen_one_thread{1};
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < nArcs; ++i) {
		const auto& H = normalizedDesignMatrices[i];
		const Eigen::DiagonalMatrix<double, Eigen::Dynamic> W(weightDiagonals[i]);
		Eigen::MatrixXd P_im = H.transpose() * W * H;
		resultNormalizedInvCovMatrices[i] = P_im + P0_matrices[i];  // == resultNormalizedInvCovMatrices[i]
	}
	auto arc_sym_residual = [&](const Eigen::MatrixXd& A){
		double num = (A - A.transpose()).norm();
		double den = A.norm();
		return (den > 0.0) ? num/den : num;
	};
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < nArcs; ++i) {
		const auto& H = normalizedDesignMatrices[i];
		const Eigen::DiagonalMatrix<double, Eigen::Dynamic> W(weightDiagonals[i]);

		Eigen::MatrixXd Pim = H.transpose() * W * H;
		Eigen::MatrixXd M = Pim + P0_matrices[i]; // per-arc info

		double res = arc_sym_residual(M);
		double minDiag = M.diagonal().minCoeff();

//		if (res > 1e-12) {
//			// Temporary symmetrization for diagnostics; comment out once fixed
//			M = 0.5 * (M + M.transpose());
//		}
		// store M back
		resultNormalizedInvCovMatrices[i] = std::move(M);
	}

	#if COV_DEBUG
	for (int i = 0; i < nArcs; ++i) {
		const auto& M = resultNormalizedInvCovMatrices[i]; // == resultNormalizedInvCovMatrices[i]
		if (!allFinite(M)) std::cerr << "NaN/Inf in resultNormalizedInvCovMatrices["<<i<<"]\n";
		double off = (M - M.transpose()).cwiseAbs().maxCoeff();
		//if (off > 1e-12) std::cerr << "Arc "<<i<<" not perfectly symmetric, max|A-AT|="<< off <<"\n";
		// Optional: quick SPD probe on the local block (small)
		Eigen::LLT<Eigen::MatrixXd> llt(M.block(0,0, 6+1, 6+1));
		if (llt.info() != Eigen::Success) std::cerr << "Arc "<<i<<" local block not SPD\n";
	}
	#endif


	std::cout<<"Normalized inverse covariance matrices computed!"<<std::endl;
	std::cout<<"Now summing up the values for the global parameters..."<<std::endl;
	Eigen::MatrixXd inverseNormalizedCovarianceMatrixSummedGP = Eigen::MatrixXd::Zero(numberOfGlobalParameters, numberOfGlobalParameters);
	for (int i = 0; i < resultNormalizedInvCovMatrices.size(); i++){
		const auto& M = resultNormalizedInvCovMatrices[i];
		const auto& GG = M.block(numberOfLocalParametersAll[i], numberOfLocalParametersAll[i], numberOfGlobalParameters, numberOfGlobalParameters);
		if (i == 0){
	        inverseNormalizedCovarianceMatrixSummedGP += GG;
		} else {
			inverseNormalizedCovarianceMatrixSummedGP += GG - P0_matrices[i].block(numberOfLocalParametersAll[i],numberOfLocalParametersAll[i],numberOfGlobalParameters,numberOfGlobalParameters);
		}
	}
	std::cout<<"Global parameters values summed up!"<<std::endl;
        std::cout<<"Now assembling the global normalized inverse covariance matrix..."<<std::endl;

	std::vector<int> prefL(nArcs + 1, 0);
	for (int i = 0; i < nArcs; ++i) prefL[i + 1] = prefL[i] + numberOfLocalParametersAll[i]-6;
	const int totalLoc  = prefL[nArcs];
	using SpMat = Eigen::SparseMatrix<double, Eigen::ColMajor>;
	std::vector<Eigen::Triplet<double>> trips;
	size_t reserve_nnz = numberOfGlobalParameters * (numberOfGlobalParameters + 1) / 2;
	for (int i = 0; i < nArcs; ++i) {
		int nLocal = numberOfLocalParametersAll[i] - stateVectorSize;
		reserve_nnz += size_t(stateVectorSize)*stateVectorSize                 // state–state lower
					 + size_t(nLocal)*nLocal           // local–local lower
					 + size_t(nLocal)*stateVectorSize // local–state (lower)
					 + size_t(numberOfGlobalParameters)*stateVectorSize                 // global–state (lower)
					 + size_t(numberOfGlobalParameters)*nLocal;             // global–local (lower)
	}
	
	trips.reserve(reserve_nnz);
	int globalOffset = totalNumberOfLocalParameters;// nArcs*(stateVectorSize+nLocal);

	auto addLower = [&](int r0, int c0,
	const Eigen::Ref<const Eigen::MatrixXd>& B,
	std::vector<Eigen::Triplet<double>>& out)
	{
	const int R = B.rows(), C = B.cols();
	for (int c = 0; c < C; ++c) {
		for (int r = 0; r < R; ++r) {
			const int Rg = r0 + r;
			const int Cg = c0 + c;
			if (Rg < Cg) continue;
			double v = B(r,c);
			if (v != 0.0) out.emplace_back(Rg,Cg, v);
		}
	}
	};

    for (int i = 0; i<nArcs;i++) {
    	const auto& M = resultNormalizedInvCovMatrices[i];
    	const int startOffset = i * stateVectorSize;
    	int nLocal = numberOfLocalParametersAll[i] - stateVectorSize;
    	const int localOffset = nArcs * stateVectorSize + i * nLocal;

        addLower(startOffset, startOffset, M.block(0, 0, stateVectorSize, stateVectorSize), trips);
    	//local-local covariance matrix
        addLower(localOffset, localOffset, M.block(stateVectorSize, stateVectorSize, nLocal, nLocal), trips);
        //arc-local covariance matrix
    	addLower(localOffset, startOffset, M.block(stateVectorSize,0,nLocal,stateVectorSize), trips);
    	//global-arc covariance matrix
        addLower(globalOffset, startOffset, M.block(stateVectorSize + nLocal,0,numberOfGlobalParameters,stateVectorSize), trips);
    	//global-local covariance matrix
    	addLower(globalOffset, localOffset, M.block(stateVectorSize + nLocal,stateVectorSize, numberOfGlobalParameters, nLocal), trips);
	}
    // global-global covariance matrix
	addLower(globalOffset, globalOffset,inverseNormalizedCovarianceMatrixSummedGP, trips);

	Eigen::SparseMatrix<double> P_global_sparse(total_size, total_size);
	P_global_sparse.setFromTriplets(trips.begin(), trips.end());
	P_global_sparse.makeCompressed();
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>, Eigen::Lower> solver;
	Eigen::SparseMatrix<double> S = P_global_sparse;
	S = S.selfadjointView<Eigen::Lower>();
	solver.compute(P_global_sparse.selfadjointView<Eigen::Lower>());
	auto D = solver.vectorD();
	double minD = D.minCoeff();
	int nonpos = (D.array() <=0).count();
	double maxD = D.maxCoeff();
	#if COV_DEBUG
	std::cerr << "Count nonpositive pivots: " << nonpos << " / " << D.size() << "\n";
	std::cerr << "minD= " << minD << "maxD= " << maxD << "\n";
	std::cerr << "P_global size: " << P_global_sparse.rows() << " x " << P_global_sparse.cols()
			  << "  nnz=" << P_global_sparse.nonZeros() << "\n";
	std::cerr << "symmetry residual: "<< std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10) << symmetryResidual(S) << " (target <= ~1e-12)\n";
	std::cerr << "min diagonal: " << minDiagonal(P_global_sparse) << " (should be > 0 for SPD)\n";
	std::cerr<< "compute info = " << int(solver.info()) << "\n" <<std::flush;
	#endif
	//std::cout<< "min diag = " << P_global_sparse.d << "\n";
	if (solver.info() != Eigen::Success) {
		throw std::runtime_error("Cholesky factorization failed (matrix may be singular or not SPD).");
	}
	
	
	// --- Factor SPD matrix with sparse Cholesky ---
	
	const int n = total_size;
	const int BlockSize = 512; // tune based on memory/cache
	Eigen::MatrixXd NormalizedCovarianceMatrixSummed(n, n);
	NormalizedCovarianceMatrixSummed.setZero();

	// Solve in column blocks: P * X = I_block
	for (int col = 0; col < n; col += BlockSize) {
		const int bs = std::min(BlockSize, n - col);
		Eigen::MatrixXd Iblock = Eigen::MatrixXd::Zero(n, bs);
		// Place identity cols
		for (int k = 0; k < bs; ++k) Iblock(col + k, k) = 1.0;

		Eigen::MatrixXd X = solver.solve(Iblock);
		if (solver.info() != Eigen::Success) {
			throw std::runtime_error("Solve failed while forming covariance.");
		}
		NormalizedCovarianceMatrixSummed.middleCols(col, bs).noalias() = X;
	}
	// --- Vectorized un-normalization: C_unnorm(i,j) = C_norm(i,j) / (s_i * s_j) ---
	Eigen::ArrayXd s = normalizationFactorsMerged.array();
	Eigen::ArrayXXd C = NormalizedCovarianceMatrixSummed.array();
	C.colwise() /= s;
	C.rowwise() /= s.transpose();
	Eigen::MatrixXd resultedUnnormalizedCovarianceMatrix = C.matrix();


	return combinedResults{resultedUnnormalizedCovarianceMatrix, P_global_sparse};
}

static void printUsage(const char* exe)
{
    std::cerr
	<< "Usage:\n"
	<< "  " << exe << "	 [--baseDir ouput --num-arcs N] \n";
}

int main(int argc, char* argv[])
{
	Eigen::initParallel();    // initialize Eigen's internal threading
    std::string baseDir;
    int start =0;
    int step =-1;
    int take =-1;
    int skip =-1;
    int numArcs = -1;
    for (int i=1; i<argc;++i){
	std::string a = argv[i];
	auto next = [&](std::string& dst){ if (++i>=argc) throw std::runtime_error("Missing value after "+a); dst = argv[i]; };

	if (a== "--baseDir") next(baseDir);
	else if (a== "--num-arcs") { std::string t; next(t); numArcs = std::atoi(t.c_str());}
	else if (a== "--start") { std::string t; next(t); start = std::atoi(t.c_str());}
	else if (a== "--step") { std::string t; next(t); step = std::atoi(t.c_str());}
	else if (a== "--take") { std::string t; next(t); take = std::atoi(t.c_str());}
	else if (a== "--skip") { std::string t; next(t); skip = std::atoi(t.c_str());}
	else { std::cerr << "Unknown arg: " << a << "\n"; printUsage(argv[0]); return 1; }
    }
    if (baseDir.empty() || numArcs <= 0)
    { std:: cerr << "Provide --base-dir AND --num-arcs\n"; printUsage(argv[0]); return 1;}
    if (start >= numArcs) throw std::runtime_error("--start >= --num-arcs");

    bool hasStep = (step>0);
    bool hasWindow = (take>0 && skip>=0);

    if (hasStep && hasWindow)
	throw std::runtime_error("Use either --step OR (--take and --skip), not both!");
    if (!hasStep && !hasWindow){
	step = 1;
	hasStep=true;
    }
    std::vector<Eigen::MatrixXd> normalizedDesignMatrices;
    std::vector<Eigen::MatrixXd> unnormalizedDesignMatrices;
    std::vector<Eigen::VectorXd> normalizationFactors;
    std::vector<Eigen::VectorXd> weightDiagonals;
    std::vector<Eigen::MatrixXd> P0_matrices;
    std::vector<double> numberOfLocalParametersAll;
    int itotalDuration = numArcs * 5;
    std::string fileTag = "accumul_InverseAprALLGlobalPars_drag_5" + std::to_string(itotalDuration);
    int used = loadArcResultsRange(baseDir,start,step,take,skip,hasStep,hasWindow,numArcs,&numberOfLocalParametersAll, &normalizedDesignMatrices,&unnormalizedDesignMatrices,&normalizationFactors,&weightDiagonals,&P0_matrices);
    if (used ==0) throw std::runtime_error("No matrices are loaded!");

    Eigen::MatrixXd resultedUnnormalizedCovarianceMatrix, P_global;

    combinedResults res = computeCovarianceMatrix(numberOfLocalParametersAll, numArcs, normalizedDesignMatrices, unnormalizedDesignMatrices, normalizationFactors, weightDiagonals,P0_matrices);
    saveMatrixBinary(baseDir + "/resultedUnnormalizedCovarianceMatrix" + fileTag + ".bin",res.resultedUnnormalizedCovarianceMatrix);


}
