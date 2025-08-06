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
