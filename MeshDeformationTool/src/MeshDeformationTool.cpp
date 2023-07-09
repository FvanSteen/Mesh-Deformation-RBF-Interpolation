#include "ReadConfigFile.h"
#include "rbfstd.h"
#include "rbfps.h"
#include "rbfds.h"
#include <iostream>
#include <ctime>
#include "MeshQuality.h"
#include "ProbParams.h"
#include "WriteResults.h"

int main(int argc, char* argv[]){

	// configuration file is given as the second command line variable, after the program name
	std::string cfgFile = argv[1];

	// initliasing the structure in which the problem parameters are saved
	struct probParams probParams;

	// Class required for reading the specified configuration and saving the required problem parameters in the probParams structure
	ReadConfigFile config(cfgFile,probParams);

	// Mesh class for reading the input mesh file and establishing mesh related information
	Mesh meshOb(probParams);

	// Class for calculating mesh quality parameters, initialised here for information on the undeformed mesh
	MeshQuality Qual(probParams, meshOb.coords);

	// Class used for setting the node types and adding control nodes for the data reduction methods
	getNodeType n(probParams, meshOb);

	// starting clock to determine CPU time of the mesh deformation routine
	std::clock_t start = std::clock();

	// Performing the regular, pseudo sliding or direct sliding RBF interpolation
	if(probParams.smode == "none"){
		rbf_std rbf(probParams, meshOb, n);
	}else if(probParams.smode == "ps"){
		rbf_ps rbf(probParams, meshOb, n);
	}else{
		rbf_ds rbf(probParams, meshOb, n);
	}

	// stopping clock after mesh deformation routine
	std::clock_t end = std::clock();
	long double time_elapsed_ms =  1.0*(end-start) / CLOCKS_PER_SEC;
	std::cout << "CPU time: " << time_elapsed_ms << " seconds" << std::endl;


	// generating mesh output file
	meshOb.writeMeshFile(probParams.directory, probParams.mesh_ifName, probParams.mesh_ofName);

	// Generate mesh quality of the deformed mesh if required
	if(probParams.generateQuality){
		Qual.getDeformedMeshQual(meshOb.coords, -1);
		std::cout << "\nMesh qualities of the deformed mesh:\nMaximum quality:\t"  << Qual.defQuals[2] << std::endl;
		std::cout << "Mean quality:\t\t" << Qual.defQuals[1] << std::endl;
		std::cout << "Minimum quality:\t" << Qual.defQuals[0] << "\n" << std::endl;
	}

}
