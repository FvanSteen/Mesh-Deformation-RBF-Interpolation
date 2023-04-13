#include "probParams.h"
#include "ReadConfigFile.h"
#include "rbfstd.h"
#include "rbfps.h"
#include "rbfds.h"
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <vector>

#include "SPDS.h"
#include "nanoflann.hpp"
#include <cstdlib>
#include <ctime>
#include <iostream>
#include "CoordTransform.h"
using namespace std;

int main()
{

//	SPDS ds;
//	ds.kdtree<double>(5,2);


//	std::cout << "DONE ALL" << std::endl;
//	std::exit(0);


	// config file containing all information required to perform the mesh deformation
	std::string configFile = "config_file.txt";

	// initliasing the structure in which the problem parameters are saved
	struct probParams probParams;

	// calling class to read the configuration file
	// Variables required to read the mesh file are public variables and later on inherited by the Mesh class
	// Variables required for performing the interpolation are saved in the probParams structure and passed on in the main rbf class.
	ReadConfigFile config(configFile,probParams);


	// lvl indicating the amount of debug messages
	int debugLvl = 3;

	// initialising class object m, reads mesh input file in constructor.
	Mesh meshOb(probParams, debugLvl);

	getNodeType n(probParams, meshOb);

	std::clock_t s = std::clock();

	if(probParams.smode == "none"){
		rbf_std rbf(probParams, meshOb, n);
	}else if(probParams.smode == "ps"){
		rbf_ps rbf(probParams, meshOb, n);
	}else{
		rbf_ds rbf(probParams, meshOb, n);
	}


	std::clock_t e = std::clock();
	long double time_elapsed_ms =  1000.0*(e-s) / CLOCKS_PER_SEC;
	std::cout << "CPU time: " << time_elapsed_ms/1000 << " ms\n";


	meshOb.writeMeshFile(probParams.mesh_ifName, probParams.mesh_ofName);


	// Initialising the rbf class with the mesh data and the problem parameters
//	rbf rbf(meshOb, probParams);

	// calling the main rbf function that performs the mesh deformation
//	rbf.RBFMain();
}
