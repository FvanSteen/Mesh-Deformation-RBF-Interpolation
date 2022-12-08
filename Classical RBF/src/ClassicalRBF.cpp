#include "rbf.h"
#include "probParams.h"
#include "ReadConfigFile.h"
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <vector>

using namespace std;


int main()
{
	// config file containing all information required to perform the mesh deformation
	std::string configFile = "config_file.txt";

	// initliasing the structure in which the problem parameters are saved
	struct probParams probParams;

	// calling class to read the configuration file
	// Variables required to read the mesh file are public variables and later on inherited by the Mesh class
	// Variables required for performing the interpolation are saved in the probParams structure and passed on in the main rbf class.
	ReadConfigFile config(configFile,probParams);

	// lvl indicating the amount of debug messages
	int debugLvl = 2;

	// initialising class object m, reads mesh input file in constructor.
	Mesh meshOb(config, debugLvl);

	// Initialising the rbf class with the mesh data and the problem parameters
	rbf rbf(meshOb, probParams);

	// calling the main rbf function that performs the mesh deformation
	rbf.RBFMain();
}
