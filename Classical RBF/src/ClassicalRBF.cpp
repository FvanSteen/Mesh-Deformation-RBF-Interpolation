#include "Mesh.h"
#include "rbf.h"
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
using namespace std;

// todo make the curved ds work for the moving periodicmode
// clean up the rbf.cpp file
// start timing the different components to see which parts take longest
// do summation instead of matrix vector multiplication



int main()
{
	int debugLvl = 2;
	const double rFactor = 2.5;

//	string ifName = "mesh_NACA0012_inv.su2";
//	string ofName = "mesh_NACA0012_inv_def.su2";
//	vector<string> intBdryTags = {"airfoil"};
//	vector<string> extBdryTags = {"farfield"};

//	const string ifName = "5x5.su2";
//	const string ofName = "5x5_def.su2";
	string ifName = "25x25.su2";
	string ofName = "25x25_def.su2";
	const vector<string> intBdryTags = {"block"}; //What if there are multiple  tags?
	const vector<string> extBdryTags = {"lower","right","upper", "left"};

//	string ifName = "5x5x5.su2";
//	string ofName = "5x5x5_def.su2";
//	const vector<string> intBdryTags = {"BLOCK"};
//	const vector<string> extBdryTags = {"FRONT", "BACK", "LEFT", "RIGHT", "LOWER", "UPPER"};
	const string slidingMode = "none";
	const bool curved = false; // only relevant in case of ds algorithm
	const string periodicMode = "none";				// none for no periodicity, periodic for making the domain periodic in a to be specified direction ,fixed for allowing periodic boundary displacement with fixed corners, moving for periodic boundary displacement with moving corners
	const vector<string> periodicBdry = {"upper","lower"};
	const string periodicDirection = "y";

	// initialising class object m, reads mesh input file in constructor.
	Mesh meshOb(ifName,ofName, intBdryTags, extBdryTags, rFactor, debugLvl, slidingMode, periodicBdry, periodicMode);


	const double xDef = -0.2, yDef = -0.32, zDef= -0.0;

	const int steps = 20;
	Eigen::VectorXd rotVec;
	Eigen::RowVectorXd rotationPnt(meshOb.nDims);
	Eigen::VectorXd displacementVector(meshOb.nDims);
	if(meshOb.nDims == 2){
		rotVec.resize(1);
		rotVec << 60;
		rotationPnt << 0.5,0.5;
		displacementVector << xDef/steps,yDef/steps;
	}else{
		rotationPnt << 0.5,0.5,0.5;
		rotVec.resize(3);
		rotVec << 0,0,60;
		displacementVector << xDef/steps,yDef/steps,zDef/steps;
	}

	rbf rbf(meshOb, displacementVector, rotVec, steps, rotationPnt, slidingMode, periodicDirection, curved);
	rbf.RBFMain();

}
