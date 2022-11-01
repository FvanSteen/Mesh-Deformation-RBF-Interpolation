#include "Mesh.h"
#include "rbf.h"
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
using Eigen::MatrixXd;
using namespace std;

int main()
{
	int debugLvl = 2;
	const double rFactor = 2.5;

//	string ifName = "mesh_NACA0012_inv.su2";
//	string ofName = "mesh_NACA0012_inv_def.su2";
//	vector<string> intBdryTags = {"airfoil"}; //What if there are multiple tags?
//	vector<string> extBdryTags = {"farfield"};
//
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
	const string slidingMode = "ps";
	const string periodicMode = "moving";
	const vector<string> periodicBdry = {"lower","upper"};
//	const vector<string> periodicBdry = {};
	// initialising class object m, reads mesh input file in constructor.
	Mesh meshOb(ifName,ofName, intBdryTags, extBdryTags, rFactor, debugLvl, slidingMode, periodicBdry, periodicMode);


	const double xDef = -0.2, yDef = -0.32, zDef= -0.0;

	const int steps = 4;
	Eigen::VectorXd rotVec;
	Eigen::RowVectorXd rotationPnt(meshOb.nDims);
	if(meshOb.nDims == 2){
		rotVec.resize(1);
		rotVec << 60;
		rotationPnt << 0.5,0.5;
	}else{
		rotationPnt << 0.5,0.5,0.5;
		rotVec.resize(3);
		rotVec << 0,0,60;
	}

	rbf rbf(meshOb, xDef, yDef, zDef, rotVec, steps, rotationPnt, slidingMode);
	rbf.RBFMain();

}
