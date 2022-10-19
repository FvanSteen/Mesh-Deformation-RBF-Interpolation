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
	int debugLvl = 3;
	const double rFactor = 2.5;

//	string ifName = "mesh_NACA0012_inv.su2";
//	string ofName = "mesh_NACA0012_inv_def.su2";
//	vector<string> intBdryTags = {"airfoil"}; //What if there are multiple tags?
//	vector<string> extBdryTags = {"farfield"};
//
//	const string ifName = "5x5.su2";
//	const string ofName = "5x5_def.su2";
//	string ifName = "25x25.su2";
//	string ofName = "25x25_def.su2";
//	const vector<string> intBdryTags = {"block"}; //What if there are multiple  tags?
//	const vector<string> extBdryTags = {"lower","right","upper", "left"};

	string ifName = "15x15x5.su2";
	string ofName = "15x15x5_def.su2";
	const vector<string> intBdryTags = {"BLOCK"};
	const vector<string> extBdryTags = {"FRONT", "BACK", "LEFT", "RIGHT", "LOWER", "UPPER"};


	const string slidingMode = "none";

	// initialising object m, reads mesh input file in constructor.
	Mesh meshOb(ifName,ofName, intBdryTags, extBdryTags, rFactor, debugLvl, slidingMode);


	const double xDef = -0.2, yDef = -0.2, zDef=0.0, rotDefDeg = 0;
	const double xRot = 0, yRot = 0, zRot=0;
	const int steps = 1;
	Eigen::RowVectorXd rotationPnt(meshOb.nDims);
	if(meshOb.nDims == 2){
		rotationPnt << 0.5,0.5;
	}else{
		rotationPnt << 0.5,0.5,0.5;
	}

	rbf rbf(meshOb, xDef, yDef, zDef, rotDefDeg, steps, rotationPnt, slidingMode);
	rbf.RBFMain();

//	const double xDef = -10, yDef = -10, rotDefDeg = -30;
//	const int steps = 20;
//	Eigen::Vector2d rotationPoint = {0.25,0.0};
//	rbf.performRbfInterpolation(xDef,yDef,rotDefDeg,steps,rotationPoint);
//	rbf.performRbfDS(xDef, yDef, rotDefDeg, steps, rotationPoint);

}
