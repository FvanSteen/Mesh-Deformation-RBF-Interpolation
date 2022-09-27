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
	double rFactor = 2.5;

//	string ifName = "mesh_NACA0012_inv.su2";
//	string ofName = "naca0012.su2";
//	vector<string> intBdryTags = {"airfoil"}; //What if there are multiple tags?
//	vector<string> extBdryTags = {"farfield"};

//	string ifName = "TestMesh.su2";
//	string ofName = "newmesh.su2";
	string ifName = "25x25mesh.su2";
	string ofName = "25x25mesh_def.su2";
	vector<string> intBdryTags = {"block"}; //What if there are multiple  tags?
	vector<string> extBdryTags = {"lower","right","upper", "left"};


	// initialising object m, reads mesh input file in constructor.
	Mesh meshOb(ifName,ofName, intBdryTags, extBdryTags, rFactor, debugLvl);


	rbf rbf(meshOb);

	const double xDef = 0, yDef = 0, rotDefDeg = 45;
	rbf.performRbfInterpolation(xDef,yDef,rotDefDeg);
}
