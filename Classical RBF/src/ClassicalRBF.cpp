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
	const string ifName = "TestMesh.su2";
	const string ofName = "TestMesh_def.su2";
//	string ifName = "25x25mesh.su2";
//	string ofName = "25x25mesh_def.su2";
	const vector<string> intBdryTags = {"block"}; //What if there are multiple  tags?
	const vector<string> extBdryTags = {"lower","right","upper", "left"};


	// initialising object m, reads mesh input file in constructor.
	Mesh meshOb(ifName,ofName, intBdryTags, extBdryTags, rFactor, debugLvl);

//	cout << "Moving nodes are: \n"<< meshOb.movingNodes << "\nSliding nodes are: \n" << meshOb.slidingNodes << endl;
//	cout << meshOb.extBdryNodes << endl;
//	cout << meshOb.intBdryNodes << endl;
	// moving nodes
	// sliding nodes
	// internal nodes



	const double xDef = -0.2, yDef = -0.2, rotDefDeg = 60;
	const int steps = 1;
	Eigen::RowVectorXd rotationPnt(2) ;
	rotationPnt << 0.5,0.5;
	rbf rbf(meshOb, xDef, yDef, rotDefDeg, steps, rotationPnt);
//	rbf.performRbfInterpolation();
//	rbf.performRbfDS();
	rbf.performRbfPS();
//	const double xDef = -10, yDef = -10, rotDefDeg = -30;
//	const int steps = 20;
//	Eigen::Vector2d rotationPoint = {0.25,0.0};
//	rbf.performRbfInterpolation(xDef,yDef,rotDefDeg,steps,rotationPoint);
//	rbf.performRbfDS(xDef, yDef, rotDefDeg, steps, rotationPoint);

}
