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

	string ifName = "TestMesh.su2";
	string ofName = "newmesh.su2";
	vector<string> intBdryTags = {"block"}; //What if there are multiple tags?
	vector<string> extBdryTags = {"lower","right","upper", "left"};


	// initialising object m, reads mesh input file in constructor.
	Mesh meshOb(ifName,ofName, intBdryTags, extBdryTags, rFactor, debugLvl);

	rbf rbf(meshOb);
	rbf.performRbfInterpolation();




//	cout << Phi << endl;

//	cout << meshOb.intBdryNodes << endl;
//	cout << dxVec.size() << endl;
//	cout << meshOb.intBdryNodes << endl;
//	cout << meshOb.bdryNodes << endl;
//
//

//	cout << dxVec << endl;
//	cout << dyVec << endl;
////	int arrDelta[4]{9,10,13,14}; // hardcoding of indices that belong to inner block
////
////	// setting displacement of inner block to 0.05;
////	for(int i=0;i<4;i++){
////		dVec(arrDelta[i]) = defX;
////	}
//
//	// solving for coefficients
//	cout << Phi.size() << '\t' << dxVec.size() << endl;
//	cout << meshOb.bdryNodes << endl;


//	cout << meshOb.extBdryNodes << endl;
//
//	cout << "size int bndry nodes: " << meshOb.intBdryNodes.size() <<endl;
//	cout << "size ext bndry nodes: " << meshOb.extBdryNodes.size() <<endl;
//	cout << meshOb.bdryNodes.size() << endl;
//	cout << Phi.rows() << '\t' << Phi.cols() << endl;



//
//
//	// setting up i-matrix phi_i,b
//	cout << meshOb.nPnts-meshOb.nBdryNodes << endl;
//	cout << meshOb.intNodes.size() << endl;

//	cout << meshOb.intNodes << endl;
//	cout << meshOb.bdryNodes << endl;


//
//	// finding displacement of internal nodes

//
//

//


//	m.writeMeshFile("TestMesh.su2", "newmesh.su2");
	cout << "Done writing mesh file" << endl;
}
