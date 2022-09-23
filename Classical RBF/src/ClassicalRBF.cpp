#include "Mesh.h"
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

	double xDomain = 1, yDomain = 1;
	auto supRad = 2.5*max(xDomain,yDomain);

//	vector<string> intBdryTags = {"airfoil"}; //What if there are multiple tags?
//	vector<string> extBdryTags = {"farfield"};
//
	vector<string> intBdryTags = {"block"}; //What if there are multiple tags?
	vector<string> extBdryTags = {"lower","right","upper", "left"};


	if(debugLvl>2){
		cout << "Domain length x: " << xDomain << endl;
		cout << "Domain length y: " << yDomain << endl;
		cout << "Support Radius r: " << supRad << endl;
		cout << "External boundary tags: ";

		for(string i : extBdryTags){
		  cout << i << " ";
		}

		cout << endl;
		cout << "Internal boundary tags: ";
		for(string i : intBdryTags){
			cout << i << " " ;
		}

		cout << endl;

	}

//	Mesh m("mesh_NACA0012_inv.su2", supRad);
	Mesh m("TestMesh.su2", supRad);
	m.findProblemChars(intBdryTags,extBdryTags);



//
////	m.obtainCoords();
//	cout << "done" << endl;
	cout << m.bdryNodes << endl;
	Eigen::MatrixXd Phi = m.interpMat(m.bdryNodes,m.bdryNodes); // makes i-matrix for finding coefficients
	cout << Phi << endl;
	Eigen::VectorXd dxVec(m.nBdryNodes); // initialise vector with known displacements
	Eigen::VectorXd dyVec(m.nBdryNodes);
	dxVec.setZero(); // set all displacements to zero
	dyVec.setZero();

//	cout << m.intBdryNodes << endl;
//
//
	double xDef = -0.15, yDef = -0.1;
//
	for(int i=0;i<4;i++){
		int j=0;
		while(j<m.nBdryNodes){
			if(m.intBdryNodes(i)==m.bdryNodes(j)){
				dxVec(j) = xDef;
				dyVec(j) = yDef;
			}
			j++;
		}
	}
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
	Eigen::VectorXd alpha_x = Phi.llt().solve(dxVec);
	Eigen::VectorXd alpha_y = Phi.llt().solve(dyVec);
//
//
//	// setting up i-matrix phi_i,b
	Eigen::MatrixXd Phi_ib = m.interpMat(m.intNodes,m.bdryNodes);
//
//	// finding displacement of internal nodes
	Eigen::VectorXd xDisp =  Phi_ib*alpha_x;
	Eigen::VectorXd yDisp =  Phi_ib*alpha_y;
//
//
	m.updateNodes(dxVec,dyVec,xDisp,yDisp);
//
//
//	m.WriteMeshFile("newmesh.su2");
	m.wmf("TestMesh.su2", "newmesh.su2");
	cout << "Done writing mesh file" << endl;
}
