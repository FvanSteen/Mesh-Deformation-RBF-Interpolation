#include "Mesh.h"
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <algorithm>
using Eigen::MatrixXd;
using namespace std;

int main()
{
	double xDomain = 1, yDomain = 1;
	auto supRad = 2.5*max(xDomain,yDomain);

	Mesh m("TestMesh.su2", supRad);

	m.findProblemChars();

	m.obtainCoords();

	Eigen::MatrixXd Phi = m.interpMat(m.bdryNodes,m.bdryNodes); // makes i-matrix for finding coefficients
	Eigen::VectorXd dVec(m.nBdryNodes); // initialise vector with known displacements
	dVec.setZero(); // set all displacements to zero


	int arrDelta[4]{9,10,13,14}; // hardcoding of indices that belong to inner block

	// setting displacement of inner block to 0.05;
	for(int i=0;i<4;i++){
		dVec(arrDelta[i]) = 0.05;
	}

	// solving for coefficients
	Eigen::VectorXd alpha = Phi.llt().solve(dVec);


	// setting up i-matrix phi_i,b
	Eigen::MatrixXd Phi_ib = m.interpMat(m.intNodes,m.bdryNodes);

	// finding displacement of internal nodes
	Eigen::VectorXd disp =  Phi_ib*alpha;


	m.updateNodes(dVec,disp);


	m.WriteMeshFile("newmesh.su2");
	cout << "Done" << endl;
}
