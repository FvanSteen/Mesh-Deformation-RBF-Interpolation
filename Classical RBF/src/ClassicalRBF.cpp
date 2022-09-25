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

	double xDomain = 40, yDomain = 40;
	auto supRad = 2.5*max(xDomain,yDomain);

//	string ifName = "mesh_NACA0012_inv.su2";
//	string ofName = "naca0012.su2";


	string ifName = "TestMesh.su2";
	string ofName = "newmesh.su2";

//	vector<string> intBdryTags = {"airfoil"}; //What if there are multiple tags?
//	vector<string> extBdryTags = {"farfield"};
//
	vector<string> intBdryTags = {"block"}; //What if there are multiple tags?
	vector<string> extBdryTags = {"lower","right","upper", "left"};

/* Move this part to the mesh object
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
*/

	Mesh m(ifName, supRad,debugLvl);
//	Mesh m("TestMesh.su2", supRad);
	m.readMeshFile(intBdryTags,extBdryTags);



//
////	m.obtainCoords();
//	cout << "done" << endl;
//	cout << m.bdryNodes << endl;
	Eigen::MatrixXd Phi = m.interpMat(m.bdryNodes,m.bdryNodes); // makes i-matrix for finding coefficients
//	cout << Phi << endl;
	Eigen::VectorXd dxVec(m.nBdryNodes); // initialise vector with known displacements
	Eigen::VectorXd dyVec(m.nBdryNodes);

	dxVec.setZero(); // set all displacements to zero
	dyVec.setZero();

//	cout << m.intBdryNodes << endl;
//	cout << dxVec.size() << endl;
//	cout << m.intBdryNodes << endl;
//	cout << m.bdryNodes << endl;
//
//
	double xDef = -0.25, yDef = -0.25;
//
	for(int i=0;i<m.intBdryNodes.size();i++){
		int j=0;
		while(j<m.nBdryNodes){
			if(m.intBdryNodes(i)==m.bdryNodes(j)){
//				cout << m.intBdryNodes(i) << '\t' << m.bdryNodes(j) << endl;

				dxVec(j) = xDef;
				dyVec(j) = yDef;
				break;
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
//	cout << Phi.size() << '\t' << dxVec.size() << endl;
//	cout << m.bdryNodes << endl;


//	cout << m.extBdryNodes << endl;
//
	cout << "size int bndry nodes: " << m.intBdryNodes.size() <<endl;
	cout << "size ext bndry nodes: " << m.extBdryNodes.size() <<endl;
	cout << m.bdryNodes.size() << endl;
	cout << Phi.rows() << '\t' << Phi.cols() << endl;
	Eigen::VectorXd alpha_x = Phi.llt().solve(dxVec);
	Eigen::VectorXd alpha_y = Phi.llt().solve(dyVec);


//
//
//	// setting up i-matrix phi_i,b
	cout << m.nPnts-m.nBdryNodes << endl;
	cout << m.intNodes.size() << endl;

//	cout << m.intNodes << endl;
//	cout << m.bdryNodes << endl;
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
	m.writeMeshFile(ifName, ofName);
//	m.writeMeshFile("TestMesh.su2", "newmesh.su2");
	cout << "Done writing mesh file" << endl;
}
