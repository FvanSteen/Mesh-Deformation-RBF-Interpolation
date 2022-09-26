#include "Mesh.h"
#include "rbf.h"
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <string>
using Eigen::MatrixXd;
using namespace std;

rbf::rbf(Mesh meshOb)
:m(meshOb)
{
	}

void rbf::performRbfInterpolation(){

	Eigen::MatrixXd Phi_bb = getPhi(m.bdryNodes,m.bdryNodes); // makes i-matrix for finding coefficients

	Eigen::VectorXd dxVec(m.bdryNodes.size()); // initialise vector with known displacements
	Eigen::VectorXd dyVec(m.bdryNodes.size());

	dxVec.setZero(); // set all displacements to zero
	dyVec.setZero();

	double xDef = -0.25, yDef = -0.25;

	for(int i=0;i<m.intBdryNodes.size();i++){
		int j=0;
		while(j<m.bdryNodes.size()){
			if(m.intBdryNodes(i)==m.bdryNodes(j)){
//				cout << meshOb.intBdryNodes(i) << '\t' << meshOb.bdryNodes(j) << endl;

				dxVec(j) = xDef;
				dyVec(j) = yDef;
				break;
			}
			j++;
		}
	}

//
//	Eigen::ArrayXi test = {1,2,3,4,5};
//	cout << test << endl;
//
//	Eigen::Array3i idx = {0,2,4};
//	cout << idx << endl;
//
//	cout << test << endl;




	Eigen::VectorXd alpha_x = Phi_bb.llt().solve(dxVec);
	Eigen::VectorXd alpha_y = Phi_bb.llt().solve(dyVec);

	Eigen::MatrixXd Phi_ib = getPhi(m.intNodes,m.bdryNodes);
	Eigen::VectorXd xDisp =  Phi_ib*alpha_x;
	Eigen::VectorXd yDisp =  Phi_ib*alpha_y;

	updateNodes(dxVec,dyVec,xDisp,yDisp);
	m.writeMeshFile(newCoords);
}


Eigen::MatrixXd rbf::getPhi(Eigen::ArrayXi idxSet1, Eigen::ArrayXi idxSet2){
	Eigen::MatrixXd Phi(idxSet1.size(), idxSet2.size());
	for(int i=0; i<idxSet1.size();i++){
		for(int j=0; j<idxSet2.size();j++){
			double dist = sqrt(pow(m.coords(idxSet1(i),0)-m.coords(idxSet2(j),0),2) + pow(m.coords(idxSet1(i),1)-m.coords(idxSet2(j),1),2));
			Phi(i,j) = rbfEval(dist);
		}
	}
	return Phi;
}

double rbf::rbfEval(double distance){
	double xi = distance/m.r;	// distance scaled by support radius
	double f_xi = pow((1-xi),4)*(4*xi+1);
	return f_xi;
}

void rbf::updateNodes(Eigen::VectorXd dxVec,Eigen::VectorXd dyVec, Eigen::VectorXd xDisp,Eigen::VectorXd yDisp){
	newCoords = m.coords;

	newCoords(m.bdryNodes,0) = m.coords(m.bdryNodes,0) + dxVec;
	newCoords(m.bdryNodes,1) = m.coords(m.bdryNodes,1) + dyVec;
	newCoords(m.intNodes,0) = m.coords(m.intNodes,0) + xDisp;
	newCoords(m.intNodes,1) = m.coords(m.intNodes,1) + yDisp;
//	return newCoords;
}
