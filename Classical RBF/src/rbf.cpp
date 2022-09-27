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
	cout << "Building matrix Phi_bb" << endl;
	Eigen::MatrixXd Phi_bb = getPhi(m.bdryNodes,m.bdryNodes); // makes i-matrix for finding coefficients

	double xDef = -0.25, yDef = -0.25;

	cout << "Building deformation vectors" << endl;
	Eigen::VectorXd dxVec = getDefVec(xDef);
	Eigen::VectorXd dyVec = getDefVec(yDef);

	cout << "Solving for interpolation coefficients" << endl;
	Eigen::VectorXd alpha_x = Phi_bb.llt().solve(dxVec);
	Eigen::VectorXd alpha_y = Phi_bb.llt().solve(dyVec);

	cout << "Building matrix Phi_ib" << endl;
	Eigen::MatrixXd Phi_ib = getPhi(m.intNodes,m.bdryNodes);

	cout << "Finding displacement internal nodes" << endl;
	Eigen::VectorXd xDisp =  Phi_ib*alpha_x;
	Eigen::VectorXd yDisp =  Phi_ib*alpha_y;

	cout << "Updating node coordinates" << endl;
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

Eigen::VectorXd rbf::getDefVec(double def){
	Eigen::VectorXd defVec = Eigen::VectorXd::Zero(m.bdryNodes.size()); 		// initialise zero vector with size of boundary nodes

	int idx = 0;
	for(int i=0; i<m.intBdryNodes.size(); i++){
		idx = idx+distance(m.bdryNodes.begin()+idx,find(m.bdryNodes.begin()+idx, m.bdryNodes.end(),m.intBdryNodes(i)));
		defVec(idx) = def;
	}
	return defVec;
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
