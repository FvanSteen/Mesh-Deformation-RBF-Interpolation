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

void rbf::performRbfInterpolation(const double& xDef, const double& yDef, const double& rotDefDeg, const int& steps, Eigen::VectorXd rotPnt){

	newCoords = m.coords;
	for(int i=0; i<steps; i++){
		cout << "Building matrix Phi_bb" << endl;
		Eigen::MatrixXd Phi_bb = getPhi(m.bdryNodes,m.bdryNodes); // makes i-matrix for finding coefficients

		cout << "Building matrix Phi_ib" << endl;
		Eigen::MatrixXd Phi_ib = getPhi(m.intNodes,m.bdryNodes);
		cout << "Building deformation vectors" << endl;

		cout << "Deformation step: " << i+1 << endl;
		Eigen::MatrixXd defMat = getDefVec(xDef/steps,yDef/steps,rotDefDeg/steps,rotPnt);

//	cout << defMat << endl;
//	Eigen::VectorXd dxVec = defMat.col(0);
//	Eigen::VectorXd dyVec = defMat(Eigen::all,1);

//	cout << dxVec << "\n\n" << dyVec << "\n\n" << endl;
		cout << "Solving for interpolation coefficients" << endl;
		Eigen::VectorXd alpha_x = Phi_bb.llt().solve(defMat.col(0));
		Eigen::VectorXd alpha_y = Phi_bb.llt().solve(defMat.col(1));

		cout << "Finding displacement internal nodes" << endl;
		Eigen::VectorXd xDisp =  Phi_ib*alpha_x;
		Eigen::VectorXd yDisp =  Phi_ib*alpha_y;

		cout << "Updating node coordinates" << endl;
		updateNodes(defMat.col(0),defMat.col(1),xDisp,yDisp);

		rotPnt(0) = rotPnt(0) + xDef/steps;
		rotPnt(1) = rotPnt(1) + yDef/steps;

	}
	m.writeMeshFile(newCoords);
}


Eigen::MatrixXd rbf::getPhi(Eigen::ArrayXi idxSet1, Eigen::ArrayXi idxSet2){
	Eigen::MatrixXd Phi(idxSet1.size(), idxSet2.size());
	for(int i=0; i<idxSet1.size();i++){
		for(int j=0; j<idxSet2.size();j++){
			double dist = sqrt(pow(newCoords(idxSet1(i),0)-newCoords(idxSet2(j),0),2) + pow(newCoords(idxSet1(i),1)-newCoords(idxSet2(j),1),2));
			Phi(i,j) = rbfEval(dist);
		}
	}
	return Phi;
}

Eigen::MatrixXd rbf::getDefVec(double xDef, double yDef, double rotDefDeg, Eigen::VectorXd rotPnt){
	Eigen::MatrixXd defMat = Eigen::MatrixXd::Zero(m.bdryNodes.size(),m.nDims); 		// initialise zero vector with size of boundary nodes
	Eigen::MatrixXd rotDef = getRotDef(rotDefDeg,rotPnt); //look at the passing by reference
	int idx = 0;
	for(int i=0; i<m.intBdryNodes.size(); i++){
		idx = idx+distance(m.bdryNodes.begin()+idx,find(m.bdryNodes.begin()+idx, m.bdryNodes.end(),m.intBdryNodes(i)));
		defMat(idx,0) = xDef + rotDef(i,0);
		defMat(idx,1) = yDef + rotDef(i,1);
	}
	return defMat;
}

Eigen::MatrixXd rbf::getRotDef(double rotDef,Eigen::VectorXd rotPnt){
	double rotRadians = rotDef*M_PI/180;
//	Eigen::Vector2d rotPoint = {0.5,0.5};
	Eigen::Matrix2d rotMat;
	rotMat << cos(rotRadians), -sin(rotRadians),sin(rotRadians), cos(rotRadians);
	Eigen::Vector2d node;
	Eigen::MatrixXd rotDefMat(m.intBdryNodes.size(),m.nDims);

	for(int i = 0; i<m.intBdryNodes.size(); i++){
			node = newCoords(m.intBdryNodes(i),Eigen::all);
			auto relDist = node-rotPnt;
			rotDefMat(i,Eigen::seqN(0,m.nDims)) = rotMat*(relDist)+rotPnt-node;
	}
	return rotDefMat;
}

double rbf::rbfEval(double distance){
	double xi = distance/m.r;	// distance scaled by support radius
	double f_xi = pow((1-xi),4)*(4*xi+1);
	return f_xi;
}

void rbf::updateNodes(Eigen::VectorXd dxVec,Eigen::VectorXd dyVec, Eigen::VectorXd xDisp,Eigen::VectorXd yDisp){
	newCoords(m.bdryNodes,0) = newCoords(m.bdryNodes,0) + dxVec;
	newCoords(m.bdryNodes,1) = newCoords(m.bdryNodes,1) + dyVec;
	newCoords(m.intNodes,0) = newCoords(m.intNodes,0) + xDisp;
	newCoords(m.intNodes,1) = newCoords(m.intNodes,1) + yDisp;
//	return newCoords;
}
