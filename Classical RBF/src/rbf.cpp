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
{}

void rbf::performRbfDS(const double& xDef, const double& yDef, const double& rotDefDeg, const int& steps, Eigen::VectorXd rotPnt){
	cout << "Performing RBF DS " << endl;
	newCoords = m.coords;
	Eigen::ArrayXXd n(m.slidingNodes.size(),2);		// all the vectors midPnts etc should be updated after each step
	Eigen::ArrayXXd t(m.slidingNodes.size(),2);
	m.getNodeVecs(m.slidingNodes, n, t);
//
//	cout << "\n Sliding Nodes:\n" << m.slidingNodes << endl;
//	cout << "\n Moving Nodes: \n" << m.movingNodes << endl;
//	cout << "\n Internal boundary nodes: \n" << m.intBdryNodes << endl;
//	cout << "\n Internal Nodes: \n" << m.intNodes << endl;

	Eigen::ArrayXi dispNodes(m.intBdryNodes.size()+m.movingNodes.size());
	dispNodes << m.intBdryNodes, m.movingNodes;

//	cout << dispNodes << endl;
	for(int i=0; i<steps; i++){

		cout << "step: " << i << endl;
		Eigen::MatrixXd defMat = getDefVecDS(xDef/steps,yDef/steps,rotDefDeg/steps,rotPnt,m.intBdryNodes);

		Eigen::MatrixXd Phi_mm = getPhi(dispNodes,dispNodes); // makes i-matrix for finding coefficients
		Eigen::MatrixXd Phi_ms = getPhi(dispNodes,m.slidingNodes);
		Eigen::MatrixXd Phi_sm = getPhi(m.slidingNodes,dispNodes);
		Eigen::MatrixXd Phi_ss = getPhi(m.slidingNodes,m.slidingNodes);

//		cout << Phi_mm.rows() << '\t' << Phi_mm.cols() << endl;
//		cout << Phi_ms.rows() << '\t' << Phi_ms.cols() << endl;
//		cout << Phi_sm.rows() << '\t' << Phi_sm.cols() << endl;
//		cout << Phi_ss.rows() << '\t' << Phi_ss.cols() << endl;

		Eigen::VectorXd defVec = Eigen::VectorXd::Zero(2*(dispNodes.size()+m.slidingNodes.size()));
		defVec(Eigen::seq(0,m.intBdryNodes.size()-1)) = defMat.col(0);
		defVec(Eigen::seq(dispNodes.size(),dispNodes.size()+m.intBdryNodes.size()-1)) = defMat.col(1);

		Eigen::MatrixXd Phi = Eigen::MatrixXd::Zero(2*(dispNodes.size()+m.slidingNodes.size()),2*(dispNodes.size()+m.slidingNodes.size()));
//		cout << "here" << endl;
	//	cout << Phi_ss << endl;
		const int a = 16; // find some way to produce const int of the various dimensions within the mesh class
		const int b = 96;

//		const int a = 200; // find some way to produce const int of the various dimensions within the mesh class
//		const int b = 50;

		Phi.block<a,a>(0,0) = Phi_mm;
		Phi.block<a,a>(a,a) = Phi_mm;
		Phi.block<a,b>(0,2*a) = Phi_ms;
		Phi.block<a,b>(a,2*a+b) = Phi_ms;
	//	cout << Phi << endl;

	//	cout << Phi_sm << endl;
	//	cout << n.col(0) << endl;
	//	cout << Phi_sm.array().colwise() * n.col(0) << endl;

		Phi.block<b,a>(2*a,0) = Phi_sm.array().colwise() * n.col(0);
		Phi.block<b,a>(2*a,a) = Phi_sm.array().colwise() * n.col(1);
		Phi.block<b,b>(2*a,2*a) = Phi_ss.array().colwise() * n.col(0);
		Phi.block<b,b>(2*a,2*a+b) = Phi_ss.array().colwise() * n.col(1);
		Eigen::VectorXd t_x = t.col(0);
		Eigen::VectorXd t_y = t.col(1);
		Phi.block<b,b>(2*a+b,2*a) = MatrixXd(t_x.asDiagonal());
		Phi.block<b,b>(2*a+b,2*a+b) = MatrixXd(t_y.asDiagonal());
	//	cout << Phi << endl;
	//	cout << Phi << endl;
	//	Eigen::VectorXd alpha = Phi.ldlt().solve(defVec);
		Eigen::VectorXd alpha = Phi.partialPivLu().solve(defVec);
//		cout << m.intNodes << endl;
		Eigen::MatrixXd Phi_im = getPhi(m.intNodes,dispNodes);
		Eigen::MatrixXd Phi_is = getPhi(m.intNodes,m.slidingNodes);

		auto dx_int = Phi_im*alpha(Eigen::seq(0,dispNodes.size()-1)) + Phi_is*alpha(Eigen::seq(2*dispNodes.size(), 2*dispNodes.size()+m.slidingNodes.size()-1));
		auto dx_slide = Phi_sm*alpha(Eigen::seq(0,dispNodes.size()-1)) + Phi_ss*alpha(Eigen::seq(2*dispNodes.size(), 2*dispNodes.size()+m.slidingNodes.size()-1));
		auto dy_int = Phi_im*alpha(Eigen::seq(dispNodes.size(),2*dispNodes.size()-1)) + Phi_is*alpha(Eigen::seq(2*dispNodes.size()+m.slidingNodes.size(), 2*(dispNodes.size()+m.slidingNodes.size())-1));
		auto dy_slide = Phi_sm*alpha(Eigen::seq(dispNodes.size(),2*dispNodes.size()-1)) + Phi_ss*alpha(Eigen::seq(2*dispNodes.size()+m.slidingNodes.size(), 2*(dispNodes.size()+m.slidingNodes.size())-1));

		newCoords(m.slidingNodes,0)+=dx_slide;
		newCoords(m.slidingNodes,1)+=dy_slide;
		newCoords(m.intNodes,0)+=dx_int;
		newCoords(m.intNodes,1)+=dy_int;
		newCoords(m.intBdryNodes,0)+=defMat.col(0);
		newCoords(m.intBdryNodes,1)+=defMat.col(1);

		rotPnt(0) = rotPnt(0) + xDef/steps;
		rotPnt(1) = rotPnt(1) + yDef/steps;
	}
	m.writeMeshFile(newCoords);

}

Eigen::MatrixXd rbf::getDefVecDS(double xDef, double yDef, double rotDefDeg, Eigen::VectorXd rotPnt,Eigen::ArrayXi intN){
	Eigen::MatrixXd defMat = Eigen::MatrixXd::Zero(intN.size(),m.nDims); 		// initialise zero vector with size of boundary nodes
	Eigen::MatrixXd rotDef = getRotDef(rotDefDeg,rotPnt); //look at the passing by reference


	defMat.array().col(0) = defMat.array().col(0) + xDef;
	defMat.array().col(1) = defMat.array().col(1) + yDef;

	defMat.array() = defMat.array() + rotDef.array();

//	std::exit(0);
//	int idx = 0;
//	for(int i=0; i<m.intBdryNodes.size(); i++){
//		idx = idx+distance(m.bdryNodes.begin()+idx,find(m.bdryNodes.begin()+idx, m.bdryNodes.end(),m.intBdryNodes(i)));
//		defMat(idx,0) = xDef + rotDef(i,0);
//		defMat(idx,1) = yDef + rotDef(i,1);
//	}
	return defMat;
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
