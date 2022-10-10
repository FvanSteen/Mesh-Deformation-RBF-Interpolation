#include "Mesh.h"
#include "rbf.h"
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <string>

rbf::rbf(Mesh& meshOb, const double xDef, const double yDef, const double rotDefDeg, const int steps, Eigen::RowVectorXd rotationPnt)
:m(meshOb), dx(xDef/steps), dy(yDef/steps), rotPnt(rotationPnt), steps(steps)
{
const double dthetaRad = rotDefDeg/steps*M_PI/180;
rotMat << cos(dthetaRad), -sin(dthetaRad),sin(dthetaRad), cos(dthetaRad);
mNodes.resize(m.intBdryNodes.size() + m.movingNodes.size());
mNodes <<  m.intBdryNodes,m.movingNodes; // idealy these are only defined in case of DS

}

void rbf::performRbfDS(){
	std::cout << "Performing RBF DS " << std::endl;

	Eigen::MatrixXd Phi_mm;
	Eigen::MatrixXd Phi_ms;
	Eigen::MatrixXd Phi_sm;
	Eigen::MatrixXd Phi_ss;


	Eigen::ArrayXXd n;		// two column array containing normal vector components
	Eigen::ArrayXXd t;		// two column array containing tangential vector components


	for (int i = 0; i < steps; i++){
		getPhi(Phi_mm, mNodes,mNodes);
		getPhi(Phi_ms, mNodes, m.slidingNodes);
		getPhi(Phi_sm, m.slidingNodes, mNodes);
		getPhi(Phi_ss, m.slidingNodes, m.slidingNodes);

		m.getExtBdryData();

		m.getNodeVecs(m.slidingNodes, n, t);
	}

//	Eigen::ArrayXXd n(m.slidingNodes.size(),2);		// all the vectors midPnts etc should be updated after each step
//	Eigen::ArrayXXd t(m.slidingNodes.size(),2);
//	m.getExtBdryData();
//
//	m.getNodeVecs(m.slidingNodes, n, t);
//
//	exit(0);
////
////	cout << "\n Sliding Nodes:\n" << m.slidingNodes << endl;
////	cout << "\n Moving Nodes: \n" << m.movingNodes << endl;
////	cout << "\n Internal boundary nodes: \n" << m.intBdryNodes << endl;
////	cout << "\n Internal Nodes: \n" << m.intNodes << endl;
//
//	Eigen::ArrayXi dispNodes(m.intBdryNodes.size()+m.movingNodes.size());
//	dispNodes << m.intBdryNodes, m.movingNodes;
//
////	cout << dispNodes << endl;
//	for(int i=0; i<steps; i++){
//
//		cout << "step: " << i << endl;
//		Eigen::MatrixXd defMat = getDefVecDS(xDef/steps,yDef/steps,rotDefDeg/steps,rotPnt,m.intBdryNodes);
//
//		Eigen::MatrixXd Phi_mm = getPhi(dispNodes,dispNodes); // makes i-matrix for finding coefficients
//		Eigen::MatrixXd Phi_ms = getPhi(dispNodes,m.slidingNodes);
//		Eigen::MatrixXd Phi_sm = getPhi(m.slidingNodes,dispNodes);
//		Eigen::MatrixXd Phi_ss = getPhi(m.slidingNodes,m.slidingNodes);
//
////		cout << Phi_mm.rows() << '\t' << Phi_mm.cols() << endl;
////		cout << Phi_ms.rows() << '\t' << Phi_ms.cols() << endl;
////		cout << Phi_sm.rows() << '\t' << Phi_sm.cols() << endl;
////		cout << Phi_ss.rows() << '\t' << Phi_ss.cols() << endl;
//
//		Eigen::VectorXd defVec = Eigen::VectorXd::Zero(2*(dispNodes.size()+m.slidingNodes.size()));
//		defVec(Eigen::seq(0,m.intBdryNodes.size()-1)) = defMat.col(0);
//		defVec(Eigen::seq(dispNodes.size(),dispNodes.size()+m.intBdryNodes.size()-1)) = defMat.col(1);
//
//		Eigen::MatrixXd Phi = Eigen::MatrixXd::Zero(2*(dispNodes.size()+m.slidingNodes.size()),2*(dispNodes.size()+m.slidingNodes.size()));
////		cout << "here" << endl;
//	//	cout << Phi_ss << endl;
//		const int a = 8; // find some way to produce const int of the various dimensions within the mesh class
//		const int b = 16;
////		const int a = 16; // find some way to produce const int of the various dimensions within the mesh class
////		const int b = 96;
//
////		const int a = 200; // find some way to produce const int of the various dimensions within the mesh class
////		const int b = 50;
//
//		Phi.block<a,a>(0,0) = Phi_mm;
//		Phi.block<a,a>(a,a) = Phi_mm;
//		Phi.block<a,b>(0,2*a) = Phi_ms;
//		Phi.block<a,b>(a,2*a+b) = Phi_ms;
//	//	cout << Phi << endl;
//
//	//	cout << Phi_sm << endl;
//	//	cout << n.col(0) << endl;
//	//	cout << Phi_sm.array().colwise() * n.col(0) << endl;
//
//		Phi.block<b,a>(2*a,0) = Phi_sm.array().colwise() * n.col(0);
//		Phi.block<b,a>(2*a,a) = Phi_sm.array().colwise() * n.col(1);
//		Phi.block<b,b>(2*a,2*a) = Phi_ss.array().colwise() * n.col(0);
//		Phi.block<b,b>(2*a,2*a+b) = Phi_ss.array().colwise() * n.col(1);
//		Eigen::VectorXd t_x = t.col(0);
//		Eigen::VectorXd t_y = t.col(1);
//		Phi.block<b,b>(2*a+b,2*a) = MatrixXd(t_x.asDiagonal());
//		Phi.block<b,b>(2*a+b,2*a+b) = MatrixXd(t_y.asDiagonal());
//	//	cout << Phi << endl;
//	//	cout << Phi << endl;
//	//	Eigen::VectorXd alpha = Phi.ldlt().solve(defVec);
//		Eigen::VectorXd alpha = Phi.partialPivLu().solve(defVec);
////		cout << m.intNodes << endl;
//		Eigen::MatrixXd Phi_im = getPhi(m.intNodes,dispNodes);
//		Eigen::MatrixXd Phi_is = getPhi(m.intNodes,m.slidingNodes);
//
//		auto dx_int = Phi_im*alpha(Eigen::seq(0,dispNodes.size()-1)) + Phi_is*alpha(Eigen::seq(2*dispNodes.size(), 2*dispNodes.size()+m.slidingNodes.size()-1));
//		auto dx_slide = Phi_sm*alpha(Eigen::seq(0,dispNodes.size()-1)) + Phi_ss*alpha(Eigen::seq(2*dispNodes.size(), 2*dispNodes.size()+m.slidingNodes.size()-1));
//		auto dy_int = Phi_im*alpha(Eigen::seq(dispNodes.size(),2*dispNodes.size()-1)) + Phi_is*alpha(Eigen::seq(2*dispNodes.size()+m.slidingNodes.size(), 2*(dispNodes.size()+m.slidingNodes.size())-1));
//		auto dy_slide = Phi_sm*alpha(Eigen::seq(dispNodes.size(),2*dispNodes.size()-1)) + Phi_ss*alpha(Eigen::seq(2*dispNodes.size()+m.slidingNodes.size(), 2*(dispNodes.size()+m.slidingNodes.size())-1));
//
//		newCoords(m.slidingNodes,0)+=dx_slide;
//		newCoords(m.slidingNodes,1)+=dy_slide;
//		newCoords(m.intNodes,0)+=dx_int;
//		newCoords(m.intNodes,1)+=dy_int;
//		newCoords(m.intBdryNodes,0)+=defMat.col(0);
//		newCoords(m.intBdryNodes,1)+=defMat.col(1);
//
//		rotPnt(0) = rotPnt(0) + xDef/steps;
//		rotPnt(1) = rotPnt(1) + yDef/steps;
//	}
//	m.writeMeshFile(newCoords);

}

Eigen::MatrixXd rbf::getDefVecDS(double xDef, double yDef, double rotDefDeg, Eigen::VectorXd rotPnt,Eigen::ArrayXi intN){
	Eigen::MatrixXd defMat = Eigen::MatrixXd::Zero(intN.size(),m.nDims); 		// initialise zero vector with size of boundary nodes
	getRotDef(); //look at the passing by reference


	defMat.array().col(0) = defMat.array().col(0) + xDef;
	defMat.array().col(1) = defMat.array().col(1) + yDef;
	//rotDef is not a variable anymore since the function getrotDef was altered
//	defMat.array() = defMat.array() + rotDef.array();
	return defMat;
}


void rbf::performRbfInterpolation(){

//	newCoords = m.coords;
	std::cout << m.bdryNodes << std::endl;
	Eigen::MatrixXd Phi_bb;
	Eigen::MatrixXd Phi_ib;

	for(int i=0; i<steps; i++){
		Eigen::VectorXd defVec = Eigen::VectorXd::Zero(m.bdryNodes.size()*m.nDims);

		std::cout << "Building matrix Phi_bb" << std::endl;
		getPhi(Phi_bb, m.bdryNodes,m.bdryNodes); // makes i-matrix for finding coefficients

		std::cout << "Building matrix Phi_ib" << std::endl;
		getPhi(Phi_ib, m.intNodes,m.bdryNodes);
		std::cout << "Building deformation vectors" << std::endl;

		std::cout << "Deformation step: " << i+1 << std::endl;

		getDefVec(defVec);

		std::cout << "Solving for interpolation coefficients" << std::endl;

		Eigen::VectorXd alpha_x = Phi_bb.llt().solve(defVec(Eigen::seqN(0,m.bdryNodes.size())));
		Eigen::VectorXd alpha_y = Phi_bb.llt().solve(defVec(Eigen::seqN(m.bdryNodes.size(),m.bdryNodes.size())));
		getDisplacement(Phi_ib, alpha_x, alpha_y, defVec);

	}
	m.writeMeshFile();
}


void rbf::getDisplacement(Eigen::MatrixXd& Phi, Eigen::VectorXd& a_x, Eigen::VectorXd& a_y, Eigen::VectorXd& defVec){
	m.coords(m.intNodes, 0) += (Phi*a_x).array();
	m.coords(m.intNodes, 1) += (Phi*a_y).array();

	m.coords(m.bdryNodes, 0) += defVec(Eigen::seqN(0,m.bdryNodes.size())).array();
	m.coords(m.bdryNodes, 1) += defVec(Eigen::seqN(m.bdryNodes.size(),m.bdryNodes.size())).array();

	rotPnt(0) += dx;
	rotPnt(1) += dy;
}

void rbf::getPhi(Eigen::MatrixXd& Phi, Eigen::ArrayXi& idxSet1, Eigen::ArrayXi& idxSet2){
	Phi.resize(idxSet1.size(), idxSet2.size());
//	Eigen::MatrixXd Phi(idxSet1.size(), idxSet2.size());
	for(int i=0; i<idxSet1.size();i++){
		for(int j=0; j<idxSet2.size();j++){
			double dist = sqrt(pow(m.coords(idxSet1(i),0)-m.coords(idxSet2(j),0),2) + pow(m.coords(idxSet1(i),1)-m.coords(idxSet2(j),1),2));
			Phi(i,j) = rbfEval(dist);
		}
	}

}

void rbf::getDefVec(Eigen::VectorXd& defVec){
	Eigen::MatrixXd intPnts(m.intBdryNodes.size(),m.nDims);
	Eigen::MatrixXd rotDef;

	intPnts = m.coords(m.intBdryNodes,Eigen::all);
	rotDef = (rotMat*(intPnts.rowwise() - rotPnt).transpose()).transpose().rowwise() +rotPnt - intPnts;

	defVec(Eigen::seqN(0,m.intBdryNodes.size())).array() += dx;
	defVec(Eigen::seqN(0,m.intBdryNodes.size())) += rotDef.col(0);
	defVec(Eigen::seqN(m.bdryNodes.size(),m.intBdryNodes.size())).array() += dy;
	defVec(Eigen::seqN(m.bdryNodes.size(),m.intBdryNodes.size())) += rotDef.col(1);
 }


void rbf::getRotDef(){
	Eigen::MatrixXd intPnts(m.intBdryNodes.size(),m.nDims);
	intPnts = m.coords(m.intBdryNodes,Eigen::all);
	Eigen::MatrixXd rotDef;
	rotDef = (rotMat*(intPnts.rowwise() - rotPnt).transpose()).transpose().rowwise() +rotPnt - intPnts ;
}

double rbf::rbfEval(double distance){
	double xi = distance/m.r;	// distance scaled by support radius
	double f_xi = pow((1-xi),4)*(4*xi+1);
	return f_xi;
}

