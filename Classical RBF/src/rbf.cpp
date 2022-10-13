#include "Mesh.h"
#include "rbf.h"
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <string>
//todo check all instances of transpose() to check for aliasing effect
rbf::rbf(Mesh& meshOb, const double xDef, const double yDef, const double rotDefDeg, const int steps, Eigen::RowVectorXd rotationPnt, const std::string& mode)
:m(meshOb), dx(xDef/steps), dy(yDef/steps), rotPnt(rotationPnt), steps(steps), mode(mode)
{
const double dthetaRad = rotDefDeg/steps*M_PI/180;
rotMat << cos(dthetaRad), -sin(dthetaRad),sin(dthetaRad), cos(dthetaRad);
// TODO suggestion: modify the assignment of mNodes to allow either DS or classic by means of an if statement
// and perhaps introduce N_in and N_ib and N_eb
mNodes.resize(m.intBdryNodes.size() + m.movingNodes.size());
mNodes <<  m.intBdryNodes,m.movingNodes; // idealy these are only defined in case of DS
N_m = mNodes.size();
N_se = m.slidingNodes.size();
std::cout << mode << std::endl;
}



void rbf::performRbfPS(){
	std::cout<< "Performing RBF PS" << std::endl;
	Eigen::ArrayXXd slidingCoords;
	slidingCoords = m.coords(m.slidingNodes, Eigen::all);

	Eigen::MatrixXd Phi_mm; //In this case only the internal boundary
	Eigen::VectorXd a_x(N_m);
	Eigen::VectorXd a_y(N_m);

	Eigen::ArrayXXd n(N_se, m.nDims);		// two column array containing normal vector components
	Eigen::ArrayXXd t(N_se, m.nDims);		// two column array containing tangential vector components
	m.getExtBdryData();
	m.getNodeVecs(n,t);
	getPhi(Phi_mm, mNodes, mNodes);

	Eigen::VectorXd defVec = Eigen::VectorXd::Zero(N_m*m.nDims);
	getDefVec(defVec);

	a_x = Phi_mm.partialPivLu().solve(defVec(Eigen::seqN(0,N_m)));
	a_y = Phi_mm.partialPivLu().solve(defVec(Eigen::seqN(N_m, N_m)));


	Eigen::MatrixXd Phi_im;
//	Eigen::ArrayXi freeNodes(m.slidingNodes.size()+m.intNodes.size());
//	freeNodes << m.slidingNodes, m.intNodes;
//	Eigen::ArrayXi freeNodes(m.slidingNodes.size());
//	freeNodes << m.slidingNodes;

	getPhi(Phi_im, m.slidingNodes, mNodes);


	slidingCoords.col(0) += (Phi_im*a_x).array();
	slidingCoords.col(1) += (Phi_im*a_y).array();


//	m.coords(mNodes,0) += (defVec(Eigen::seqN(0,N_m))).array();
//	m.coords(mNodes,1) += (defVec(Eigen::seqN(N_m,N_m))).array();
//	std::cout << m.coords << std::endl;
	Eigen::ArrayXXd finalDef(m.slidingNodes.size(),m.nDims);
	for(int i=0; i<m.slidingNodes.size(); i++){
//		Eigen::VectorXd dx = m.coords.row(m.slidingNodes(i)) - initCoords.row(m.slidingNodes(i));
		Eigen::VectorXd dx = slidingCoords.row(i) - m.coords.row(m.slidingNodes(i));
		Eigen::VectorXd n_local = n.row(i);
		finalDef.row(i) = dx-dx.dot(n_local)*n_local;
	}
	std::cout << finalDef <<std::endl;
	m.coords(m.slidingNodes,Eigen::all) += finalDef;
//	m.coords(m.slidingNodes,Eigen::all) = slidingCoords;

	m.writeMeshFile();

	std::cout << "Done" << std::endl;
}

void rbf::performRbfDS(){
	std::cout << "Performing RBF DS " << std::endl;

	Eigen::MatrixXd Phi_mm;
	Eigen::MatrixXd Phi_ms;
	Eigen::MatrixXd Phi_sm;
	Eigen::MatrixXd Phi_ss;
	Eigen::MatrixXd Phi_im;
	Eigen::MatrixXd Phi_is;


	Eigen::ArrayXXd n(N_se, m.nDims);		// two column array containing normal vector components
	Eigen::ArrayXXd t(N_se, m.nDims);		// two column array containing tangential vector components

	Eigen::MatrixXd Phi(2*(N_m+N_se),2*(N_m+N_se));

	Eigen::VectorXd alpha(2*(N_m+N_se));


	for (int i = 0; i < steps; i++){
		std::cout << "Deformation step: " << i+1 << std::endl;
		getPhi(Phi_mm, mNodes,mNodes);
		getPhi(Phi_ms, mNodes, m.slidingNodes);
		getPhi(Phi_sm, m.slidingNodes, mNodes);
		getPhi(Phi_ss, m.slidingNodes, m.slidingNodes);

		// obtaining midpoints and corresponding normals and tangentials
		m.getExtBdryData();

		// obtaining normals and tangentials at the sliding nodes as a weighted average of the midpoints
		m.getNodeVecs(n,t);
		Eigen::VectorXd defVec = Eigen::VectorXd::Zero((N_m+N_se)*m.nDims);
		getDefVec(defVec);

		getPhiDS(Phi,Phi_mm,Phi_ms, Phi_sm, Phi_ss, n, t);

		alpha = Phi.partialPivLu().solve(defVec);
		getPhi(Phi_im, m.intNodes, mNodes); // reduced matrix since a part has a zero deformation
		getPhi(Phi_is, m.intNodes, m.slidingNodes);

		getDisplacementDS(Phi_im, Phi_is, Phi_sm, Phi_ss, alpha, defVec);
	}

	m.writeMeshFile();
}


void rbf::getPhiDS(Eigen::MatrixXd& Phi,Eigen::MatrixXd& Phi_mm,Eigen::MatrixXd& Phi_ms, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_ss, Eigen::ArrayXXd& n, Eigen::ArrayXXd& t){
	// blocks related to the known displacements
	Phi.block(0, 0, N_m, N_m) = Phi_mm;
	Phi.block(N_m, N_m, N_m, N_m) = Phi_mm;
	Phi.block(0, 2*N_m, N_m, N_se) = Phi_ms;
	Phi.block(N_m, 2*N_m + N_se, N_m, N_se) = Phi_ms;

	// blocks related to the zero normal displacement condition
	Phi.block(2*N_m, 0, N_se, N_m) = Phi_sm.array().colwise() * n.col(0);
	Phi.block(2*N_m, N_m, N_se, N_m) = Phi_sm.array().colwise() * n.col(1);
	Phi.block(2*N_m, 2*N_m, N_se, N_se) = Phi_ss.array().colwise() * n.col(0);
	Phi.block(2*N_m, 2*N_m + N_se, N_se, N_se) = Phi_ss.array().colwise() * n.col(1);

	//blocks related to the zero tangential contribution condition
	Phi.block(2*N_m + N_se, 2*N_m, N_se, N_se) = Eigen::MatrixXd(t.col(0).matrix().asDiagonal());
	Phi.block(2*N_m + N_se, 2*N_m + N_se, N_se, N_se) = Eigen::MatrixXd(t.col(1).matrix().asDiagonal());
}

void rbf::getDisplacementDS(Eigen::MatrixXd& Phi_im, Eigen::MatrixXd& Phi_is,Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_ss, Eigen::VectorXd& alpha, Eigen::VectorXd& defVec){

	m.coords(m.intNodes, 0) += (Phi_im*alpha(Eigen::seqN(0,N_m)) + Phi_is*alpha(Eigen::seqN(2*N_m, N_se))).array();
	m.coords(m.intNodes, 1) += (Phi_im*alpha(Eigen::seqN(N_m,N_m)) + Phi_is*alpha(Eigen::seqN(2*N_m+N_se, N_se))).array();

	m.coords(m.slidingNodes, 0) += (Phi_sm*alpha(Eigen::seqN(0,N_m)) + Phi_ss*alpha(Eigen::seqN(2*N_m, N_se))).array();
	m.coords(m.slidingNodes, 1) += (Phi_sm*alpha(Eigen::seqN(N_m,N_m)) + Phi_ss*alpha(Eigen::seqN(2*N_m+N_se, N_se))).array();

	m.coords(m.intBdryNodes, 0) += (defVec(Eigen::seqN(0,m.intBdryNodes.size()))).array();
	m.coords(m.intBdryNodes, 1) += (defVec(Eigen::seqN(N_m,m.intBdryNodes.size()))).array();

	rotPnt(0) += dx;
	rotPnt(1) += dy;
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
		// todo these should be initialized outside the for loop
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

	// todo change expression in the classic rbf so that these statements work
	// in case of the classical RBF the defVec is constructed such that it first contains the xdef of the int bdry then
	// zero xdef of the ext bdry. Followed by the ydef of the int and ext bdries.
	defVec(Eigen::seqN(defVec.size()/m.nDims,m.intBdryNodes.size())).array() += dy;
	defVec(Eigen::seqN(defVec.size()/m.nDims,m.intBdryNodes.size())) += rotDef.col(1);

	// for the DS RBF the vector has length nDim*N_m. First it contains the int bdry xdef, static zero xdef of the ext bdry,
	// int bdry ydef and static ydef of the ext bdry.
//	defVec(Eigen::seqN(N_m,m.intBdryNodes.size())).array() += dy;
//	defVec(Eigen::seqN(N_m,m.intBdryNodes.size())) += rotDef.col(1);

	// for the pseudo only the xdef and ydef of the int bdry nodes are relevant for the first deformation step
	// This was without the zero displacement on the corners
//	defVec(Eigen::seqN(m.intBdryNodes.size(),m.intBdryNodes.size())).array() += dy;
//	defVec(Eigen::seqN(m.intBdryNodes.size(),m.intBdryNodes.size())) += rotDef.col(1);
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

