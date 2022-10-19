#include "Mesh.h"
#include "rbf.h"
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <string>

rbf::rbf(Mesh& meshOb, const double xDef, const double yDef,const double zDef, const double rotDefDeg, const int steps, Eigen::RowVectorXd& rotationPnt, const std::string& slidingMode)
:m(meshOb), dx(xDef/steps), dy(yDef/steps), dz(zDef/steps), rotPnt(rotationPnt), steps(steps), mode(slidingMode)
{
	const double dthetaRad = rotDefDeg/steps*M_PI/180;
	rotMat << cos(dthetaRad), -sin(dthetaRad),sin(dthetaRad), cos(dthetaRad);

	std::cout << "Current sliding mode: " << mode << std::endl;
	if(mode=="none"){
		mNodes.resize(m.N_ib+m.N_eb);
		mNodes << m.intBdryNodes, m.extBdryNodes;
	}else if(mode=="ds"){
		mNodes.resize(m.N_ib+m.N_es);
		mNodes << m.intBdryNodes, m.extStaticNodes;

	}else if(mode=="ps"){
		mNodes.resize(m.N_ib+m.N_eb);
		mNodes << m.intBdryNodes, m.extStaticNodes, m.slidingNodes;
		mNodesPro.resize(m.N_ib+m.N_es);
		mNodesPro << m.intBdryNodes, m.extStaticNodes;
		N_mPro = m.N_ib+m.N_es;
	}
	N_m = mNodes.size();

	// todo could be passed as argument straight from ClassicalRBF.cpp
	dVec.resize(m.nDims);
	if(m.nDims==2){
		dVec << dx,dy;
	}
	else if(m.nDims==3){
		dVec << dx,dy,dz;
	}
}

void rbf::RBFMain(){
	if(mode=="none"){
		RBF_standard();
	}
	else if(mode=="ds"){
		RBF_DS();
	}
	else if(mode=="ps"){
		RBF_PS();
	}
}

void rbf::RBF_PS(){
	std::cout<< "Performing RBF PS" << std::endl;


	Eigen::MatrixXd Phi_mmPro, Phi_sm, Phi_mm, Phi_im; 	//In this case only the internal boundary and static ext bdry
	Eigen::VectorXd defVecPro, defVec;
	Eigen::ArrayXXd delta(m.N_se, m.nDims), n(m.N_se, m.nDims), finalDef(m.N_se,m.nDims);

	for(int i = 0; i < steps; i++){
		std::cout << "Deformation step: " << i+1 << std::endl;

		getPhi(Phi_mmPro, mNodesPro, mNodesPro);
		getPhi(Phi_sm, m.slidingNodes, mNodesPro);
		getPhi(Phi_mm, mNodes, mNodes);
		getPhi(Phi_im, m.intNodes,mNodes);

		defVecPro = Eigen::VectorXd::Zero(N_mPro*m.nDims);
		getDefVec(defVecPro, N_mPro);

		performRBF_PS(Phi_mmPro, Phi_sm, Phi_mm, Phi_im, defVecPro, delta, n, finalDef, defVec);
	}

	m.writeMeshFile();
}

void rbf::performRBF_PS(Eigen::MatrixXd& Phi_mmPro, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_mm, Eigen::MatrixXd& Phi_im, Eigen::VectorXd& defVecPro,Eigen::ArrayXXd& delta,Eigen::ArrayXXd& n, Eigen::ArrayXXd& finalDef, Eigen::VectorXd& defVec){

	for(int dim = 0; dim < m.nDims; dim++){
		delta.col(dim) = (Phi_sm*(Phi_mmPro.partialPivLu().solve(defVecPro(Eigen::seqN(dim*N_mPro,N_mPro))))).array();
	}

	m.getExtBdryData();
	m.getNormals(n);

	for(int i=0; i<m.N_se; i++){
		finalDef.row(i) = (delta.row(i)).matrix().transpose() - (delta.row(i).matrix()).dot((n.row(i)).matrix().transpose())*(n.row(i)).matrix().transpose();
	}


	defVec = Eigen::VectorXd::Zero(N_m*m.nDims);
	for(int dim = 0; dim< m.nDims; dim++){
		defVec(Eigen::seqN(dim*N_m,N_mPro)) = defVecPro(Eigen::seqN(dim*N_mPro,N_mPro));
		defVec(Eigen::seqN(dim*N_m+N_mPro,m.N_se)) = finalDef.col(dim);
	}

	performRBF(Phi_mm,Phi_im,defVec);
}

void rbf::RBF_DS(){
	std::cout << "Performing RBF DS " << std::endl;

	Eigen::MatrixXd Phi_mm, Phi_ms, Phi_sm, Phi_ss, Phi_im, Phi_is, Phi;
	Eigen::ArrayXXd n(m.N_se, m.nDims), t(m.N_se, m.nDims);		// two column array containing normal vector components
	Eigen::VectorXd defVec, alpha(m.nDims*(N_m+m.N_se));

	for (int i = 0; i < steps; i++){
		std::cout << "Deformation step: " << i+1 << std::endl;
		getPhi(Phi_mm, mNodes,mNodes);
		getPhi(Phi_ms, mNodes, m.slidingNodes);
		getPhi(Phi_sm, m.slidingNodes, mNodes);
		getPhi(Phi_ss, m.slidingNodes, m.slidingNodes);
		getPhi(Phi_im, m.intNodes, mNodes);
		getPhi(Phi_is, m.intNodes, m.slidingNodes);


		m.getExtBdryData();
		m.getNodeVecs(n,t);

		defVec = Eigen::VectorXd::Zero((N_m+m.N_se)*m.nDims);
		getDefVec(defVec,N_m);

		getPhiDS(Phi,Phi_mm,Phi_ms, Phi_sm, Phi_ss, n, t);

		performRBF_DS(Phi, Phi_im, Phi_is, Phi_sm, Phi_ss,defVec, alpha);
	}

	m.writeMeshFile();
}

void rbf::performRBF_DS(Eigen::MatrixXd& Phi, Eigen::MatrixXd& Phi_im, Eigen::MatrixXd& Phi_is, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_ss, Eigen::VectorXd& defVec, Eigen::VectorXd& alpha){

	alpha = Phi.partialPivLu().solve(defVec);

	for (int dim = 0; dim < m.nDims; dim++){
		m.coords(m.intNodes, dim) += (Phi_im*alpha(Eigen::seqN(dim*(N_m+m.N_se),N_m)) + Phi_is*alpha(Eigen::seqN(dim*(N_m+m.N_se)+N_m, m.N_se))).array();
		m.coords(m.slidingNodes, dim) += (Phi_sm*alpha(Eigen::seqN(dim*(N_m+m.N_se),N_m)) + Phi_ss*alpha(Eigen::seqN(dim*(N_m+m.N_se)+N_m, m.N_se))).array();
		m.coords(m.intBdryNodes, dim) += (defVec(Eigen::seqN(dim*N_m,m.N_ib))).array();
		rotPnt(dim) += dVec(dim);
	}
}


void rbf::getPhiDS(Eigen::MatrixXd& Phi,Eigen::MatrixXd& Phi_mm,Eigen::MatrixXd& Phi_ms, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_ss, Eigen::ArrayXXd& n, Eigen::ArrayXXd& t){

	Phi = Eigen::MatrixXd::Zero(m.nDims*(N_m+m.N_se),m.nDims*(N_m+m.N_se));

	for(int dim = 0; dim< m.nDims; dim++){
		// blocks related to the known displacements
		Phi.block(dim*N_m, dim*(N_m+m.N_se), N_m, N_m) = Phi_mm;
		Phi.block(dim*N_m, dim*(N_m+m.N_se)+N_m, N_m, m.N_se) = Phi_ms;

		// blocks related to the zero normal displacement condition
		Phi.block(2*N_m, dim*(N_m+m.N_se), m.N_se, N_m) = Phi_sm.array().colwise() * n.col(dim);
		Phi.block(2*N_m, dim*(N_m+m.N_se)+N_m, m.N_se, m.N_se) = Phi_ss.array().colwise() * n.col(dim);

		//blocks related to the zero tangential contribution condition
		Phi.block(2*N_m + m.N_se, dim*(N_m+m.N_se)+N_m, m.N_se, m.N_se) = Eigen::MatrixXd(t.col(dim).matrix().asDiagonal());
	}
}



void rbf::RBF_standard(){
	Eigen::MatrixXd Phi_mm, Phi_im;
	Eigen::VectorXd defVec;
	for(int i=0; i<steps; i++){
		std::cout << "Deformation step: " << i+1 << std::endl;
		defVec = Eigen::VectorXd::Zero(N_m*m.nDims);
		std::cout << "Obtaining Phi_mm" << std::endl;
		getPhi(Phi_mm, mNodes,mNodes);
		std::cout << "Obtaining Phi_im" << std::endl;
		getPhi(Phi_im, m.intNodes,mNodes);
		std::cout << "Obtaining deformation vector" << std::endl;
		getDefVec(defVec,N_m);
		std::cout << "Performing RBF" << std::endl;
		performRBF(Phi_mm, Phi_im, defVec);
	}
	m.writeMeshFile();
}

void rbf::performRBF(Eigen::MatrixXd& Phi_mm, Eigen::MatrixXd& Phi_im, Eigen::VectorXd& defVec){
	for(int dim = 0; dim < m.nDims; dim++){
		std::cout << "Solving for dimension: " << dim << std::endl;
		m.coords(m.intNodes,dim) += (Phi_im*(Phi_mm.llt().solve(defVec(Eigen::seqN(dim*N_m,N_m))))).array();
		m.coords(mNodes,dim) += defVec(Eigen::seqN(dim*N_m,N_m)).array();
		rotPnt(dim) += dVec(dim);
	}
}


void rbf::getPhi(Eigen::MatrixXd& Phi, Eigen::ArrayXi& idxSet1, Eigen::ArrayXi& idxSet2){
	Phi.resize(idxSet1.size(), idxSet2.size());
	double dist;
	for(int i=0; i<idxSet1.size();i++){
		for(int j=0; j<idxSet2.size();j++){
			if(m.nDims == 2){
				dist = sqrt(pow(m.coords(idxSet1(i),0)-m.coords(idxSet2(j),0),2) + pow(m.coords(idxSet1(i),1)-m.coords(idxSet2(j),1),2));
			}
			else if(m.nDims == 3){
				dist = sqrt(pow(m.coords(idxSet1(i),0)-m.coords(idxSet2(j),0),2) + pow(m.coords(idxSet1(i),1)-m.coords(idxSet2(j),1),2) + pow(m.coords(idxSet1(i),2)-m.coords(idxSet2(j),2),2));
			}
			Phi(i,j) = pow((1-(dist/m.r)),4)*(4*(dist/m.r)+1);
		}
	}

}

void rbf::getDefVec(Eigen::VectorXd& defVec, int& N){
	Eigen::MatrixXd intPnts(m.N_ib,m.nDims);
	Eigen::MatrixXd rotDef;// not suitable for 3D probably

	intPnts = m.coords(m.intBdryNodes,Eigen::all);
//	rotDef = (rotMat*(intPnts.rowwise() - rotPnt).transpose()).transpose().rowwise() +rotPnt - intPnts;

	for(int dim = 0; dim < m.nDims; dim++){
		defVec(Eigen::seqN(dim*N, m.N_ib)).array() += dVec(dim);
//		defVec(Eigen::seqN(dim*N, m.N_ib)) += rotDef.col(dim);
	}
 }

//void rbf::rbfEval(double distance){
////	double xi = distance/m.r;	// distance scaled by support radius
//	f_xi = pow((1-(distance/m.r)),4)*(4*(distance/m.r)+1);
////	return f_xi;
//}
