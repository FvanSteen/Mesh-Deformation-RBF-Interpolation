#include "Mesh.h"
#include "rbf.h"
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <string>

rbf::rbf(Mesh& meshOb, const double xDef, const double yDef,const double zDef, Eigen::VectorXd& rVec, const int steps, Eigen::RowVectorXd& rotationPnt, const std::string& slidingMode)
:m(meshOb), dx(xDef/steps), dy(yDef/steps), dz(zDef/steps), rotPnt(rotationPnt), steps(steps), smode(slidingMode), rotVec(rVec)
{
	// todo specifying the periodic direction(s)
	// todo smode is already in the mesh class -> replace all instances with m.smode
	std::cout << "Current sliding mode: " << smode << std::endl;
	if(smode=="none"){
		mNodes.resize(m.N_ib+m.N_eb);
		mNodes << m.intBdryNodes, m.extBdryNodes;
	}else if(smode=="ds"){
		mNodes.resize(m.N_ib+m.N_es);
		mNodes << m.intBdryNodes, m.extStaticNodes;

	}else if(smode=="ps"){
		mNodes.resize(m.N_ib + m.N_es + m.N_se);
		mNodes << m.intBdryNodes, m.extStaticNodes, m.slidingEdgeNodes;

		if(m.pmode == "fixed"){
			mNodesPro.resize(m.N_ib+m.N_es);
			mNodesPro << m.intBdryNodes, m.extStaticNodes;
			N_mPro = m.N_ib+m.N_es;
		}
		else if(m.pmode == "moving"){
			mNodesPro.resize(m.N_ib);
			mNodesPro << m.intBdryNodes;
			N_mPro = m.N_ib;
			sNodes.resize(m.N_se + m.N_es);
			sNodes << m.slidingEdgeNodes, m.extStaticNodes;
			N_s = sNodes.size();
		}

		iNodes.resize(m.N_i + m.N_p);
		iNodes << m.intNodes, m.periodicNodes;
	}
	N_m = mNodes.size();

	// todo an vector with displacements could be passed as argument straight from ClassicalRBF.cpp
	dVec.resize(m.nDims);
	if(m.nDims==2){
		dVec << dx,dy;
	}
	else if(m.nDims==3){
		dVec << dx,dy,dz;
	}

	getRotationalMat();

}

void rbf::RBFMain(){
	if(smode=="none"){
		RBF_standard();
	}
	else if(smode=="ds"){
		if(m.nDims ==3){
			RBF_DS_3D();
		}else if(m.nDims ==2){
			RBF_DS();
		}
	}
	else if(smode=="ps"){
		RBF_PS();
	}
}

void rbf::RBF_PS(){
	std::cout<< "Performing RBF PS" << std::endl;


	Eigen::MatrixXd Phi_mmPro, Phi_sm, Phi_mm, Phi_im; 	//In this case only the internal boundary and static ext bdry
	Eigen::VectorXd defVecPro, defVec;
//	Eigen::ArrayXXd delta(m.N_se, m.nDims), finalDef(m.N_se,m.nDims);
	Eigen::ArrayXXd delta(N_s, m.nDims), finalDef(m.N_se,m.nDims);

	for(int i = 0; i < steps; i++){
		std::cout << "Deformation step: " << i+1 << std::endl;

		getPhi(Phi_mmPro, mNodesPro, mNodesPro);
//		getPhi(Phi_sm, m.slidingEdgeNodes, mNodesPro);
		getPhi(Phi_sm, sNodes, mNodesPro);
		getPhi(Phi_mm, mNodes, mNodes);
//		getPhi(Phi_im, m.intNodes,mNodes);
		getPhi(Phi_im, iNodes, mNodes);

		defVecPro = Eigen::VectorXd::Zero(N_mPro*m.nDims);
		getDefVec(defVecPro, N_mPro);

		performRBF_PS(Phi_mmPro, Phi_sm, Phi_mm, Phi_im, defVecPro, delta, finalDef, defVec);
//		if(i==1){
//			m.writeMeshFile();
//			std::exit(0);
//		}
	}

	m.writeMeshFile();
}

void rbf::performRBF_PS(Eigen::MatrixXd& Phi_mmPro, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_mm, Eigen::MatrixXd& Phi_im, Eigen::VectorXd& defVecPro,Eigen::ArrayXXd& delta, Eigen::ArrayXXd& finalDef, Eigen::VectorXd& defVec){

	for(int dim = 0; dim < m.nDims; dim++){
		delta.col(dim) = (Phi_sm*(Phi_mmPro.llt().solve(defVecPro(Eigen::seqN(dim*N_mPro,N_mPro))))).array();
	}


//	m.coords(m.slidingEdgeNodes,Eigen::all) = m.coords(m.slidingEdgeNodes,Eigen::all) + delta;
	m.coords(sNodes,Eigen::all) = m.coords(sNodes,Eigen::all) + delta;
	for(int x =0; x < m.nDims;x++ ){
		m.coords(mNodesPro,x) = m.coords(mNodesPro,x) + defVecPro(Eigen::seqN(x*N_mPro,N_mPro)).array();
	}
//	std::cout << m.coords << std::endl;

	m.getVecs();
	for(int i=0; i<m.N_se; i++){
		finalDef.row(i) = (delta.row(i)).matrix().transpose() - (delta.row(i).matrix()).dot((m.n.row(i)).matrix().transpose())*(m.n.row(i)).matrix().transpose();
	}
	std::cout << finalDef.rows() << std::endl;
//	std::cout << finalDef << std::endl;
//	m.coords(m.slidingEdgeNodes,Eigen::all) = m.coords(m.slidingEdgeNodes,Eigen::all) + finalDef;
	m.writeMeshFile();
	std::exit(0);

	defVec = Eigen::VectorXd::Zero(N_m*m.nDims);
	for(int dim = 0; dim< m.nDims; dim++){
		defVec(Eigen::seqN(dim*N_m,N_mPro)) = defVecPro(Eigen::seqN(dim*N_mPro,N_mPro));
		defVec(Eigen::seqN(dim*N_m+N_mPro,m.N_se)) = finalDef.col(dim);
	}

	performRBF(Phi_mm,Phi_im,defVec);
}

void rbf::RBF_DS_3D(){
	std::cout << "Performing RBF DS in three dimensions" << std::endl;

	Eigen::MatrixXd Phi_mm, Phi_ms, Phi_sm, Phi_ss, Phi_im, Phi_is, Phi;
//	Eigen::ArrayXXd n(m.N_se, m.nDims), t(m.N_se, m.nDims);		// two column array containing normal vector components
	Eigen::VectorXd defVec, alpha(m.nDims*(N_m+m.N_se+m.N_ss));
//	Eigen::ArrayXXd t_se(m.N_se,m.nDims),n1_se(m.N_se,m.nDims), n2_se(m.N_se,m.nDims), n_ss(m.N_ss, m.nDims),t1_ss(m.N_ss, m.nDims),t2_ss(m.N_ss, m.nDims);
	Eigen::MatrixXd Phi_me, Phi_ee, Phi_es, Phi_em, Phi_se, Phi_ie;
	for (int i = 0; i < steps; i++){
		std::cout << "Deformation step: " << i+1 << std::endl;
		getPhi(Phi_mm, mNodes,mNodes);
		getPhi(Phi_me, mNodes, m.slidingEdgeNodes);
		getPhi(Phi_ms, mNodes, m.slidingSurfNodes);

		getPhi(Phi_em, m.slidingEdgeNodes, mNodes);
		getPhi(Phi_ee, m.slidingEdgeNodes, m.slidingEdgeNodes);
		getPhi(Phi_es, m.slidingEdgeNodes, m.slidingSurfNodes);

		getPhi(Phi_sm, m.slidingSurfNodes, mNodes);
		getPhi(Phi_se, m.slidingSurfNodes, m.slidingEdgeNodes);
		getPhi(Phi_ss, m.slidingSurfNodes, m.slidingSurfNodes);
		getPhi(Phi_im, m.intNodes, mNodes);
		getPhi(Phi_ie, m.intNodes, m.slidingEdgeNodes);
		getPhi(Phi_is, m.intNodes, m.slidingSurfNodes);




		std::cout << "Reached the point at which normal and tangential vectors should be found" << std::endl;
		m.getVecs();
//		m.getVecs3D(t_se, n1_se, n2_se, n_ss,t1_ss,t2_ss); // function that obtains all the normal and tan vectors


		defVec = Eigen::VectorXd::Zero((N_m+m.N_se+m.N_ss)*m.nDims);
		getDefVec(defVec,N_m);
//		std::cout << defVec << std::endl;

		getPhiDS_3D(Phi,Phi_mm, Phi_me, Phi_ms, Phi_em, Phi_ee, Phi_es, Phi_sm, Phi_se, Phi_ss);
//		std::cout << Phi << std::endl;
//		std::cout << "DONE" << std::endl;
//		std::exit(0);
		performRBF_DS_3D(Phi, Phi_im, Phi_ie, Phi_is, Phi_em, Phi_ee, Phi_es, Phi_sm, Phi_se, Phi_ss,defVec, alpha);
	}

	m.writeMeshFile();
}

void rbf::performRBF_DS_3D(Eigen::MatrixXd& Phi, Eigen::MatrixXd& Phi_im, Eigen::MatrixXd& Phi_ie, Eigen::MatrixXd& Phi_is, Eigen::MatrixXd& Phi_em, Eigen::MatrixXd& Phi_ee, Eigen::MatrixXd& Phi_es, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_se, Eigen::MatrixXd& Phi_ss, Eigen::VectorXd& defVec, Eigen::VectorXd& alpha){

	alpha = Phi.partialPivLu().solve(defVec);
	std::cout << "SOLUTION WAS FOUND" << std::endl;

	for (int dim = 0; dim < m.nDims; dim++){

		m.coords(m.intNodes, dim) += (Phi_im*alpha(Eigen::seqN(dim*(N_m+m.N_se+m.N_ss),N_m)) + Phi_ie*alpha(Eigen::seqN(dim*(N_m+m.N_se+m.N_ss)+N_m, m.N_se)) + Phi_is*alpha(Eigen::seqN(dim*(N_m+m.N_se+m.N_ss)+N_m+m.N_se, m.N_ss)) ).array();
		m.coords(m.slidingEdgeNodes, dim) += (Phi_em*alpha(Eigen::seqN(dim*(N_m+m.N_se+m.N_ss),N_m)) + Phi_ee*alpha(Eigen::seqN(dim*(N_m+m.N_se+m.N_ss)+N_m, m.N_se)) + Phi_es*alpha(Eigen::seqN(dim*(N_m+m.N_se+m.N_ss)+N_m+m.N_se, m.N_ss)) ).array();
		m.coords(m.slidingSurfNodes, dim) += (Phi_sm*alpha(Eigen::seqN(dim*(N_m+m.N_se+m.N_ss),N_m)) + Phi_se*alpha(Eigen::seqN(dim*(N_m+m.N_se+m.N_ss)+N_m, m.N_se)) + Phi_ss*alpha(Eigen::seqN(dim*(N_m+m.N_se+m.N_ss)+N_m+m.N_se, m.N_ss)) ).array();
		m.coords(m.intBdryNodes, dim) += (defVec(Eigen::seqN(dim*N_m,m.N_ib))).array();
		rotPnt(dim) += dVec(dim);
	}
}

void rbf::getPhiDS_3D(Eigen::MatrixXd& Phi, Eigen::MatrixXd& Phi_mm, Eigen::MatrixXd&  Phi_me, Eigen::MatrixXd&  Phi_ms, Eigen::MatrixXd& Phi_em, Eigen::MatrixXd& Phi_ee, Eigen::MatrixXd& Phi_es, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_se, Eigen::MatrixXd& Phi_ss){

	Phi = Eigen::MatrixXd::Zero(m.nDims*(N_m+m.N_se+m.N_ss),m.nDims*(N_m+m.N_se+m.N_ss));

	for(int dim = 0; dim< m.nDims; dim++){
		// blocks related to the known displacements
		Phi.block(dim*N_m, dim*(N_m+m.N_se+m.N_ss), N_m, N_m) = Phi_mm;
		Phi.block(dim*N_m, dim*(N_m+m.N_se+m.N_ss)+N_m, N_m, m.N_se) = Phi_me;
		Phi.block(dim*N_m, dim*(N_m+m.N_se+m.N_ss)+N_m+m.N_se, N_m, m.N_ss) = Phi_ms;

		//blocks realteð to the first zero normal displacement condition of the sliding edge nodes
		Phi.block(3*N_m, dim*(N_m+m.N_se+m.N_ss),					m.N_se, N_m ) = Phi_em.array().colwise() * m.n1_se.col(dim);
		Phi.block(3*N_m, dim*(N_m+m.N_se+m.N_ss) + N_m,				m.N_se, m.N_se ) = Phi_ee.array().colwise() * m.n1_se.col(dim);
		Phi.block(3*N_m, dim*(N_m+m.N_se+m.N_ss) + N_m + m.N_se,	m.N_se,	m.N_ss ) = Phi_es.array().colwise() * m.n1_se.col(dim);
		//blocks realteð to the second zero normal displacement condition of the sliding edge nodes
		Phi.block(3*N_m+m.N_se, dim*(N_m+m.N_se+m.N_ss),					m.N_se, N_m ) = Phi_em.array().colwise() * m.n2_se.col(dim);
		Phi.block(3*N_m+m.N_se, dim*(N_m+m.N_se+m.N_ss) + N_m,				m.N_se, m.N_se ) = Phi_ee.array().colwise() * m.n2_se.col(dim);
		Phi.block(3*N_m+m.N_se, dim*(N_m+m.N_se+m.N_ss) + N_m + m.N_se,		m.N_se,	m.N_ss ) = Phi_es.array().colwise() * m.n2_se.col(dim);
		// blocks related to the zero tangential contribution of the sliding edge nodes
		Phi.block(3*N_m+2*m.N_se, dim*(N_m+m.N_se+m.N_ss) + N_m,	m.N_se, m.N_se) = Eigen::MatrixXd(m.t_se.col(dim).matrix().asDiagonal());
		// blocks related to the zero normal displacement condition of the sliding surface nodes
		Phi.block(3*N_m+3*m.N_se, dim*(N_m+m.N_se+m.N_ss),		m.N_ss,N_m) = Phi_sm.array().colwise() * m.n_ss.col(dim);
		Phi.block(3*N_m+3*m.N_se, dim*(N_m+m.N_se+m.N_ss)+N_m,	m.N_ss,m.N_se) = Phi_se.array().colwise() * m.n_ss.col(dim);
		Phi.block(3*N_m+3*m.N_se, dim*(N_m+m.N_se+m.N_ss)+N_m+m.N_se,	m.N_ss,m.N_ss) = Phi_ss.array().colwise() * m.n_ss.col(dim);

		// blocks related to the zero tangential contribution of the sliding surface nodes.
		Phi.block(3*N_m+3*m.N_se+m.N_ss, dim*(N_m+m.N_se+m.N_ss)+N_m+m.N_se, m.N_ss,m.N_ss) = Eigen::MatrixXd(m.t1_ss.col(dim).matrix().asDiagonal());
		Phi.block(3*N_m+3*m.N_se+2*m.N_ss, dim*(N_m+m.N_se+m.N_ss)+N_m+m.N_se, m.N_ss,m.N_ss) = Eigen::MatrixXd(m.t2_ss.col(dim).matrix().asDiagonal());
	}


}


void rbf::RBF_DS(){
	std::cout << "Performing RBF DS " << std::endl;

	Eigen::MatrixXd Phi_mm, Phi_ms, Phi_sm, Phi_ss, Phi_im, Phi_is, Phi;
//	Eigen::ArrayXXd n(m.N_se, m.nDims), t(m.N_se, m.nDims);		// two column array containing normal vector components
	Eigen::VectorXd defVec, alpha(m.nDims*(N_m+m.N_se));

	for (int i = 0; i < steps; i++){
		std::cout << "Deformation step: " << i+1 << std::endl;
		getPhi(Phi_mm, mNodes,mNodes);
		getPhi(Phi_ms, mNodes, m.slidingEdgeNodes);
		getPhi(Phi_sm, m.slidingEdgeNodes, mNodes);
		getPhi(Phi_ss, m.slidingEdgeNodes, m.slidingEdgeNodes);
		getPhi(Phi_im, m.intNodes, mNodes);
		getPhi(Phi_is, m.intNodes, m.slidingEdgeNodes);



		m.getVecs();


//		m.getExtBdryData();
//		m.getNodeVecs(n,t);

		defVec = Eigen::VectorXd::Zero((N_m+m.N_se)*m.nDims);
		getDefVec(defVec,N_m);

		getPhiDS(Phi,Phi_mm,Phi_ms, Phi_sm, Phi_ss, m.n, m.t);

		performRBF_DS(Phi, Phi_im, Phi_is, Phi_sm, Phi_ss,defVec, alpha);
	}

	m.writeMeshFile();
}

void rbf::performRBF_DS(Eigen::MatrixXd& Phi, Eigen::MatrixXd& Phi_im, Eigen::MatrixXd& Phi_is, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_ss, Eigen::VectorXd& defVec, Eigen::VectorXd& alpha){

	alpha = Phi.householderQr().solve(defVec);

	for (int dim = 0; dim < m.nDims; dim++){
		m.coords(m.intNodes, dim) += (Phi_im*alpha(Eigen::seqN(dim*(N_m+m.N_se),N_m)) + Phi_is*alpha(Eigen::seqN(dim*(N_m+m.N_se)+N_m, m.N_se))).array();
		m.coords(m.slidingEdgeNodes, dim) += (Phi_sm*alpha(Eigen::seqN(dim*(N_m+m.N_se),N_m)) + Phi_ss*alpha(Eigen::seqN(dim*(N_m+m.N_se)+N_m, m.N_se))).array();
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
		std::cout << defVec << std::endl;
		std::cout << "Performing RBF" << std::endl;
		performRBF(Phi_mm, Phi_im, defVec);
	}
	m.writeMeshFile();
}

void rbf::performRBF(Eigen::MatrixXd& Phi_mm, Eigen::MatrixXd& Phi_im, Eigen::VectorXd& defVec){
	for(int dim = 0; dim < m.nDims; dim++){
		std::cout << "Solving for dimension: " << dim << std::endl;
//		m.coords(m.intNodes,dim) += (Phi_im*(Phi_mm.householderQr().solve(defVec(Eigen::seqN(dim*N_m,N_m))))).array();
		m.coords(iNodes,dim) += (Phi_im*(Phi_mm.householderQr().solve(defVec(Eigen::seqN(dim*N_m,N_m))))).array();
		m.coords(mNodes,dim) += defVec(Eigen::seqN(dim*N_m,N_m)).array();
		rotPnt(dim) += dVec(dim);
	}
//	m.writeMeshFile();
//	std::exit(0);

}


void rbf::getPhi(Eigen::MatrixXd& Phi, Eigen::ArrayXi& idxSet1, Eigen::ArrayXi& idxSet2){
	Phi.resize(idxSet1.size(), idxSet2.size());
	double dist;
	for(int i=0; i<idxSet1.size();i++){
		for(int j=0; j<idxSet2.size();j++){
			if(m.nDims == 2){
				// Euclidian distance

//				dist = sqrt(pow(m.coords(idxSet1(i),0)-m.coords(idxSet2(j),0),2) + pow(m.coords(idxSet1(i),1)-m.coords(idxSet2(j),1),2));

				dist = sqrt(pow(m.coords(idxSet1(i),0)-m.coords(idxSet2(j),0),2) + pow(1/M_PI*sin( (m.coords(idxSet1(i),1)-m.coords(idxSet2(j),1))*M_PI/1),2));
//				std::cout << dist << std::endl;
//				std::cout << 'x' << '\t' << m.coords(idxSet1(i),0)-m.coords(idxSet2(j),0) << std::endl;
//				std::cout << 'y' << '\t' << m.coords(idxSet1(i),1)-m.coords(idxSet2(j),1) << std::endl;

			}
			else if(m.nDims == 3){
				dist = sqrt(pow(m.coords(idxSet1(i),0)-m.coords(idxSet2(j),0),2) + pow(m.coords(idxSet1(i),1)-m.coords(idxSet2(j),1),2) + pow(m.coords(idxSet1(i),2)-m.coords(idxSet2(j),2),2));
			}
			Phi(i,j) = pow((1-(dist/m.r)),4)*(4*(dist/m.r)+1);
		}
	}
//	std::exit(0);
}

void rbf::getDefVec(Eigen::VectorXd& defVec, int& N){
	Eigen::MatrixXd intPnts(m.N_ib,m.nDims);
	Eigen::MatrixXd rotDef;

	intPnts = m.coords(m.intBdryNodes,Eigen::all);
//	std::cout << intPnts.rowwise() - rotPnt << std::endl;
//	std::cout << (rotMatZ*(intPnts.rowwise() - rotPnt).transpose()).transpose() << std::endl;

	if(m.nDims == 2){
		rotDef = (rotMat*(intPnts.rowwise() - rotPnt).transpose()).transpose().rowwise() +rotPnt - intPnts;
	}
	else if(m.nDims == 3){
		rotDef = Eigen::MatrixXd::Zero(m.N_ib,m.nDims);

		if(rotVec[0] != 0){
			rotDef+= (rotMatX*(intPnts.rowwise() - rotPnt).transpose()).transpose().rowwise() +rotPnt - intPnts;
		}
		if(rotVec[1] != 0){
			rotDef+= (rotMatY*(intPnts.rowwise() - rotPnt).transpose()).transpose().rowwise() +rotPnt - intPnts;
		}
		if(rotVec[2] != 0){
			rotDef+= (rotMatZ*(intPnts.rowwise() - rotPnt).transpose()).transpose().rowwise() +rotPnt - intPnts;
		}
	}

	for(int dim = 0; dim < m.nDims; dim++){
		defVec(Eigen::seqN(dim*N, m.N_ib)).array() += dVec(dim);
		defVec(Eigen::seqN(dim*N, m.N_ib)) += rotDef.col(dim);
	}
 }

void rbf::getRotationalMat(){
	if(m.nDims == 2){
		const double theta = rotVec[0]/steps*M_PI/180;
		rotMat << 	cos(theta), -sin(theta),
					sin(theta),	cos(theta);
	}
	if(m.nDims == 3){
		const double x_theta = rotVec[0]/steps*M_PI/180;
		const double y_theta = rotVec[1]/steps*M_PI/180;
		const double z_theta = rotVec[2]/steps*M_PI/180;
		rotMatX << 	1,	0,	0,
					0,	cos(x_theta),	-sin(x_theta),
					0,	sin(x_theta),	cos(x_theta);
		rotMatY << 	cos(y_theta),	0,	sin(y_theta),
					0,	1,	0,
					-sin(y_theta),	0,	cos(y_theta);
		rotMatZ <<	cos(z_theta),	-sin(z_theta),	 0,
					sin(z_theta),	cos(z_theta),	 0,
					0,	0,	1;

	}
}

//void rbf::rbfEval(double distance){
////	double xi = distance/m.r;	// distance scaled by support radius
//	f_xi = pow((1-(distance/m.r)),4)*(4*(distance/m.r)+1);
////	return f_xi;
//}
