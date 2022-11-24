#include "getNodeType.h"
#include "rbf.h"

#include "rbfstd.h"
#include "rbfps.h"
//#include "rbfds.h"
#include "greedy.h"
//#include "selectRbf.h"
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <Eigen/Dense>

rbf::rbf(Mesh& meshObject, probParams& probParamsObject)
:m(meshObject), params(probParamsObject)
//:mPtr(meshPtr), dVec(displacementVector), rotPnt(rotationPnt), steps(steps), smode(slidingMode), rotVec(rVec), curved(curved), perDir(periodicDirection), dataRed(dataReduction)
{
	// todo specifying the periodic direction(s)
	// todo smode is already in the mesh class -> replace all instances with m.smode

//	Greedy go;

}

void rbf::RBFMain(){

//	selectRbf<rbf_std> rbfObject(m, dVec, rotPnt, rotVec, steps, smode, curved, perDir);
//	std::exit(0);
	// called through the initliaser of the the rbf std class

//	rbfGenFunc genFuns(m,params);


//	std::cout << ptr << std::endl;


//
	getNodeType n(m);


	// make a structure with the various node types and numbers
	// use that as function


//todo can be replaced by a switch
	if(params.sMode=="none"){
//		 initialising the rbf standard class and performing the interpolation
		rbf_std rbf(m, params);
		rbf.perform_rbf(n);

	}
//	else if(smode=="ds"){
//		if(m.nDims ==3){
////			RBF_DS_3D();
//		}else if(m.nDims ==2){
//			rbf_ds rbf(m, dVec, rotPnt, rotVec, steps, smode, curved, perDir);
//			rbf.perform_rbf(n.iNodes,n.mNodes,n.mNodesStd,n.sNodes,n.N_i,n.N_m,n.N_mStd,n.N_s);
////			RBF_DS();
//		}
//	}
	else if(params.sMode=="ps"){
////		RBF_PS();
//		rbf_ps *ptr;
		rbf_ps rbf(m,params);

		if(params.dataRed){
			std::cout << "perform Greedy" << std::endl;
//			greedy go(n,rbf);
			rbf.perform_rbf(n);
		}else{
			rbf.perform_rbf(n);
		}
////		rbf.perform_rbf();
//		ptr = &rbf;
//
//		greedy go(ptr);
//
//
////		std::cout << ptr << std::endl;
////		greedy go(rbf);
//
//
	}


}




//void rbf::RBF_DS_3D(){
//	std::cout << "Performing RBF DS in three dimensions" << std::endl;
//
//	Eigen::MatrixXd Phi_mm, Phi_ms, Phi_sm, Phi_ss, Phi_im, Phi_is, Phi;
////	Eigen::ArrayXXd n(m.N_se, m.nDims), t(m.N_se, m.nDims);		// two column array containing normal vector components
//	Eigen::VectorXd defVec, alpha(m.nDims*(N_m+m.N_se+m.N_ss));
////	Eigen::ArrayXXd t_se(m.N_se,m.nDims),n1_se(m.N_se,m.nDims), n2_se(m.N_se,m.nDims), n_ss(m.N_ss, m.nDims),t1_ss(m.N_ss, m.nDims),t2_ss(m.N_ss, m.nDims);
//	Eigen::MatrixXd Phi_me, Phi_ee, Phi_es, Phi_em, Phi_se, Phi_ie;
//	for (int i = 0; i < steps; i++){
//		std::cout << "Deformation step: " << i+1 << std::endl;
//		getPhi(Phi_mm, mNodes,mNodes);
//		getPhi(Phi_me, mNodes, m.slidingEdgeNodes);
//		getPhi(Phi_ms, mNodes, m.slidingSurfNodes);
//
//		getPhi(Phi_em, m.slidingEdgeNodes, mNodes);
//		getPhi(Phi_ee, m.slidingEdgeNodes, m.slidingEdgeNodes);
//		getPhi(Phi_es, m.slidingEdgeNodes, m.slidingSurfNodes);
//
//		getPhi(Phi_sm, m.slidingSurfNodes, mNodes);
//		getPhi(Phi_se, m.slidingSurfNodes, m.slidingEdgeNodes);
//		getPhi(Phi_ss, m.slidingSurfNodes, m.slidingSurfNodes);
//		getPhi(Phi_im, m.intNodes, mNodes);
//		getPhi(Phi_ie, m.intNodes, m.slidingEdgeNodes);
//		getPhi(Phi_is, m.intNodes, m.slidingSurfNodes);
//
//
//
//
//		std::cout << "Reached the point at which normal and tangential vectors should be found" << std::endl;
//		m.getVecs();
////		m.getVecs3D(t_se, n1_se, n2_se, n_ss,t1_ss,t2_ss); // function that obtains all the normal and tan vectors
//
//
//		defVec = Eigen::VectorXd::Zero((N_m+m.N_se+m.N_ss)*m.nDims);
//		getDefVec(defVec,N_m,(N_m+m.N_se+m.N_ss));
////		std::cout << defVec << std::endl;
//
//		getPhiDS_3D(Phi,Phi_mm, Phi_me, Phi_ms, Phi_em, Phi_ee, Phi_es, Phi_sm, Phi_se, Phi_ss);
////		std::cout << Phi << std::endl;
////		std::cout << "DONE" << std::endl;
////		std::exit(0);
//		performRBF_DS_3D(Phi, Phi_im, Phi_ie, Phi_is, Phi_em, Phi_ee, Phi_es, Phi_sm, Phi_se, Phi_ss,defVec, alpha);
//	}

//	m.writeMeshFile();
//}

//void rbf::performRBF_DS_3D(Eigen::MatrixXd& Phi, Eigen::MatrixXd& Phi_im, Eigen::MatrixXd& Phi_ie, Eigen::MatrixXd& Phi_is, Eigen::MatrixXd& Phi_em, Eigen::MatrixXd& Phi_ee, Eigen::MatrixXd& Phi_es, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_se, Eigen::MatrixXd& Phi_ss, Eigen::VectorXd& defVec, Eigen::VectorXd& alpha){
//
//	alpha = Phi.partialPivLu().solve(defVec);
//	std::cout << "SOLUTION WAS FOUND" << std::endl;
//
//	for (int dim = 0; dim < m.nDims; dim++){
//
//		m.coords(m.intNodes, dim) += (Phi_im*alpha(Eigen::seqN(dim*(N_m+m.N_se+m.N_ss),N_m)) + Phi_ie*alpha(Eigen::seqN(dim*(N_m+m.N_se+m.N_ss)+N_m, m.N_se)) + Phi_is*alpha(Eigen::seqN(dim*(N_m+m.N_se+m.N_ss)+N_m+m.N_se, m.N_ss)) ).array();
//		m.coords(m.slidingEdgeNodes, dim) += (Phi_em*alpha(Eigen::seqN(dim*(N_m+m.N_se+m.N_ss),N_m)) + Phi_ee*alpha(Eigen::seqN(dim*(N_m+m.N_se+m.N_ss)+N_m, m.N_se)) + Phi_es*alpha(Eigen::seqN(dim*(N_m+m.N_se+m.N_ss)+N_m+m.N_se, m.N_ss)) ).array();
//		m.coords(m.slidingSurfNodes, dim) += (Phi_sm*alpha(Eigen::seqN(dim*(N_m+m.N_se+m.N_ss),N_m)) + Phi_se*alpha(Eigen::seqN(dim*(N_m+m.N_se+m.N_ss)+N_m, m.N_se)) + Phi_ss*alpha(Eigen::seqN(dim*(N_m+m.N_se+m.N_ss)+N_m+m.N_se, m.N_ss)) ).array();
//		m.coords(m.intBdryNodes, dim) += (defVec(Eigen::seqN(dim*N_m,m.N_ib))).array();
//		rotPnt(dim) += dVec(dim);
//	}
//}

//void rbf::getPhiDS_3D(Eigen::MatrixXd& Phi, Eigen::MatrixXd& Phi_mm, Eigen::MatrixXd&  Phi_me, Eigen::MatrixXd&  Phi_ms, Eigen::MatrixXd& Phi_em, Eigen::MatrixXd& Phi_ee, Eigen::MatrixXd& Phi_es, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_se, Eigen::MatrixXd& Phi_ss){
//
//	Phi = Eigen::MatrixXd::Zero(m.nDims*(N_m+m.N_se+m.N_ss),m.nDims*(N_m+m.N_se+m.N_ss));
//
//	for(int dim = 0; dim< m.nDims; dim++){
//		// blocks related to the known displacements
//		Phi.block(dim*N_m, dim*(N_m+m.N_se+m.N_ss), N_m, N_m) = Phi_mm;
//		Phi.block(dim*N_m, dim*(N_m+m.N_se+m.N_ss)+N_m, N_m, m.N_se) = Phi_me;
//		Phi.block(dim*N_m, dim*(N_m+m.N_se+m.N_ss)+N_m+m.N_se, N_m, m.N_ss) = Phi_ms;
//
//		//blocks realteð to the first zero normal displacement condition of the sliding edge nodes
//		Phi.block(3*N_m, dim*(N_m+m.N_se+m.N_ss),					m.N_se, N_m ) = Phi_em.array().colwise() * m.n1_se.col(dim);
//		Phi.block(3*N_m, dim*(N_m+m.N_se+m.N_ss) + N_m,				m.N_se, m.N_se ) = Phi_ee.array().colwise() * m.n1_se.col(dim);
//		Phi.block(3*N_m, dim*(N_m+m.N_se+m.N_ss) + N_m + m.N_se,	m.N_se,	m.N_ss ) = Phi_es.array().colwise() * m.n1_se.col(dim);
//		//blocks realteð to the second zero normal displacement condition of the sliding edge nodes
//		Phi.block(3*N_m+m.N_se, dim*(N_m+m.N_se+m.N_ss),					m.N_se, N_m ) = Phi_em.array().colwise() * m.n2_se.col(dim);
//		Phi.block(3*N_m+m.N_se, dim*(N_m+m.N_se+m.N_ss) + N_m,				m.N_se, m.N_se ) = Phi_ee.array().colwise() * m.n2_se.col(dim);
//		Phi.block(3*N_m+m.N_se, dim*(N_m+m.N_se+m.N_ss) + N_m + m.N_se,		m.N_se,	m.N_ss ) = Phi_es.array().colwise() * m.n2_se.col(dim);
//		// blocks related to the zero tangential contribution of the sliding edge nodes
//		Phi.block(3*N_m+2*m.N_se, dim*(N_m+m.N_se+m.N_ss) + N_m,	m.N_se, m.N_se) = Eigen::MatrixXd(m.t_se.col(dim).matrix().asDiagonal());
//		// blocks related to the zero normal displacement condition of the sliding surface nodes
//		Phi.block(3*N_m+3*m.N_se, dim*(N_m+m.N_se+m.N_ss),		m.N_ss,N_m) = Phi_sm.array().colwise() * m.n_ss.col(dim);
//		Phi.block(3*N_m+3*m.N_se, dim*(N_m+m.N_se+m.N_ss)+N_m,	m.N_ss,m.N_se) = Phi_se.array().colwise() * m.n_ss.col(dim);
//		Phi.block(3*N_m+3*m.N_se, dim*(N_m+m.N_se+m.N_ss)+N_m+m.N_se,	m.N_ss,m.N_ss) = Phi_ss.array().colwise() * m.n_ss.col(dim);
//
//		// blocks related to the zero tangential contribution of the sliding surface nodes.
//		Phi.block(3*N_m+3*m.N_se+m.N_ss, dim*(N_m+m.N_se+m.N_ss)+N_m+m.N_se, m.N_ss,m.N_ss) = Eigen::MatrixXd(m.t1_ss.col(dim).matrix().asDiagonal());
//		Phi.block(3*N_m+3*m.N_se+2*m.N_ss, dim*(N_m+m.N_se+m.N_ss)+N_m+m.N_se, m.N_ss,m.N_ss) = Eigen::MatrixXd(m.t2_ss.col(dim).matrix().asDiagonal());
//	}


//}











