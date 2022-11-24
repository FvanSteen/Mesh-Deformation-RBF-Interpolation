///*
// * rbfds.cpp
// *
// *  Created on: 17 nov. 2022
// *      Author: floyd
// */
//
//#include "rbfds.h"
//#include <iostream>
//#include <Eigen/Dense>
//#include <chrono>
//rbf_ds::rbf_ds(Mesh& meshObject, Eigen::VectorXd& dVec, Eigen::RowVectorXd& rotPnt, Eigen::VectorXd& rotVec, const int& steps, const std::string& smode, const bool& curved, const std::string& pDir)
//:rbf_std(meshObject,dVec, rotPnt, rotVec, steps, smode, curved, pDir)
//{
//
//
//	std::cout << "Initialised the ds class" << std::endl;
//
//
////	iNodes.resize(m.N_i+m.N_p);
////	iNodes << m.intNodes, m.periodicNodes;
////
////	if(m.pmode == "moving"){
////		mNodes.resize(m.N_ib);
////		mNodes << m.intBdryNodes;
////		sNodes.resize(m.N_se+m.N_es);
////		sNodes << m.slidingEdgeNodes, m.extStaticNodes;
////	}else{
////		mNodes.resize(m.N_ib+m.N_es);
////		mNodes << m.intBdryNodes, m.extStaticNodes;
////		sNodes.resize(m.N_se);
////		sNodes << m.slidingEdgeNodes;
////	}
////	if(curved){
////		// todo rename to make clearer
////		mNodesStd.resize(m.N_ib+m.N_es+m.N_se);
////		mNodesStd << m.intBdryNodes,m.extStaticNodes, m.slidingEdgeNodes;
////		N_mStd = mNodesStd.size();
////	}
////
////	N_i = iNodes.size();
////	N_m = mNodes.size();
////	N_s = sNodes.size();
////	N_mStd = mNodesStd.size();
////
////	// todo
////
////	perform_rbf_ds();
//}
//
//void rbf_ds::perform_rbf(Eigen::ArrayXi& iNodes, Eigen::ArrayXi& mNodes, Eigen::ArrayXi& mNodesStd, Eigen::ArrayXi& sNodes, int& N_i, int& N_m, int& N_mStd, int& N_s){
//
//
//	projection* p;
//	if(curved){
//		projection proObject;
//		p = &proObject;
//	}
//
//	auto start = std::chrono::high_resolution_clock::now();
//	std::cout << "Performing RBF DS " << std::endl;
//
//	Eigen::MatrixXd Phi_mm, Phi_ms, Phi_sm, Phi_ss, Phi_im, Phi_is, Phi;
//	Eigen::VectorXd defVec, alpha(m.nDims*(N_m+N_s));
//
//	for (int i = 0; i < steps; i++){
//		auto start2 = std::chrono::high_resolution_clock::now();
//		std::cout << "Deformation step: " << i+1 << std::endl;
//
//		auto starti = std::chrono::high_resolution_clock::now();
//		getPhi(Phi_mm, mNodes,mNodes);
//		getPhi(Phi_ms, mNodes, sNodes);
//		getPhi(Phi_sm, sNodes, mNodes);
//		getPhi(Phi_ss, sNodes, sNodes);
//
//		auto stopi = std::chrono::high_resolution_clock::now();
//		auto durationi = std::chrono::duration_cast<std::chrono::microseconds>(stopi-starti);
//		std::cout << "Obtaining phi with s and m: \t"<<  durationi.count()/1e6 << " seconds"<< std::endl;
//
//		starti = std::chrono::high_resolution_clock::now();
//		getPhi(Phi_im, iNodes, mNodes);
//		getPhi(Phi_is, iNodes, sNodes);
//		stopi = std::chrono::high_resolution_clock::now();
//		durationi = std::chrono::duration_cast<std::chrono::microseconds>(stopi-starti);
//		std::cout << "obtaining phi with i: \t"<<  durationi.count()/1e6 << " seconds"<< std::endl;
//
//
//
//		if(curved || i==0){
////		if(i==0){
//			// getVecs obtains average vector at the nodes
//			starti = std::chrono::high_resolution_clock::now();
//			m.getVecs();
//			stopi = std::chrono::high_resolution_clock::now();
//			durationi = std::chrono::duration_cast<std::chrono::microseconds>(stopi-starti);
//			std::cout << "obtaining normals nodes: \t"<<  durationi.count()/1e6 << " seconds"<< std::endl;
//
//
//			// getMidPnts obtains vectors at midpoint of boundary segments
//
//			starti = std::chrono::high_resolution_clock::now();
//			m.getMidPnts();
//			stopi = std::chrono::high_resolution_clock::now();
//			durationi = std::chrono::duration_cast<std::chrono::microseconds>(stopi-starti);
//			std::cout << "obtaining normals segments: \t"<<  durationi.count()/1e6 << " seconds"<< std::endl;
//
//
//		}
//
////		defVec = Eigen::VectorXd::Zero((N_m+N_s)*m.nDims);
//		starti = std::chrono::high_resolution_clock::now();
//		defVec = Eigen::VectorXd::Zero((N_m+N_s)*m.nDims);
//		getDefVec(defVec,N_m);
//		stopi = std::chrono::high_resolution_clock::now();
//		durationi = std::chrono::duration_cast<std::chrono::microseconds>(stopi-starti);
//		std::cout << "obtaining deformation vector: \t"<<  durationi.count()/1e6 << " seconds"<< std::endl;
//
//
//		starti = std::chrono::high_resolution_clock::now();
//		getPhiDS(Phi,Phi_mm,Phi_ms, Phi_sm, Phi_ss, m.n, m.t,N_m,N_s);
//		stopi = std::chrono::high_resolution_clock::now();
//		durationi = std::chrono::duration_cast<std::chrono::microseconds>(stopi-starti);
//		std::cout << "assembly of Phi: \t"<<  durationi.count()/1e6 << " seconds"<< std::endl;
//
//		starti = std::chrono::high_resolution_clock::now();
//		// todo check which items can be omitted
//		performRBF_DS(Phi, Phi_im, Phi_is, Phi_sm, Phi_ss,defVec, alpha, p,iNodes,mNodes,mNodesStd, sNodes, N_i, N_m, N_mStd, N_s);
//		stopi = std::chrono::high_resolution_clock::now();
//		durationi = std::chrono::duration_cast<std::chrono::microseconds>(stopi-starti);
//		std::cout << "obtaining solution and updating nodes: \t"<<  durationi.count()/1e6 << " seconds"<< std::endl;
//
//
//		auto stop2 = std::chrono::high_resolution_clock::now();
//		auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(stop2-start2);
//		std::cout <<  "Runtime whole step duration: \t"<<  duration2.count()/1e6 << " seconds"<< std::endl;
//	}
//	auto stop = std::chrono::high_resolution_clock::now();
//	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
//	std::cout <<  "Runtime duration: \t"<<  duration.count()/1e6 << " seconds"<< std::endl;
//	m.writeMeshFile();
//}
//
//void rbf_ds::performRBF_DS(Eigen::MatrixXd& Phi, Eigen::MatrixXd& Phi_im, Eigen::MatrixXd& Phi_is, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_ss, Eigen::VectorXd& defVec, Eigen::VectorXd& alpha,projection* proPnt, Eigen::ArrayXi& iNodes, Eigen::ArrayXi& mNodes, Eigen::ArrayXi& mNodesStd, Eigen::ArrayXi& sNodes, int& N_i, int& N_m, int& N_mStd, int& N_s){
//	auto starti = std::chrono::high_resolution_clock::now();
//	alpha = Phi.fullPivLu().solve(defVec);
//	auto stopi = std::chrono::high_resolution_clock::now();
//	auto durationi = std::chrono::duration_cast<std::chrono::microseconds>(stopi-starti);
//	std::cout << "solving for alpha: \t"<<  durationi.count()/1e6 << " seconds"<< std::endl;
//	if(curved){
//		Eigen::ArrayXXd delta(N_s, m.nDims), finalDef(N_s,m.nDims);
//
//		// find displacement
//		for (int dim = 0; dim < m.nDims; dim++){
//			delta.col(dim) = (Phi_sm*alpha(Eigen::seqN(dim*(N_m+N_s),N_m)) + Phi_ss*alpha(Eigen::seqN(dim*(N_m+N_s)+N_m, N_s))).array();
//
//		}
////		m.coords(sNodes,Eigen::seq(0,1)) += delta;
//
//
//
//		// calling project function to find the final deformation after the projection
//		auto starti = std::chrono::high_resolution_clock::now();
//		proPnt->project(m,sNodes,delta,finalDef,pVec);
//
//		auto stopi = std::chrono::high_resolution_clock::now();
//		auto durationi = std::chrono::duration_cast<std::chrono::microseconds>(stopi-starti);
//		std::cout << "projection: \t"<<  durationi.count()/1e6 << " seconds"<< std::endl;
////		m.coords(sNodes,Eigen::seq(0,1)) += finalDef;
////		m.writeMeshFile();
////
//
////		std::exit(0);
////		defVec = Eigen::VectorXd::Zero(N_m2*m.nDims);
//		defVec = Eigen::VectorXd::Zero(N_mStd*m.nDims);
//		getDefVec(defVec,N_mStd);
////		std::cout << defVec << std::endl;
//
////		std::cout << defVec << std::endl;
//		for(int dim = 0; dim< m.nDims; dim++){
////			defVec(Eigen::seqN(dim*(N_m2)+m.N_ib+m.N_es,N_s)) = finalDef.col(dim);
//			defVec(Eigen::seqN(dim*(N_mStd)+N_m,N_s)) = finalDef.col(dim);
//			//todo difference in last two lines, for fixed/moving
//		}
//
//
//
//
////		std::cout << finalDef << std::endl;
////		std::cout << defVec << std::endl;
////		std::exit(0);
//		starti = std::chrono::high_resolution_clock::now();
//		Eigen::MatrixXd Phi_mm2, Phi_im2;
//
//		getPhi(Phi_mm2, mNodesStd,mNodesStd);
//		getPhi(Phi_im2,iNodes,mNodesStd);
//		Eigen::ArrayXXd d;
//		performRBF(Phi_mm2,Phi_im2,defVec,mNodesStd,iNodes,N_mStd,d);
//		stopi = std::chrono::high_resolution_clock::now();
//		durationi = std::chrono::duration_cast<std::chrono::microseconds>(stopi-starti);
//		std::cout << "Doing second classical rbf: \t"<<  durationi.count()/1e6 << " seconds"<< std::endl;
//
//	}
//	else{
//		auto start2 = std::chrono::high_resolution_clock::now();
//		for (int dim = 0; dim < m.nDims; dim++){
//			m.coords(iNodes, dim) += (Phi_im*alpha(Eigen::seqN(dim*(N_m+N_s),N_m)) + Phi_is*alpha(Eigen::seqN(dim*(N_m+N_s)+N_m, N_s))).array();
//			m.coords(sNodes, dim) += (Phi_sm*alpha(Eigen::seqN(dim*(N_m+N_s),N_m)) + Phi_ss*alpha(Eigen::seqN(dim*(N_m+N_s)+N_m, N_s))).array();
//			m.coords(m.intBdryNodes, dim) += (defVec(Eigen::seqN(dim*N_m,m.N_ib))).array();
//			rotPnt(dim) += dVec(dim);
//		}
//		auto stop2 = std::chrono::high_resolution_clock::now();
//		auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(stop2-start2);
//		std::cout << "Updating Nodes: \t"<<  duration2.count()/1e6 << " seconds"<< std::endl;
//	}
//}
//
//
//
//void rbf_ds::getPhiDS(Eigen::MatrixXd& Phi,Eigen::MatrixXd& Phi_mm,Eigen::MatrixXd& Phi_ms, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_ss, Eigen::ArrayXXd& n, Eigen::ArrayXXd& t,int& N_m, int& N_s){
//
//	Phi = Eigen::MatrixXd::Zero(m.nDims*(N_m+N_s),m.nDims*(N_m+N_s));
//
//
//	if(m.pmode == "moving"){
//		n.conservativeResize(N_s,m.nDims);
//		t.conservativeResize(N_s,m.nDims);
//		for(int i = 0; i<m.N_es; i++){
//			n.row(m.N_se+i) = pnVec;
//			t.row(m.N_se+i) = pVec;
//		}
//	}
//
//
//	for(int dim = 0; dim< m.nDims; dim++){
//		// blocks related to the known displacements
//		Phi.block(dim*N_m, dim*(N_m+N_s), N_m, N_m) = Phi_mm;
//		Phi.block(dim*N_m, dim*(N_m+N_s)+N_m, N_m, N_s) = Phi_ms;
//
//		// blocks related to the zero normal displacement condition
//		Phi.block(2*N_m, dim*(N_m+N_s), N_s, N_m) = Phi_sm.array().colwise() * n.col(dim);
//		Phi.block(2*N_m, dim*(N_m+N_s)+N_m, N_s, N_s) = Phi_ss.array().colwise() * n.col(dim);
//
//		//blocks related to the zero tangential contribution condition
//		Phi.block(2*N_m + N_s, dim*(N_m+N_s)+N_m, N_s, N_s) = Eigen::MatrixXd(t.col(dim).matrix().asDiagonal());
//	}
//
//
//}
//
