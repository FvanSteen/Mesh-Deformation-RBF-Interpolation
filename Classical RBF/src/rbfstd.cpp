#include "rbfstd.h"

#include <Eigen/Dense>
#include <chrono>
#include <iostream>
//rbf_std::rbf_std(Mesh& meshObject, Eigen::VectorXd& dVec, Eigen::RowVectorXd& rotPnt, Eigen::VectorXd& rotVec, const int& steps, const std::string& smode, const bool& curved, const std::string& pDir)
rbf_std::rbf_std(Mesh& meshObject, struct probParams& probParamsObject)
:rbfGenFunc(meshObject, probParamsObject)
{

//	std::cout << n.iNodes << std::endl;
//	std::exit(0);
//
//	int N_i = m.N_i+m.N_p;
//	int N_m = m.N_ib+ m.N_es + m.N_se;
//	Eigen::ArrayXi iNodes(N_i);
//	Eigen::ArrayXi mNodes(N_m);
//	iNodes << m.intNodes, m.periodicNodes;
//	mNodes << m.intBdryNodes, m.extStaticNodes, m.slidingEdgeNodes;

//	perform_rbf(iNodes,mNodes,N_m);

}


void rbf_std::perform_rbf(getNodeType& n){

	std::cout << "Performing standard rbf interpolation" << std::endl;

	auto start = std::chrono::high_resolution_clock::now();
	Eigen::MatrixXd Phi_mm, Phi_im;
	Eigen::VectorXd defVec;

	for(int i=0; i<params.steps; i++){
		auto start2 = std::chrono::high_resolution_clock::now();
		std::cout << "Deformation step: " << i+1 << std::endl;

		std::cout << "Obtaining Phi_mm" << std::endl;
		auto starti = std::chrono::high_resolution_clock::now();
		getPhi(Phi_mm, n.mNodes,n.mNodes);

		auto stopi = std::chrono::high_resolution_clock::now();
		auto durationi = std::chrono::duration_cast<std::chrono::microseconds>(stopi-starti);
		std::cout << "Runtime duration obtaining Phi_mm: \t"<<  durationi.count()/1e6 << " seconds"<< std::endl;

		std::cout << "Obtaining Phi_im" << std::endl;
		starti = std::chrono::high_resolution_clock::now();
		getPhi(Phi_im, n.iNodes,n.mNodes);
		stopi = std::chrono::high_resolution_clock::now();
		durationi = std::chrono::duration_cast<std::chrono::microseconds>(stopi-starti);
		std::cout << "Runtime duration obtaining Phi_im: \t"<<  durationi.count()/1e6 << " seconds"<< std::endl;

		std::cout << "Obtaining deformation vector" << std::endl;
		starti = std::chrono::high_resolution_clock::now();
		defVec = Eigen::VectorXd::Zero(n.N_m*m.nDims);
		getDefVec(defVec,n.N_m,n.ibNodes);


		stopi = std::chrono::high_resolution_clock::now();
		durationi = std::chrono::duration_cast<std::chrono::microseconds>(stopi-starti);
		std::cout << "Runtime duration obtaining defVec: \t"<<  durationi.count()/1e6 << " seconds"<< std::endl;

		std::cout << "Performing RBF" << std::endl;
		starti = std::chrono::high_resolution_clock::now();
//		Eigen::ArrayXXd d;
		performRBF(Phi_mm, Phi_im, defVec,n.mNodes,n.iNodes, n.N_m);
		stopi = std::chrono::high_resolution_clock::now();
		durationi = std::chrono::duration_cast<std::chrono::microseconds>(stopi-starti);
		std::cout << "Runtime duration obtaining solution and updating nodes: \t"<<  durationi.count()/1e6 << " seconds"<< std::endl;




		auto stop2 = std::chrono::high_resolution_clock::now();
		auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(stop2-start2);
		std::cout << "Runtime duration whole step: \t"<<  duration2.count()/1e6 << " seconds"<< std::endl;

	}
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);

	std::cout << "Runtime duration: \t"<<  duration.count()/1e6 << " seconds"<< std::endl;
	m.writeMeshFile();
}


void rbf_std::performRBF(Eigen::MatrixXd& Phi_mm, Eigen::MatrixXd& Phi_im, Eigen::VectorXd& defVec, Eigen::ArrayXi& movingNodes, Eigen::ArrayXi& internalNodes,int& N){
	if(params.dataRed){
		d.resize(internalNodes.size(),m.nDims);
	}
//	std::cout << m.coords(movingNodes,Eigen::all) << std::endl;
//	double side1,side2;
////	std::cout << (m.coords.row(movingNodes(1))-m.coords.row(movingNodes(0)))<< std::endl;
////	std::cout << (m.coords.row(movingNodes(2))-m.coords.row(movingNodes(0)))<< std::endl;
//	side1 = (m.coords.row(movingNodes(1))-m.coords.row(movingNodes(0))).matrix().norm();
//	side2 = (m.coords.row(movingNodes(2))-m.coords.row(movingNodes(0))).matrix().norm();
////	std::cout << side1 << std::endl;
////	std::cout << side2 << std::endl;
//	std::cout << side1*side2 << std::endl;
//	std::cout << std::endl;
	for(int dim = 0; dim < m.nDims; dim++){
		std::cout << "Solving for dimension: " << dim << std::endl;
		if(params.dataRed){
			d.col(dim) = (Phi_im*(Phi_mm.fullPivHouseholderQr().solve(defVec(Eigen::seqN(dim*N,N))))).array();
		}else{
			m.coords(internalNodes,dim) += (Phi_im*(Phi_mm.fullPivHouseholderQr().solve(defVec(Eigen::seqN(dim*N,N))))).array();
			m.coords(movingNodes,dim) += defVec(Eigen::seqN(dim*N,N)).array();
			params.rotPnt(dim) += params.dVec(dim);
		}
//		m.coords(internalNodes,dim) += (Phi_im*(Phi_mm.fullPivLu().solve(defVec(Eigen::seqN(dim*N,N))))).array();
//		m.coords(internalNodes,dim) += d.col(dim);
//		m.coords(movingNodes,dim) += defVec(Eigen::seqN(dim*N,N)).array();
//		m.coords(movingNodes,dim) += (Phi_mm*(Phi_mm.fullPivLu().solve(defVec(Eigen::seqN(dim*N,N))))).array();

	}


//	 std::cout << gf.mPtr->coords(gf.mPtr->intBdryNodes,Eigen::all) << std::endl;

//	std::cout << (m.coords.row(movingNodes(1))-m.coords.row(movingNodes(0))).matrix().norm() << std::endl;
//	side1 = (m.coords.row(movingNodes(1))-m.coords.row(movingNodes(0))).matrix().norm();
//	side2 = (m.coords.row(movingNodes(2))-m.coords.row(movingNodes(0))).matrix().norm();
//	std::cout << side1 << std::endl;
//	std::cout << side2 << std::endl;
//	std::cout << side1*side2 << std::endl;
//	std::exit(0);
//	std::cout << m.coords(movingNodes,Eigen::all) << std::endl;

//	std::cout << "\n rotation point \n" << std::endl;
//	std::cout << gf.rotPnt << std::endl;
//
//	std::cout << "\n delta vector \n" << std::endl;
//	std::cout << gf.dVec << std::endl;
//	std::cout << "\n rotMat: \n" << std::endl;
//	std::cout << gf.rotMat << std::endl;
//	std::exit(0);

//	m.writeMeshFile();
//	std::exit(0);
}

void rbf_std::updateNodes(Eigen::VectorXd& defVec, Eigen::ArrayXi& movingNodes, Eigen::ArrayXi& internalNodes,int& N){
	for(int dim = 0; dim < m.nDims; dim++){
		m.coords(internalNodes,dim) += d.col(dim);//(Phi_im*(Phi_mm.fullPivHouseholderQr().solve(defVec(Eigen::seqN(dim*N,N))))).array();
		m.coords(movingNodes,dim) += defVec(Eigen::seqN(dim*N,N)).array();
		params.rotPnt(dim) += params.dVec(dim);
	}
}
