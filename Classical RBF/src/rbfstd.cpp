#include "rbfstd.h"
#include "greedy.h"
#include <Eigen/Dense>
#include <chrono>
#include <iostream>
//rbf_std::rbf_std(Mesh& meshObject, Eigen::VectorXd& dVec, Eigen::RowVectorXd& rotPnt, Eigen::VectorXd& rotVec, const int& steps, const std::string& smode, const bool& curved, const std::string& pDir)
rbf_std::rbf_std(Mesh& meshObject, struct probParams& probParamsObject)
:rbfGenFunc(meshObject, probParamsObject)
{
}


void rbf_std::perform_rbf(getNodeType& n){
	std::cout << "Performing standard rbf interpolation" << std::endl;

//	n.assignNodeTypes();

//	std::cout << "moving: \n" << *n.mPtr << std::endl;
//	std::cout << "int: \n" << *n.iPtr << std::endl;

	auto start = std::chrono::high_resolution_clock::now();
	Eigen::MatrixXd Phi_mm, Phi_im, Phi_imGreedy;
	Eigen::VectorXd defVec,defVecStd;


	// node containing max error, iter for nr of greedy iterations
	int iter;
	Eigen::ArrayXi maxErrorNodes;
	if(params.dataRed){
//		maxErrorNode = m.intBdryNodes(0);
		n.addControlNode(m.intBdryNodes(0));
	}
	// max error.
	double error;

	greedy go;


	for(int i=0; i<params.steps; i++){
		std::cout << "Deformation step: " << i+1 << " out of "<< params.steps << std::endl;
		error = 1;
		iter = 0;
		while(error>params.tol){

			if(iter!=0){
				for(int node = 0; node < maxErrorNodes.size(); node++){
					n.addControlNode(maxErrorNodes(node));
				}
			}

			std::cout << "Obtaining Phi_mm" << std::endl;
			getPhi(Phi_mm, *n.mPtr, *n.mPtr); // could simply pass the pointer here

			std::cout << "Obtaining Phi_im" << std::endl;
			getPhi(Phi_im, *n.iPtr, *n.mPtr);

			if(i==0 || params.dataRed){
				std::cout << "Obtaining deformation vector" << std::endl;
				defVec = Eigen::VectorXd::Zero(n.N_m*m.nDims);
				getDefVec(defVec,n.N_m,params.steps,*n.mPtr);
			}


			std::cout << "Performing RBF" << std::endl;
			performRBF(Phi_mm, Phi_im, defVec,*n.mPtr,*n.iPtr, n.N_m);

			if(params.dataRed){

				//next statement should also take into account the periodic nodes
				if(m.N_i == n.N_i){
					std::cout << "error zet to zero" << std::endl;
					error = 0;
				}else{
					go.getError(n,m,d,error,maxErrorNodes, params.smode, mIndex, displacement,pnVec);
				}
				std::cout << "error: \t"<< error <<" at node: \t" << maxErrorNodes(0)<< std::endl;

			}else{
				error = 0;
			}
			iter++;
		}
		auto stop = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
		std::cout <<  "Runtime duration: \t"<<  duration.count()/1e6 << " seconds"<< std::endl;

//		std::cout << iter-1 << std::endl;
//		m.coords(*n.iPtr, Eigen::all) +=d;
//		m.writeMeshFile();
//		std::exit(0);

		if(params.dataRed){
			updateNodes(Phi_imGreedy,n, defVec);
			go.correction(m,n, params.gamma);
			std::cout << "number of control nodes: " << n.N_m << std::endl;
		}

	}
//	std::cout << "number of control nodes: " << n.N_m << std::endl;
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);

	std::cout << "Runtime duration: \t"<<  duration.count()/1e6 << " seconds"<< std::endl;
	m.writeMeshFile();
}


void rbf_std::performRBF(Eigen::MatrixXd& Phi_mm, Eigen::MatrixXd& Phi_im, Eigen::VectorXd& defVec, Eigen::ArrayXi& movingNodes, Eigen::ArrayXi& internalNodes,int& N){
	alpha.resize(N*m.nDims);

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
//	std::cout << movingNodes << std::endl;
	for(int dim = 0; dim < m.nDims; dim++){

//		std::cout << "Solving for dimension: " << dim << std::endl;
		alpha(Eigen::seqN(dim*N,N)) = Phi_mm.fullPivHouseholderQr().solve(defVec(Eigen::seqN(dim*N,N)));
		if(params.dataRed){
			d.col(dim) = Phi_im*alpha(Eigen::seqN(dim*N,N));
//			std::cout << (Phi_mm.fullPivHouseholderQr().solve(defVec(Eigen::seqN(dim*N,N)))) << std::endl;
		}else{
			m.coords(internalNodes,dim) += (Phi_im*alpha(Eigen::seqN(dim*N,N))).array();
			m.coords(movingNodes,dim) += defVec(Eigen::seqN(dim*N,N)).array();
//			m.coords(internalNodes,dim) += (Phi_im*(Phi_mm.fullPivHouseholderQr().solve(defVec(Eigen::seqN(dim*N,N))))).array();
//			m.coords(movingNodes,dim) += defVec(Eigen::seqN(dim*N,N)).array();
//			params.rotPnt(dim) += params.dVec(dim);
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

void rbf_std::updateNodes(Eigen::MatrixXd& Phi_imGreedy, getNodeType& n, Eigen::VectorXd& defVec){
//	for(int dim = 0; dim < m.nDims; dim++){
//		m.coords(internalNodes,dim) += d.col(dim);//(Phi_im*(Phi_mm.fullPivHouseholderQr().solve(defVec(Eigen::seqN(dim*N,N))))).array();
//		m.coords(movingNodes,dim) += defVec(Eigen::seqN(dim*N,N)).array();
////		params.rotPnt(dim) += params.dVec(dim);
//	}

	int N_m;
	Eigen::ArrayXi* ptr;

	if(m.smode == "none"){
		ptr = n.mPtr;
		N_m = n.N_m;
	}else{
		ptr = n.mStdPtr;
		N_m = n.N_mStd;
	}

	getPhi(Phi_imGreedy,*n.iPtrGrdy,*ptr);

//	std::cout << Phi_imG reedy*alpha(Eigen::seqN(0,n.N_m)) << std::endl;
//	std::cout << Phi_imGreedy*alpha(Eigen::seqN(n.N_m,n.N_m)) << std::endl;

//	std::cout << d << std::endl;
//	std::cout << defVec << std::endl;
	m.coords(*n.iPtr,Eigen::all) += d;

//	m.coords(*n.mPtr,0) += defVec(Eigen::seqN(0,n.N_m)).array();
//	m.coords(*n.mPtr,1) += defVec(Eigen::seqN(n.N_m,n.N_m)).array();
//	std::cout << *n.mStdPtr << std::endl;

	m.coords(*n.iPtrGrdy,0) +=  (Phi_imGreedy*alpha(Eigen::seqN(0,N_m))).array();

	m.coords(*n.iPtrGrdy,1) +=  (Phi_imGreedy*alpha(Eigen::seqN(N_m,N_m))).array();


}

void rbf_std::getExactDef(getNodeType& n, Eigen::VectorXd& exactDeformation){

	int N = m.N_ib-n.N_ib;
	exactDeformation = Eigen::VectorXd::Zero(m.nDims*N);
//	std::cout << n.iNodes(Eigen::seq(rbf.m.N_i+rbf.m.N_p, rbf.m.N_i+rbf.m.N_p+rbf.m.N_ib-n.N_ib-1)) << std::endl;

	Eigen::ArrayXi ibNodes = n.iNodes(Eigen::seq(m.N_i+m.N_p, m.N_i+m.N_p+m.N_ib-n.N_ib-1));
//	getDefVec(exactDeformation, N, params.steps);
//	std::cout << exactDeformation << std::endl;
}
