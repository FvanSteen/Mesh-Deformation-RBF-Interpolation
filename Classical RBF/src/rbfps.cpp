#include "rbfps.h"

#include <iostream>
#include <chrono>
#include <Eigen/Dense>

//rbf_ps::rbf_ps(Mesh& meshObject, Eigen::VectorXd& dVec, Eigen::RowVectorXd& rotPnt, Eigen::VectorXd& rotVec, const int& steps, const std::string& smode, const bool& curved, const std::string& pDir)
rbf_ps::rbf_ps(Mesh& meshObject, struct probParams& probParamsObject)
:rbf_std(meshObject, probParamsObject)
//:rbf_std(meshObject,dVec, rotPnt, rotVec, steps, smode, curved, pDir)
{
	std::cout << "in the ps class" << std::endl;

	projection proObject;
	p = &proObject;


//	perform_rbf(iNodes,mNodes,mNodesPro,sNodes,N_i,N_m,N_mPro,N_s);

//	iNodes.resize(m.N_i+m.N_p);
//	iNodes << m.intNodes, m.periodicNodes;

//	mNodes.resize(m.N_ib + m.N_es + m.N_se);
//	mNodes << m.intBdryNodes, m.extStaticNodes, m.slidingEdgeNodes;
//
//
//
//	if(m.pmode == "moving"){
//		mNodesPro.resize(m.N_ib);
//		mNodesPro << m.intBdryNodes;
//
//		sNodes.resize(m.N_se + m.N_es);
//		sNodes << m.slidingEdgeNodes, m.extStaticNodes;
//
//	}else{
//		mNodesPro.resize(m.N_ib+m.N_es);
//		mNodesPro << m.intBdryNodes, m.extStaticNodes;
//
//		sNodes.resize(m.N_se); // sNodes is always the sliding nodes in 2D
//		sNodes << m.slidingEdgeNodes;
//	}
//	N_i = iNodes.size();
//	N_m = mNodes.size();
//	N_s = sNodes.size();
//	N_mPro = mNodesPro.size();


}

void rbf_ps::perform_rbf(getNodeType& n){
	std::cout<< "Performing RBF PS" << std::endl;
	auto start = std::chrono::high_resolution_clock::now();

	Eigen::MatrixXd Phi_mm, Phi_sm, Phi_mmStd, Phi_im; 	//In this case only the internal boundary and static ext bdry
	Eigen::VectorXd defVec, defVecStd;
	Eigen::ArrayXXd delta, finalDef;

	if(params.dataRed){
		n.GreedyInit();
		n.greedyNodes(m.intBdryNodes(0),params.sMode);

	}

	int maxErrorNode;
	greedy go;

	int iter;
	for(int i = 0; i < params.steps; i++){
		std::cout << "Deformation step: " << i+1 << " out of "<< params.steps << std::endl;
		double error = 1;
		iter = 0;
		while(error > params.tolerance){
			if(iter!=0){
				n.greedyNodes(maxErrorNode,params.sMode);
//				std::cout << n.N_mStd << std::endl;
			}
			delta.resize(m.N_se, m.nDims);
			finalDef.resize(m.N_se,m.nDims);

			getPhi(Phi_mm, m.mNodes, m.mNodes);

			getPhi(Phi_sm, m.seNodes, m.mNodes);
			getPhi(Phi_mmStd, m.mNodesStd, m.mNodesStd);

			getPhi(Phi_im, m.iNodes, m.mNodesStd);

			std::cout << "got Phi's" << std::endl;

			defVec = Eigen::VectorXd::Zero(m.N_m*m.nDims);
			std::cout << "initialised deformation vector " << std::endl;
			getDefVec(defVec, m.N_m,m.intBdryNodes);
			std::cout << "got deformation vector " << std::endl;

			if(params.curved || i==0){
				m.getMidPnts();

				m.getVecs();
			}


			std::cout << "normals are obtained" << std::endl;
			performRBF_PS(Phi_mm, Phi_sm, Phi_mmStd, Phi_im, defVec, delta, finalDef, defVecStd, n);

			std::cout << "performed rbf PS" << std::endl;

			if(params.dataRed){
				// exact deformation of internal boundary nodes not included in the interpolation
				Eigen::VectorXd exactDef;
				getExactDef(n, exactDef);

				//next statement should also take into account the periodic nodes
				if(m.N_i == n.N_i){
					std::cout << "error zet to zero" << std::endl;
					error = 0;
				}else{
					go.getError(n,m,d,exactDef,error,maxErrorNode,params.sMode);
				}
				std::cout << "error: \t"<< error <<" at node: \t" << maxErrorNode<< std::endl;
			}else{
				error = 0;
			}
			iter++;

		}
		std::cout << iter << std::endl;
		if(params.dataRed){
			updateNodes(defVecStd,n.mNodesStd,n.iNodes,n.N_mStd);
//			std::cout << n.mNodesStd << std::endl;
		}


	}
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
	std::cout <<  "Runtime duration: \t"<<  duration.count()/1e6 << " seconds"<< std::endl;
	m.writeMeshFile();
}

void rbf_ps::performRBF_PS(Eigen::MatrixXd& Phi_mm, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_mmStd, Eigen::MatrixXd& Phi_im, Eigen::VectorXd& defVec,Eigen::ArrayXXd& delta, Eigen::ArrayXXd& finalDef, Eigen::VectorXd& defVecStd,getNodeType& n){
	auto starti = std::chrono::high_resolution_clock::now();



	for(int dim = 0; dim < m.nDims; dim++){

		delta.col(dim) = (Phi_sm*(Phi_mm.fullPivLu().solve(defVec(Eigen::seqN(dim*m.N_m,m.N_m))))).array();

	}


	auto stopi = std::chrono::high_resolution_clock::now();
	auto durationi = std::chrono::duration_cast<std::chrono::microseconds>(stopi-starti);
	std::cout <<  "obtaining solution sliding nodes: \t"<<  durationi.count()/1e6 << " seconds"<< std::endl;

	std::cout << m.seNodes << std::endl;
	std::cout << delta << std::endl;


	// todo make pVec an array instead of vector
//	m.getVecs();
//	for(int i=0; i<m.N_se; i++){
//		finalDef.row(i) = (delta.row(i)).matrix().transpose() - (delta.row(i).matrix()).dot((m.n.row(i)).matrix().transpose())*(m.n.row(i)).matrix().transpose();
////				std::cout << "finalDef \t" << finalDef.row(i) << std::endl;
////				std::cout << "delta \t" << delta.row(i) << std::endl;
//	}
	// updating the midpoints on the external boundary of the mesh
//	std::cout << delta << std::endl;
	starti = std::chrono::high_resolution_clock::now();
	p->project(m, m.seNodes, delta, finalDef, pVec);

//	std::cout << finalDef << std::endl;
//	std::exit(0);
	stopi = std::chrono::high_resolution_clock::now();
	durationi = std::chrono::duration_cast<std::chrono::microseconds>(stopi-starti);
	std::cout <<  "projection: \t"<<  durationi.count()/1e6 << " seconds"<< std::endl;
//	std::cout << finalDef << std::endl;


	//todo make an if statement for fixed vertices
//	if(m.pmode == "moving"){
//		for(int i=0; i < (m.N_m - m.N_ib); i++){
//			finalDef.row(m.N_se+i) = pVec.transpose().array()*delta.row(m.N_se+i);
//		std::cout << "finalDef \t" << finalDef.row(m.N_se+i) << std::endl;
//		std::cout << "delta \t" << delta.row(m.N_se+i) << std::endl;
//			m.coords.row(m.extStaticNodes(i)) +=  pVec.transpose().array()*delta.row(m.N_se+i);
//		}
//	}



//	std::cout << finalDef << std::endl;
//	std::exit(0);

//	m.coords(m.extStaticNodes,1) += delta(Eigen::seq(Eigen::last+1-4,Eigen::last), 1);
//	m.coords(sNodes,Eigen::all) = m.coords(sNodes,Eigen::all) + finalDef;
//	m.writeMeshFile();
	std::cout << "obtaining defVecStd" << std::endl;
//	std::cout << n.N_mStd << std::endl;
//	std::cout << n.N_m << std::endl;
//	std::cout << n.N_s << std::endl;
	starti = std::chrono::high_resolution_clock::now();
	defVecStd = Eigen::VectorXd::Zero(m.N_mStd*m.nDims);
//	std::cout << N_m << '\t' << N_mPro <<  '\t' << m.N_se << std::endl;
	for(int dim = 0; dim< m.nDims; dim++){
		defVecStd(dim*m.N_mStd+m.ibIndices) = defVec(dim*m.N_m+m.ibIndices);
//		defVec(Eigen::seqN(dim*N_m+N_mPro,m.N_se)) = finalDef.col(dim);
		defVecStd(Eigen::seqN(dim*m.N_mStd+m.N_m,m.N_se)) = finalDef.col(dim);
	}

//	std::cout << defVec << std::endl;
//	std::cout << '\n' << finalDef << std::endl;
//	std::cout << '\n' << defVecStd << std::endl;


//	std::cout << "obtained defVecStd" << std::endl;
//	std::cout << "finalDef: \n" << finalDef << "\n" << std::endl;
//	std::cout << "defVec: \n" << defVec << std::endl;


	//for greedy:
//	Eigen::ArrayXXd d;

//	rbf_std stand;

	performRBF(Phi_mmStd,Phi_im,defVecStd,m.mNodesStd,m.iNodes,m.N_mStd);
	stopi = std::chrono::high_resolution_clock::now();
	durationi = std::chrono::duration_cast<std::chrono::microseconds>(stopi-starti);
	std::cout <<  "performing classical rbf: \t"<<  durationi.count()/1e6 << " seconds"<< std::endl;
}



