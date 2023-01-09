#include "rbfps.h"

#include <iostream>
#include <chrono>
#include <Eigen/Dense>

//rbf_ps::rbf_ps(Mesh& meshObject, Eigen::VectorXd& dVec, Eigen::RowVectorXd& rotPnt, Eigen::VectorXd& rotVec, const int& steps, const std::string& smode, const bool& curved, const std::string& pDir)
rbf_ps::rbf_ps(Mesh& meshObject, struct probParams& probParamsObject)
:rbf_std(meshObject, probParamsObject)
{
	std::cout << "in the ps class" << std::endl;

	projection proObject;
	p = &proObject;
}

void rbf_ps::perform_rbf(getNodeType& n){
	std::cout<< "Performing RBF PS" << std::endl;
	auto start = std::chrono::high_resolution_clock::now();

	Eigen::MatrixXd Phi_mm, Phi_sm, Phi_mmStd, Phi_im,Phi_imGrdy; 	//In this case only the internal boundary and static ext bdry
	Eigen::VectorXd defVec, defVecStd;
	Eigen::ArrayXXd delta, finalDef;


//	int maxErrorNode;
	Eigen::ArrayXi maxErrorNodes;

	if(params.dataRed){
		n.addControlNode(m.intBdryNodes(0));
	}

	greedy go;
	int iter;
	double error;
	for(int i = 0; i < params.steps; i++){
		std::cout << "Deformation step: " << i+1 << " out of "<< params.steps << std::endl;
		error = 1;
		iter = 0;
		while(error > params.tol){


			if(iter!=0){
				for(int node = 0; node < maxErrorNodes.size(); node++){
					n.addControlNode(maxErrorNodes(node));
				}

//				std::cout << "Moving nodes: \n" << *n.mPtr << "\n internal nodes (same): \n" << *n.iPtr << "\n sliding edge Nodes: \n" << *n.sePtr << "\n mnodesStd: \n" << *n.mStdPtr << std::endl;
			}
			delta.resize(n.N_se, m.nDims);
			finalDef.resize(n.N_se,m.nDims);

			getPhi(Phi_mm, *n.mPtr, *n.mPtr);

			getPhi(Phi_sm, *n.sePtr, *n.mPtr);
			getPhi(Phi_mmStd, *n.mStdPtr, *n.mStdPtr);

			getPhi(Phi_im, *n.iPtr, *n.mStdPtr);



			if(i==0 || params.dataRed){
				defVec = Eigen::VectorXd::Zero(n.N_m*m.nDims);
				getDefVec(defVec,n.N_m,params.steps,*n.mPtr);
			}


			if(params.curved || i==0){
				m.getMidPnts();
				m.getVecs();
			}

			performRBF_PS(Phi_mm, Phi_sm, Phi_mmStd, Phi_im, defVec, delta, finalDef, defVecStd, n);


			if(params.dataRed){
				go.getError(n,m,d,error, maxErrorNodes, params.smode, mIndex, displacement,pnVec);
				std::cout << "error: \t"<< error <<" at node: \t" << maxErrorNodes(0)<< std::endl;
			}else{
				error = 0;
			}
//			if(iter == 5){
//				m.coords(*n.iPtr, Eigen::all) +=d;
//				m.writeMeshFile();
//				std::exit(0);
//			}


			iter++;
		}
//		auto stop = std::chrono::high_resolution_clock::now();
//		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
//		std::cout <<  "Runtime duration: \t"<<  duration.count()/1e6 << " seconds"<< std::endl;

//		std::cout << iter-1 << std::endl;
//		m.coords(*n.iPtr, Eigen::all) +=d;
//		m.writeMeshFile();
//		std::exit(0);
//		std::cout << iter << std::endl;
		if(params.dataRed){
			updateNodes(Phi_imGrdy, n, defVec);
//			m.writeMeshFile();
//			std::exit(0);
//			std::cout << "DOING AN UPDATE" << std::endl;
			go.correction(m,n, params.gamma);
		}
//		m.writeMeshFile();
//		std::exit(0);


	}
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
	std::cout <<  "Runtime duration: \t"<<  duration.count()/1e6 << " seconds"<< std::endl;
	m.writeMeshFile();
}

void rbf_ps::performRBF_PS(Eigen::MatrixXd& Phi_mm, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_mmStd, Eigen::MatrixXd& Phi_im, Eigen::VectorXd& defVec,Eigen::ArrayXXd& delta, Eigen::ArrayXXd& finalDef, Eigen::VectorXd& defVecStd,getNodeType& n){
	auto starti = std::chrono::high_resolution_clock::now();



	for(int dim = 0; dim < m.nDims; dim++){
		delta.col(dim) = (Phi_sm*(Phi_mm.fullPivLu().solve(defVec(Eigen::seqN(dim*n.N_m,n.N_m))))).array();
	}


	auto stopi = std::chrono::high_resolution_clock::now();
	auto durationi = std::chrono::duration_cast<std::chrono::microseconds>(stopi-starti);
//	std::cout <<  "obtaining solution sliding nodes: \t"<<  durationi.count()/1e6 << " seconds"<< std::endl;
//	m.coords(*n.sePtr,Eigen::all) += delta;
//	m.writeMeshFile();
//	std::exit(0);

//
//	std::cout << *n.sePtr << std::endl;
//	std::cout << delta << std::endl;


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

	/* This is the old projection function
	p->project(m, *n.sePtr, delta, finalDef, pVec);
	*/
	p->projectIter(m, *n.sePtr, delta, finalDef);
//	std::cout << finalDef << std::endl;
//	m.coords(*n.sePtr,Eigen::all) += finalDef;
//	m.writeMeshFile();
//	std::exit(0);
//	std::exit(0);
	stopi = std::chrono::high_resolution_clock::now();
	durationi = std::chrono::duration_cast<std::chrono::microseconds>(stopi-starti);
//	std::cout <<  "projection: \t"<<  durationi.count()/1e6 << " seconds"<< std::endl;
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
//	std::cout << "obtaining defVecStd" << std::endl;
//	std::cout << n.N_mStd << std::endl;
//	std::cout << n.N_m << std::endl;
//	std::cout << n.N_s << std::endl;
	starti = std::chrono::high_resolution_clock::now();
	defVecStd = Eigen::VectorXd::Zero(n.N_mStd*m.nDims);
//	std::cout << N_m << '\t' << N_mPro <<  '\t' << m.N_se << std::endl;

	for(int dim = 0; dim< m.nDims; dim++){
		defVecStd(Eigen::seqN(dim*n.N_mStd,n.N_m)) = defVec(Eigen::seqN(dim*n.N_m,n.N_m));

//		defVec(Eigen::seqN(dim*N_m+N_mPro,m.N_se)) = finalDef.col(dim);
		defVecStd(Eigen::seqN(dim*n.N_mStd+n.N_m,n.N_se)) = finalDef.col(dim);

	}

//	std::cout << defVec << std::endl;
//	std::cout << '\n' << finalDef << std::endl;
//	std::cout << '\n' << defVecStd << std::endl;
//	std::exit(0);


//	std::cout << "obtained defVecStd" << std::endl;
//	std::cout << "finalDef: \n" << finalDef << "\n" << std::endl;
//	std::cout << "defVec: \n" << defVec << std::endl;


	//for greedy:
//	Eigen::ArrayXXd d;

//	rbf_std stand;

	performRBF(Phi_mmStd,Phi_im,defVecStd,*n.mStdPtr,*n.iPtr,n.N_mStd);
	stopi = std::chrono::high_resolution_clock::now();
	durationi = std::chrono::duration_cast<std::chrono::microseconds>(stopi-starti);
//	std::cout <<  "performing classical rbf: \t"<<  durationi.count()/1e6 << " seconds"<< std::endl;
}



