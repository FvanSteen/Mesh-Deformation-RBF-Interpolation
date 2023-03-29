#include "rbfps.h"

#include <iostream>
#include <chrono>
#include <Eigen/Dense>


//rbf_ps::rbf_ps(Mesh& meshObject, Eigen::VectorXd& dVec, Eigen::RowVectorXd& rotPnt, Eigen::VectorXd& rotVec, const int& steps, const std::string& smode, const bool& curved, const std::string& pDir)
rbf_ps::rbf_ps( struct probParams& probParamsObject, Mesh& meshObject, getNodeType& n)
:rbfGenFunc(meshObject, probParamsObject)
{
	perform_rbf(n);
}

void rbf_ps::perform_rbf(getNodeType& n){
	std::cout<< "Performing RBF PS" << std::endl;

	std::clock_t s = std::clock();

	Eigen::VectorXd defVec, defVec_b;
	Eigen::ArrayXXd delta, finalDef;


	Eigen::ArrayXi maxErrorNodes;
	greedy go(m, params, exactDisp, movingIndices, alpha, d, periodicVec);

	int iter, lvl;
	double maxError;
	bool iterating;


	for(int i = 0; i < params.steps; i++){
		std::cout << "Deformation step: " << i+1 << " out of "<< params.steps << std::endl;
		maxError = 1;
		iter = 0;
		lvl = 0;

		if((params.dataRed && i==0) || params.multiLvl ){
			go.setInitMaxErrorNodes();
		}

		iterating = true;
		if(params.curved || i==0){
			m.getMidPnts(params);
			m.getVecs();
		}

		while(iterating){

			if(params.dataRed){
//				for(int node = 0; node < maxErrorNodes.size(); node++){
//					n.addControlNode(maxErrorNodes(node), params.smode, m);
//				}
			}

			getPhis(n);

			if(lvl > 0){
				getDefVec(defVec_b, n, go.errorPrevLvl, n.N_b);
			}else if(i==0 || params.dataRed){
				getDefVec(defVec, n.N_c, n.cPtr);
			}

			if(lvl!=0){
				performRBF(Phis.Phi_bb, Phis.Phi_ib, defVec_b, n.bPtr, n.iPtr, n.N_b);
			}else{
				performRBF_PS(PhiPtr, defVec, delta, finalDef, defVec_b, n);
			}


			if(params.dataRed){

//				std::cout << "getting error" << std::endl;
				go.getError(n,d, lvl);
				std::cout << "error: \t"<< maxError <<" at node: \t" << maxErrorNodes(0)<< std::endl;

				if(maxError < params.tol){
					iterating = false;
					maxErrorNodes.resize(0);
				}

			}else{
				iterating = false;
			}

			if(params.multiLvl && (maxError/go.maxErrorPrevLvl < params.tolCrit || iterating == false)){

				go.setLevelParams(n,lvl, d, alpha, defVec_b, n.bPtr, n.N_b);

				std::cout << "LEVEL: " << lvl << " HAS BEEN DONE" << std::endl;
				lvl++;


				n.assignNodeTypesGrdy(m);

				if(maxError < params.tol){
					iterating = false;
					go.getAlphaVector();
				}

			}

			iter++;

		}


		if(params.dataRed){
//			m.coords(*n.iPtr, Eigen::all) += (d-go.error);
//			m.writeMeshFile(params.mesh_ifName, params.mesh_ofName);
//			std::exit(0);
//			std::cout << "Control nodes: \n" << *n.bPtr << std::endl;
//			std::cout << "internal nodes: \n" << *n.iPtr << std::endl;
//			for(int i = 0; i < m.nDims;i++){
//				m.coords(*n.bPtr, i) += defVec_b(Eigen::seqN(i*n.N_b, n.N_b)).array();
//			}
//			for(int i = 0; i < n.N_s; i++){
//				std::cout << (*n.sPtr)(i) << ",";
//			}
//			std::cout << std::endl;
//			m.coords(*n.iPtr, Eigen::all) += (d - go.error);
//			m.writeMeshFile(params.mesh_ifName, params.mesh_ofName);
//			std::cout << "DONE " << std::endl;
//			std::exit(0);
			updateNodes(n,defVec_b, go.d_step, go.alpha_step, go.ctrlPtr);

			go.correction(m,n, params.gamma, params.multiLvl);


		}

		std::cout << "number of control nodes: " << n.N_b << std::endl;

	}
	std::clock_t e = std::clock();
	long double time_elapsed_ms =  1000.0*(e-s) / CLOCKS_PER_SEC;
	std::cout << "CPU time: " << time_elapsed_ms/1000 << " ms\n";
}

void rbf_ps::performRBF_PS(PhiStruct* PhiPtr, Eigen::VectorXd& defVec,Eigen::ArrayXXd& delta, Eigen::ArrayXXd& finalDef, Eigen::VectorXd& defVec_b,getNodeType& n){
//todo implementation of doing the edge first and then the surfaces
//	if(m.nDims == 3){
//		delta.resize(n.N_se,m.nDims);
//		finalDef.resize(n.N_se, m.nDims);
//
//		for(int dim = 0; dim < m.nDims; dim++){
//			delta.col(dim) = (Phi_sc(Eigen::seqN(0,),)*(Phi_cc.fullPivLu().solve(defVec(Eigen::seqN(dim*n.N_c,n.N_c))))).array();
//		}
//	}
	delta.resize(n.N_s, m.nDims);
	finalDef.resize(n.N_s,m.nDims);

	for(int dim = 0; dim < m.nDims; dim++){
		delta.col(dim) = (PhiPtr->Phi_sc*(PhiPtr->Phi_cc.fullPivLu().solve(defVec(Eigen::seqN(dim*n.N_c,n.N_c))))).array();
	}

//	for(int dim = 0; dim < m.nDims; dim++){
//		m.coords(*n.cPtr, dim) += (defVec(Eigen::seqN(dim*n.N_c, n.N_c))).array();
//	}



	p.project(m, n, delta, finalDef, periodicVec);




//	m.coords(*n.sPtr,Eigen::all) += finalDef;
//	m.writeMeshFile(params.mesh_ifName, params.mesh_ofName);
//	std::exit(0);

	getDefVec(defVec_b, defVec, n, finalDef);
	std::cout << "obtained defVec for the std rbf" << std::endl;
	performRBF(PhiPtr->Phi_bb,PhiPtr->Phi_ib,defVec_b,n.bPtr,n.iPtr,n.N_b);
	std::cout << "performed rbf" << std::endl;
}



