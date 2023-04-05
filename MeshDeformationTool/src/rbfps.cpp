#include "rbfps.h"

#include <iostream>
#include <chrono>
#include <Eigen/Dense>
#include "greedy.h"

//rbf_ps::rbf_ps(Mesh& meshObject, Eigen::VectorXd& dVec, Eigen::RowVectorXd& rotPnt, Eigen::VectorXd& rotVec, const int& steps, const std::string& smode, const bool& curved, const std::string& pDir)
rbf_ps::rbf_ps( struct probParams& probParamsObject, Mesh& meshObject, getNodeType& n)
:rbfGenFunc(meshObject, probParamsObject)
{
	if(params.dataRed){
		greedy g(m, params, exactDisp, movingIndices, alpha, d);
		perform_rbf(n,g);
	}else{
		perform_rbf(n);
	}

}

void rbf_ps::perform_rbf(getNodeType& n){
	std::cout<< "Performing RBF PS without data reduction" << std::endl;

	std::clock_t s = std::clock();

	Eigen::VectorXd defVec, defVec_b;
	Eigen::ArrayXXd delta, finalDef;

	for(int i = 0; i < params.steps; i++){
		std::cout << "Deformation step: " << i+1 << " out of "<< params.steps << std::endl;

		if(params.curved || i==0){
			m.getMidPnts(params);
			m.getVecs(params);
		}

		getPhis(n, 0);

		if(i == 0){
			getDefVec(defVec, n.N_m, n.mPtr);
		}
		performRBF_PS(PhiPtr, defVec, delta, finalDef, defVec_b, n);

	}

	std::clock_t e = std::clock();
	long double time_elapsed_ms =  1000.0*(e-s) / CLOCKS_PER_SEC;
	std::cout << "CPU time: " << time_elapsed_ms/1000 << " ms\n";

}

void rbf_ps::perform_rbf(getNodeType& n, greedy& g){
	std::cout<< "Performing RBF PS with data reduction" << std::endl;

	std::clock_t s = std::clock();

	Eigen::VectorXd defVec, defVec_b;
	Eigen::ArrayXXd delta, finalDef;


	int iter, lvl;
	bool iterating = true;


	for(int i = 0; i < params.steps; i++){
		std::cout << "Deformation step: " << i+1 << " out of "<< params.steps << std::endl;

		iter = 0;
		lvl = 0;

		if(params.curved || i==0){
			m.getMidPnts(params);
			m.getVecs(params);
		}

		while(iterating){

			n.addControlNodes(g.maxErrorNodes, params.smode, m);


//			std::cout << "\n" <<  n.addedNodes.type << "\n" << std::endl;
			getPhis(n, iter);
//			if(Phis.Phi_cc.cols() >= 8){
//				std::cout << Phis.Phi_cc << std::endl;
//				std::exit(0);
//			}
			if(lvl > 0){
				getDefVec(defVec_b, n, g.errorPrevLvl, n.N_c);
			}else{
				getDefVec(defVec, n.N_m, n.mPtr);
			}

			if(lvl!=0){
				performRBF(Phis.Phi_cc, Phis.Phi_ic, defVec_b, n.cPtr, n.iPtr, n.N_c);
			}else{
				performRBF_PS(PhiPtr, defVec, delta, finalDef, defVec_b, n);
			}


			g.getError(n,d, lvl);
			std::cout << "error: \t"<< g.maxError <<" at node: \t" << g.maxErrorNodes(0)<< std::endl;

			if(g.maxError < params.tol){
				iterating = false;
				if(params.multiLvl == false){
					g.maxErrorNodes.resize(0);
				}
			}

			if(params.multiLvl && (g.maxError/g.maxErrorPrevLvl < params.tolCrit || iterating == false)){

				g.setLevelParams(n,lvl, d, alpha, defVec_b, n.cPtr, n.N_c);

				std::cout << "LEVEL: " << lvl << " HAS BEEN DONE" << std::endl;
				lvl++;
				iter = -1;

				n.assignNodeTypesGrdy(m);

				if(iterating == false){
					g.getAlphaVector();
					g.setInitMaxErrorNodes();
				}

			}
			iter++;

		}

		updateNodes(n,defVec_b, g.d_step, g.alpha_step, g.ctrlPtr);
		g.correction(m,n, params.gamma, params.multiLvl);
		iterating = true;

		std::cout << "number of control nodes: " << n.N_c << std::endl;

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
	delta.resize(n.N_se, m.nDims);
	finalDef.resize(n.N_se,m.nDims);

	for(int dim = 0; dim < m.nDims; dim++){
		delta.col(dim) = (PhiPtr->Phi_em*(PhiPtr->Phi_mm.fullPivLu().solve(defVec(Eigen::seqN(dim*n.N_m,n.N_m))))).array();
	}


//	for(int dim = 0; dim < m.nDims; dim++){
//		m.coords(*n.cPtr, dim) += (defVec(Eigen::seqN(dim*n.N_c, n.N_c))).array();
//	}


	p.project(m, n, delta, finalDef, m.periodicVec);



//	m.coords(*n.sPtr,Eigen::all) += finalDef;
//	m.writeMeshFile(params.mesh_ifName, params.mesh_ofName);
//	std::exit(0);

	getDefVec(defVec_b, defVec, n, finalDef);
	std::cout << "obtained defVec for the std rbf" << std::endl;
//	std::cout << defVec_b << "\n" << std::endl;
//	std::cout << PhiPtr->Phi_cc << "\n\n" << PhiPtr->Phi_ic << std::endl;
	performRBF(PhiPtr->Phi_cc,PhiPtr->Phi_ic,defVec_b,n.cPtr,n.iPtr,n.N_c);
	std::cout << "performed rbf" << std::endl;
}



