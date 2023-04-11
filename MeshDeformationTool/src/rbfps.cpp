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

	for(int i = 0; i < params.steps; i++){
		std::cout << "Deformation step: " << i+1 << " out of "<< params.steps << std::endl;

		if(params.curved || i==0){
			m.getMidPnts(params);
			m.getVecs(params);
		}

		getPhis(n, 0);

		if(i == 0){
			getDefVec(defVec_m, n.N_m, n.mPtr);
		}

		// find the finaldef of the edge nodes
		pseudo_sliding_edge(PhiPtr, n);

		// make a second deformation vector that also contains the displacement of the edges
		if(m.nDims  == 3){
			getDefVec(defVec_me, defVec_m, n, proj_disp_edge, n.N_m+n.N_se, n.N_m);

			pseudo_sliding_surf(PhiPtr, n);

			getDefVec(defVec_all, defVec_me, n, proj_disp_all, n.N_c, n.N_m+n.N_se);

		}else{
			getDefVec(defVec_all, defVec_m, n, proj_disp_edge, n.N_c, n.N_m); // check if right?
		}

		performRBF(PhiPtr->Phi_cc,PhiPtr->Phi_ic,defVec_all,n.cPtr,n.iPtr,n.N_c);

	}
}

void rbf_ps::perform_rbf(getNodeType& n, greedy& g){
	std::cout<< "Performing RBF PS with data reduction" << std::endl;

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


//			std::cout << "\n\n errornodes:\n" <<  g.maxErrorNodes << "\n" << std::endl;
//			std::cout << "\n added node types:\n" << n.addedNodes.type << "\n" << std::endl;
			getPhis(n, iter);

			if(lvl != 0){
				getDefVec(defVec_all, n, g.errorPrevLvl, n.N_c);
				performRBF(Phis.Phi_cc, Phis.Phi_ic, defVec_all, n.cPtr, n.iPtr, n.N_c);
			}else{
				getDefVec(defVec_m, n.N_m, n.mPtr);
				pseudo_sliding_edge(PhiPtr, n);

				if(m.nDims == 2){
					getDefVec(defVec_all, defVec_m, n, proj_disp_edge, n.N_c, n.N_m);
					performRBF(PhiPtr->Phi_cc,PhiPtr->Phi_ic,defVec_all,n.cPtr,n.iPtr,n.N_c);
				}else if (m.nDims == 3){
					getDefVec(defVec_me, defVec_m, n, proj_disp_edge, n.N_m+n.N_se, n.N_m);
					performRBF(PhiPtr->Phi_meme,PhiPtr->Phi_ic_reduced,defVec_me,n.cPtr,n.iPtr_reduced,n.N_m+n.N_se);
				}
			}


			g.getError(n,d, lvl);
			std::cout << "error: \t"<< g.maxError <<" at node: \t" << g.maxErrorNodes(0)<< std::endl;

			if(g.maxError < params.tol){
				iterating = false;
				if(params.multiLvl == false){
					g.maxErrorNodes.resize(0);
				}
			}

			if(params.multiLvl  && m.nDims == 2  && (g.maxError/g.maxErrorPrevLvl < params.tolCrit || iterating == false)){

				g.setLevelParams(n,lvl, d, alpha, defVec_all, n.cPtr, n.N_c);
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


		int iter_surf = 0;

		if(m.nDims == 3){
			iterating = true;

			pseudo_sliding_surf(PhiPtr,n);
			getDefVec(defVec_all, defVec_me, n, proj_disp_all, n.N_c, n.N_m+n.N_se);
			performRBF(PhiPtr->Phi_cc,PhiPtr->Phi_ic,defVec_all,n.cPtr,n.iPtr,n.N_c);

			g.getError(n,d,lvl);
		}

		while(iterating){


			n.addControlNodes(g.maxErrorNodes, params.smode, m);
//			std::cout << "\n added node types:\n" << n.addedNodes.type << "\n" << std::endl;

			getPhis(n, iter_surf);

			if(lvl != 0){
				getDefVec(defVec_all, n, g.errorPrevLvl, n.N_c);

			}else{

				if(n.addedNodes.type.minCoeff() == 0){
					getDefVec(defVec_m, n.N_m, n.mPtr);
				}


				if(n.addedNodes.type.minCoeff() <= 1){
					pseudo_sliding_edge(PhiPtr, n);
					getDefVec(defVec_me, defVec_m, n, proj_disp_edge, n.N_m+n.N_se, n.N_m);
				}


				pseudo_sliding_surf(PhiPtr, n);
				getDefVec(defVec_all, defVec_me, n, proj_disp_all, n.N_c, n.N_m+n.N_se);

			}
			performRBF(PhiPtr->Phi_cc,PhiPtr->Phi_ic,defVec_all,n.cPtr,n.iPtr,n.N_c);

			g.getError(n,d,lvl);
			std::cout << "error: \t"<< g.maxError <<" at node: \t" << g.maxErrorNodes(0)<< std::endl;

			if(g.maxError < params.tol){
				iterating = false;
				g.maxErrorNodes.resize(0);
			}


			if(params.multiLvl && (g.maxError/g.maxErrorPrevLvl < params.tolCrit || iterating == false)){

				g.setLevelParams(n,lvl, d, alpha, defVec_all, n.cPtr, n.N_c);
				std::cout << "LEVEL: " << lvl << " HAS BEEN DONE" << std::endl;
				lvl++;
				iter_surf = -1;

				n.assignNodeTypesGrdy(m);

				if(iterating == false){
					g.getAlphaVector();
					g.setInitMaxErrorNodes();
				}

			}
			iter_surf++;
		}


		updateNodes(n,defVec_all, g.d_step, g.alpha_step, g.ctrlPtr);
		g.correction(m,n, params.gamma, params.multiLvl);
		iterating = true;

		std::cout << "number of control nodes: " << n.N_c << std::endl;

	}
}



void rbf_ps::pseudo_sliding_surf(PhiStruct* PhiPtr, getNodeType& n){
	free_disp_all.resize(n.N_ss, m.nDims);
	proj_disp_all.resize(n.N_ss, m.nDims);

	for(int dim = 0; dim < m.nDims; dim++){
		free_disp_all(Eigen::seqN(0,n.N_ss),dim) = (PhiPtr->Phi_sme*(PhiPtr->Phi_meme.fullPivLu().solve(defVec_me(Eigen::seqN(dim*(n.N_m+n.N_se),n.N_m+n.N_se))))).array();
	}

	p.projectSurf(m, n.ssPtr, free_disp_all, proj_disp_all, m.periodicVec, 0, n.N_ss, 1);
}


void rbf_ps::pseudo_sliding_edge(PhiStruct* PhiPtr, getNodeType& n){
	free_disp_edge.resize(n.N_se, m.nDims);
	proj_disp_edge.resize(n.N_se,m.nDims);
	for(int dim = 0; dim < m.nDims; dim++){
		free_disp_edge(Eigen::seqN(0,n.N_se),dim) = (PhiPtr->Phi_em*(PhiPtr->Phi_mm.fullPivLu().solve(defVec_m(Eigen::seqN(dim*n.N_m,n.N_m))))).array();
	}


	p.projectEdge(m, n.sePtr, free_disp_edge, proj_disp_edge, m.periodicVec, 0, n.N_se, 1);

}



