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
		greedy g(m, params, exactDispPtr, movingIndices, alpha, d);
		perform_rbf(n,g);
	}else{
		perform_rbf(n);
	}

}

void rbf_ps::perform_rbf(getNodeType& n){
	std::cout<< "Performing RBF PS without data reduction" << std::endl;

	if(params.ptype){
		transform.cart_to_polar_cylindrical(m.coords, m.coords_polar_cylindrical);
	}
	for(int i = 0; i < params.steps; i++){
		std::cout << "Deformation step: " << i+1 << " out of "<< params.steps << std::endl;

		if(params.curved || i==0){
			m.getMidPnts(params);
		}

		getPhis(n, 0);

		if(i == 0){
			getDefVec(defVec_m, n.N_m, n.mPtr);
		}

		// find the finaldef of the edge nodes
		pseudo_sliding_edge(PhiPtr, n);

//		m.coords(*n.sePtr, Eigen::all) += proj_disp_edge;
//		for(int dim = 0; dim <m.nDims; dim++){
//			m.coords(*n.mPtr, dim) += defVec_m(Eigen::seqN(dim*n.N_m,n.N_m)).array();
//		}
//		m.writeMeshFile(params.mesh_ifName, params.mesh_ofName);
//		std::exit(0);

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

	if(params.ptype){
		transform.polar_cylindrical_to_cart(m.coords_polar_cylindrical,m.coords);
	}
}

void rbf_ps::perform_rbf(getNodeType& n, greedy& g){
	std::cout<< "Performing RBF PS with data reduction" << std::endl;


	int iter, lvl;
	bool iterating = true;

	for(int i = 0; i < params.steps; i++){
		std::cout << "Deformation step: " << i+1 << " out of "<< params.steps << std::endl;
		if(params.ptype){
			transform.cart_to_polar_cylindrical(m.coords, m.coords_polar_cylindrical);
		}

		iter = 0;
		lvl = 0;

		if(params.curved || i==0){
			m.getMidPnts(params);
//			m.getVecs(params);
		}

		while(iterating){
			std::cout << "iter: " << iter << std::endl;
			n.addControlNodes(g.maxErrorNodes, params.smode, m);
//			std::cout << "control nodes:\n";
//			for(auto x : *n.cPtr){
//				std::cout << x << ", ";
//			}
//
//			std::cout << std::endl;

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


			g.getError(n, d, lvl);
			std::cout << "error: \t"<< g.maxError <<" at node: \t" << g.maxErrorNodes(0)<< std::endl;

//			for(auto x : *n.iPtr_reduced){
//				std::cout << x << ", ";
//			}
//			std::cout << std::endl;



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

//				if(lvl == 1){
//					g.getAlphaVector();
//					std::cout << *g.d_step << std::endl;
//					std::cout << (*g.d_step).rows() << std::endl;
//					std::cout << "\n\n" << g.errorPrevLvl << std::endl;
////					updateNodes(n,defVec_all, g.d_step, g.alpha_step, g.ctrlPtr);
//					m.coords_polar_cylindrical(*n.iPtr,Eigen::all) += (*g.d_step - g.errorPrevLvl);
////					for(int dim = 0; dim <m.nDims; dim++){
////						m.coords_polar_cylindrical
////					}
//					if(params.ptype){
//						transform.polar_cylindrical_to_cart(m.coords_polar_cylindrical, m.coords);
//					}
//	//				g.correction(m,n, params.gamma, params.multiLvl);
//					m.writeMeshFile(params.mesh_ifName, params.mesh_ofName);
//					std::exit(0);
//				}

			}
			iter++;

		}



//		for(auto x : m.seNodes){
//			std::cout << x << ", ";
//		}
//		std::cout << std::endl;
//
//		for(auto x : *n.cPtr){
//			std::cout << x << ", ";
//		}
//		std::cout << std::endl;
//
//		(*m.ptrCoords)(*n.iPtr_reduced,Eigen::all) += d;
////				(*m.ptrCoords)(*n.sePtr,Eigen::all) += free_disp_edge;
//		for(int i = 0; i < m.nDims; i++){
//			(*m.ptrCoords)(*n.cPtr, i) += defVec_me(Eigen::seqN(i*n.N_c, n.N_c)).array();
//		}
//
//
//		transform.polar_cylindrical_to_cart(m.coords_polar_cylindrical, m.coords);
//		m.writeMeshFile(params.mesh_ifName, params.mesh_ofName);
//		std::exit(0);


		std::cout << "MOVING TO THE SURFACE NODES\n";





		int iter_surf = 0;

		if(m.nDims == 3){
			iterating = true;
			pseudo_sliding_surf(PhiPtr,n);
			getDefVec(defVec_all, defVec_me, n, proj_disp_all, n.N_c, n.N_m+n.N_se);
			performRBF(PhiPtr->Phi_cc,PhiPtr->Phi_ic,defVec_all,n.cPtr,n.iPtr,n.N_c);

			g.getError(n,d,lvl);




		}

		while(iterating){

			std::cout << "iter_surf: " << iter_surf << std::endl;

			n.addControlNodes(g.maxErrorNodes, params.smode, m);
//			std::cout << "\n added node types:\n" << n.addedNodes.type << "\n" << std::endl;
			std::cout << "control nodes:\n";
			for(auto x : *n.cPtr){
				std::cout << x << ", ";
			}
			std::cout << std::endl;
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


			/*if(iter_surf==6){
				std::cout << "Moving control nodes:\n";
				for(auto x : *n.mPtr){
					std::cout << x << ", ";
				}
				std::cout << std::endl;

				std::cout << "Sliding edge control nodes:\n";
				for(auto x : *n.sePtr){
					std::cout << x << ", ";
				}
				std::cout << std::endl;

				std::cout << "Sliding surface control nodes:\n";
				for(auto x : *n.ssPtr){
					std::cout << x << ", ";
				}
				std::cout << std::endl;

				for(int i =0;i<m.nDims; i++){
					m.coords_polar_cylindrical(*n.cPtr,i) += defVec_all(Eigen::seqN(i*n.N_c, n.N_c)).array();
				}

				transform.polar_cylindrical_to_cart(m.coords_polar_cylindrical, m.coords);
	//			m.coords(m.ssNodes,Eigen::all) += (d(Eigen::lastN(m.N_ss),Eigen::all);
				m.coords(*n.iPtr, Eigen::all) += (d);

				m.writeMeshFile(params.mesh_ifName, params.mesh_ofName);
				std::cout << "DONE\n";

				std::exit(0);

			}*/


//			if(iter_surf == 10){
//
//				std::cout << "sliding edge nodes\n";
//				for(auto x : *n.sePtr){
//					std::cout << x << ", ";
//				}
//				std::cout << std::endl;
//				std::cout << "moving nodes\n";
//				for(auto x : *n.mPtr){
//					std::cout << x << ", ";
//				}
//				std::cout << std::endl;
//				std::cout << "sliding surface nodes\n";
//				for(auto x : *n.ssPtr){
//					std::cout << x << ", ";
//				}
//				std::cout << std::endl;
//				for(int i = 0; i < m.nDims; i++){
//					(*m.ptrCoords)(*n.mPtr, i) += defVec_m(Eigen::seqN(i*n.N_m, n.N_m)).array();
//				}
//
//				m.coords_polar_cylindrical(*n.sePtr, Eigen::all) += proj_disp_edge;
//				m.coords_polar_cylindrical(*n.ssPtr, Eigen::all) += proj_disp_all;
//				transform.polar_cylindrical_to_cart(m.coords_polar_cylindrical, m.coords);
//
////				m.coords(*n.ssPtr, Eigen::all) += free_disp_all;
//				m.coords(*n.iPtr, Eigen::all) += (d-g.error);
//				m.writeMeshFile(params.mesh_ifName, params.mesh_ofName);
//				std::cout << "DONE\n";
//				std::exit(0);
//			}

//			if(iter_surf == 0){
//				for(auto x : *n.sePtr){
//					std::cout << x << ", ";
//				}
//				std::cout << std::endl;
//
//				for(auto x : *n.mPtr){
//					std::cout << x << ", ";
//				}
//				std::cout << std::endl;
//
//				for(auto x : *n.ssPtr){
//					std::cout << x << ", ";
//				}
//				std::cout << std::endl;
//
//				for(auto x : *n.cPtr){
//					std::cout << x << ", ";
//				}
//				std::cout << std::endl;
//
////				(*m.ptrCoords)(*n.iPtr,Eigen::all) += d;
//				(*m.ptrCoords)(*n.sePtr,Eigen::all) += proj_disp_edge;
//				(*m.ptrCoords)(*n.ssPtr,Eigen::all) += proj_disp_all;
//				for(int i = 0; i < m.nDims; i++){
//					(*m.ptrCoords)(*n.mPtr, i) += defVec_m(Eigen::seqN(i*n.N_m, n.N_m)).array();
//				}


//				transform.polar_cylindrical_to_cart(m.coords_polar_cylindrical, m.coords);
//				m.writeMeshFile(params.mesh_ifName, params.mesh_ofName);
//				std::exit(0);
//			}


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
//		std::string tempFile = "stator_per_ffd_mod_uncorrected.su2";
//		m.writeMeshFile(params.mesh_ifName, tempFile);
		g.correction(m,n, params.gamma, params.multiLvl);
		std::cout << "sliding edge nodes\n";
		for(auto x : *n.sePtr){
			std::cout << x << ", ";
		}
		std::cout << std::endl;
		std::cout << "moving nodes\n";
		for(auto x : *n.mPtr){
			std::cout << x << ", ";
		}
		std::cout << std::endl;
		std::cout << "sliding surface nodes\n";
		for(auto x : *n.ssPtr){
			std::cout << x << ", ";
		}

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

	if(params.ptype){
		transform.disp_to_cart(free_disp_all, *n.ssPtr, n.N_ss, m);
	}

	p.projectSurf(m, n.ssPtr, free_disp_all, proj_disp_all, 0, n.N_ss, 1, params.ptype);

	if(params.ptype){
		transform.disp_to_polar_cylindrical(proj_disp_all, *n.ssPtr, n.N_ss, m);
	}
}


void rbf_ps::pseudo_sliding_edge(PhiStruct* PhiPtr, getNodeType& n){
	free_disp_edge.resize(n.N_se, m.nDims);
	proj_disp_edge.resize(n.N_se,m.nDims);
	for(int dim = 0; dim < m.nDims; dim++){
		free_disp_edge(Eigen::seqN(0,n.N_se),dim) = (PhiPtr->Phi_em*(PhiPtr->Phi_mm.fullPivLu().solve(defVec_m(Eigen::seqN(dim*n.N_m,n.N_m))))).array();
	}

	if(params.ptype){
		transform.disp_to_cart(free_disp_edge, *n.sePtr, n.N_se, m);
	}

	p.projectEdge(m, n.sePtr, free_disp_edge, proj_disp_edge, 0, n.N_se, 1, params.ptype);

	if(params.ptype){
		transform.disp_to_polar_cylindrical(proj_disp_edge, *n.sePtr, n.N_se, m);
	}
}



