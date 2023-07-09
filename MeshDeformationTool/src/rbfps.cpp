#include "rbfps.h"

#include <iostream>
#include <chrono>
#include <Eigen/Dense>
#include "greedy.h"


//rbf_ps::rbf_ps(Mesh& meshObject, Eigen::VectorXd& dVec, Eigen::RowVectorXd& rotPnt, Eigen::VectorXd& rotVec, const int& steps, const std::string& smode, const bool& curved, const std::string& pDir)
rbf_ps::rbf_ps( struct probParams& probParamsObject, Mesh& meshObject, getNodeType& n)
:rbfGenFunc(meshObject, probParamsObject)
{

	// if data reduction is applied then initialise the class with greedy function and perform the rbf interpolation
	if(params.dataRed){
		// creating a convergence history file
		w.createConvHistFile(params.directory);
		greedy g(m, params, exactDispPtr, movingIndices, alpha, d);
		perform_rbf(n,g);
	}else{
		perform_rbf(n);
	}

}

void rbf_ps::perform_rbf(getNodeType& n){
	std::cout<< "\nPerforming Pseudo sliding RBF interpolation\n" << std::endl;


	// transform to polar/ cylindrical coordinates if required
	if(params.ptype){
		transform.cart_to_polar_cylindrical(m.coords, m.coords_polar_cylindrical);
	}

	// for-loop going through the deformation steps
	for(int i = 0; i < params.steps; i++){
		std::cout << "Deformation step: " << i+1 << " out of "<< params.steps << std::endl;

		// obtaining the midpoints and vectors of the boundary elements
		if(params.curved || i==0){
			m.getMidPnts(params);
		}

		// obtaining the interpolation matrices
		getPhis(n, 0);

		// obtaining the deformation vector
		if(i == 0){
			getDefVec(defVec_m, n.N_m, n.mPtr);
		}

		// perform the sliding of the edge nodes
		pseudo_sliding_edge(PhiPtr, n);


		if(m.nDims  == 3){
			// setup new deformation vector that includes the found projection of the edge nodes
			getDefVec(defVec_me, defVec_m, n, proj_disp_edge, n.N_m+n.N_se, n.N_m);

			// perform the sliding of the surface nodes
			pseudo_sliding_surf(PhiPtr, n);

			// setup new deformation vector that includes the found projectino of the surface nodes
			getDefVec(defVec_all, defVec_me, n, proj_disp_all, n.N_c, n.N_m+n.N_se);

		}else{
			// setup new deformation vector that includes the found projectino of the edge nodes
			getDefVec(defVec_all, defVec_m, n, proj_disp_edge, n.N_c, n.N_m);
		}

		// perform regular RBF interpolation with the found projected position in the deformation vector
		performRBF(PhiPtr->Phi_cc,PhiPtr->Phi_ic,defVec_all,n.cPtr,n.iPtr,n.N_c);

		// transform to Cartesian coordinates
		if(params.ptype){
			transform.polar_cylindrical_to_cart(m.coords_polar_cylindrical,m.coords);
		}

	}
}

void rbf_ps::perform_rbf(getNodeType& n, greedy& g){
	std::cout<< "\nPerforming Pseudo sliding RBF interpolation with data reduction\n" << std::endl;

	// clcoks for keeping track of CPU time
	std::clock_t s = std::clock();
	std::clock_t e;


	// integers for number of iterations and levels
	int iter, lvl;

	// iterating bool set to True
	bool iterating = true;

	// for-loop going through the deformation steps
	for(int i = 0; i < params.steps; i++){
		std::cout << "Deformation step: " << i+1 << " out of "<< params.steps << std::endl;

		// transformation to polar/ cylindrical coordinates
		if(params.ptype){
			transform.cart_to_polar_cylindrical(m.coords, m.coords_polar_cylindrical);
		}

		iter = 0;
		lvl = 0;


		// obtaining midpoints and vectors of the boundary elements
		if(params.curved || i==0){
			m.getMidPnts(params);
		}


		// iterating loop
		while(iterating){
			// adding control nodes
			n.addControlNodes(g.maxErrorNodes, params.smode, m);

			// obtaining interpolation matrices
			getPhis(n, iter);


			if(lvl != 0){

				// obtain deformation vector
				getDefVec(defVec_all, n, g.errorPrevLvl, n.N_c);
				// perform RBF interpolation
				performRBF(Phis.Phi_cc, Phis.Phi_ic, defVec_all, n.cPtr, n.iPtr, n.N_c);
			}else{
				// get deformation vector
				getDefVec(defVec_m, n.N_m, n.mPtr);

				// perform sliding of the edge nodes
				pseudo_sliding_edge(PhiPtr, n);

				if(m.nDims == 2){
					// setup deformation vector that includes the projection of the edge nodes
					getDefVec(defVec_all, defVec_m, n, proj_disp_edge, n.N_c, n.N_m);
					// perform RBF interpolation
					performRBF(PhiPtr->Phi_cc,PhiPtr->Phi_ic,defVec_all,n.cPtr,n.iPtr,n.N_c);

				}else if (m.nDims == 3){
					// setup deformation vector that includes the projection of the surface nodes
					getDefVec(defVec_me, defVec_m, n, proj_disp_edge, n.N_m+n.N_se, n.N_m);
					// perform RBF interpolation
					performRBF(PhiPtr->Phi_meme,PhiPtr->Phi_ic_reduced,defVec_me,n.cPtr,n.iPtr_reduced,n.N_m+n.N_se);


				}
			}

			// obtain error
			g.getError(n, d, lvl);
			std::cout << "error: \t"<< g.maxError <<" at node: \t" << g.maxErrorNodes(0)<< std::endl;

			// writing intermediate results of the convergence history
			e = std::clock();
			long double time_elapsed_ms =  1000.0*(e-s) / CLOCKS_PER_SEC;
			std::cout << "CPU time: " << time_elapsed_ms/1000 << " ms\n";
			w.setIntResults(params.directory, i, lvl, g.maxError, time_elapsed_ms, n.N_c);


			// check if error tolerance is reached
			if(g.maxError < params.tol){
				iterating = false;
				if(params.multiLvl == false){
					g.maxErrorNodes.resize(0);
				}
			}

			// check for multi-level criteria
			if(params.multiLvl  && m.nDims == 2  && ((params.mCrit == "size" && n.N_c == params.lvlSize) ||  (params.mCrit == "tol" &&  g.maxError/g.maxErrorPrevLvl < params.tolCrit) || iterating == false)){

				// saving lvl parameters
				g.setLevelParams(n,lvl, d, alpha, defVec_all, n.cPtr, n.N_c);
				std::cout << "Level: " << lvl << " has been performed" <<std::endl;


				lvl++;
				iter = -1;

				// reset the node types for the next level
				n.assignNodeTypesGrdy(m);


				// if the tolerance is reached
				if(iterating == false){
					g.getAlphaVector();
					g.setInitMaxErrorNodes();
				}
			}
			iter++;
		}


		// After the edges have reached the greedy tolerance the surface nodes are considered
		// Single iteration is done to find the maxErrorNodes of the surface
		int iter_surf = 0;

		if(m.nDims == 3){
			iterating = true;

			// performing sliding of the surface nodes
			pseudo_sliding_surf(PhiPtr,n);

			// setup deformation vector that includes the projection of the surface nodes
			getDefVec(defVec_all, defVec_me, n, proj_disp_all, n.N_c, n.N_m+n.N_se);

			// performing regular RBF
			performRBF(PhiPtr->Phi_cc,PhiPtr->Phi_ic,defVec_all,n.cPtr,n.iPtr,n.N_c);

			// obtaining error
			g.getError(n,d,lvl);
			std::cout << "error: \t"<< g.maxError <<" at node: \t" << g.maxErrorNodes(0)<< std::endl;



		}

		// iterating until tolerance is reached
		while(iterating){

			// adding control nodes
			n.addControlNodes(g.maxErrorNodes, params.smode, m);

			// obtain interpolation matrices
			getPhis(n, iter_surf);


			if(lvl != 0){
				// obtain deformation vector
				getDefVec(defVec_all, n, g.errorPrevLvl, n.N_c);

			}else{
				if(n.addedNodes.type.minCoeff() == 0){
					// if moving nodes are added then deformation vector has to be updated
					getDefVec(defVec_m, n.N_m, n.mPtr);
				}


				// if moving or edge nodes are added as control nodes then the projection of the edges has to be re-performed
				if(n.addedNodes.type.minCoeff() <= 1){
					// perform sliding of the edge
					pseudo_sliding_edge(PhiPtr, n);

					// setup deformation vector that includes projection of the edge nodes
					getDefVec(defVec_me, defVec_m, n, proj_disp_edge, n.N_m+n.N_se, n.N_m);
				}

				// perform sliding of the surface
				pseudo_sliding_surf(PhiPtr, n);

				// setup deformation vector that includes the projection of the surfaces
				getDefVec(defVec_all, defVec_me, n, proj_disp_all, n.N_c, n.N_m+n.N_se);
			}

			// perform regular RBF interpolation
			performRBF(PhiPtr->Phi_cc,PhiPtr->Phi_ic,defVec_all,n.cPtr,n.iPtr,n.N_c);

			// Obtaining error
			g.getError(n,d,lvl);
			std::cout << "error: \t"<< g.maxError <<" at node: \t" << g.maxErrorNodes(0)<< std::endl;


			// check if tolerance is reached
			if(g.maxError < params.tol){
				iterating = false;
				g.maxErrorNodes.resize(0);
			}


			// multi-level criteria check
			if(params.multiLvl  && m.nDims == 2  && ((params.mCrit == "size" && n.N_c == params.lvlSize) ||  (params.mCrit == "tol" &&  g.maxError/g.maxErrorPrevLvl < params.tolCrit) || iterating == false)){

				// saving lvl parameters
				g.setLevelParams(n,lvl, d, alpha, defVec_all, n.cPtr, n.N_c);
				std::cout << "Level: " << lvl << " has been done" << std::endl;

				lvl++;
				iter_surf = -1;

				// reset of node types for next level
				n.assignNodeTypesGrdy(m);

				// if tolerance is reached
				if(iterating == false){
					g.getAlphaVector();
					g.setInitMaxErrorNodes();
				}

			}
			iter_surf++;
		}


		// performing correction and updating the node coordinates
		g.correction(m,n, params.gamma, params.multiLvl);
		updateNodes(n,defVec_all, g.d_step, g.alpha_step, g.ctrlPtr);

		iterating = true;
	}
}

/* pseudo_sliding_surf
 *
 * function that performs the sliding of the surface nodes
 */

void rbf_ps::pseudo_sliding_surf(PhiStruct* PhiPtr, getNodeType& n){

	// resizing of the displacement arrays
	free_disp_all.resize(n.N_ss, m.nDims);
	proj_disp_all.resize(n.N_ss, m.nDims);


	// determine free displacement of the surface nodes
	for(int dim = 0; dim < m.nDims; dim++){
		free_disp_all(Eigen::seqN(0,n.N_ss),dim) = (PhiPtr->Phi_sme*(PhiPtr->Phi_meme.fullPivLu().solve(defVec_me(Eigen::seqN(dim*(n.N_m+n.N_se),n.N_m+n.N_se))))).array();
	}

	// transform to cartesian coordinates, since projection is done in Cartesian coordinates
	if(params.ptype){
		transform.disp_to_cart(free_disp_all, *n.ssPtr, n.N_ss, m);
	}

	// projection of the surface noeds
	p.projectSurf(m, n.ssPtr, free_disp_all, proj_disp_all, 0, n.N_ss, 1, params.ptype);

	// transform to polar/cylindrical coordinates
	if(params.ptype){
		transform.disp_to_polar_cylindrical(proj_disp_all, *n.ssPtr, n.N_ss, m);
	}
}

/* pseudo_sliding_surf
 *
 * function that performs the sliding of the edge nodes
 */

void rbf_ps::pseudo_sliding_edge(PhiStruct* PhiPtr, getNodeType& n){

	// resizing of the displacement arrays
	free_disp_edge.resize(n.N_se, m.nDims);
	proj_disp_edge.resize(n.N_se,m.nDims);


	// determine free displacement of the edge nodes
	for(int dim = 0; dim < m.nDims; dim++){
		free_disp_edge(Eigen::seqN(0,n.N_se),dim) = (PhiPtr->Phi_em*(PhiPtr->Phi_mm.fullPivLu().solve(defVec_m(Eigen::seqN(dim*n.N_m,n.N_m))))).array();
	}

	// transform to cartesian coordinates, since projection is done in Cartesian coordinates
	if(params.ptype){
		transform.disp_to_cart(free_disp_edge, *n.sePtr, n.N_se, m);
	}

	// projection of the edge nodes
	p.projectEdge(m, n.sePtr, free_disp_edge, proj_disp_edge, 0, n.N_se, 1, params.ptype);

	// transform to polar/cylindrical coordinates
	if(params.ptype){
		transform.disp_to_polar_cylindrical(proj_disp_edge, *n.sePtr, n.N_se, m);
	}
}



