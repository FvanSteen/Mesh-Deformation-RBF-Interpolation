#include "rbfstd.h"
#include <Eigen/Dense>
#include <chrono>
#include <iostream>
#include "greedy.h"
#include <ctime>


rbf_std::rbf_std(struct probParams& probParamsObject, Mesh& meshObject, getNodeType& n)
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


void rbf_std::perform_rbf(getNodeType& n){
	std::cout << "\nPerforming regular RBF interpolation\n" << std::endl;
	// transform to polar/ cylindrical coordinates if required
	if(params.ptype){
		transform.cart_to_polar_cylindrical(m.coords, m.coords_polar_cylindrical);
	}

	// for-loop going through the deformation steps
	for(int i=0; i<params.steps; i++){
		std::cout << "Deformation step: " << i+1 << " out of "<< params.steps << std::endl;

		// obtaining the interpolation matrices
		getPhis(n, 0);

		// if its the first step then the deformation vector has to be established
		if(i == 0){
			getDefVec(defVec, n.N_c, n.cPtr);
		}

		// performing the RBF interpolation
		performRBF(Phis.Phi_cc, Phis.Phi_ic, defVec, n.cPtr, n.iPtr, n.N_m);
	}

	// transform polar/cylindrical coordinates back to Cartesian coordinates
	if(params.ptype){
		transform.polar_cylindrical_to_cart(m.coords_polar_cylindrical, m.coords );
	}
}

void rbf_std::perform_rbf(getNodeType& n, greedy& g){
	std::cout << "\nPerforming regular RBF interpolation with data reduction\n" << std::endl;

	// clocks for keeping track of CPU time
	std::clock_t s = std::clock();
	std::clock_t e;


	// integers for number of iterations and levels
	int iter, lvl;

	// iterating bool set to True
	bool iterating = true;


	// for-loop going through the deformation steps
	for(int i=0; i < params.steps; i++){
		std::cout << "Deformation step: " << i+1 << " out of "<< params.steps << std::endl;

		// transformation to polar/ cylindrical coordinates
		if(params.ptype){
			transform.cart_to_polar_cylindrical(m.coords, m.coords_polar_cylindrical);
		}

		iter = 0;
		lvl = 0;

		while(iterating){
			// adding control nodes
			n.addControlNodes(g.maxErrorNodes, params.smode, m);

			// obtaining interpolation matrices
			getPhis(n, iter);

			// obtaining deformation vectors
			if(lvl > 0){
				getDefVec(defVec, n, g.errorPrevLvl, n.N_m);
			}else{
				getDefVec(defVec, n.N_c, n.cPtr);
			}


			// performing the RBF interpolation
			performRBF(Phis.Phi_cc, Phis.Phi_ic, defVec, n.cPtr, n.iPtr, n.N_m);

			// obtaining error
			g.getError(n, d, lvl);
			std::cout << "error: \t"<< g.maxError <<" at node: \t" << g.maxErrorNodes(0)<< std::endl;


			// writing intermediate results of the convergence history
			e = std::clock();
			long double time_elapsed_ms =  1000.0*(e-s) / CLOCKS_PER_SEC;
			std::cout << "CPU time: " << time_elapsed_ms/1000 << " ms\n";
			w.setIntResults(params.directory, i, lvl, g.maxError, time_elapsed_ms, n.N_c);


			// check if the error tolerance is reached
			if(g.maxError < params.tol){
				iterating = false;
				if(params.multiLvl == false){
					g.maxErrorNodes.resize(0);
				}
			}

			// check for the multi-level criteria
			if(params.multiLvl && ( (params.mCrit == "size" && n.N_c == params.lvlSize) ||  (params.mCrit == "tol" &&  g.maxError/g.maxErrorPrevLvl < params.tolCrit) || iterating == false)){

				// saving lvl parameters
				g.setLevelParams( n, lvl, d, alpha, defVec, n.cPtr, n.N_c);
				std::cout << "Level: " << lvl << " has been performed" <<std::endl;


				lvl++;
				iter = -1;

				// reset the nodetypes for the next level
				n.assignNodeTypesGrdy(m);

				// if the tolerance is reached,
				if(iterating == false){
					g.getAlphaVector();
					g.setInitMaxErrorNodes();
				}
			}
			iter++;
		}

		// performing the correction and updating the node coordinates
		g.correction(m,n,params.gamma, params.multiLvl);
		updateNodes(n, defVec, g.d_step, g.alpha_step, g.ctrlPtr);

		iterating = true;

	}
}







