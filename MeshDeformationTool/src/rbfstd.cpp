#include "rbfstd.h"
#include <Eigen/Dense>
#include <chrono>
#include <iostream>
#include "greedy.h"
#include "WriteResults.h"
#include <ctime>
//rbf_std::rbf_std(Mesh& meshObject, Eigen::VectorXd& dVec, Eigen::RowVectorXd& rotPnt, Eigen::VectorXd& rotVec, const int& steps, const std::string& smode, const bool& curved, const std::string& pDir)
rbf_std::rbf_std(struct probParams& probParamsObject, Mesh& meshObject, getNodeType& n)
:rbfGenFunc(meshObject, probParamsObject)
{
	if(params.dataRed){
		greedy g(m, params, exactDispPtr, movingIndices, alpha, d);
		perform_rbf(n,g);
	}else{
		perform_rbf(n);
	}
}


void rbf_std::perform_rbf(getNodeType& n){
	if(params.ptype){
		transform.cart_to_polar_cylindrical(m.coords, m.coords_polar_cylindrical);
	}

	for(int i=0; i<params.steps; i++){
		std::cout << "Deformation step: " << i+1 << " out of "<< params.steps << std::endl;
		getPhis(n, 0);
		if(i == 0){
			getDefVec(defVec, n.N_c, n.cPtr);

		}
		performRBF(Phis.Phi_cc, Phis.Phi_ic, defVec, n.cPtr, n.iPtr, n.N_m);
	}


	std::cout << "number of control nodes: " << n.N_m << std::endl;
	if(params.ptype){
		transform.polar_cylindrical_to_cart(m.coords_polar_cylindrical, m.coords );
	}
}

void rbf_std::perform_rbf(getNodeType& n, greedy& g){
	std::cout << "Performing standard rbf interpolation w/ greedy" << std::endl;


	WriteResults w;
	w.createConvHistFile(params.convHistFile);




	int iter, lvl;
	bool iterating = true;

	if(params.ptype){
		transform.cart_to_polar_cylindrical(m.coords, m.coords_polar_cylindrical);
	}


	for(int i=0; i < params.steps; i++){
		std::cout << "Deformation step: " << i+1 << " out of "<< params.steps << std::endl;

		iter = 0;
		lvl = 0;

		while(iterating){


			n.addControlNodes(g.maxErrorNodes, params.smode, m);

			getPhis(n, iter);

			if(lvl > 0){
				getDefVec(defVec, n, g.errorPrevLvl, n.N_m);
			}else{
				getDefVec(defVec, n.N_c, n.cPtr);
			}


			performRBF(Phis.Phi_cc, Phis.Phi_ic, defVec, n.cPtr, n.iPtr, n.N_m);


			g.getError(n, d, lvl);
			std::cout << "error: \t"<< g.maxError <<" at node: \t" << g.maxErrorNodes(0)<< std::endl;


//				auto stop = std::chrono::high_resolution_clock::now();
//				auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
//				w.setIntResults(i, lvl, maxError, abs(go.error).mean(), duration.count()/1e6, params.convHistFile, n.N_c);

			if(g.maxError < params.tol){
				iterating = false;
				if(params.multiLvl == false){
					g.maxErrorNodes.resize(0);
				}
			}

			if(params.multiLvl && (g.maxError/g.maxErrorPrevLvl < params.tolCrit || iterating == false)){
				std::cout << "multi lvl criterium reached\n";
				g.setLevelParams( n, lvl, d, alpha, defVec, n.cPtr, n.N_c);


//				w.setIntResults(i, lvl, maxError, abs(go.error).mean(), duration.count()/1e6, params.convHistFile, n.N_m);
				std::cout << "LEVEL: " << lvl << " HAS BEEN DONE" <<std::endl;


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



		// related to the data reduction
		updateNodes(n, defVec, g.d_step, g.alpha_step, g.ctrlPtr);
		g.correction(m,n,params.gamma, params.multiLvl);

		iterating = true;

	}
	std::cout << "Number of different control nodes: " << g.ctrlNodesAll.size() << std::endl;
	std::cout << "number of control nodes: " << n.N_c << std::endl;

}







