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
		greedy g(m, params, exactDisp, movingIndices, alpha, d);
		perform_rbf(n,g);
	}else{
		perform_rbf(n);
	}
}


void rbf_std::perform_rbf(getNodeType& n){
	std::clock_t s = std::clock();
	Eigen::VectorXd defVec;

	for(int i=0; i<params.steps; i++){
		std::cout << "Deformation step: " << i+1 << " out of "<< params.steps << std::endl;
		getPhis(n, 0);
		if(i == 0){
			getDefVec(defVec, n.N_m, n.mPtr);
		}
		performRBF(Phis.Phi_cc, Phis.Phi_ic, defVec, n.mPtr, n.iPtr, n.N_m);
	}

	std::cout << "number of control nodes: " << n.N_m << std::endl;
	std::clock_t e = std::clock();
	long double time_elapsed_ms =  1000.0*(e-s) / CLOCKS_PER_SEC;
	std::cout << "CPU time: " << time_elapsed_ms/1000 << " ms\n";

}

void rbf_std::perform_rbf(getNodeType& n, greedy& g){
	std::cout << "Performing standard rbf interpolation w/ greedy" << std::endl;


	WriteResults w;
	w.createConvHistFile(params.convHistFile);

	std::clock_t s = std::clock();

	Eigen::VectorXd defVec;


	int iter, lvl;
	bool iterating = true;



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
				getDefVec(defVec, n.N_m, n.mPtr);
			}


			performRBF(Phis.Phi_cc, Phis.Phi_ic, defVec, n.mPtr, n.iPtr, n.N_m);


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

				g.setLevelParams( n, lvl, d, alpha, defVec, n.mPtr, n.N_m);


//				w.setIntResults(i, lvl, maxError, abs(go.error).mean(), duration.count()/1e6, params.convHistFile, n.N_m);
				std::cout << "LEVEL: " << lvl << " HAS BEEN DONE" <<std::endl;


				lvl++;
				iter = -1;
				// reset the nodetypes for the next level
				n.assignNodeTypesGrdy(m);

				// if the tolerance is reached,
				if(iterating == false){
					g.getAlphaVector();
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
	std::cout << "number of control nodes: " << n.N_m << std::endl;


	std::clock_t e = std::clock();
	long double time_elapsed_ms =  1000.0*(e-s) / CLOCKS_PER_SEC;
	std::cout << "CPU time: " << time_elapsed_ms/1000 << " ms\n";

}







