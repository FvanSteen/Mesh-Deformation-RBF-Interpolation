#include "rbfstd.h"
#include <Eigen/Dense>
#include <chrono>
#include <iostream>
#include "greedy.h"
#include "WriteResults.h"
//rbf_std::rbf_std(Mesh& meshObject, Eigen::VectorXd& dVec, Eigen::RowVectorXd& rotPnt, Eigen::VectorXd& rotVec, const int& steps, const std::string& smode, const bool& curved, const std::string& pDir)
rbf_std::rbf_std(struct probParams& probParamsObject, Mesh& meshObject, getNodeType& n)
:rbfGenFunc(meshObject, probParamsObject)
{
	perform_rbf(n);
}


void rbf_std::perform_rbf(getNodeType& n){
	std::cout << "Performing standard rbf interpolation" << std::endl;

	WriteResults w;
	w.createConvHistFile(params.convHistFile);

	projection p(pVec);

	auto start = std::chrono::high_resolution_clock::now();\

	// matrices and vector required to solve the rbf interpolation
	Eigen::MatrixXd Phi_cc, Phi_ic, Phi_icGrdy;
	Eigen::VectorXd defVec;

	double maxError;

	int iter, lvl;
	bool iterating;


	Eigen::ArrayXi maxErrorNodes;
	greedy go;

//	int lvlSize = 16;

	Eigen::VectorXd* alpha_step;
	Eigen::ArrayXXd* d_step;
	Eigen::ArrayXi* ctrlPtr;
	if(params.dataRed){
		if(params.multiLvl){
			alpha_step = &go.alphaGrdy;
			d_step = &go.delta;
			ctrlPtr = &go.ctrlNodesAll;
		}else{
			alpha_step = &alpha;
			d_step = &d;
		}
	}

	m.r = 10;
	for(int i=0; i<params.steps; i++){
		std::cout << "Deformation step: " << i+1 << " out of "<< params.steps << std::endl;
		maxError = 1;
		iter = 0;
		lvl = 0;


		if((params.dataRed && i==0) || params.multiLvl ){
			go.setInitMaxErrorNodes(m, m.coords, exactDisp, movingIndices, maxErrorNodes, params.doubleEdge);
			std::cout << "initial selected Nodes:\n" << maxErrorNodes << std::endl;
		}

		iterating = true;

		while(iterating){

			if(params.dataRed){
				for(int node = 0; node < maxErrorNodes.size(); node++){
//					if(n.N_c < lvlSize){
						n.addControlNode(maxErrorNodes(node), params.smode, m);
//					}
				}
			}

			getPhis(Phi_cc, Phi_ic, n.cPtr, n.iPtr);


			if(lvl > 0){
				getDefVec(defVec, n, go.errorPrevLvl, n.N_c);
			}else if(i==0 || params.dataRed){
				getDefVec(defVec, n.N_c, n.cPtr);
			}



			performRBF(Phi_cc, Phi_ic, defVec,n.cPtr, n.iPtr, n.N_c);


			if(params.dataRed){
				go.getError(m,n,d,maxError,maxErrorNodes, movingIndices, exactDisp,pnVec,p, params.multiLvl, lvl, params.doubleEdge);
				std::cout << "error: \t"<< maxError <<" at node: \t" << maxErrorNodes(0)<< std::endl;


				auto stop = std::chrono::high_resolution_clock::now();
				auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
				w.setIntResults(i, lvl, maxError, abs(go.error).mean(), duration.count()/1e6, params.convHistFile, n.N_c);


				if(maxError < params.tol){ // first part is added
					iterating = false;
					maxErrorNodes.resize(0);
				}

			}else{

				iterating = false;
			}
//			if(lvl == 12 && n.N_m == 8){
//				std::cout << *n.mPtr << std::endl;
////				std::cout << go.error.matrix().rowwise().norm().transpose() << std::endl;
//				std::cout << "Level: " << lvl << '\t' << "Nc: " << n.N_m << std::endl;
////				go.getAlphaVector();
////				updateNodes(Phi_imGreedy,n, defVec, d_step, alpha_step,ctrlPtr);
//				go.setLevelParams(m,n,lvl,n.N_m, d, alpha, maxError);
//				m.coords(*n.iPtr,Eigen::all) += go.delta;
//				m.writeMeshFile();
//				std::exit(0);
//			}



//			if(params.multiLvl && n.N_c == lvlSize){
//				std::cout << "levelsize: " << lvlSize << std::endl;
//				if(maxError > 0.5*go.maxErrorPrevLvl && lvl !=0){
//					lvlSize += 16;
//				}else{
			if(params.multiLvl && (maxError/go.maxErrorPrevLvl < params.tolCrit || iterating == false)){

//				std::cout << go.maxErrorPrevLvl << '\t' << maxError << '\t' << maxError/go.maxErrorPrevLvl << std::endl;
//				std::cout << "Nm: " << n.N_m << std::endl;
//				for(int a = 0; a< n.N_m; a++){
//					std::cout << (*n.mPtr)(a) << ", ";
//				}
//				std::cout << std::endl;
//				m.coords(*n.iPtr,Eigen::all) += d;


//				if(maxError > 0.5*go.maxErrorPrevLvl && lvl !=0){
////					std::cout << "level: " << lvl << "\nerror: " << maxError << "\nPrevious error: " << go.maxErrorPrevLvl << std::endl;
//					params.lvlSize +=4;// params.lvlSize*2;
////					std::cout << params.lvlSize << '\t' << params.lvlSizeInit <<std::endl;
//
//				}else{


					go.setLevelParams(m,n, lvl, d, alpha, maxError, defVec, n.cPtr, n.N_c);

	//				if(lvl == 4){
	//					m.coords(*n.iPtr, Eigen::all) += *d_step;
	//					m.writeMeshFile();
	//					std::exit(0);
	//				}
					auto stop = std::chrono::high_resolution_clock::now();
					auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);

	//				w.setIntResults(i, lvl, maxError, abs(go.error).mean(), duration.count()/1e6, params.convHistFile, n.N_m);
	//				std::cout << "LEVEL: " << lvl << " HAS BEEN DONE" <<std::endl;

		//				std::cout << "mean error: " << abs(go.error).mean() << std::endl;

					lvl++;

					n.assignNodeTypesGrdy(m);

					if(maxError < params.tol){
						iterating = false;
						go.getAlphaVector();
					}
//					lvlSize = 16;

//					std::cout << "Maximum deformation: "<<  maxError << std::endl;

//					m.r = 2.5*maxError/0.0780524;
//					m.r = 100*maxError/0.707107;
//					std::cout<< m.r << std::endl;
//				}
			}

			iter++;

		}


		if(params.dataRed){

//			if(params.multiLvl == 0){
//
//				auto stop = std::chrono::high_resolution_clock::now();
//				auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
//				w.setIntResults(i, lvl, maxError, abs(go.error).mean(), duration.count()/1e6, params.convHistFile);
//			}


			updateNodes(Phi_icGrdy,n,defVec, d_step, alpha_step,ctrlPtr);
			go.correction(m,n,params.gamma, params.multiLvl);
		}

	}
	std::cout << "Number of different control nodes: " << go.ctrlNodesAll.size() << std::endl;
	std::cout << "number of control nodes: " << n.N_c << std::endl;
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);

	std::cout << "Runtime duration: \t"<<  duration.count()/1e6 << " seconds"<< std::endl;
}







