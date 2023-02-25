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

	projection* p;

	projection proObject;
	p = &proObject;

	auto start = std::chrono::high_resolution_clock::now();
	Eigen::MatrixXd Phi_mm, Phi_im, Phi_imGreedy;
	Eigen::VectorXd defVec;

	double maxError;

	int iter, lvl;
	bool iterating;

	Eigen::ArrayXi maxErrorNodes;

	greedy go;



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


	for(int i=0; i<params.steps; i++){
		std::cout << "Deformation step: " << i+1 << " out of "<< params.steps << std::endl;
		maxError = 1;
		iter = 0;
		lvl = 0;


		if((params.dataRed && i==0) || params.multiLvl ){
			go.setInitMaxErrorNodes(m, m.coords, exactDisp, movingIndices, maxErrorNodes);
			std::cout << "initial selected Nodes:\n" << maxErrorNodes << std::endl;
		}

		iterating = true;
//		if(params.multiLvl && i != 0){
//			n.addControlNode(m.intBdryNodes(0));
//			n.addControlNode(m.intBdryNodes(m.intBdryNodes.size()-1));
//		}
		while(iterating){

			if(params.dataRed){
				int node = 0;
				while(node < maxErrorNodes.size()){
					n.addControlNode(maxErrorNodes(node), params.smode, m);
					node++;
				}
			}

			getPhis(Phi_mm, Phi_im, n.mPtr, n.iPtr);

			if(i==0 || params.dataRed){
				getDefVec(defVec, n, lvl, go.errorPrevLvl);
			}

			performRBF(Phi_mm, Phi_im, defVec,*n.mPtr,*n.iPtr, n.N_m);



//			std::cout << d << std::endl;
			if(params.dataRed){
				go.getError(m,n,d,maxError,maxErrorNodes, movingIndices, exactDisp,pnVec,p, params.multiLvl, lvl);
				std::cout << "error: \t"<< maxError <<" at node: \t" << maxErrorNodes(0)<< std::endl;


				auto stop = std::chrono::high_resolution_clock::now();
				auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
				w.setIntResults(i, lvl, maxError, abs(go.error).mean(), duration.count()/1e6, params.convHistFile, n.N_m);


				if(maxError < params.tol){
					iterating = false;
				}

			}else{

//				maxError = 0;
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



//			if(params.multiLvl && n.N_m == params.lvlSize){
			if(params.multiLvl && (maxError/go.maxErrorPrevLvl < 0.5 || iterating == false)){

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


				go.setLevelParams(m,n,lvl,params.lvlSize, d, alpha, maxError);
//				if(lvl == 4){
//					m.coords(*n.iPtr, Eigen::all) += *d_step;
//					m.writeMeshFile();
//					std::exit(0);
//				}
				auto stop = std::chrono::high_resolution_clock::now();
				auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);

//				w.setIntResults(i, lvl, maxError, abs(go.error).mean(), duration.count()/1e6, params.convHistFile, n.N_m);
//					std::cout << "LEVEL: " << lvl << " HAS BEEN DONE" << "\tTime: "<< duration.count()/1e6 << std::endl;

	//				std::cout << "mean error: " << abs(go.error).mean() << std::endl;

				lvl++;

				n.assignNodeTypesGrdy(m);

				if(maxError < params.tol){
					iterating = false;
					go.getAlphaVector();
				}
//					params.lvlSize = params.lvlSizeInit;

//					std::cout << "Maximum deformation: "<<  maxError << std::endl;

//					m.r = 2.5*maxError/0.0780524;
//					m.r = 100*maxError/0.707107;
//					std::cout<< m.r << std::endl;
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


			updateNodes(Phi_imGreedy,n,defVec, d_step, alpha_step,ctrlPtr);
			go.correction(m,n,params.gamma);
		}

	}
	std::cout << "Number of different control nodes: " << go.ctrlNodesAll.size() << std::endl;
	std::cout << "number of control nodes: " << n.N_m << std::endl;
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);

	std::cout << "Runtime duration: \t"<<  duration.count()/1e6 << " seconds"<< std::endl;
	m.writeMeshFile(params.mesh_ifName, params.mesh_ofName);
}


void rbf_std::performRBF(Eigen::MatrixXd& Phi_mm, Eigen::MatrixXd& Phi_im, Eigen::VectorXd& defVec, Eigen::ArrayXi& movingNodes, Eigen::ArrayXi& internalNodes,int& N){
	alpha.resize(N*m.nDims);

	if(params.dataRed){
		d.resize(internalNodes.size(),m.nDims);
	}

	for(int dim = 0; dim < m.nDims; dim++){
//		std::cout << "Solving for dimension: " << dim << std::endl;
		auto start = std::chrono::high_resolution_clock::now();
		alpha(Eigen::seqN(dim*N,N)) = Phi_mm.fullPivHouseholderQr().solve(defVec(Eigen::seqN(dim*N,N))); //fullPivHouseholderQr()

		auto stop = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);

//		std::cout << "Time to solve for alpha: \t"<<  duration.count()/1e6 << " seconds"<< std::endl;
		if(params.dataRed){
			d.col(dim) = Phi_im*alpha(Eigen::seqN(dim*N,N));
		}else{
			m.coords(internalNodes,dim) += (Phi_im*alpha(Eigen::seqN(dim*N,N))).array();
			m.coords(movingNodes,dim) += defVec(Eigen::seqN(dim*N,N)).array();
		}
	}

//	std::cout << "performRBF is done" << std::endl;
}

void rbf_std::updateNodes(Eigen::MatrixXd& Phi_imGreedy, getNodeType& n, Eigen::VectorXd& defVec, Eigen::ArrayXXd* d_step, Eigen::VectorXd* alpha_step, Eigen::ArrayXi* ctrlPtr){

	int N_m;
	Eigen::ArrayXi* ptr;

	if(params.multiLvl){
		ptr = ctrlPtr;
		N_m = (*ctrlPtr).size();

//		m.coords(*n.iPtrGrdy,Eigen::all) += deltaInternal;
	}else{
		if(params.smode == "none"){
			ptr = n.mPtr;
			N_m = n.N_m;
		}else{
			ptr = n.mStdPtr;
			N_m = n.N_mStd;
		}
	}

//	std::cout << *alpha_step << std::endl;
	auto start = std::chrono::high_resolution_clock::now();

	getPhi(Phi_imGreedy,n.iPtrGrdy,ptr);


	m.coords(*n.iPtr, Eigen::all) += *d_step;



	for(int dim = 0; dim < m.nDims; dim++){
		m.coords(*n.iPtrGrdy,dim) +=  (Phi_imGreedy*(*alpha_step)(Eigen::seqN(dim*N_m,N_m))).array();
	}

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);

	std::cout << "Time for updating internal nodes: \t"<<  duration.count()/1e6 << " seconds"<< std::endl;


//	auto stop = std::chrono::high_resolution_clock::now();
//	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
//	std::cout << "time: \t"<<  duration.count()/1e6 << " seconds"<< std::endl;

}


