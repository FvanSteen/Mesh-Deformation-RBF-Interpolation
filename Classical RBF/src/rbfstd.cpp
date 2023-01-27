#include "rbfstd.h"
#include <Eigen/Dense>
#include <chrono>
#include <iostream>
#include "greedy.h"
//rbf_std::rbf_std(Mesh& meshObject, Eigen::VectorXd& dVec, Eigen::RowVectorXd& rotPnt, Eigen::VectorXd& rotVec, const int& steps, const std::string& smode, const bool& curved, const std::string& pDir)
rbf_std::rbf_std(Mesh& meshObject, struct probParams& probParamsObject)
:rbfGenFunc(meshObject, probParamsObject)
{
}


void rbf_std::perform_rbf(getNodeType& n){
	std::cout << "Performing standard rbf interpolation" << std::endl;
	projection* p;

	projection proObject;
	p = &proObject;

	auto start = std::chrono::high_resolution_clock::now();
	Eigen::MatrixXd Phi_mm, Phi_im, Phi_imGreedy;
	Eigen::VectorXd defVec,defVecStd;

	double maxError;

	int iter, lvl;
	bool iterating;

	Eigen::ArrayXi maxErrorNodes;


	if(params.dataRed){
		n.addControlNode(m.intBdryNodes(0));
		n.addControlNode(m.intBdryNodes(m.intBdryNodes.size()-1));
	}

	greedy go;

	for(int i=0; i<params.steps; i++){
		std::cout << "Deformation step: " << i+1 << " out of "<< params.steps << std::endl;
		maxError = 1;
		iter = 0;
		lvl = 0;
		iterating = true;
		if(params.multiLvl && i != 0){
			n.addControlNode(m.intBdryNodes(0));
			n.addControlNode(m.intBdryNodes(m.intBdryNodes.size()-1));
		}
		while(iterating){

			if(iter!=0){
				for(int node = 0; node < maxErrorNodes.size(); node++){
					n.addControlNode(maxErrorNodes(node));
				}
			}

//			std::cout << "Obtaining Phi_mm" << std::endl;
			getPhi(Phi_mm, *n.mPtr, *n.mPtr);

//			std::cout << "Obtaining Phi_im" << std::endl;
			getPhi(Phi_im, *n.iPtr, *n.mPtr);

			if(i==0 || params.dataRed){
				getDefVec(defVec, n, lvl, go.errorPrevLvl);
			}



			performRBF(Phi_mm, Phi_im, defVec,*n.mPtr,*n.iPtr, n.N_m);

//			std::cout << d << std::endl;
			if(params.dataRed){
				go.getError(m,n,d,maxError,maxErrorNodes, movingIndices, exactDisp,pnVec,p, params.multiLvl, lvl);
				std::cout << "error: \t"<< maxError <<" at node: \t" << maxErrorNodes(0)<< std::endl;


				if(params.multiLvl == false && maxError < params.tol){

					iterating = false;
				}

			}else{

//				maxError = 0;
				iterating = false;
			}

			if(params.multiLvl && n.N_m >= params.lvlSize){

				getPhi(Phi_imGreedy, *n.iPtrGrdy, *n.mPtr);
				go.setLevelParams(m,n,lvl,params.lvlSize, d, alpha, Phi_imGreedy);

				std::cout << "LEVEL: " << lvl << " HAS BEEN DONE" << std::endl;
				lvl++;

				n.assignNodeTypesGreedy();

				if(maxError < params.tol){
					iterating = false;
				}

			}

			iter++;

		}


		if(params.dataRed){
			updateNodes(Phi_imGreedy,n, defVec,go.delta, go.deltaInternal);
			go.correction(m,n,params.gamma);
		}

	}
//	std::cout << "number of control nodes: " << n.N_m << std::endl;
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);

	std::cout << "Runtime duration: \t"<<  duration.count()/1e6 << " seconds"<< std::endl;
	m.writeMeshFile();
}


void rbf_std::performRBF(Eigen::MatrixXd& Phi_mm, Eigen::MatrixXd& Phi_im, Eigen::VectorXd& defVec, Eigen::ArrayXi& movingNodes, Eigen::ArrayXi& internalNodes,int& N){
	alpha.resize(N*m.nDims);

	if(params.dataRed){
		d.resize(internalNodes.size(),m.nDims);
	}

	for(int dim = 0; dim < m.nDims; dim++){
		std::cout << "Solving for dimension: " << dim << std::endl;
		alpha(Eigen::seqN(dim*N,N)) = Phi_mm.fullPivHouseholderQr().solve(defVec(Eigen::seqN(dim*N,N)));
		if(params.dataRed){
			d.col(dim) = Phi_im*alpha(Eigen::seqN(dim*N,N));

		}else{
			m.coords(internalNodes,dim) += (Phi_im*alpha(Eigen::seqN(dim*N,N))).array();
			m.coords(movingNodes,dim) += defVec(Eigen::seqN(dim*N,N)).array();

		}
	}
}

void rbf_std::updateNodes(Eigen::MatrixXd& Phi_imGreedy, getNodeType& n, Eigen::VectorXd& defVec, Eigen::ArrayXXd& delta, Eigen::ArrayXXd& deltaInternal){

	if(params.multiLvl){
		m.coords(*n.iPtr, Eigen::all) += delta;
		m.coords(*n.iPtrGrdy,Eigen::all) += deltaInternal;
	}else{
		int N_m;
		Eigen::ArrayXi* ptr;

		if(m.smode == "none"){
			ptr = n.mPtr;
			N_m = n.N_m;
		}else{
			ptr = n.mStdPtr;
			N_m = n.N_mStd;
		}

		getPhi(Phi_imGreedy,*n.iPtrGrdy,*ptr);

		m.coords(*n.iPtr,Eigen::all) += d;

		for(int dim = 0; dim < m.nDims; dim++){
			m.coords(*n.iPtrGrdy,dim) +=  (Phi_imGreedy*alpha(Eigen::seqN(dim*N_m,N_m))).array();
		}
	}
}


