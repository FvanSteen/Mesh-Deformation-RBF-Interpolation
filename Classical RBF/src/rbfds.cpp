#include "rbfds.h"
#include <iostream>
#include <Eigen/Dense>
#include <chrono>

rbf_ds::rbf_ds(Mesh& meshObject, struct probParams& probParamsObject)
:rbf_std(meshObject, probParamsObject)
{
	std::cout << "Initialised the ds class" << std::endl;
}

void rbf_ds::perform_rbf(getNodeType& n){
	std::cout << "Performing RBF DS " << std::endl;
	projection* p;

	projection proObject;
	p = &proObject;


	auto start = std::chrono::high_resolution_clock::now();


	Eigen::MatrixXd Phi_mm, Phi_ms, Phi_sm, Phi_ss, Phi_im, Phi_is, Phi,Phi_imGrdy;
	Eigen::VectorXd defVec;

//	int maxErrorNode;
	Eigen::ArrayXi maxErrorNodes;
	if(params.dataRed){
//		maxErrorNode = ;
		n.addControlNode(m.intBdryNodes(0));
	}
	std::cout << "'internal' nodes: \n"<< *n.iPtr << "\n moving nodes: \n" << *n.mPtr << "\n sliding nodes: \n" << *n.sePtr << "\nmSTD nodes: \n" << *n.mStdPtr << std::endl;
	greedy go;
	int iter;
	double error;
	int lvl = 0;
	for (int i = 0; i < params.steps; i++){
		std::cout << "Deformation step: " << i+1 << " out of "<< params.steps << std::endl;
		error = 1;
		iter = 0;
		while(error > params.tol){

			if(iter!=0){
				for(int node = 0; node < maxErrorNodes.size(); node++){
					n.addControlNode(maxErrorNodes(node));
				}
			}
			std::cout << "Moving nodes: \n " << *n.mPtr << std::endl;
			std::cout << "sliding nodes: \n " << *n.sePtr << std::endl;

			getPhi(Phi_mm, *n.mPtr,*n.mPtr);
			getPhi(Phi_ms, *n.mPtr, *n.sePtr);
			getPhi(Phi_sm, *n.sePtr, *n.mPtr);
			getPhi(Phi_ss, *n.sePtr, *n.sePtr);

			getPhi(Phi_im, *n.iPtr, *n.mPtr);
			getPhi(Phi_is, *n.iPtr, *n.sePtr);


			if(params.curved || i==0){
		//		if(i==0){
				// getVecs obtains average vector at the nodes

				m.getVecs();

				// getMidPnts obtains vectors at midpoint of boundary segments


				m.getMidPnts();
			}

		//		defVec = Eigen::VectorXd::Zero((N_m+N_s)*m.nDims);



//			getDefVec(defVec,n.N_m,params.steps,*n.mPtr);
			getDefVec(defVec, n, lvl, go.errorPrevLvl);

			getPhiDS(Phi,Phi_mm,Phi_ms, Phi_sm, Phi_ss, m.n, m.t,n.N_m,n.N_se, *n.sePtr);


			// todo check which items can be omitted
			performRBF_DS(n, Phi, Phi_im, Phi_is, Phi_sm, Phi_ss, defVec, p,*n.iPtr, *n.iPtrGrdy,*n.mPtr,*n.mStdPtr, *n.sePtr, n.N_i, n.N_m, n.N_mStd, n.N_se);

//			std::exit(0);
			if(params.dataRed){

				go.getError(m,n, d,error,maxErrorNodes,movingIndices, exactDisp,pnVec,p, params.multiLvl, lvl);
//				std::cout << "error: \t"<< error <<" at node: \t" << maxErrorNode<< std::endl;

			}else{
				error = 0;
			}
//			if(iter == 4){
//				m.coords(*n.iPtr, Eigen::all) += d;
//				m.writeMeshFile();
//				std::exit(0);
//			}
			iter++;


		}

		if(params.dataRed){
//			std::cout << *n.mStdPtr << std::endl;
			updateNodes(Phi_imGrdy, n, defVec,go.delta, go.deltaInternal);
//			std::cout << m.coords << std::endl;
			std::cout << "DOING AN UPDATE" << std::endl;
			go.correction( m,n,params.gamma);
		}
//		m.coords(*n.iPtr, Eigen::all) += d;
//		m.writeMeshFile();
//		std::exit(0);
	}
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
	std::cout <<  "Runtime duration: \t"<<  duration.count()/1e6 << " seconds"<< std::endl;
	m.writeMeshFile();
}

void rbf_ds::performRBF_DS(getNodeType& n, Eigen::MatrixXd& Phi, Eigen::MatrixXd& Phi_im, Eigen::MatrixXd& Phi_is, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_ss, Eigen::VectorXd& defVec, projection* proPnt, Eigen::ArrayXi& iNodes,Eigen::ArrayXi& iNodesGrdy, Eigen::ArrayXi& mNodes, Eigen::ArrayXi& mNodesStd, Eigen::ArrayXi& sNodes, int& N_i, int& N_m, int& N_mStd, int& N_s){

	alpha = Phi.fullPivLu().solve(defVec);


	if(params.dataRed){
		d.resize(iNodes.size(),m.nDims);
	}


	if(params.curved){
		Eigen::ArrayXXd delta(N_s, m.nDims), finalDef(N_s,m.nDims);

		// find displacement
		for (int dim = 0; dim < m.nDims; dim++){
			delta.col(dim) = (Phi_sm*alpha(Eigen::seqN(dim*(N_m+N_s),N_m)) + Phi_ss*alpha(Eigen::seqN(dim*(N_m+N_s)+N_m, N_s))).array();

		}
//		m.coords(sNodes,Eigen::seq(0,1)) += delta;
//		std::cout << "delta; \n" << delta << std::endl;


		// calling project function to find the final deformation after the projection

		proPnt->project(m,sNodes,delta,finalDef,pVec);
//		std::cout << finalDef << std::endl;



//		m.coords(sNodes,Eigen::seq(0,1)) += finalDef;
//		m.writeMeshFile();
//

//		std::exit(0);
//		defVec = Eigen::VectorXd::Zero(N_m2*m.nDims);
		Eigen::VectorXd defVecStd;
		defVecStd = Eigen::VectorXd::Zero(N_mStd*m.nDims);
//		getDefVec(defVec,N_mStd,params.steps);
//		std::cout << defVec << std::endl;

		for(int dim = 0; dim< m.nDims; dim++){
			defVecStd(Eigen::seqN(dim*n.N_mStd,n.N_m)) = defVec(Eigen::seqN(dim*n.N_m,n.N_m));
			defVecStd(Eigen::seqN(dim*n.N_mStd+n.N_m,n.N_se)) = finalDef.col(dim);
//			std::cout << "\n" << defVec<< '\n'<< std::endl;
	//		defVec(Eigen::seqN(dim*N_m+N_mPro,m.N_se)) = finalDef.col(dim);

		}
//		for(int dim = 0; dim< m.nDims; dim++){
//
////			std::cout << "\n" << defVec<< '\n'<< std::endl;
//		}

//		std::cout << defVecStd << std::endl;




//		std::cout << finalDef << std::endl;
//		std::cout << defVec << std::endl;
//		std::exit(0);

		Eigen::MatrixXd Phi_mm2, Phi_im2;

		//todo give n as argument in the function :(
		getPhi(Phi_mm2, *n.mStdPtr,*n.mStdPtr);
		getPhi(Phi_im2,*n.iPtr,*n.mStdPtr);

		performRBF(Phi_mm2,Phi_im2,defVecStd,mNodesStd,iNodes,N_mStd);
//		std::cout << d << std::endl;
//		std::cout << "std rbf has been performed " << std::endl;
//		if(N_s == 1){
//			std::exit(0);
//		}
	}
	else{
		for (int dim = 0; dim < m.nDims; dim++){
			if(params.dataRed){
				d.col(dim) = Phi_im*alpha(Eigen::seqN(dim*(N_m+N_s),N_m)) + Phi_is*alpha(Eigen::seqN(dim*(N_m+N_s)+N_m, N_s));
			}
			else{
				m.coords(iNodes, dim) += (Phi_im*alpha(Eigen::seqN(dim*(N_m+N_s),N_m)) + Phi_is*alpha(Eigen::seqN(dim*(N_m+N_s)+N_m, N_s))).array();

				m.coords(sNodes, dim) += (Phi_sm*alpha(Eigen::seqN(dim*(N_m+N_s),N_m)) + Phi_ss*alpha(Eigen::seqN(dim*(N_m+N_s)+N_m, N_s))).array();

				m.coords(mNodes, dim) += (defVec(Eigen::seqN(dim*N_m,N_m))).array();

			}
		}

	}
//	std::cout << "DISPLACEMENT \n\n" << d << std::endl;
}



void rbf_ds::getPhiDS(Eigen::MatrixXd& Phi,Eigen::MatrixXd& Phi_mm,Eigen::MatrixXd& Phi_ms, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_ss, Eigen::ArrayXXd& n, Eigen::ArrayXXd& t,int& N_m, int& N_s, Eigen::ArrayXi& sNodes){

	Phi = Eigen::MatrixXd::Zero(m.nDims*(N_m+N_s),m.nDims*(N_m+N_s));


	if(m.pmode == "moving"){
		n.conservativeResize(N_s,m.nDims);
		t.conservativeResize(N_s,m.nDims);
		for(int i = N_s - m.staticNodes.size(); i< N_s; i++){
			n.row(i) = pnVec;
			t.row(i) = pVec;
		}
	}
	Eigen::ArrayXi indices;
	getIdxSlidingNodes(sNodes,indices);
//	std::cout << "indices \n" << indices << std::endl;
	for(int dim = 0; dim< m.nDims; dim++){
		// blocks related to the known displacements
		Phi.block(dim*N_m, dim*(N_m+N_s), N_m, N_m) = Phi_mm;
		Phi.block(dim*N_m, dim*(N_m+N_s)+N_m, N_m, N_s) = Phi_ms;


		Phi.block(2*N_m, dim*(N_m+N_s), N_s, N_m) = Phi_sm.array().colwise() * n(indices, dim);
		Phi.block(2*N_m, dim*(N_m+N_s)+N_m, N_s, N_s) = Phi_ss.array().colwise() * n(indices, dim);


		//blocks related to the zero tangential contribution condition
		Eigen::VectorXd diag(sNodes.size());
		diag = t(indices,dim);

//		Phi.block(2*N_m + N_s, dim*(N_m+N_s)+N_m, N_s, N_s) = Eigen::MatrixXd(t.col(dim).matrix().asDiagonal());
		Phi.block(2*N_m + N_s, dim*(N_m+N_s)+N_m, N_s, N_s) = diag.asDiagonal();

	}
//	std::cout << "PHI \n\n" << Phi << std::endl;
//	std::exit(0);

}

void rbf_ds::getIdxSlidingNodes(Eigen::ArrayXi& sNodes, Eigen::ArrayXi& idx){
	//todo make some if statement in case it is not in the datareduction mode

	idx.resize(sNodes.size());
	for(int i = 0; i < sNodes.size(); i++){
		idx(i) = std::distance(std::begin(m.seNodes),std::find(std::begin(m.seNodes),std::end(m.seNodes),sNodes(i)));
	}

}


