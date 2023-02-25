#include "rbfds.h"
#include <iostream>
#include <Eigen/Dense>
#include <chrono>

rbf_ds::rbf_ds(struct probParams& probParamsObject, Mesh& meshObject, getNodeType& n)
:rbf_std(probParamsObject, meshObject, n)
{
	std::cout << "Initialised the ds class" << std::endl;
}

void rbf_ds::perform_rbf(getNodeType& n){
	std::cout << "Performing RBF DS " << std::endl;
	projection* p;

	projection proObject;
	p = &proObject;


	auto start = std::chrono::high_resolution_clock::now();


	Eigen::MatrixXd Phi_mm, Phi_ms, Phi_sm, Phi_ss, Phi_im, Phi_is, Phi,Phi_imGrdy,Phi_mmStd,Phi_imStd;
	//for 3D
	Eigen::MatrixXd Phi_me, Phi_em, Phi_es, Phi_ee,Phi_se,Phi_ie;
	Eigen::VectorXd defVec;

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

	int iter, lvl;
	double maxError;
	bool iterating;

	for (int i = 0; i < params.steps; i++){
		std::cout << "Deformation step: " << i+1 << " out of "<< params.steps << std::endl;
		maxError = 1;
		iter = 0;
		lvl = 0;

		if((params.dataRed && i==0) || params.multiLvl ){
			go.setInitMaxErrorNodes(m, m.coords, exactDisp, movingIndices, maxErrorNodes);
		}

		if(params.curved || i==0){
			m.getVecs();
			m.getMidPnts(params);
		}

		iterating = true;

		while(iterating){

			if(params.dataRed){
				for(int node = 0; node < maxErrorNodes.size(); node++){
					n.addControlNode(maxErrorNodes(node), params.smode, m);
				}
			}
//			std::cout << "Moving nodes: \n " << *n.mPtr << std::endl;
//			std::cout << "sliding nodes: \n " << *n.sPtr << std::endl;
			if (m.nDims == 2){
				getPhi(Phi_mm, n.mPtr,n.mPtr);
				getPhi(Phi_ms, n.mPtr, n.sPtr);
				getPhi(Phi_sm, n.sPtr, n.mPtr);
				getPhi(Phi_ss, n.sPtr, n.sPtr);

				getPhi(Phi_im, n.iPtr, n.mPtr);
				getPhi(Phi_is, n.iPtr, n.sPtr);
			}else if(m.nDims == 3){
//				Eigen::MatrixXd Phi_mm, Phi_ms, Phi_sm, Phi_ss, Phi_im, Phi_is, Phi;
//			//	Eigen::ArrayXXd n(m.N_se, m.nDims), t(m.N_se, m.nDims);		// two column array containing normal vector components
//				Eigen::VectorXd defVec, alpha(m.nDims*(N_m+m.N_se+m.N_ss));
//			//	Eigen::ArrayXXd t_se(m.N_se,m.nDims),n1_se(m.N_se,m.nDims), n2_se(m.N_se,m.nDims), n_ss(m.N_ss, m.nDims),t1_ss(m.N_ss, m.nDims),t2_ss(m.N_ss, m.nDims);
//				Eigen::MatrixXd Phi_me, Phi_ee, Phi_es, Phi_em, Phi_se, Phi_ie;

				getPhi(Phi_mm,n.mPtr,n.mPtr);
				getPhi(Phi_me, n.mPtr, n.sePtr);
				getPhi(Phi_ms, n.mPtr, n.ssPtr);

				getPhi(Phi_em, n.sePtr, n.mPtr);
				getPhi(Phi_ee, n.sePtr, n.sePtr);
				getPhi(Phi_es, n.sePtr, n.ssPtr);

				getPhi(Phi_sm, n.ssPtr, n.mPtr);
				getPhi(Phi_se, n.ssPtr, n.sePtr);
				getPhi(Phi_ss, n.ssPtr, n.ssPtr);

				getPhi(Phi_im, n.iPtr, n.mPtr);
				getPhi(Phi_ie, n.iPtr, n.sePtr);
				getPhi(Phi_is, n.iPtr, n.ssPtr);
			}

			if(i==0 || params.dataRed){
				getDefVec(defVec, n, lvl, go.errorPrevLvl);
			}




		//		defVec = Eigen::VectorXd::Zero((N_m+N_s)*m.nDims);



//			getDefVec(defVec,n.N_m,params.steps,*n.mPtr);
//			getDefVec(defVec, n, lvl, go.errorPrevLvl);
			if(lvl!=0){
				getPhi(Phi_mmStd, n.mStdPtr,n.mStdPtr);
				getPhi(Phi_imStd, n.iPtr,n.mStdPtr);
				performRBF(Phi_mmStd, Phi_imStd, defVec, *n.mStdPtr, *n.iPtr, n.N_mStd);
			}else{
				if(m.nDims == 2){
					getPhiDS(Phi,Phi_mm,Phi_ms, Phi_sm, Phi_ss, m.n, m.t,n.N_m,n.N_s, *n.sPtr);
					// todo check which items can be omitted
					performRBF_DS(n, Phi, Phi_im, Phi_is, Phi_sm, Phi_ss, defVec, p,*n.iPtr, *n.iPtrGrdy,*n.mPtr,*n.mStdPtr, *n.sPtr, n.N_i, n.N_m, n.N_mStd, n.N_s);
				}else if (m.nDims == 3){

					getPhiDS_3D(Phi,Phi_mm, Phi_me, Phi_ms, Phi_em, Phi_ee, Phi_es, Phi_sm, Phi_se, Phi_ss, n);

					performRBF_DS_3D( Phi,  Phi_im,  Phi_ie,  Phi_is,  Phi_em,  Phi_ee, Phi_es, Phi_sm, Phi_se, Phi_ss,defVec,  n);
				}


			}


//			std::exit(0);
			if(params.dataRed){

				go.getError(m,n, d,maxError,maxErrorNodes,movingIndices, exactDisp,pnVec,p, params.multiLvl, lvl);
				std::cout << "error: \t"<< maxError <<" at node: \t" << maxErrorNodes(0) << std::endl;

				if(maxError < params.tol){
					iterating = false;
				}
			}else{
				iterating = false;
			}
//			if(iter == 4){
//				m.coords(*n.iPtr, Eigen::all) += d;
//				m.writeMeshFile();
//				std::exit(0);
//			}

			if(params.multiLvl && (maxError/go.maxErrorPrevLvl < 0.1 || iterating == false)){


//				std::cout << *n.mStdPtr << std::endl;

				go.setLevelParams(m,n,lvl,params.lvlSize, d, alpha, maxError);

//				std::cout << "LEVEL: " << lvl << " HAS BEEN DONE" << std::endl;
				lvl++;

				n.assignNodeTypesGrdy(m);

				if(maxError < params.tol){
					iterating = false;
					go.getAlphaVector();
				}

			}

			iter++;


		}

		if(params.dataRed){
//			updateNodes(Phi_imGrdy, n, defVec,go.delta, go.deltaInternal, go.alphaSum);
			updateNodes(Phi_imGrdy,n,defVec, d_step, alpha_step,ctrlPtr);
			std::cout << "DOING AN UPDATE" << std::endl;
			go.correction( m,n,params.gamma);
		}

	}
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
	std::cout <<  "Runtime duration: \t"<<  duration.count()/1e6 << " seconds"<< std::endl;
	m.writeMeshFile(params.mesh_ifName, params.mesh_ofName);
}

void rbf_ds::performRBF_DS(getNodeType& n, Eigen::MatrixXd& Phi, Eigen::MatrixXd& Phi_im, Eigen::MatrixXd& Phi_is, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_ss, Eigen::VectorXd& defVec, projection* proPnt, Eigen::ArrayXi& iNodes,Eigen::ArrayXi& iNodesGrdy, Eigen::ArrayXi& mNodes, Eigen::ArrayXi& mNodesStd, Eigen::ArrayXi& sNodes, int& N_i, int& N_m, int& N_mStd, int& N_s){

	alpha = Phi.fullPivHouseholderQr().solve(defVec);


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
			defVecStd(Eigen::seqN(dim*n.N_mStd+n.N_m,n.N_s)) = finalDef.col(dim);
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
		getPhi(Phi_mm2, n.mStdPtr,n.mStdPtr);
		getPhi(Phi_im2,n.iPtr,n.mStdPtr);

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
//				std::cout << "solving for dimensions: " << dim << std::endl;
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


	if(params.pmode == "moving"){
		n.conservativeResize(N_s,m.nDims);
		t.conservativeResize(N_s,m.nDims);
		for(int i = N_s - m.verticesNodes.size(); i< N_s; i++){
			n.row(i) = pnVec;
			t.row(i) = pVec;
		}
	}
	Eigen::ArrayXi indices;
	int type = 0;
	getIdxSlidingNodes(sNodes,indices, type);


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

void rbf_ds::getIdxSlidingNodes(Eigen::ArrayXi& sNodes, Eigen::ArrayXi& idx,int type){
	//todo make some if statement in case it is not in the datareduction mode

	Eigen::ArrayXi* ptr;
	if(type == 0){
		ptr = &m.seNodes;
	}else if(type == 1){
		ptr = &m.ssNodes;
	}

	idx.resize(sNodes.size());
	for(int i = 0; i < sNodes.size(); i++){
		idx(i) = std::distance(std::begin(*ptr),std::find(std::begin(*ptr),std::end(*ptr),sNodes(i)));
	}

}


void rbf_ds::getPhiDS_3D(Eigen::MatrixXd& Phi, Eigen::MatrixXd& Phi_mm, Eigen::MatrixXd&  Phi_me, Eigen::MatrixXd&  Phi_ms, Eigen::MatrixXd& Phi_em, Eigen::MatrixXd& Phi_ee, Eigen::MatrixXd& Phi_es, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_se, Eigen::MatrixXd& Phi_ss, getNodeType& n){

	Phi = Eigen::MatrixXd::Zero(m.nDims*(n.N_mStd),m.nDims*(n.N_mStd));

	for(int dim = 0; dim< m.nDims; dim++){
		// blocks related to the known displacements
		Phi.block(dim*n.N_m, dim*(n.N_mStd), n.N_m, n.N_m) = Phi_mm;
		Phi.block(dim*n.N_m, dim*(n.N_mStd)+n.N_m, n.N_m, n.N_se) = Phi_me;
		Phi.block(dim*n.N_m, dim*(n.N_mStd)+n.N_m+ n.N_se, n.N_m, n.N_ss) = Phi_ms;

		int type = 0;
		Eigen::ArrayXi seIndices;
		getIdxSlidingNodes(*n.sePtr, seIndices, type);

		//blocks realteð to the first zero normal displacement condition of the sliding edge nodes
		Phi.block(3*n.N_m, dim*(n.N_mStd),					n.N_se, n.N_m ) = Phi_em.array().colwise() * m.n1_se(seIndices, dim);  //m.n1_se.col(dim);
		Phi.block(3*n.N_m, dim*(n.N_mStd) + n.N_m,			n.N_se, n.N_se ) = Phi_ee.array().colwise() * m.n1_se(seIndices,dim);
		Phi.block(3*n.N_m, dim*(n.N_mStd) + n.N_m + n.N_se,	n.N_se,	n.N_ss ) = Phi_es.array().colwise() * m.n1_se(seIndices, dim);

		//blocks realteð to the second zero normal displacement condition of the sliding edge nodes
		Phi.block(3*n.N_m+n.N_se, dim*(n.N_mStd),					n.N_se, n.N_m ) = Phi_em.array().colwise() * m.n2_se(seIndices, dim);
		Phi.block(3*n.N_m+n.N_se, dim*(n.N_mStd) + n.N_m,			n.N_se, n.N_se ) = Phi_ee.array().colwise() * m.n2_se(seIndices, dim);
		Phi.block(3*n.N_m+n.N_se, dim*(n.N_mStd) + n.N_m + n.N_se,	n.N_se,	n.N_ss ) = Phi_es.array().colwise() * m.n2_se(seIndices, dim);

		// blocks related to the zero tangential contribution of the sliding edge nodes
		Eigen::ArrayXd diag;
		diag = m.t_se(seIndices,dim);

		Phi.block(3*n.N_m+2*n.N_se, dim*(n.N_mStd) + n.N_m,	n.N_se, n.N_se) = Eigen::MatrixXd(diag.matrix().asDiagonal());


		type = 1;
		Eigen::ArrayXi ssIndices;
		getIdxSlidingNodes(*n.ssPtr, ssIndices, type);

		// blocks related to the zero normal displacement condition of the sliding surface nodes
		Phi.block(3*n.N_m+3*n.N_se, dim*(n.N_mStd),		n.N_ss,n.N_m) = Phi_sm.array().colwise() * m.n_ss(ssIndices,dim);
		Phi.block(3*n.N_m+3*n.N_se, dim*(n.N_mStd)+n.N_m,	n.N_ss,n.N_se) = Phi_se.array().colwise() * m.n_ss(ssIndices,dim);
		Phi.block(3*n.N_m+3*n.N_se, dim*(n.N_mStd)+n.N_m+n.N_se,	n.N_ss,n.N_ss) = Phi_ss.array().colwise() * m.n_ss(ssIndices, dim);

		Eigen::ArrayXd diag2;
		diag2 = m.t1_ss(ssIndices,dim);
		// blocks related to the zero tangential contribution of the sliding surface nodes.
		Phi.block(3*n.N_m+3*n.N_se+n.N_ss, dim*(n.N_mStd)+n.N_m+n.N_se, n.N_ss,n.N_ss) = Eigen::MatrixXd(diag2.matrix().asDiagonal());
		diag2 = m.t2_ss(ssIndices,dim);
		Phi.block(3*n.N_m+3*n.N_se+2*n.N_ss, dim*(n.N_mStd)+n.N_m+n.N_se, n.N_ss,n.N_ss) = Eigen::MatrixXd(diag2.matrix().asDiagonal());

	}


}

void rbf_ds::performRBF_DS_3D(Eigen::MatrixXd& Phi, Eigen::MatrixXd& Phi_im, Eigen::MatrixXd& Phi_ie, Eigen::MatrixXd& Phi_is, Eigen::MatrixXd& Phi_em, Eigen::MatrixXd& Phi_ee, Eigen::MatrixXd& Phi_es, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_se, Eigen::MatrixXd& Phi_ss, Eigen::VectorXd& defVec, getNodeType& n){
	std::cout << "deze" << std::endl;
	alpha = Phi.partialPivLu().solve(defVec);
	std::cout << "SOLUTION WAS FOUND" << std::endl;

	if(params.dataRed){
		d.resize(n.N_i,m.nDims);
	}

	for (int dim = 0; dim < m.nDims; dim++){
		std::cout << "dim: " << dim << std::endl;
		if(params.dataRed){
			d.col(dim) = (Phi_im*alpha(Eigen::seqN(dim*(n.N_mStd),n.N_m)) + Phi_ie*alpha(Eigen::seqN(dim*(n.N_mStd)+n.N_m, n.N_se)) + Phi_is*alpha(Eigen::seqN(dim*(n.N_mStd)+n.N_m+n.N_se, n.N_ss)) ).array();
		}else{

			m.coords(*n.iPtr, dim) += (Phi_im*alpha(Eigen::seqN(dim*(n.N_mStd),n.N_m)) + Phi_ie*alpha(Eigen::seqN(dim*(n.N_mStd)+n.N_m, n.N_se)) + Phi_is*alpha(Eigen::seqN(dim*(n.N_mStd)+n.N_m+n.N_se, n.N_ss)) ).array();
			m.coords(*n.sePtr, dim) += (Phi_em*alpha(Eigen::seqN(dim*(n.N_mStd),n.N_m)) + Phi_ee*alpha(Eigen::seqN(dim*(n.N_mStd)+n.N_m, n.N_se)) + Phi_es*alpha(Eigen::seqN(dim*(n.N_mStd)+n.N_m+n.N_se, n.N_ss)) ).array();
			m.coords(*n.ssPtr, dim) += (Phi_sm*alpha(Eigen::seqN(dim*(n.N_mStd),n.N_m)) + Phi_se*alpha(Eigen::seqN(dim*(n.N_mStd)+n.N_m, n.N_se)) + Phi_ss*alpha(Eigen::seqN(dim*(n.N_mStd)+n.N_m+n.N_se, n.N_ss)) ).array();
			m.coords(*n.mPtr, dim) += (defVec(Eigen::seqN(dim*n.N_m,n.N_m))).array();
		}

	}

}
