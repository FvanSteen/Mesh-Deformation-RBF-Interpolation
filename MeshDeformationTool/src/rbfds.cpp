#include "rbfds.h"
#include <iostream>
#include <Eigen/Dense>
#include <chrono>

rbf_ds::rbf_ds(struct probParams& probParamsObject, Mesh& meshObject, getNodeType& n)
:rbfGenFunc(meshObject,probParamsObject)
{
	std::cout << "Initialised the ds class" << std::endl;
	perform_rbf(n);
}

void rbf_ds::perform_rbf(getNodeType& n){
	std::cout << "Performing RBF DS " << std::endl;

	projection p(pVec);


	auto start = std::chrono::high_resolution_clock::now();


	Eigen::MatrixXd Phi_cc, Phi_cs, Phi_sc, Phi_ss, Phi_ic, Phi_is, Phi,Phi_ibGrdy,Phi_mmStd,Phi_imStd;
	//for 3D
	Eigen::MatrixXd Phi_me, Phi_em, Phi_es, Phi_ee,Phi_se,Phi_ie;
	Eigen::VectorXd defVec, defVec_b;

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
			go.setInitMaxErrorNodes(m, m.coords, exactDisp, movingIndices, maxErrorNodes, params.doubleEdge);
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



			if (m.nDims == 2){
				getPhis(Phi_cc, Phi_cs, Phi_sc, Phi_ss, Phi_ic, Phi_is, n);
			}else if(m.nDims == 3){

//				getPhi(Phi_mm,n.cPtr,n.cPtr);
//				getPhi(Phi_me, n.cPtr, n.sePtr);
//				getPhi(Phi_ms, n.cPtr, n.ssPtr);
//
//				getPhi(Phi_em, n.sePtr, n.cPtr);
//				getPhi(Phi_ee, n.sePtr, n.sePtr);
//				getPhi(Phi_es, n.sePtr, n.ssPtr);
//
//				getPhi(Phi_sm, n.ssPtr, n.cPtr);
//				getPhi(Phi_se, n.ssPtr, n.sePtr);
//				getPhi(Phi_ss, n.ssPtr, n.ssPtr);
//
//				getPhi(Phi_im, n.iPtr, n.cPtr);
//				getPhi(Phi_ie, n.iPtr, n.sePtr);
//				getPhi(Phi_is, n.iPtr, n.ssPtr);
			}

			if(i==0 || params.dataRed){
				getDefVec(defVec, n.N_b, n.cPtr);
			}



		//		defVec = Eigen::VectorXd::Zero((N_m+N_s)*m.nDims);



//			getDefVec(defVec,n.N_m,params.steps,*n.mPtr);
//			getDefVec(defVec, n, lvl, go.errorPrevLvl);
			if(lvl!=0){
				getPhi(Phi_mmStd, n.bPtr,n.bPtr);
				getPhi(Phi_imStd, n.iPtr,n.bPtr);
				performRBF(Phi_mmStd, Phi_imStd, defVec, n.bPtr, n.iPtr, n.N_mStd);
			}else{
				if(m.nDims == 2){
					getPhiDS(Phi,Phi_cc,Phi_cs, Phi_sc, Phi_ss, m.n, m.t,n);

					// todo check which items can be omitted
					performRBF_DS(n, Phi, Phi_cc, Phi_cs, Phi_ic, Phi_is, Phi_sc, Phi_ss, defVec, defVec_b, p,*n.iPtr, *n.iPtr,*n.cPtr,*n.bPtr, *n.sPtr, n.N_i, n.N_m, n.N_mStd, n.N_s);
				}else if (m.nDims == 3){

//					getPhiDS_3D(Phi,Phi_mm, Phi_me, Phi_ms, Phi_em, Phi_ee, Phi_es, Phi_sm, Phi_se, Phi_ss, n);

//					performRBF_DS_3D( Phi,  Phi_im,  Phi_ie,  Phi_is,  Phi_em,  Phi_ee, Phi_es, Phi_sm, Phi_se, Phi_ss,defVec,  n);
				}


			}


//			std::exit(0);
			if(params.dataRed){

				go.getError(m,n, d,maxError,maxErrorNodes,movingIndices, exactDisp,pnVec,p, params.multiLvl, lvl, params.doubleEdge);
				std::cout << "error: \t"<< maxError <<" at node: \t" << maxErrorNodes(0) << std::endl;

				if(maxError < params.tol){
					iterating = false;
					maxErrorNodes.resize(0);
				}
			}else{
				iterating = false;
			}
//			if(iter == 2){
//				m.coords(*n.iPtr, Eigen::all) += (d);
//				m.writeMeshFile(params.mesh_ifName, params.mesh_ofName);
//				std::cout << *n.cPtr << "\n\n\n" << *n.bPtr << std::endl;
//				std::exit(0);
//			}

			if(params.multiLvl && (maxError/go.maxErrorPrevLvl < 0.1 || iterating == false)){


//				std::cout << *n.mStdPtr << std::endl;

				go.setLevelParams(m,n,lvl, d, alpha, maxError, defVec,n.cPtr, n.N_c);

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


			if(params.curved == false){
				setDefVec_b(defVec, defVec_b, n, Phi_sc, Phi_ss);
			}
			updateNodes(Phi_ibGrdy,n,defVec_b, d_step, alpha_step,ctrlPtr);
//			m.writeMeshFile(params.mesh_ifName, params.mesh_ofName);
//			std::exit(0);
			std::cout << "DOING AN UPDATE" << std::endl;
			go.correction( m,n,params.gamma, params.multiLvl);
		}

	}
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
	std::cout <<  "Runtime duration: \t"<<  duration.count()/1e6 << " seconds"<< std::endl;
}

void rbf_ds::setDefVec_b(Eigen::VectorXd& defVec, Eigen::VectorXd& defVec_b, getNodeType& n, Eigen::MatrixXd& Phi_sc, Eigen::MatrixXd& Phi_ss){
	Eigen::ArrayXXd ds(n.N_s,m.nDims);
	for(int dim = 0; dim < m.nDims; dim++){
		ds.col(dim) = (Phi_sc*alpha(Eigen::seqN(dim*(n.N_b),n.N_c)) + Phi_ss*alpha(Eigen::seqN(dim*(n.N_b)+n.N_c, n.N_s))).array();
	}
	getDefVec(defVec_b,defVec,n,ds);
}

void rbf_ds::performRBF_DS(getNodeType& n, Eigen::MatrixXd& Phi,Eigen::MatrixXd& Phi_cc, Eigen::MatrixXd& Phi_cs, Eigen::MatrixXd& Phi_ic, Eigen::MatrixXd& Phi_is, Eigen::MatrixXd& Phi_sc, Eigen::MatrixXd& Phi_ss, Eigen::VectorXd& defVec, Eigen::VectorXd& defVec_b, projection& p, Eigen::ArrayXi& iNodes,Eigen::ArrayXi& iNodesGrdy, Eigen::ArrayXi& mNodes, Eigen::ArrayXi& mNodesStd, Eigen::ArrayXi& sNodes, int& N_i, int& N_m, int& N_mStd, int& N_s){

	alpha = Phi.fullPivHouseholderQr().solve(defVec);

	if(params.dataRed){
		d.resize(n.N_i,m.nDims);
	}


	if(params.curved){
		Eigen::ArrayXXd delta(n.N_s, m.nDims), finalDef(n.N_s,m.nDims);

		for (int dim = 0; dim < m.nDims; dim++){
			delta.col(dim) = (Phi_sc*alpha(Eigen::seqN(dim*(n.N_b),n.N_c)) + Phi_ss*alpha(Eigen::seqN(dim*(n.N_b)+n.N_c, n.N_s))).array();
		}

		p.projectIter(m,*n.sPtr,delta,finalDef,n.N_s);

		getDefVec(defVec_b, defVec, n, finalDef);


		Eigen::MatrixXd Phi_bb(n.N_b,n.N_b), Phi_ib(n.N_i, n.N_b);
		Phi_bb.block(0,0, n.N_c, n.N_c) = Phi_cc;
		Phi_bb.block(0,n.N_c, n.N_c, n.N_s) = Phi_cs;
		Phi_bb.block(n.N_c, n.N_c, n.N_s, n.N_s) = Phi_ss;
		Phi_bb.block(n.N_c, 0, n.N_s, n.N_c) = Phi_sc;

		Phi_ib.block(0,0,n.N_i, n.N_c) = Phi_ic;
		Phi_ib.block(0,n.N_c,n.N_i, n.N_s) = Phi_is;

		performRBF(Phi_bb,Phi_ib,defVec_b,n.bPtr, n.iPtr, n.N_b);

	}
	else{
		for (int dim = 0; dim < m.nDims; dim++){

			if(params.dataRed){
				d.col(dim) = Phi_ic*alpha(Eigen::seqN(dim*(n.N_b),n.N_c)) + Phi_is*alpha(Eigen::seqN(dim*(n.N_b)+n.N_c, n.N_s));
			}
			else{
//				std::cout << "solving for dimensions: " << dim << std::endl;
				m.coords(*n.iPtr, dim) += (Phi_ic*alpha(Eigen::seqN(dim*(n.N_b),n.N_c)) + Phi_is*alpha(Eigen::seqN(dim*(n.N_b)+n.N_c, n.N_s))).array();

				m.coords(*n.sPtr, dim) += (Phi_sc*alpha(Eigen::seqN(dim*(n.N_b),n.N_c)) + Phi_ss*alpha(Eigen::seqN(dim*(n.N_b)+n.N_c, n.N_s))).array();

				m.coords(*n.cPtr, dim) += (defVec(Eigen::seqN(dim*n.N_c,n.N_c))).array();

			}
		}

	}
}



void rbf_ds::getPhiDS(Eigen::MatrixXd& Phi,Eigen::MatrixXd& Phi_mm,Eigen::MatrixXd& Phi_ms, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_ss, Eigen::ArrayXXd& nVec, Eigen::ArrayXXd& tVec, getNodeType& n){

	Phi = Eigen::MatrixXd::Zero(m.nDims*(n.N_b),m.nDims*(n.N_b));

	if(params.pmode == "moving"){
		nVec.conservativeResize(m.N_se,m.nDims);
		tVec.conservativeResize(m.N_se,m.nDims);
		for(int i = m.N_se - m.periodicVerticesNodes.size(); i< n.N_s; i++){
			nVec.row(i) = pnVec;
			tVec.row(i) = pVec;
		}
	}

	Eigen::ArrayXi indices;
	int type = 0;
	getIdxSlidingNodes(n.sPtr,indices, type);


	for(int dim = 0; dim< m.nDims; dim++){
		// blocks related to the known displacements
		Phi.block(dim*n.N_c, dim*(n.N_b), n.N_c, n.N_c) = Phi_mm;
		Phi.block(dim*n.N_c, dim*(n.N_b)+n.N_c, n.N_c, n.N_s) = Phi_ms;


		Phi.block(2*n.N_c, dim*(n.N_b), n.N_s, n.N_c) = Phi_sm.array().colwise() * nVec(indices, dim);
		Phi.block(2*n.N_c, dim*(n.N_b)+n.N_c, n.N_s, n.N_s) = Phi_ss.array().colwise() * nVec(indices, dim);


		//blocks related to the zero tangential contribution condition
		Eigen::VectorXd diag(n.N_s);
		diag = tVec(indices,dim);

//		Phi.block(2*N_m + N_s, dim*(N_m+N_s)+N_m, N_s, N_s) = Eigen::MatrixXd(t.col(dim).matrix().asDiagonal());
		Phi.block(2*n.N_c + n.N_s, dim*(n.N_b)+n.N_c, n.N_s, n.N_s) = diag.asDiagonal();

	}
//	std::cout << "PHI \n\n" << Phi << std::endl;
//	std::exit(0);

}

void rbf_ds::getIdxSlidingNodes(Eigen::ArrayXi* sPtr, Eigen::ArrayXi& idx,int type){
	//todo make some if statement in case it is not in the datareduction mode

	Eigen::ArrayXi* ptr;
	if(type == 0){
		ptr = &m.seNodes;
	}else if(type == 1){
		ptr = &m.ssNodes;
	}

	idx.resize((*sPtr).size());
	for(int i = 0; i < (*sPtr).size(); i++){
		idx(i) = std::distance(std::begin(*ptr),std::find(std::begin(*ptr),std::end(*ptr),(*sPtr)(i)));
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
		getIdxSlidingNodes(n.sePtr, seIndices, type);

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
		getIdxSlidingNodes(n.ssPtr, ssIndices, type);

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
			m.coords(*n.cPtr, dim) += (defVec(Eigen::seqN(dim*n.N_m,n.N_m))).array();
		}

	}

}
