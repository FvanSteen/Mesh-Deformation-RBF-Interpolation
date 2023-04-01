#include "rbfds.h"
#include <iostream>
#include <Eigen/Dense>
#include <chrono>


rbf_ds::rbf_ds(struct probParams& probParamsObject, Mesh& meshObject, getNodeType& n)
:rbfGenFunc(meshObject,probParamsObject)
{
	std::cout << "Initialised the ds class" << std::endl;
	if(params.dataRed){
		greedy g(m, params, exactDisp, movingIndices, alpha, d);
		perform_rbf(n,g);
	}else{
		perform_rbf(n);
	}
}

void rbf_ds::perform_rbf(getNodeType& n){
	std::cout << "Performing RBF DS without data reduction" << std::endl;
	std::clock_t s = std::clock();

	Eigen::MatrixXd Phi;
	Eigen::VectorXd defVec, defVec_b;

	for (int i = 0; i < params.steps; i++){
		std::cout << "Deformation step: " << i+1 << " out of "<< params.steps << std::endl;

		if(params.curved || i==0){
			m.getVecs(params);
			m.getMidPnts(params);
		}

		if(i==0){
			getDefVec(defVec, n.N_c, n.mPtr);
		}

		getPhis(n, 0);

		getPhiDS(Phi, PhiPtr ,n);

		performRBF_DS(n, Phi, PhiPtr, defVec, defVec_b);
	}

	std::clock_t e = std::clock();
	long double time_elapsed_ms =  1000.0*(e-s) / CLOCKS_PER_SEC;
	std::cout << "CPU time: " << time_elapsed_ms/1000 << " ms\n";

}

void rbf_ds::perform_rbf(getNodeType& n, greedy& g){
	std::cout << "Performing RBF DS " << std::endl;


	std::clock_t s = std::clock();
	Eigen::MatrixXd Phi;
	Eigen::VectorXd defVec, defVec_b;



	int iter, lvl;
	bool iterating = true;

	for (int i = 0; i < params.steps; i++){
		std::cout << "Deformation step: " << i+1 << " out of "<< params.steps << std::endl;

		iter = 0;
		lvl = 0;

		if(params.curved || i==0){
			m.getVecs(params);
			m.getMidPnts(params);
		}

		while(iterating){
			n.addControlNodes(g.maxErrorNodes, params.smode, m);

			// obtaining the interpolation matrices
			getPhis(n, iter);

			// obtaining the deformation vector

			getDefVec(defVec, n.N_c, n.mPtr);

			if(lvl!=0){
				getPhi(PhiPtr->Phi_cc, n.cPtr,n.cPtr);
				getPhi(PhiPtr->Phi_ic, n.iPtr,n.cPtr);
				performRBF(PhiPtr->Phi_cc, PhiPtr->Phi_ic, defVec, n.cPtr, n.iPtr, n.N_c);
			}else{
				getPhiDS(Phi, PhiPtr ,n);
				performRBF_DS(n, Phi, PhiPtr, defVec, defVec_b);
			}



			g.getError(n, d, lvl);
			std::cout << "error: \t"<< g.maxError <<" at node: \t" << g.maxErrorNodes(0)<< std::endl;

			if(g.maxError < params.tol){
				iterating = false;
				if(params.multiLvl == false){
					g.maxErrorNodes.resize(0);
				}
			}

			if(params.multiLvl && (g.maxError/g.maxErrorPrevLvl < params.tolCrit || iterating == false)){
				// todo direct sliding with multi level
				std::cout << "LEVEL: " << lvl << " HAS BEEN DONE" << std::endl;
				g.setLevelParams(n,lvl, d, alpha, defVec_b,n.cPtr, n.N_c);

				lvl++;
				iter = -1;

				n.assignNodeTypesGrdy(m);

				if(iterating == false){
					g.getAlphaVector();
					g.setInitMaxErrorNodes();
				}
			}

			iter++;
		}


		if(params.curved == false){
			setDefVec_b(defVec, defVec_b, n, PhiPtr);
		}

		updateNodes(n, defVec_b, g.d_step, g.alpha_step, g.ctrlPtr);
		g.correction( m,n,params.gamma, params.multiLvl);
		iterating = true;
	}

	std::clock_t e = std::clock();
	long double time_elapsed_ms =  1000.0*(e-s) / CLOCKS_PER_SEC;
	std::cout << "CPU time: " << time_elapsed_ms/1000 << " ms\n";
}

void rbf_ds::setDefVec_b(Eigen::VectorXd& defVec, Eigen::VectorXd& defVec_b, getNodeType& n, PhiStruct* PhiPtr ){
	Eigen::ArrayXXd ds(n.N_se+n.N_ss,m.nDims);
	if(m.nDims == 2){
		for(int dim = 0; dim < m.nDims; dim++){
			ds.col(dim) = (PhiPtr->Phi_em*alpha(Eigen::seqN(dim*(n.N_c),n.N_m)) + PhiPtr->Phi_ee*alpha(Eigen::seqN(dim*(n.N_c)+n.N_m, n.N_se))).array();
		}
	} else if(m.nDims == 3){
		for(int dim = 0; dim < m.nDims; dim++){
			ds(Eigen::seqN(0, n.N_se),dim) = (PhiPtr->Phi_em*alpha(Eigen::seqN(dim*(n.N_c),n.N_m)) + PhiPtr->Phi_ee*alpha(Eigen::seqN(dim*(n.N_c)+n.N_m, n.N_se)) + PhiPtr->Phi_es*alpha(Eigen::seqN(dim*(n.N_c)+n.N_m+n.N_se, n.N_ss))).array();
			ds(Eigen::seqN(n.N_se, n.N_ss),dim) = (PhiPtr->Phi_sm*alpha(Eigen::seqN(dim*(n.N_c),n.N_m)) + PhiPtr->Phi_se*alpha(Eigen::seqN(dim*(n.N_c)+n.N_m, n.N_se)) + PhiPtr->Phi_ss*alpha(Eigen::seqN(dim*(n.N_c)+n.N_m+n.N_se, n.N_ss))).array();
		}
	}
	getDefVec(defVec_b,defVec,n,ds);
}

void rbf_ds::performRBF_DS(getNodeType& n, Eigen::MatrixXd& Phi,  PhiStruct* PhiPtr , Eigen::VectorXd& defVec, Eigen::VectorXd& defVec_b){
	alpha = Phi.colPivHouseholderQr().solve(defVec);

	if(params.dataRed){
		d.resize(n.N_i,m.nDims);
	}

	if(m.nDims == 2){
		if(params.curved){
			Eigen::ArrayXXd delta(n.N_se, m.nDims), finalDef(n.N_se,m.nDims);

			for (int dim = 0; dim < m.nDims; dim++){
				delta.col(dim) = (PhiPtr->Phi_em*alpha(Eigen::seqN(dim*(n.N_c),n.N_m)) + PhiPtr->Phi_ee*alpha(Eigen::seqN(dim*(n.N_c)+n.N_m, n.N_se))).array();
			}


			p.project(m,n,delta, finalDef, m.periodicVec);
			getDefVec(defVec_b, defVec, n, finalDef);

			performRBF(PhiPtr->Phi_cc,PhiPtr->Phi_ic,defVec_b,n.cPtr, n.iPtr, n.N_c);

		}
		else{
			for (int dim = 0; dim < m.nDims; dim++){

				if(params.dataRed){

					d.col(dim) = PhiPtr->Phi_ic*alpha(Eigen::seqN(dim*n.N_c,n.N_c));
				}
				else{
					m.coords(*n.iPtr, dim) += (PhiPtr->Phi_ic*alpha(Eigen::seqN(dim*n.N_c, n.N_c))).array();

					m.coords(*n.sePtr, dim) += (PhiPtr->Phi_em*alpha(Eigen::seqN(dim*(n.N_c),n.N_m)) + PhiPtr->Phi_ee*alpha(Eigen::seqN(dim*(n.N_c)+n.N_m, n.N_se))).array();

					m.coords(*n.mPtr, dim) += (defVec(Eigen::seqN(dim*n.N_m,n.N_m))).array();

				}
			}
		}
	}else if(m.nDims == 3){
		if(params.curved){
			Eigen::ArrayXXd delta(n.N_se+n.N_ss, m.nDims), finalDef(n.N_se+n.N_ss,m.nDims);


			for (int dim = 0; dim < m.nDims; dim++){
				delta(Eigen::seqN(0,n.N_se),dim) = (PhiPtr->Phi_em*alpha(Eigen::seqN(dim*(n.N_c),n.N_m)) + PhiPtr->Phi_ee*alpha(Eigen::seqN(dim*(n.N_c)+n.N_m, n.N_se)) + PhiPtr->Phi_es*alpha(Eigen::seqN(dim*(n.N_c)+n.N_m+n.N_se, n.N_ss))).array();
				delta(Eigen::seqN(n.N_se,n.N_ss),dim) = (PhiPtr->Phi_sm*alpha(Eigen::seqN(dim*(n.N_c),n.N_m)) + PhiPtr->Phi_se*alpha(Eigen::seqN(dim*(n.N_c)+n.N_m, n.N_se)) + PhiPtr->Phi_ss*alpha(Eigen::seqN(dim*(n.N_c)+n.N_m+n.N_se, n.N_ss))).array();
			}



			p.project(m, n, delta, finalDef, m.periodicVec);

			getDefVec(defVec_b, defVec, n, finalDef);

			performRBF(PhiPtr->Phi_cc,PhiPtr->Phi_ic,defVec_b,n.cPtr, n.iPtr, n.N_c);

		}else{
			for (int dim = 0; dim < m.nDims; dim++){
				if(params.dataRed){
					d.col(dim) = PhiPtr->Phi_ic*alpha(Eigen::seqN(dim*n.N_c,n.N_c));

				}else{

					m.coords(*n.iPtr, dim) += (PhiPtr->Phi_ic*alpha(Eigen::seqN(dim*n.N_c, n.N_c))).array();
					m.coords(*n.sePtr, dim) += (PhiPtr->Phi_em*alpha(Eigen::seqN(dim*(n.N_c),n.N_m)) + PhiPtr->Phi_ee*alpha(Eigen::seqN(dim*(n.N_c)+n.N_m, n.N_se)) + PhiPtr->Phi_es*alpha(Eigen::seqN(dim*(n.N_c)+n.N_m+n.N_se, n.N_ss)) ).array();
					m.coords(*n.ssPtr, dim) += (PhiPtr->Phi_sm*alpha(Eigen::seqN(dim*(n.N_c),n.N_m)) + PhiPtr->Phi_se*alpha(Eigen::seqN(dim*(n.N_c)+n.N_m, n.N_se)) + PhiPtr->Phi_ss*alpha(Eigen::seqN(dim*(n.N_c)+n.N_m+n.N_se, n.N_ss)) ).array();
					m.coords(*n.mPtr, dim) += (defVec(Eigen::seqN(dim*n.N_m,n.N_m))).array();
				}

			}
		}
	}
}



void rbf_ds::getPhiDS(Eigen::MatrixXd& Phi, PhiStruct* PhiPtr, getNodeType& n){
	//todo perfect the next parts.
	Phi = Eigen::MatrixXd::Zero(m.nDims*(n.N_c),m.nDims*(n.N_c));
	if(m.nDims == 2){

		Eigen::ArrayXi indices;

		getIdxSlidingNodes(n.sePtr,indices, m.seNodes);

		for(int dim = 0; dim< m.nDims; dim++){
			// blocks related to the known displacements
			Phi.block(dim*n.N_m, dim*(n.N_c), n.N_m, n.N_m) = PhiPtr->Phi_mm;
			Phi.block(dim*n.N_m, dim*(n.N_c)+n.N_m, n.N_m, n.N_se) = PhiPtr->Phi_me;

			Phi.block(2*n.N_m, dim*(n.N_c), n.N_se, n.N_m) = PhiPtr->Phi_em.array().colwise() * m.n(indices, dim);
			Phi.block(2*n.N_m, dim*(n.N_c)+n.N_m, n.N_se, n.N_se) = PhiPtr->Phi_ee.array().colwise() * m.n(indices, dim);

			//blocks related to the zero tangential contribution condition
			Eigen::VectorXd diag(n.N_se);
			diag = m.t(indices,dim);

			Phi.block(2*n.N_m + n.N_se, dim*(n.N_c)+n.N_m, n.N_se, n.N_se) = diag.asDiagonal();

		}

	}else if(m.nDims ==3){

		for(int dim = 0; dim< m.nDims; dim++){
			// blocks related to the known displacements
			Phi.block(dim*n.N_m, dim*(n.N_c), n.N_m, n.N_m) = PhiPtr->Phi_mm;
			Phi.block(dim*n.N_m, dim*(n.N_c)+n.N_m, n.N_m, n.N_se) = PhiPtr->Phi_me;
			Phi.block(dim*n.N_m, dim*(n.N_c)+n.N_m+ n.N_se, n.N_m, n.N_ss) = PhiPtr->Phi_ms;


			Eigen::ArrayXi seIndices;
			getIdxSlidingNodes(n.sePtr, seIndices, m.seNodes);

			//blocks realted to the first zero normal displacement condition of the sliding edge nodes
			Phi.block(3*n.N_m, dim*(n.N_c),					n.N_se, n.N_m ) = PhiPtr->Phi_em.array().colwise() * m.n1_se(seIndices, dim);  //m.n1_se.col(dim);
			Phi.block(3*n.N_m, dim*(n.N_c) + n.N_m,			n.N_se, n.N_se ) = PhiPtr->Phi_ee.array().colwise() * m.n1_se(seIndices,dim);
			Phi.block(3*n.N_m, dim*(n.N_c) + n.N_m + n.N_se,	n.N_se,	n.N_ss ) = PhiPtr->Phi_es.array().colwise() * m.n1_se(seIndices, dim);

			//blocks realteð to the second zero normal displacement condition of the sliding edge nodes
			Phi.block(3*n.N_m+n.N_se, dim*(n.N_c),					n.N_se, n.N_m ) = PhiPtr->Phi_em.array().colwise() * m.n2_se(seIndices, dim);
			Phi.block(3*n.N_m+n.N_se, dim*(n.N_c) + n.N_m,			n.N_se, n.N_se ) = PhiPtr->Phi_ee.array().colwise() * m.n2_se(seIndices, dim);
			Phi.block(3*n.N_m+n.N_se, dim*(n.N_c) + n.N_m + n.N_se,	n.N_se,	n.N_ss ) = PhiPtr->Phi_es.array().colwise() * m.n2_se(seIndices, dim);

			// blocks related to the zero tangential contribution of the sliding edge nodes
			Eigen::ArrayXd diag;
			diag = m.t_se(seIndices,dim);

			Phi.block(3*n.N_m+2*n.N_se, dim*(n.N_c) + n.N_m,	n.N_se, n.N_se) = Eigen::MatrixXd(diag.matrix().asDiagonal());

			Eigen::ArrayXi ssIndices;
			getIdxSlidingNodes(n.ssPtr, ssIndices, m.ssNodes);

			// blocks related to the zero normal displacement condition of the sliding surface nodes
			Phi.block(3*n.N_m+3*n.N_se, dim*(n.N_c),		n.N_ss,n.N_m) = PhiPtr->Phi_sm.array().colwise() * m.n_ss(ssIndices,dim);
			Phi.block(3*n.N_m+3*n.N_se, dim*(n.N_c)+n.N_m,	n.N_ss,n.N_se) = PhiPtr->Phi_se.array().colwise() * m.n_ss(ssIndices,dim);
			Phi.block(3*n.N_m+3*n.N_se, dim*(n.N_c)+n.N_m+n.N_se,	n.N_ss,n.N_ss) = PhiPtr->Phi_ss.array().colwise() * m.n_ss(ssIndices, dim);

			Eigen::ArrayXd diag2;
			diag2 = m.t1_ss(ssIndices,dim);
			// blocks related to the zero tangential contribution of the sliding surface nodes.
			Phi.block(3*n.N_m+3*n.N_se+n.N_ss, dim*(n.N_c)+n.N_m+n.N_se, n.N_ss,n.N_ss) = Eigen::MatrixXd(diag2.matrix().asDiagonal());
			diag2 = m.t2_ss(ssIndices,dim);
			Phi.block(3*n.N_m+3*n.N_se+2*n.N_ss, dim*(n.N_c)+n.N_m+n.N_se, n.N_ss,n.N_ss) = Eigen::MatrixXd(diag2.matrix().asDiagonal());

		}
	}

}

void rbf_ds::getIdxSlidingNodes(Eigen::ArrayXi* sPtr, Eigen::ArrayXi& idx, Eigen::ArrayXi& sNodesInit){
	idx.resize((*sPtr).size());
	if(params.dataRed){
		for(int i = 0; i < (*sPtr).size(); i++){
			idx(i) = std::distance(std::begin(sNodesInit),std::find(std::begin(sNodesInit),std::end(sNodesInit),(*sPtr)(i)));
		}
	}else{
		idx = Eigen::ArrayXi::LinSpaced(sNodesInit.size(), 0, sNodesInit.size()-1);
	}
}


