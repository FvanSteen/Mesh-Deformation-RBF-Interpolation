#include "rbfds.h"
#include <iostream>
#include <Eigen/Dense>
#include <chrono>
#include "SPDS.h"

rbf_ds::rbf_ds(struct probParams& probParamsObject, Mesh& meshObject, getNodeType& n)
:rbfGenFunc(meshObject,probParamsObject)
{
	std::cout << "Initialised the ds class" << std::endl;
	perform_rbf(n);
}

void rbf_ds::perform_rbf(getNodeType& n){
	std::cout << "Performing RBF DS " << std::endl;

//	auto start = std::chrono::high_resolution_clock::now();
	std::clock_t s = std::clock();
	Eigen::MatrixXd Phi;
	Eigen::VectorXd defVec, defVec_b;

	Eigen::ArrayXi maxErrorNodes;
	greedy go(m, params, exactDisp, movingIndices, alpha, d, periodicVec);


	int iter, lvl;
	double maxError;
	bool iterating;

	for (int i = 0; i < params.steps; i++){
		std::cout << "Deformation step: " << i+1 << " out of "<< params.steps << std::endl;
		maxError = 1;
		iter = 0;
		lvl = 0;

		if((params.dataRed && i==0) || params.multiLvl ){
			go.setInitMaxErrorNodes(m);
		}

		if(params.curved || i==0){
			m.getVecs();
			m.getMidPnts(params);
		}


		if(params.pmode == "moving"){
			if(m.nDims == 2){
				m.n.conservativeResize(m.N_se,m.nDims);
				m.t.conservativeResize(m.N_se,m.nDims);
				for(int i = m.N_se - m.periodicVerticesNodes.size(); i< m.N_se; i++){
					m.n.row(i) = periodicNormalVec1;
					m.t.row(i) = periodicVec;
				}
			}else if(m.nDims == 3){
				m.n1_se.conservativeResize(m.N_se,m.nDims);
				m.n2_se.conservativeResize(m.N_se,m.nDims);
				m.t_se.conservativeResize(m.N_se, m.nDims);
				for(int i = m.N_se - m.periodicVerticesNodes.size(); i< m.N_se; i++){
					m.n1_se.row(i) = periodicNormalVec1;
					m.n2_se.row(i) = periodicNormalVec2;
					m.t_se.row(i) = periodicVec;
				}
			}

		}
		iterating = true;

		while(iterating){

			if(params.dataRed){
				for(int node = 0; node < maxErrorNodes.size(); node++){
					n.addControlNode(maxErrorNodes(node), params.smode, m);
				}
			}


			// obtaining the interpolation matrices
			getPhis(n);

			// obtaining the deformation vector
			if(i==0 || params.dataRed){
				getDefVec(defVec, n.N_b, n.cPtr);
			}

			if(lvl!=0){
				getPhi(PhiPtr->Phi_bb, n.bPtr,n.bPtr);
				getPhi(PhiPtr->Phi_ib, n.iPtr,n.bPtr);
				performRBF(PhiPtr->Phi_bb, PhiPtr->Phi_ib, defVec, n.bPtr, n.iPtr, n.N_b);
			}else{
				getPhiDS(Phi, PhiPtr ,n);

				performRBF_DS(n, Phi, PhiPtr, defVec, defVec_b);

			}


			if(params.dataRed){
				go.getError(m,n, d, lvl);
				std::cout << "error: \t"<< maxError <<" at node: \t" << maxErrorNodes(0)<< std::endl;
				if(maxError < params.tol){
					iterating = false;
					maxErrorNodes.resize(0);
				}
			}else{
				iterating = false;
			}


			if(params.multiLvl && (maxError/go.maxErrorPrevLvl < 0.1 || iterating == false)){
				go.setLevelParams(m,n,lvl, d, alpha, defVec,n.cPtr, n.N_c);
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
				setDefVec_b(defVec, defVec_b, n, PhiPtr);
			}

			updateNodes(n, defVec_b, go.d_step, go.alpha_step, go.ctrlPtr);
			std::cout << "DOING AN UPDATE" << std::endl;
			go.correction( m,n,params.gamma, params.multiLvl);
		}
	}
//	auto stop = std::chrono::high_resolution_clock::now();
//	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
//	std::cout <<  "Runtime duration: \t"<<  duration.count()/1e6 << " seconds"<< std::endl;
	std::clock_t e = std::clock();
//	auto stop = std::chrono::high_resolution_clock::now();
//	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
	long double time_elapsed_ms =  1000.0*(e-s) / CLOCKS_PER_SEC;
	std::cout << "CPU time: " << time_elapsed_ms/1000 << " ms\n";
}

void rbf_ds::setDefVec_b(Eigen::VectorXd& defVec, Eigen::VectorXd& defVec_b, getNodeType& n, PhiStruct* PhiPtr ){
	Eigen::ArrayXXd ds(n.N_s,m.nDims);
	if(m.nDims == 2){
		for(int dim = 0; dim < m.nDims; dim++){
			ds.col(dim) = (PhiPtr->Phi_sc*alpha(Eigen::seqN(dim*(n.N_b),n.N_c)) + PhiPtr->Phi_ss*alpha(Eigen::seqN(dim*(n.N_b)+n.N_c, n.N_s))).array();
		}
	} else if(m.nDims == 3){
		for(int dim = 0; dim < m.nDims; dim++){
			ds(Eigen::seqN(0, n.N_se),dim) = (PhiPtr->Phi_ec*alpha(Eigen::seqN(dim*(n.N_b),n.N_c)) + PhiPtr->Phi_ee*alpha(Eigen::seqN(dim*(n.N_b)+n.N_c, n.N_se)) + PhiPtr->Phi_es*alpha(Eigen::seqN(dim*(n.N_b)+n.N_c+n.N_se, n.N_ss))).array();
			ds(Eigen::seqN(n.N_se, n.N_ss),dim) = (PhiPtr->Phi_sc*alpha(Eigen::seqN(dim*(n.N_b),n.N_c)) + PhiPtr->Phi_se*alpha(Eigen::seqN(dim*(n.N_b)+n.N_c, n.N_se)) + PhiPtr->Phi_ss*alpha(Eigen::seqN(dim*(n.N_b)+n.N_c+n.N_se, n.N_ss))).array();
		}
	}
	getDefVec(defVec_b,defVec,n,ds);
}

void rbf_ds::performRBF_DS(getNodeType& n, Eigen::MatrixXd& Phi,  PhiStruct* PhiPtr , Eigen::VectorXd& defVec, Eigen::VectorXd& defVec_b){

	alpha = Phi.householderQr().solve(defVec);

	if(params.dataRed){
		d.resize(n.N_i,m.nDims);
	}

	if(m.nDims == 2){
		if(params.curved){
			Eigen::ArrayXXd delta(n.N_s, m.nDims), finalDef(n.N_s,m.nDims);

			for (int dim = 0; dim < m.nDims; dim++){
				delta.col(dim) = (PhiPtr->Phi_sc*alpha(Eigen::seqN(dim*(n.N_b),n.N_c)) + PhiPtr->Phi_ss*alpha(Eigen::seqN(dim*(n.N_b)+n.N_c, n.N_s))).array();
			}

//			p.projectIter(m,*n.sPtr,delta,finalDef,n.N_s);
			SPDS SPDSobj;
			SPDSobj.project(m,n,delta, finalDef, periodicVec);
			getDefVec(defVec_b, defVec, n, finalDef);

			performRBF(PhiPtr->Phi_bb,PhiPtr->Phi_ib,defVec_b,n.bPtr, n.iPtr, n.N_b);

		}
		else{
			for (int dim = 0; dim < m.nDims; dim++){

				if(params.dataRed){
					d.col(dim) = PhiPtr->Phi_ic*alpha(Eigen::seqN(dim*(n.N_b),n.N_c)) + PhiPtr->Phi_is*alpha(Eigen::seqN(dim*(n.N_b)+n.N_c, n.N_s));
				}
				else{
	//				std::cout << "solving for dimensions: " << dim << std::endl;
					m.coords(*n.iPtr, dim) += (PhiPtr->Phi_ic*alpha(Eigen::seqN(dim*(n.N_b),n.N_c)) + PhiPtr->Phi_is*alpha(Eigen::seqN(dim*(n.N_b)+n.N_c, n.N_s))).array();

					m.coords(*n.sPtr, dim) += (PhiPtr->Phi_sc*alpha(Eigen::seqN(dim*(n.N_b),n.N_c)) + PhiPtr->Phi_ss*alpha(Eigen::seqN(dim*(n.N_b)+n.N_c, n.N_s))).array();

					m.coords(*n.cPtr, dim) += (defVec(Eigen::seqN(dim*n.N_c,n.N_c))).array();

				}
			}
		}
	}else if(m.nDims == 3){
		if(params.curved){
			Eigen::ArrayXXd delta(n.N_s, m.nDims), finalDef(n.N_s,m.nDims);


			for (int dim = 0; dim < m.nDims; dim++){
				delta(Eigen::seqN(0,n.N_se),dim) = (PhiPtr->Phi_ec*alpha(Eigen::seqN(dim*(n.N_b),n.N_c)) + PhiPtr->Phi_ee*alpha(Eigen::seqN(dim*(n.N_b)+n.N_c, n.N_se)) + PhiPtr->Phi_es*alpha(Eigen::seqN(dim*(n.N_b)+n.N_c+n.N_se, n.N_ss))).array();
				delta(Eigen::seqN(n.N_se,n.N_ss),dim) = (PhiPtr->Phi_sc*alpha(Eigen::seqN(dim*(n.N_b),n.N_c)) + PhiPtr->Phi_se*alpha(Eigen::seqN(dim*(n.N_b)+n.N_c, n.N_se)) + PhiPtr->Phi_ss*alpha(Eigen::seqN(dim*(n.N_b)+n.N_c+n.N_se, n.N_ss))).array();
			}


			SPDS SPDSobj;
			SPDSobj.project(m, n, delta, finalDef, periodicVec);

			getDefVec(defVec_b, defVec, n, finalDef);


			Eigen::MatrixXd Phi_bb(n.N_b,n.N_b), Phi_ib(n.N_i, n.N_b);
			Phi_bb.block(0,0, n.N_c, n.N_c) = PhiPtr->Phi_cc;
			Phi_bb.block(0,n.N_c, n.N_c, n.N_se) = PhiPtr->Phi_ce;
			Phi_bb.block(0,n.N_c+n.N_se, n.N_c, n.N_ss) = PhiPtr->Phi_cs;

			Phi_bb.block(n.N_c, 0, n.N_se, n.N_c) = PhiPtr->Phi_ec;
			Phi_bb.block(n.N_c, n.N_c, n.N_se, n.N_se) = PhiPtr->Phi_ee;
			Phi_bb.block(n.N_c, n.N_c+n.N_se, n.N_se, n.N_ss) = PhiPtr->Phi_es;

			Phi_bb.block(n.N_c+n.N_se , 0, n.N_ss, n.N_c) = PhiPtr->Phi_sc;
			Phi_bb.block(n.N_c+n.N_se, n.N_c, n.N_ss, n.N_se) = PhiPtr->Phi_se;
			Phi_bb.block(n.N_c+n.N_se, n.N_c+n.N_se, n.N_ss, n.N_ss) = PhiPtr->Phi_ss;

			Phi_ib.block(0,0,n.N_i, n.N_c) = PhiPtr->Phi_ic;
			Phi_ib.block(0,n.N_c,n.N_i, n.N_se) = PhiPtr->Phi_ie;
			Phi_ib.block(0,n.N_c+n.N_se,n.N_i, n.N_ss) = PhiPtr->Phi_is;

			performRBF(Phi_bb,Phi_ib,defVec_b,n.bPtr, n.iPtr, n.N_b);



		}else{
			for (int dim = 0; dim < m.nDims; dim++){
				if(params.dataRed){
					d.col(dim) = (PhiPtr->Phi_ic*alpha(Eigen::seqN(dim*(n.N_b),n.N_c)) + PhiPtr->Phi_ie*alpha(Eigen::seqN(dim*(n.N_b)+n.N_c, n.N_se)) + PhiPtr->Phi_is*alpha(Eigen::seqN(dim*(n.N_b)+n.N_c+n.N_se, n.N_ss)) ).array();
				}else{

					m.coords(*n.iPtr, dim) += (PhiPtr->Phi_ic*alpha(Eigen::seqN(dim*(n.N_b),n.N_c)) + PhiPtr->Phi_ie*alpha(Eigen::seqN(dim*(n.N_b)+n.N_c, n.N_se)) + PhiPtr->Phi_is*alpha(Eigen::seqN(dim*(n.N_b)+n.N_c+n.N_se, n.N_ss)) ).array();
					m.coords(*n.sePtr, dim) += (PhiPtr->Phi_ec*alpha(Eigen::seqN(dim*(n.N_b),n.N_c)) + PhiPtr->Phi_ee*alpha(Eigen::seqN(dim*(n.N_b)+n.N_c, n.N_se)) + PhiPtr->Phi_es*alpha(Eigen::seqN(dim*(n.N_b)+n.N_c+n.N_se, n.N_ss)) ).array();
					m.coords(*n.ssPtr, dim) += (PhiPtr->Phi_sc*alpha(Eigen::seqN(dim*(n.N_b),n.N_c)) + PhiPtr->Phi_se*alpha(Eigen::seqN(dim*(n.N_b)+n.N_c, n.N_se)) + PhiPtr->Phi_ss*alpha(Eigen::seqN(dim*(n.N_b)+n.N_c+n.N_se, n.N_ss)) ).array();
					m.coords(*n.cPtr, dim) += (defVec(Eigen::seqN(dim*n.N_c,n.N_c))).array();
				}

			}
		}
	}
}



void rbf_ds::getPhiDS(Eigen::MatrixXd& Phi, PhiStruct* PhiPtr, getNodeType& n){
	//todo perfect the next parts. N_s should perhaps be distinct in N_se and N_ss in 2D to make the function the same
	Phi = Eigen::MatrixXd::Zero(m.nDims*(n.N_b),m.nDims*(n.N_b));
	if(m.nDims == 2){

		Eigen::ArrayXi indices;
		int type = 0;
		getIdxSlidingNodes(n.sPtr,indices, type);


		for(int dim = 0; dim< m.nDims; dim++){
			// blocks related to the known displacements
			Phi.block(dim*n.N_c, dim*(n.N_b), n.N_c, n.N_c) = PhiPtr->Phi_cc;
			Phi.block(dim*n.N_c, dim*(n.N_b)+n.N_c, n.N_c, n.N_s) = PhiPtr->Phi_cs;


			Phi.block(2*n.N_c, dim*(n.N_b), n.N_s, n.N_c) = PhiPtr->Phi_sc.array().colwise() * m.n(indices, dim);
			Phi.block(2*n.N_c, dim*(n.N_b)+n.N_c, n.N_s, n.N_s) = PhiPtr->Phi_ss.array().colwise() * m.n(indices, dim);


			//blocks related to the zero tangential contribution condition
			Eigen::VectorXd diag(n.N_s);
			diag = m.t(indices,dim);

	//		Phi.block(2*N_m + N_s, dim*(N_m+N_s)+N_m, N_s, N_s) = Eigen::MatrixXd(t.col(dim).matrix().asDiagonal());
			Phi.block(2*n.N_c + n.N_s, dim*(n.N_b)+n.N_c, n.N_s, n.N_s) = diag.asDiagonal();

		}

	}else if(m.nDims ==3){

		for(int dim = 0; dim< m.nDims; dim++){
			// blocks related to the known displacements
			Phi.block(dim*n.N_c, dim*(n.N_b), n.N_c, n.N_c) = PhiPtr->Phi_cc;
			Phi.block(dim*n.N_c, dim*(n.N_b)+n.N_c, n.N_c, n.N_se) = PhiPtr->Phi_ce;
			Phi.block(dim*n.N_c, dim*(n.N_b)+n.N_c+ n.N_se, n.N_c, n.N_ss) = PhiPtr->Phi_cs;

			int type = 0;
			Eigen::ArrayXi seIndices;
			getIdxSlidingNodes(n.sePtr, seIndices, type);

			//blocks realted to the first zero normal displacement condition of the sliding edge nodes
			Phi.block(3*n.N_c, dim*(n.N_b),					n.N_se, n.N_c ) = PhiPtr->Phi_ec.array().colwise() * m.n1_se(seIndices, dim);  //m.n1_se.col(dim);
			Phi.block(3*n.N_c, dim*(n.N_b) + n.N_c,			n.N_se, n.N_se ) = PhiPtr->Phi_ee.array().colwise() * m.n1_se(seIndices,dim);
			Phi.block(3*n.N_c, dim*(n.N_b) + n.N_c + n.N_se,	n.N_se,	n.N_ss ) = PhiPtr->Phi_es.array().colwise() * m.n1_se(seIndices, dim);

			//blocks realteð to the second zero normal displacement condition of the sliding edge nodes
			Phi.block(3*n.N_c+n.N_se, dim*(n.N_b),					n.N_se, n.N_c ) = PhiPtr->Phi_ec.array().colwise() * m.n2_se(seIndices, dim);
			Phi.block(3*n.N_c+n.N_se, dim*(n.N_b) + n.N_c,			n.N_se, n.N_se ) = PhiPtr->Phi_ee.array().colwise() * m.n2_se(seIndices, dim);
			Phi.block(3*n.N_c+n.N_se, dim*(n.N_b) + n.N_c + n.N_se,	n.N_se,	n.N_ss ) = PhiPtr->Phi_es.array().colwise() * m.n2_se(seIndices, dim);

			// blocks related to the zero tangential contribution of the sliding edge nodes
			Eigen::ArrayXd diag;
			diag = m.t_se(seIndices,dim);

			Phi.block(3*n.N_c+2*n.N_se, dim*(n.N_b) + n.N_c,	n.N_se, n.N_se) = Eigen::MatrixXd(diag.matrix().asDiagonal());


			type = 1;
			Eigen::ArrayXi ssIndices;
			getIdxSlidingNodes(n.ssPtr, ssIndices, type);

			// blocks related to the zero normal displacement condition of the sliding surface nodes
			Phi.block(3*n.N_c+3*n.N_se, dim*(n.N_b),		n.N_ss,n.N_c) = PhiPtr->Phi_sc.array().colwise() * m.n_ss(ssIndices,dim);
			Phi.block(3*n.N_c+3*n.N_se, dim*(n.N_b)+n.N_c,	n.N_ss,n.N_se) = PhiPtr->Phi_se.array().colwise() * m.n_ss(ssIndices,dim);
			Phi.block(3*n.N_c+3*n.N_se, dim*(n.N_b)+n.N_c+n.N_se,	n.N_ss,n.N_ss) = PhiPtr->Phi_ss.array().colwise() * m.n_ss(ssIndices, dim);

			Eigen::ArrayXd diag2;
			diag2 = m.t1_ss(ssIndices,dim);
			// blocks related to the zero tangential contribution of the sliding surface nodes.
			Phi.block(3*n.N_c+3*n.N_se+n.N_ss, dim*(n.N_b)+n.N_c+n.N_se, n.N_ss,n.N_ss) = Eigen::MatrixXd(diag2.matrix().asDiagonal());
			diag2 = m.t2_ss(ssIndices,dim);
			Phi.block(3*n.N_c+3*n.N_se+2*n.N_ss, dim*(n.N_b)+n.N_c+n.N_se, n.N_ss,n.N_ss) = Eigen::MatrixXd(diag2.matrix().asDiagonal());

		}
	}

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


