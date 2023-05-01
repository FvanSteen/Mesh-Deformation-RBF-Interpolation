#include "rbfds.h"
#include <iostream>
#include <Eigen/Dense>
#include <chrono>


rbf_ds::rbf_ds(struct probParams& probParamsObject, Mesh& meshObject, getNodeType& n)
:rbfGenFunc(meshObject,probParamsObject)
{
	std::cout << "Initialised the ds class" << std::endl;
	if(params.dataRed){
		greedy g(m, params, disp, movingIndices, alpha, d);
		perform_rbf(n,g);
	}else{
		perform_rbf(n);
	}
}

void rbf_ds::perform_rbf(getNodeType& n){
	std::cout << "Performing RBF DS without data reduction" << std::endl;
	if(params.ptype){
		transform.cart_to_polar_cylindrical(m.coords, m.coords_polar_cylindrical);
	}

	for (int i = 0; i < params.steps; i++){
		std::cout << "Deformation step: " << i+1 << " out of "<< params.steps << std::endl;

		if(params.curved || i==0){
			m.getVecs(params);
			m.getMidPnts(params);
		}

		if(i==0){
			getDefVec(defVec_ds, n.N_c, n.mPtr);
		}

		getPhis(n, 0);

		getPhiDS(n, PhiPtr);

		performRBF_DS(n, PhiPtr);
	}

	if(params.ptype){
		transform.polar_cylindrical_to_cart(m.coords_polar_cylindrical,m.coords);
	}

}

void rbf_ds::perform_rbf(getNodeType& n, greedy& g){
	std::cout << "Performing RBF DS " << std::endl;

	int iter, lvl;
	bool iterating = true;



	for (int i = 0; i < params.steps; i++){
		std::cout << "Deformation step: " << i+1 << " out of "<< params.steps << std::endl;
		if(params.ptype){
			transform.cart_to_polar_cylindrical(m.coords, m.coords_polar_cylindrical);
		}
		iter = 0;
		lvl = 0;

		if(params.curved || i==0){
			m.getVecs(params);
			m.getMidPnts(params);
		}

		while(iterating){

			n.addControlNodes(g.maxErrorNodes, params.smode, m);
			std::cout << "control nodes:\n";
			for(auto x : *n.cPtr){
				std::cout << x << ", ";
			}
			std::cout << std::endl;
			// obtaining the interpolation matrices
			getPhis(n, iter);

			// obtaining the deformation vector
			if(lvl!=0){
				getDefVec(defVec_all, n, g.errorPrevLvl, n.N_c);
				performRBF(PhiPtr->Phi_cc, PhiPtr->Phi_ic, defVec_all, n.cPtr, n.iPtr, n.N_c);
			}else{
				getDefVec(defVec_ds, n.N_c, n.mPtr);
				getPhiDS(n, PhiPtr);
				performRBF_DS(n, PhiPtr);
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
				std::cout << "LEVEL: " << lvl << " HAS BEEN DONE" << std::endl;

				if(lvl == 0){
					setDefVec_all(n, PhiPtr);
				}

				g.setLevelParams(n,lvl, d, alpha, defVec_all,n.cPtr, n.N_c);

				lvl++;
				iter = -1;


				n.assignNodeTypesGrdy(m);

				if(iterating == false){
					g.getAlphaVector();
					g.setInitMaxErrorNodes();
				}
//				if(lvl == 3){
//					g.getAlphaVector();
//					updateNodes(n, defVec_all, g.d_step, g.alpha_step, g.ctrlPtr);
//
//					std::exit(0);
//
//				}

			}

			iter++;
		}


		if(params.curved == false && params.multiLvl == false){
			setDefVec_all(n, PhiPtr);
		}
		std::cout << "updating nodes\n";
		updateNodes(n, defVec_all, g.d_step, g.alpha_step, g.ctrlPtr);
//		transform.polar_cylindrical_to_cart(m.coords_polar_cylindrical, m.coords);
//		m.writeMeshFile(params.mesh_ifName, params.mesh_ofName);
//		std::exit(0);
		std::cout << "performing correction\n";
		g.correction( m,n,params.gamma, params.multiLvl);
		iterating = true;
	}
}

void rbf_ds::setDefVec_all(getNodeType& n, PhiStruct* PhiPtr){
	Eigen::ArrayXXd ds(n.N_se+n.N_ss,m.nDims);

	for(int dim = 0; dim < m.nDims; dim++){
		ds(Eigen::seqN(0, n.N_se),dim) = (PhiPtr->Phi_ec*alpha(Eigen::seqN(dim*(n.N_c),n.N_c)) ).array();
		ds(Eigen::seqN(n.N_se, n.N_ss),dim) = (PhiPtr->Phi_sc*alpha(Eigen::seqN(dim*(n.N_c),n.N_c)) ).array();

	}

	getDefVec(defVec_all,defVec_ds,n,ds, n.N_c, n.N_m);
}

void rbf_ds::performRBF_DS(getNodeType& n, PhiStruct* PhiPtr){
	alpha = Phi.colPivHouseholderQr().solve(defVec_ds);

	if(params.dataRed){
		d.resize(n.N_i,m.nDims);
	}

	if(params.curved){
		Eigen::ArrayXXd delta(n.N_se+n.N_ss, m.nDims), finalDef(n.N_se+n.N_ss,m.nDims);

		for (int dim = 0; dim < m.nDims; dim++){
			delta(Eigen::seqN(0,n.N_se),dim) = (PhiPtr->Phi_ec*alpha(Eigen::seqN(dim*(n.N_c),n.N_c))).array();
			delta(Eigen::seqN(n.N_se,n.N_ss),dim) = (PhiPtr->Phi_sc*alpha(Eigen::seqN(dim*(n.N_c),n.N_c)) ).array();
		}

		p.project(m, n, delta, finalDef);

		getDefVec(defVec_all, defVec_ds, n, finalDef, n.N_c, n.N_m);

		performRBF(PhiPtr->Phi_cc,PhiPtr->Phi_ic,defVec_all,n.cPtr, n.iPtr, n.N_c);

	}else{
		for (int dim = 0; dim < m.nDims; dim++){
			if(params.dataRed){
				d.col(dim) = PhiPtr->Phi_ic*alpha(Eigen::seqN(dim*n.N_c,n.N_c));
			}else{
				(*m.ptrCoords)(*n.iPtr, dim) += (PhiPtr->Phi_ic*alpha(Eigen::seqN(dim*n.N_c, n.N_c))).array();
				(*m.ptrCoords)(*n.sePtr, dim) += (PhiPtr->Phi_ec*alpha(Eigen::seqN(dim*(n.N_c),n.N_c))).array();
				(*m.ptrCoords)(*n.ssPtr, dim) += (PhiPtr->Phi_sc*alpha(Eigen::seqN(dim*(n.N_c),n.N_c)) ).array();
				(*m.ptrCoords)(*n.mPtr, dim) += (defVec_ds(Eigen::seqN(dim*n.N_m,n.N_m))).array();
			}
		}
	}


}



void rbf_ds::getPhiDS(getNodeType& n, PhiStruct* PhiPtr){

	// setting the overal matrix to correct size
	Phi = Eigen::MatrixXd::Zero(m.nDims*(n.N_c),m.nDims*(n.N_c));

	Eigen::ArrayXi seIndices;
	getIdxSlidingNodes(n.sePtr, seIndices, m.seNodes);

	Eigen::ArrayXi ssIndices;
	getIdxSlidingNodes(n.ssPtr, ssIndices, m.ssNodes);

	Eigen::VectorXd diagonalEdge(n.N_se), diagonalSurf(n.N_ss);
	for(int dim = 0; dim < m.nDims; dim++){
		// block related to the known displacements
		Phi.block(dim*n.N_m, dim*n.N_c, n.N_m, n.N_c) = PhiPtr->Phi_mc;


		// blocks related to the zero normal displacement conditions of the sliding edge nodes
		// first zero normal condition
		Phi.block(m.nDims*n.N_m, dim*n.N_c,	n.N_se,	n.N_c) = PhiPtr->Phi_ec.array().colwise() * m.n1_se(seIndices, dim);

		// second zero normal condition (only in 3D)
		if(m.nDims == 3){
			Phi.block(m.nDims*n.N_m+n.N_se, dim*n.N_c,	n.N_se,	n.N_c ) = PhiPtr->Phi_ec.array().colwise() * m.n2_se(seIndices, dim);
		}

		// block related to the zero tangential contribution of the sliding edge nodes
		diagonalEdge = m.t_se(seIndices,dim);
		Phi.block(m.nDims*n.N_m+(m.nDims-1)*n.N_se, dim*(n.N_c) + n.N_m, n.N_se, n.N_se) = diagonalEdge.asDiagonal();


		// block related to the zero normal displacement condition of the sliding surface nodes
		Phi.block(m.nDims*n.N_m+m.nDims*n.N_se, dim*n.N_c,	n.N_ss,n.N_c) = PhiPtr->Phi_sc.array().colwise() * m.n_ss(ssIndices, dim);

		// blocks related to the zero tangential contribution of the sliding surface nodes.
		// first tangential condition
		diagonalSurf = m.t1_ss(ssIndices,dim);
		Phi.block(m.nDims*n.N_m+m.nDims*n.N_se+n.N_ss, dim*(n.N_c)+n.N_m+n.N_se, n.N_ss,n.N_ss) = diagonalSurf.asDiagonal();

		// second tangential condition
		diagonalSurf= m.t2_ss(ssIndices,dim);
		Phi.block(m.nDims*n.N_m+m.nDims*n.N_se+(m.nDims-1)*n.N_ss, dim*(n.N_c)+n.N_m+n.N_se, n.N_ss,n.N_ss) = diagonalSurf.asDiagonal();

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


