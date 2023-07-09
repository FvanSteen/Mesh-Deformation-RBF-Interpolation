#include "rbfds.h"
#include <iostream>
#include <Eigen/Dense>
#include <chrono>

rbf_ds::rbf_ds(struct probParams& probParamsObject, Mesh& meshObject, getNodeType& n)
:rbfGenFunc(meshObject,probParamsObject)
{
	// if data reduction is applied then initialise the class with greedy function and perform the rbf interpolation
	if(params.dataRed){
		// creating a convergence history file
		w.createConvHistFile(params.directory);
		greedy g(m, params, exactDispPtr, movingIndices, alpha, d);
		perform_rbf(n,g);
	}else{
		perform_rbf(n);
	}
}

void rbf_ds::perform_rbf(getNodeType& n){
	std::cout << "Performing direct sliding RBF interpolation" << std::endl;

	// transform to polar/ cylindrical coordinates if required
	if(params.ptype){
		transform.cart_to_polar_cylindrical(m.coords, m.coords_polar_cylindrical);
	}


	// for-loop going through the deformation steps
	for (int i = 0; i < params.steps; i++){
		std::cout << "Deformation step: " << i+1 << " out of "<< params.steps << std::endl;

		if(params.curved || i==0){
			// obtaining midpoints and corresponding vectors of boundary elements
			m.getMidPnts(params);
			// obtaining vectors at the sliding nodes
			m.getVecs(params);
		}


		// obtaining the deformation vector
		if(i==0){
			getDefVec(defVec_ds, n.N_c, n.mPtr);
		}

		// obtian deformation vector
		getPhis(n, 0);

		// setup interpolation matrix for the direct sliding
		getPhiDS(n, PhiPtr);

		// perform the direct sliding RBF interpolation
		performRBF_DS(n, PhiPtr);

	}

	// transform to Cartesian coordinates
	if(params.ptype){
		transform.polar_cylindrical_to_cart(m.coords_polar_cylindrical,m.coords);
	}

}

void rbf_ds::perform_rbf(getNodeType& n, greedy& g){
	std::cout << "Performing direct sliding RBF interpolation with data reduction" << std::endl;

	// clocks for keeping track of CPU time
	std::clock_t s = std::clock();
	std::clock_t e = std::clock();

	// integers for number of iterations and levels
	int iter, lvl;

	// iterating bool set to True
	bool iterating = true;

	// for-loop going through the deformation steps
	for (int i = 0; i < params.steps; i++){
		std::cout << "Deformation step: " << i+1 << " out of "<< params.steps << std::endl;

		// transformation to polar/ cylindrical coordinates
		if(params.ptype){
			transform.cart_to_polar_cylindrical(m.coords, m.coords_polar_cylindrical);
		}

		iter = 0;
		lvl = 0;

		// in case of curved boundaries the midpoints and normal/ tangential vectors are updated
		// or in case of the first step the midpoints and vectors are determined
		if(params.curved || i==0){
			m.getMidPnts(params);
			m.getVecs(params);
		}


		while(iterating){

			// adding control nodes
			n.addControlNodes(g.maxErrorNodes, params.smode, m);

			// obtaining the interpolation matrices
			getPhis(n, iter);


			if(lvl!=0){
				// obtaining the deformation vector
				getDefVec(defVec_all, n, g.errorPrevLvl, n.N_c);

				// performing regular RBF interpolation
				performRBF(PhiPtr->Phi_cc, PhiPtr->Phi_ic, defVec_all, n.cPtr, n.iPtr, n.N_c);
			}else{
				// obtaining deformation vector
				getDefVec(defVec_ds, n.N_c, n.mPtr);

				// setup the direct sliding interpolation matrix
				getPhiDS(n, PhiPtr);

				// perform the direct sliding rbf interpolation
				performRBF_DS(n, PhiPtr);
			}

			// obtaining error
			g.getError(n, d, lvl);
			std::cout << "error: \t"<< g.maxError <<" at node: \t" << g.maxErrorNodes(0)<< std::endl;

			// writing intermediate results of the convergence history
			e = std::clock();
			long double time_elapsed_ms =  1000.0*(e-s) / CLOCKS_PER_SEC;
			std::cout << "CPU time: " << time_elapsed_ms/1000 << " ms\n";
			w.setIntResults(params.directory, i, lvl, g.maxError, time_elapsed_ms,  n.N_c);


			// check if the error tolerance is reached
			if(g.maxError < params.tol){
				iterating = false;
				if(params.multiLvl == false){
					g.maxErrorNodes.resize(0);
				}
			}

			// check for the multi-level criteria
			if(params.multiLvl && (g.maxError/g.maxErrorPrevLvl < params.tolCrit || iterating == false)){
				std::cout << "Level: " << lvl << " has been done" << std::endl;

				// establish the deformation vector
				if(lvl == 0){
					setDefVec_all(n, PhiPtr);
				}

				// saving level parameters
				g.setLevelParams(n,lvl, d, alpha, defVec_all,n.cPtr, n.N_c);

				lvl++;
				iter = -1;

				// reset node types for the next level
				n.assignNodeTypesGrdy(m);


				// if the tolerance is reached
				if(iterating == false){
					g.getAlphaVector();
					g.setInitMaxErrorNodes();
				}
			}
			iter++;
		}

		// setup the deformation vector required for updating the nodes
		if(params.curved == false && params.multiLvl == false){
			setDefVec_all(n, PhiPtr);
		}

		// performing the correction and updating the node coordinates
		g.correction( m,n,params.gamma, params.multiLvl);
		updateNodes(n, defVec_all, g.d_step, g.alpha_step, g.ctrlPtr);

		iterating = true;

	}
}


/* setDefVec_all
 *
 * Determine the displacement of the edge and surface nodes and place them in a deformation vector
 *
 */

void rbf_ds::setDefVec_all(getNodeType& n, PhiStruct* PhiPtr){
	Eigen::ArrayXXd ds(n.N_se+n.N_ss,m.nDims);

	for(int dim = 0; dim < m.nDims; dim++){
		ds(Eigen::seqN(0, n.N_se),dim) = (PhiPtr->Phi_ec*alpha(Eigen::seqN(dim*(n.N_c),n.N_c)) ).array();
		ds(Eigen::seqN(n.N_se, n.N_ss),dim) = (PhiPtr->Phi_sc*alpha(Eigen::seqN(dim*(n.N_c),n.N_c)) ).array();

	}
	getDefVec(defVec_all,defVec_ds,n,ds, n.N_c, n.N_m);
}


/* performRBF_DS
 *
 * Performing the direct sliding RBF interpolation
 *
 */

void rbf_ds::performRBF_DS(getNodeType& n, PhiStruct* PhiPtr){

	// determine the interpolation coefficients
	alpha = Phi.colPivHouseholderQr().solve(defVec_ds);

	// resize the displacement array
	if(params.dataRed){
		d.resize(n.N_i,m.nDims);
	}

	// in case of curved boundaries a projection is performed to ensure that nodes are on the boundary
	if(params.curved){

		Eigen::ArrayXXd delta(n.N_se+n.N_ss, m.nDims), finalDef(n.N_se+n.N_ss,m.nDims);

		// displacement resulting from the interpolation
		for (int dim = 0; dim < m.nDims; dim++){
			delta(Eigen::seqN(0,n.N_se),dim) = (PhiPtr->Phi_ec*alpha(Eigen::seqN(dim*(n.N_c),n.N_c))).array();
			delta(Eigen::seqN(n.N_se,n.N_ss),dim) = (PhiPtr->Phi_sc*alpha(Eigen::seqN(dim*(n.N_c),n.N_c)) ).array();
		}

		// projection of the sliding nodes
		p.project(m, n, delta, finalDef, params.ptype);

		// setup deformation vector that includes the new position of the sliding nodes
		getDefVec(defVec_all, defVec_ds, n, finalDef, n.N_c, n.N_m);


		// perform regular RBF interpolation with the new deformation vector
		performRBF(PhiPtr->Phi_cc,PhiPtr->Phi_ic,defVec_all,n.cPtr, n.iPtr, n.N_c);

	}
	// Perform the direct sliding RBF interpolation
	else{
		for (int dim = 0; dim < m.nDims; dim++){

			// determine resulting displacement of the boundary nodes
			if(params.dataRed){
				d.col(dim) = PhiPtr->Phi_ic*alpha(Eigen::seqN(dim*n.N_c,n.N_c));
			}
			// Directly apply the resulting deformation to the coordinates
			else{
				(*m.ptrCoords)(*n.iPtr, dim) += (PhiPtr->Phi_ic*alpha(Eigen::seqN(dim*n.N_c, n.N_c))).array();
				(*m.ptrCoords)(*n.sePtr, dim) += (PhiPtr->Phi_ec*alpha(Eigen::seqN(dim*(n.N_c),n.N_c))).array();
				(*m.ptrCoords)(*n.ssPtr, dim) += (PhiPtr->Phi_sc*alpha(Eigen::seqN(dim*(n.N_c),n.N_c)) ).array();
				(*m.ptrCoords)(*n.mPtr, dim) += (defVec_ds(Eigen::seqN(dim*n.N_m,n.N_m))).array();
			}
		}
	}


}



/* getPhiDS
 *
 * Function that sets up the interpolation matrix for the direct sliding RBF interpolation.
 * See section 2.4.1 and equation 2.22 of the manuscript for details
 */

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

/* getIdxSlidingNodes
 *
 * Function that finds the indices of the sliding nodes that are selected as control nodes
 *
 */

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


