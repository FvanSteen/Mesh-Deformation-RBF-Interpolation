#include "rbfGenFunc.h"
#include <iostream>
#include <Eigen/Dense>
#include <chrono>
#include <fstream>
#include <sstream>
#include "CoordTransform.h"

rbfGenFunc::rbfGenFunc(Mesh& meshObject, struct probParams& probParamsObject)
:m(meshObject), params(probParamsObject)
{

	// resizing the arrays containing the moving node indices and their displacement
	movingIndices.resize(m.N_nonzeroDisp);
	exactDisp.resize(m.N_nonzeroDisp,m.nDims);

	// reading the file containing the prescribed displacement of the moving nodes
	readDisplacementFile();


	// setting pointers to the displacement in Cartesian coordinates and indices
	exactDispPtr = &exactDisp;
	dispIdx = &movingIndices;

	// in case of a rotational periodic domain the prescribed deformation is transformed into polar/cylindrical coordinates
	if(params.ptype){
		// resize
		exactDisp_polar_cylindrical.resize(m.N_nonzeroDisp, m.nDims);

		// transformation
		transform.vector_cart_to_polar_cylindrical(exactDisp, exactDisp_polar_cylindrical, movingIndices, m.coords);

		// divide by number of deformation steps
		exactDisp_polar_cylindrical = exactDisp_polar_cylindrical/params.steps;

		// also a pointer for the displacement in polar/cylindrical coordinates
		disp = &exactDisp_polar_cylindrical;

	}else{
		disp = &exactDisp;
	}

	// divide by number of deformation steps
	exactDisp = exactDisp/params.steps;

	// Pointer to a structure containing the interpolation matrices
	PhiPtr = &Phis;
}


/* getPhis
 *
 * Obtains all the interpolation matrices the different submatrices like Phi_mm are part of the matrix Phi_cc.
 * Only Phi_cc and Phi_ic have to be determined and the other relevant matrices can be found as blocks in those matrices.
 */


void rbfGenFunc::getPhis(getNodeType& n, int iter){

	// in case of data reduction and the iteration is not equal to zero
	if(params.dataRed && iter != 0){
		// add rows and columns to the relevant matrices based on the added control nodes
		getPhisReduced(n);
	}else{
		// Calculate the entire interpolation matrix
		getPhisFull(n);
	}

	// Finding the interpolation matrices for the pseudo and direct sliding methods.
	// these matrices are part of the Phi_cc or Phi_ic matrices.
	if(params.smode == "ps"){
		Phis.Phi_mm.resize(n.N_m, n.N_m);
		Phis.Phi_mm = Phis.Phi_cc.block(0,0,n.N_m, n.N_m);

		Phis.Phi_em.resize(n.N_se, n.N_m);
		Phis.Phi_em = Phis.Phi_cc.block(n.N_m, 0, n.N_se, n.N_m);

		if(m.nDims == 3){
			Phis.Phi_meme.resize(n.N_m + n.N_se,n.N_m + n.N_se);
			Phis.Phi_meme = Phis.Phi_cc.block(0,0,n.N_m+n.N_se, n.N_m+n.N_se);

			Phis.Phi_sme.resize(n.N_ss, n.N_m+n.N_se);
			Phis.Phi_sme = Phis.Phi_cc.block(n.N_m+n.N_se, 0, n.N_ss, n.N_m+n.N_se);

			Phis.Phi_ic_reduced = Phis.Phi_ic.block(0,0,n.N_i- m.N_ss + n.N_ss, n.N_m+n.N_se);
		}

	}else if(params.smode == "ds"){

		Phis.Phi_mc.resize(n.N_m, n.N_c);
		Phis.Phi_mc = Phis.Phi_cc.block(0,0,n.N_m, n.N_c);

		Phis.Phi_ec.resize(n.N_se, n.N_c);
		Phis.Phi_ec = Phis.Phi_cc.block(n.N_m,0,n.N_se, n.N_c);

		Phis.Phi_sc.resize(n.N_ss, n.N_c);
		Phis.Phi_sc = Phis.Phi_cc.block(n.N_m+n.N_se,0,n.N_ss, n.N_c);
	}
}


/*
 * getPhisFull
 *
 * Determines the entire interpolation matrices
 */

void rbfGenFunc::getPhisFull(getNodeType& n){

	// interpolation matrix of the control nodes w.r.t. the control nodes
	getPhi(Phis.Phi_cc, n.cPtr, n.cPtr);
	// interpolation matrix of the internal nodes w.r.t. the control nodes
	getPhi(Phis.Phi_ic, n.iPtr, n.cPtr);

}

/* getPhisReduced
 *
 * Only add rows and columns based on the added control nodes
 */

void rbfGenFunc::getPhisReduced(getNodeType& n){
	// 0 adds row, 1 adds column, 2 adds both

	// removing rows from the Phi_ic matrix
	getReducedPhi(Phis.Phi_ic, n);

	// adjusting sizes of the interpolation matrices
	adjustPhi(Phis.Phi_cc, n, 2);
	adjustPhi(Phis.Phi_ic, n, 1);

	// adding the new rows and/ or columns to the matrices
	for(int i = 0; i < n.addedNodes.idx.size(); i++){
		switch (n.addedNodes.type[i]){

		// in case of moving node
			case 0:
				getPhi(Phis.Phi_cc, n.cPtr, n.cPtr, n.addedNodes.idx[i], 2);
				getPhi(Phis.Phi_ic, n.iPtr, n.cPtr, n.addedNodes.idx[i], 1);
				break;

		// in case of edge node
			case 1:
				getPhi(Phis.Phi_cc, n.cPtr, n.cPtr, n.addedNodes.idx[i]+n.N_m, 2);
				getPhi(Phis.Phi_ic, n.iPtr, n.cPtr, n.addedNodes.idx[i]+n.N_m, 1);
				break;

		// in case of surface node
			case 2:
				getPhi(Phis.Phi_cc, n.cPtr, n.cPtr, n.addedNodes.idx[i]+n.N_m+n.N_se, 2);
				getPhi(Phis.Phi_ic, n.iPtr, n.cPtr, n.addedNodes.idx[i]+n.N_m+n.N_se, 1);
				break;
		}

	}
}

/* getReducedPhi
 *
 * removing the rows of the added control nodes from the Phi_ic matrix
 *
 */

void rbfGenFunc::getReducedPhi(Eigen::MatrixXd& Phi, getNodeType& n){
	for(int i = 0; i < n.addedNodes.idx_i.size(); i++){
		Phi.middleRows(n.addedNodes.idx_i[i],Phi.rows()-n.addedNodes.idx_i[i]-1) = Phi.bottomRows(Phi.rows()-n.addedNodes.idx_i[i]-1);
	}
	Phi.conservativeResize(n.N_i,Phi.cols());
}


/* adjustPhi
 *
 * Adjusting the size of the relevant interpolation matrices based on the added control nodes
 * and shifting the rows/ columns accordingly
 */

void rbfGenFunc::adjustPhi(Eigen::MatrixXd& Phi, getNodeType& n,  int type){

	int idx = -1;
	switch(type){
		case 0:
			Phi.conservativeResize(Phi.rows()+n.addedNodes.idx.size(), Phi.cols());
			break;
		case 1:
			Phi.conservativeResize(Phi.rows(), Phi.cols()+n.addedNodes.idx.size());
			break;
		case 2:
			Phi.conservativeResize(Phi.rows()+n.addedNodes.idx.size(), Phi.cols()+n.addedNodes.idx.size());
			break;
	}


	for(int i = 0; i < n.addedNodes.idx.size(); i++){
		int idx_shift = 0;
		if(n.addedNodes.idx.size() > 1 && n.addedNodes.type[1] < n.addedNodes.type[0]){
			idx_shift = -1;
		}

		switch(n.addedNodes.type[i]){
		case 0:
			idx = n.addedNodes.idx[i];
			break;
		case 1:
			idx = n.addedNodes.idx[i] + n.N_m + idx_shift;
			break;
		case 2:
			idx = n.addedNodes.idx[i] + n.N_m + n.N_se + idx_shift;
			break;
		}

		if(idx < 0){
			idx = 0;
		}

		if( type > 0){
			int col = Phi.cols()-1;
			while(col > idx){
				Phi.col(col) = Phi.col(col-1);
				col--;
			}
		}
		if( type != 1){
			int row = Phi.rows()-1;
			while(row > idx){
				Phi.row(row) = Phi.row(row-1);
				row--;
			}
		}
	}

}


/* getPhi
 *
 * Adding rows or columns to the interpolation matrices based on the added control nodes
 *
 */

void rbfGenFunc::getPhi(Eigen::MatrixXd& Phi, Eigen::ArrayXi* idxSet1, Eigen::ArrayXi* idxSet2, int idx, int type){

	double distance;

	//depending on type add row, column or both
	switch (type) {

	// adding a row at row idx
		case 0:
			for(int col = 0; col < Phi.cols(); col++){
				distance = getDistance((*idxSet1)[idx],(*idxSet2)[col]);
				Phi(idx,col) = rbfEval(distance);
			}
			break;

	// adding column at column idx
		case 1:
			for(int row = 0; row < Phi.rows(); row++){
				distance = getDistance((*idxSet1)[row], (*idxSet2)[idx]);
				Phi(row,idx) = rbfEval(distance);
			}
			break;

	// adding both a row and column at idx
		case 2:
			for(int col = 0; col < Phi.cols(); col++){
				distance = getDistance((*idxSet1)[idx],(*idxSet2)[col]);
				Phi(idx,col) = rbfEval(distance);
			}

			Phi.col(idx) = Phi.row(idx);
			break;
	}
}

/* getPhi
 *
 * Determining all the elements of the interpolation matrix
 *
 */


void rbfGenFunc::getPhi(Eigen::MatrixXd& Phi, Eigen::ArrayXi* idxSet1, Eigen::ArrayXi* idxSet2){
	Phi.resize((*idxSet1).size(), (*idxSet2).size());

	double distance;
	for(int i=0; i<(*idxSet1).size();i++){
		for(int j=0; j<(*idxSet2).size();j++){

			// determine distance between nodes
			distance = getDistance((*idxSet1)[i],(*idxSet2)[j]);

			// set element to RBF evaluation value
			Phi(i,j) = rbfEval(distance);
		}
	}
}

/* getDistance
 *
 * Function for calculating the distance. Either Euclidean or for a periodic RBF
 * (Some duplicate code that could be improved)
 */



double rbfGenFunc::getDistance(int node1, int node2){
	double dist = 0.0;

	// in case of a rotationally periodic domain obtain the distance in polar/cylindrical coordinates
	// In polar/ cylindrical coordinates only periodicity in the theta direction is supported
	if(params.ptype){
		if(m.nDims == 2){
			if(params.pmode != "none"){
				dist = sqrt(pow(m.coords_polar_cylindrical(node1,0),2) + pow(m.coords_polar_cylindrical(node2,0),2) -2*m.coords_polar_cylindrical(node1,0)*m.coords_polar_cylindrical(node2,0)*cos(m.periodic_length/M_PI*sin( (m.coords_polar_cylindrical(node2,1)-m.coords_polar_cylindrical(node1,1))*M_PI/m.periodic_length) ));
			}else{
				dist = sqrt(pow(m.coords_polar_cylindrical(node1,0),2) + pow(m.coords_polar_cylindrical(node2,0),2) -2*m.coords_polar_cylindrical(node1,0)*m.coords_polar_cylindrical(node2,0)*cos( (m.coords_polar_cylindrical(node2,1)-m.coords_polar_cylindrical(node1,1))));
			}

		}else if(m.nDims == 3){
			if(params.pmode != "none"){
				dist = sqrt(pow(m.coords_polar_cylindrical(node1,0),2) + pow(m.coords_polar_cylindrical(node2,0),2) -2*m.coords_polar_cylindrical(node1,0)*m.coords_polar_cylindrical(node2,0)*cos(m.periodic_length/M_PI*sin( (m.coords_polar_cylindrical(node2,1)-m.coords_polar_cylindrical(node1,1))*M_PI/m.periodic_length)) + pow(m.coords_polar_cylindrical(node2,2) - m.coords_polar_cylindrical(node1,2),2) );
			}else{
				dist = sqrt(pow(m.coords_polar_cylindrical(node1,0),2) + pow(m.coords_polar_cylindrical(node2,0),2) -2*m.coords_polar_cylindrical(node1,0)*m.coords_polar_cylindrical(node2,0)*cos((m.coords_polar_cylindrical(node2,1)-m.coords_polar_cylindrical(node1,1))) + pow(m.coords_polar_cylindrical(node2,2) - m.coords_polar_cylindrical(node1,2),2) );
			}
		}
	}

	// Translational periodic domains:
	else{
		// if 2D
		if(m.nDims == 2){

			// if periodic
			if(params.pmode != "none"){
				dist = 0.0;

				// loopingthrough the dimensions to check if coordinate is periodic, if so then add the distance with the sine function.
				// otherwise use Euclidean distance
				for(int dim = 0; dim < m.nDims; dim++){
					if(dim == params.pDir){
						// sine function
						dist += pow(m.periodic_length/M_PI*sin( (m.coords(node1,dim)-m.coords(node2,dim))*M_PI/m.periodic_length),2); // Eq. 2.11 from the manuscript
					}
					else{
						// Euclidean distance
						dist += pow(m.coords(node1,dim)-m.coords(node2,dim),2);
					}
				}

				// take the square root
				dist = sqrt(dist);
			}
			// determine the distance for all coordinates at once
			else{
				dist = sqrt(pow(m.coords(node1,0)-m.coords(node2,0),2) + pow(m.coords(node1,1)-m.coords(node2,1),2) );
			}
		}

		// if 3D
		else if(m.nDims == 3){


			// in case of periodicity
			if(params.pmode != "none"){
				dist = 0;

				// looping through dimension to check for periodic coordinate
				for(int dim = 0; dim < m.nDims; dim++){

					if(dim == params.pDir){
						// add contribution of the sine function
						dist += pow(m.periodic_length/M_PI*sin( (m.coords(node1,dim)-m.coords(node2,dim))*M_PI/m.periodic_length),2);
					}else{
						// Euclidean distance contribution
						dist += pow(m.coords(node1,dim)-m.coords(node2,dim),2);
					}
				}

				// take square root
				dist = sqrt(dist);
			}

			// determine the distance for all coordinates at once
			else{
				dist = sqrt(pow(m.coords(node1,0)-m.coords(node2,0),2) + pow(m.coords(node1,1)-m.coords(node2,1),2) + pow(m.coords(node1,2)-m.coords(node2,2),2));
			}
		}
	}
	return dist;
}


/* getDefVec
 *
 * Determining the deformation vector by looping through the control nodes and if they are among the indices in the deformation file
 * then their displacement is included in the vector
 */


void rbfGenFunc::getDefVec(Eigen::VectorXd& defVec, int N_m, Eigen::ArrayXi* mPtr){

	defVec = Eigen::VectorXd::Zero(N_m*m.nDims);
	int idx;
	int size = (*mPtr).size();

	//loop through the control nodes
	for(int i = 0; i < size; i++){

		// Checking whether they are among the moving nodes
		idx = std::distance(std::begin(*dispIdx), std::find(std::begin(*dispIdx), std::end(*dispIdx),(*mPtr)(i)));
		if(idx!= (*dispIdx).size()){

			// setting displacements in the deformation vector
			for(int dim = 0; dim < m.nDims; dim++){
				defVec(dim*size+i) = (*disp)(idx,dim);
			}
		}
	}
}

/* getDefVec
 *
 * Merging the deformation vectors of the moving nodes defVec with the found projected displacment in finalDef, this results in defVec_b
 */

void rbfGenFunc::getDefVec(Eigen::VectorXd& defVec_b, Eigen::VectorXd& defVec, getNodeType& n,Eigen::ArrayXXd& finalDef, int N, int N_init){
	defVec_b = Eigen::VectorXd::Zero(N*m.nDims);

	for(int dim = 0; dim< m.nDims; dim++){
		defVec_b(Eigen::seqN(dim*N,N_init)) = defVec(Eigen::seqN(dim*N_init,N_init));
		defVec_b(Eigen::seqN(dim*N+N_init,N-N_init)) = finalDef.col(dim);
	}
}


/* readDisplacementFile
 *
 * Function that reads the displacement file and save the prescribed deformation and its indices
 *
 */

void rbfGenFunc::readDisplacementFile(){
	std::cout << "Reading the displacement file :" << params.dispFile << std::endl;

	std::string line;
	std::ifstream file(params.directory + "\\" + params.dispFile);

	int lineNo = 0;
	if(file.is_open()){
		while(getline(file, line)){
			std::stringstream ss(line);
			if(m.nDims == 2){
				ss >> movingIndices(lineNo) >> exactDisp(lineNo,0) >> exactDisp(lineNo,1);
			}else if(m.nDims == 3){
				ss >> movingIndices(lineNo) >> exactDisp(lineNo,0) >> exactDisp(lineNo,1) >> exactDisp(lineNo,2);
			}
			lineNo++;
		}
	}else{ std::cout << "Unable to open the displacement file" << std::endl;
		std::exit(0);
	}
}



/* rbfEval
 *
 * Evaluates the radial basis function
 *
 */


double rbfGenFunc::rbfEval(double distance){
	double xi = distance/m.r;	// distance scaled by support radius
	double f_xi = 0;
	if(xi < 1){
		f_xi = pow((1-(distance/m.r)),4)*(4*(distance/m.r)+1);
	}
	return f_xi;
}


/*getDefVec
 *
 * Establishes the deformation vector in case of the multi-level greedy algorithm from the second level onwards
 * The deformation is equal to the errors of the previous level
 *
 */

void rbfGenFunc::getDefVec(Eigen::VectorXd& defVec, getNodeType& n, Eigen::ArrayXXd& errors, int N_m){

	defVec.resize(m.nDims*N_m);
	for(int dim = 0; dim< m.nDims; dim++){
		defVec(Eigen::seqN(N_m*dim,N_m)) = -errors(n.cNodesIdx, dim);
	}
}

/* performRBF
 *
 * function that performs the RBF interpolation. It solves for the interpolation coefficients and find the displacement of the internal nodes
 *
 * in case of data reduction the displacement is determined. Otherwise, the displacement is directly applied
 * Here the colPivHouseholderQr() method of the Eigen library is used for solving the interpolation system
 * Other methods are available. See https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html for more details
 *
 */

void rbfGenFunc::performRBF(Eigen::MatrixXd& Phi_cc, Eigen::MatrixXd& Phi_ic, Eigen::VectorXd& defVec, Eigen::ArrayXi* cNodes, Eigen::ArrayXi* iNodes, int N){

	// resizing of the interpolation coefficients array
	alpha.resize(N*m.nDims);

	// resizing of the displacement array
	if(params.dataRed){
		d.resize((*iNodes).size(),m.nDims);
	}

	// solving for alpha
	for(int dim = 0; dim < m.nDims; dim++){

		alpha(Eigen::seqN(dim*N,N)) = Phi_cc.colPivHouseholderQr().solve(defVec(Eigen::seqN(dim*N,N)));

		// find resulting displacement or apply displacement directly to the coordinates
		if(params.dataRed){
			d.col(dim) = Phi_ic*alpha(Eigen::seqN(dim*N,N));
		}else{
			(*m.ptrCoords)(*iNodes,dim) += (Phi_ic*alpha(Eigen::seqN(dim*N,N))).array();
			(*m.ptrCoords)(*cNodes,dim) += defVec(Eigen::seqN(dim*N,N)).array();
		}
	}
}

/* updateNodes
 *
 * Updating the nodes once the greedy algorithm has reached the desired tolerance
 *
 */

void rbfGenFunc::updateNodes(getNodeType& n, Eigen::VectorXd& defVec, Eigen::ArrayXXd* d_step, Eigen::VectorXd* alpha_step, Eigen::ArrayXi* ctrlPtr){
	Eigen::MatrixXd Phi_icGrdy;

	int N_m;
	Eigen::ArrayXi* ptr;

	if(params.multiLvl){
		// in case of multilevel ptr is the go ptr that points to all the ctrl points that are selected so far
		ptr = ctrlPtr;
		N_m = (*ctrlPtr).size();
	}else{
		if(params.smode == "none"){
			ptr = n.mPtr;
			N_m = n.N_m;
		}else{
			ptr = n.cPtr;
			N_m = n.N_c;
		}
	}

	// determine interpolation matric of the actual internal nodes w.r.t. the selected control nodes

	getPhi(Phi_icGrdy,n.iPtrGrdy,ptr);
	for(int dim = 0; dim < m.nDims; dim++){

		// apply deformation of the control nodes
		if(params.multiLvl == false){
			(*m.ptrCoords)(*ptr,dim) += (defVec(Eigen::seqN(dim*N_m,N_m))).array();
		}

		// apply deformation of the internal nodes
		(*m.ptrCoords)(*n.iPtrGrdy,dim) +=  (Phi_icGrdy*(*alpha_step)(Eigen::seqN(dim*N_m,N_m))).array();
	}

	// if rotational periodic transform back to cartesian coordinates
	if(params.ptype){
		transform.polar_cylindrical_to_cart(m.coords_polar_cylindrical, m.coords);
	}

	// apply the deformation of the unselected boundary nodes
	m.coords(*n.iPtr, Eigen::all) += *d_step;
}
