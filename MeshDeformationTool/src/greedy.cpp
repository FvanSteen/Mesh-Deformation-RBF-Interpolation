#include "greedy.h"

#include <iostream>
#include <Eigen/Dense>
#include <iterator>
#include <math.h>
#include "SPDS.h"
#include "nanoflann.hpp"


greedy::greedy(Mesh& m, probParams& params, Eigen::ArrayXXd* disp, Eigen::ArrayXi& movingIndices, Eigen::VectorXd& alpha, Eigen::ArrayXXd& d)
{

	// setting pointers

	mPtr = &m;						// to the mesh class
	exctDispPtr = disp;				// to the displacement
	mIdxPtr = &movingIndices;		// to the moving indices
	p = &params;					// to the problem parameters


	// set pointers to the interpolation coefficients, displacement and control nodes for the multi level greedy algorithm
	if((*p).multiLvl){
		alpha_step = &alphaGrdy;
		d_step = &delta;
		ctrlPtr = &ctrlNodesAll;
	}else{
		alpha_step = &alpha;
		d_step = &d;
	}

	// Setting the initial control nodes based on the max displacement and a double edged selection
	setInitMaxErrorNodes();
}


/* getError
 *
 * This function calculates the error at the boundary nodes
 * The error calculation is done in Cartesian coordinates since the projection algorithm relies on finding the nearest midpoint
 * in terms of the Euclidean distance
 */

void greedy::getError(getNodeType& n, Eigen::ArrayXXd& d, int lvl){

	// resizing of the error and error angle arrays
	error.resize(d.rows(),(*mPtr).nDims);
	errorAngle.resize(d.rows(), (*mPtr).nDims-1);

	// if the domain is rotationally periodic then the displacement is transformed to Cartesian coordinates
	if((*p).ptype){
		if( (*p).smode == "ps" &&  d.rows() == (*n.iPtr_reduced).size())
			transform.disp_to_cart(d, *n.iPtr_reduced, (*n.iPtr_reduced).size(), *mPtr);
		else
			transform.disp_to_cart(d, *n.iPtr, n.N_i, *mPtr);
	}


	// obtaining error based on whether a single or multi-level greedy algorithm is used
	if((*p).multiLvl && lvl>0){
		getErrorMultiLvl(n,d);
	}else{
		getErrorSingleLvl(n,d);
	}



    // if all boundary nodes are selected as control nodes then the array with displacements has size zero and no maxErrorNodes can be selected
	if(d.rows() == 0){
		maxError = 0;
		maxErrorNodes.resize(1);
		maxErrorNodes << -1;
	}
	// Finding node with the maximum error magnitude and second node is found with the double edged selection if required
	else{

		// Select node with max error magnitude
		int idxMax;
		Eigen::ArrayXd errorSquaredNorm;
		errorSquaredNorm = error.rowwise().squaredNorm();
		errorSquaredNorm.maxCoeff(&idxMax);

		// and the maximum error magnitude
		maxError = error.row(idxMax).matrix().norm();

		// double edged selection
		if((*p).doubleEdge){

			// establish error angles
			getErrorAngle();

			// find node with largest error with an angle of 90 deg or more w.r.t. the max error node
			int idxMaxAngle = getDoubleEdgeError(idxMax, n.N_i, error);

			// if no node is found with the double edged selection then only the max error node is selected
			if(idxMaxAngle == -1){
				maxErrorNodes.resize(1);
				maxErrorNodes << (*n.iPtr)(idxMax);
			}
			// else both nodes are included
			else{
				maxErrorNodes.resize(2);
				maxErrorNodes << (*n.iPtr)(idxMax), (*n.iPtr)(idxMaxAngle);
			}
		}
		// no double edged selection
		else{
			maxErrorNodes.resize(1);
			maxErrorNodes << (*n.iPtr)(idxMax);
		}

		// If the addition of 2 control nodes will mean that the levelsize exceeds the set level size then only a single node will be included
		if( (*p).mCrit == "size" && n.N_c + maxErrorNodes.size() > (*p).lvlSize){
			maxErrorNodes.conservativeResize(1);
		}
	}
}

/* getErrorMovingNodes
 *
 * Finds the error of the moving nodes by considering the difference between the found displacement and the prescribed displacement
 *
 */


void greedy::getErrorMovingNodes(Eigen::ArrayXi* nodes, Eigen::ArrayXXd& d, size_t N){
	int idx_m;

	for(size_t i = 0; i < N; i++){
		idx_m = std::distance(std::begin(*mIdxPtr), std::find(std::begin(*mIdxPtr),std::end(*mIdxPtr),(*nodes)(i)));

		// if among the nodes with nonzero displacement
		if(idx_m !=  (*mIdxPtr).size()){
			error.row(i) = d.row(i) - (*exctDispPtr).row(idx_m);
		}
		// else the node is fixed and has a zero prescribed displacement
		else{
			error.row(i) = d.row(i);
		}
	}
}


/* getErrorAngle
 *
 * Function for finding the angle of the error. In case of 2D a single angle is determined and for 3D, 2 angles are found
 *
 */


void greedy::getErrorAngle(){
	for(int i =0 ; i < errorAngle.rows(); i++){
		errorAngle(i,0) = atan2(error(i,1),error(i,0));
		if((*mPtr).nDims == 3){
			errorAngle(i,1) = atan2(  sqrt(pow(error(i,0),2) + pow(error(i,1),2)), error(i,2));
		}
	}
}

/* getErrorSingleLvl
 *
 * Function for finding the error in case of the single level greedy algorithm
 *
 */

void greedy::getErrorSingleLvl(getNodeType& n, Eigen::ArrayXXd& d){

	// ending index of the different node types
	size_t m_end, se_end, ss_end;
	m_end = (*mPtr).N_m-n.N_m;
	se_end = m_end + (*mPtr).N_se-n.N_se;

	// get the error of the moving nodes
	getErrorMovingNodes(n.iPtr,  d,  m_end);


	// get error of the sliding edge nodes
	if(m_end != size_t(d.rows())){
		SPDSobj.projectEdge(*mPtr, n.iPtr, d, error, m_end, se_end, 0, (*p).ptype);
	}

	// get the error of the sliding surface nodes
	if(se_end != size_t(d.rows())){
		ss_end = se_end + (*mPtr).N_ss-n.N_ss;
		SPDSobj.projectSurf(*mPtr, n.iPtr, d, error, se_end, ss_end, 0, (*p).ptype);
	}


}

/* getErrorMultiLvl
 *
 * Find the error in case of the multi-level greedy algorithm
 * The error is the difference between the found displacement and the error of the previous level,
 * which is the deformation applied to the current level
 */

void greedy::getErrorMultiLvl(getNodeType& n, Eigen::ArrayXXd& d){
	error = d + errorPrevLvl(n.iNodesIdx, Eigen::all);
}



/* getDoubleEdgeError
 *
 * determining the second selected control node with a double edged control node selection routine
 *
 */
int greedy::getDoubleEdgeError(int idxMax, int N_i, Eigen::ArrayXXd& error){


	Eigen::ArrayXd maxErrorAngle = errorAngle.row(idxMax);
	// obtaining the error direction w.r.t. the direction of the max magnitude error
	errorAngle.rowwise() -= maxErrorAngle.transpose();

	// array containing the errors of the nodes that have an angle > 90 deg w.r.t. the max magnitude error direction
	Eigen::ArrayXd largeAngleError(error.rows()-1);
	// corresponding indices
	Eigen::ArrayXi largeAngleIdx(error.rows()-1);
	// Size of these arrays is unknown beforehand, but has a maximum of N_i-1 elements



	// counting the number of error with an relative angle > 90 def
	int cnt = 0;
	// reference angle of 90 deg.
	double refAngle = M_PI/2;
	// second reference angle of 360-90 = 270 deg.
	double refAngle2 = 3*refAngle;

	// for each boundary node
	for(int i=0; i<error.rows(); i++){
		// if the relative angle is between 90 and 270 degrees the squared norm is included in the array and its index is saved.

		// 2D just a single angle is considered
		if (errorAngle.cols() == 1){
			if(abs(errorAngle(i,0)) >= refAngle && abs(errorAngle(i,0)) <= refAngle2){
				largeAngleError(cnt) = error.row(i).matrix().squaredNorm();
				largeAngleIdx(cnt) = i;
				cnt++;
			}
		}

		// For 3D the two angles need to be considered in order to have
		else{
			if(abs(errorAngle(i,0)) >= refAngle && abs(errorAngle(i,0)) <= refAngle2 && abs(errorAngle(i,1)) >= refAngle && abs(errorAngle(i,1)) <= refAngle2){
				largeAngleError(cnt) = error.row(i).matrix().squaredNorm();
				largeAngleIdx(cnt) = i;
				cnt++;
			}
		}
	}

	// if in 3D no nodes are found then only the first angle is considered
	if (cnt == 0 && errorAngle.cols() == 2){
		for(int i=0; i<error.rows(); i++){
			if(abs(errorAngle(i,0)) >= refAngle && abs(errorAngle(i,0)) <= refAngle2){
				largeAngleError(cnt) = error.row(i).matrix().squaredNorm();
				largeAngleIdx(cnt) = i;
				cnt++;
			}
		}
	}

	int idxMaxLargeAngle;

	if(cnt!= 0){
		// resizing the largeAngleError array to only include the values included in the for-loop
		largeAngleError.conservativeResize(cnt);

		// finding the index where the errors is maximum
		largeAngleError.maxCoeff(&idxMaxLargeAngle);

		// setting the integer equal to the corresponding index of all boundary nodes considered
		idxMaxLargeAngle = largeAngleIdx(idxMaxLargeAngle);
	}
	// if no node is found then set index to -1
	else{
		std::cout << "no error found with a large angle" << std::endl;
		idxMaxLargeAngle = -1;
	}

	return idxMaxLargeAngle;


}


/* correction
 * Applies the error as deformation to the boundary to ensure that all boundary points have zero error
 *
 */

void greedy::correction(Mesh& m, getNodeType& n, double& gamma, bool& multiLvl){


	// get correct error pointer for single/multi level greedy
	Eigen::ArrayXXd* errorPtr;
	if(multiLvl){
		errorPtr = &errorPrevLvl;
	}else{
		errorPtr = &error;
	}

	// if there is a nonzero error array
	if((*errorPtr).rows() != 0){

		// integer for storing the index where the error is largest
		int idxMax;
		// finding the node with the maximum error
		Eigen::ArrayXd errorSquaredNorm;
		errorSquaredNorm = (*errorPtr).rowwise().squaredNorm();
		errorSquaredNorm.maxCoeff(&idxMax);

		// returning the largest error
		double maxError = error.row(idxMax).matrix().norm();

		// applying the correction using a nearest neighbour search for all internal nodes
		SPDSobj.kdt_NNSearch(*n.iPtr, *n.iPtrGrdy,  m.coords, m.nDims, gamma, maxError, errorPtr);
	}

	m.coords(*n.iPtr , Eigen::all) -= *errorPtr;

}


/* setLevelParams
 *
 * After satisfying the level criterium the level parameters are saved with this function
 * All used control nodes are stored along with the interpolation coefficients and displacement of the level
 * This information is later used in the update of the node coordinates
 */


void greedy::setLevelParams(getNodeType& n, int lvl, Eigen::ArrayXXd& d, Eigen::VectorXd& alpha,Eigen::VectorXd& defVec, Eigen::ArrayXi* cPtr, int N_c){

	// setting maxErrorPrevLvl to the current error
	maxErrorPrevLvl = maxError;

	// initializing variables that keep track of the control node information
	if(lvl == 0){
		delta = Eigen::ArrayXXd::Zero((*mPtr).N_m+(*mPtr).N_se+(*mPtr).N_ss, (*mPtr).nDims);
		alphaTotal.resize(0,0);
		ctrlNodesAll.resize(0);
	}

	// array with the new control nodes and interpolation coefficients
	Eigen::ArrayXi newCtrlNodes(N_c);
	Eigen::ArrayXXd newAlpha(N_c, (*mPtr).nDims);


	int cnt = 0;
	int idx, dim;
	for(int i = 0; i < N_c; i++){

		// check for each control if its already included
		idx = std::distance(std::begin(ctrlNodesAll), std::find(std::begin(ctrlNodesAll), std::end(ctrlNodesAll), (*cPtr)(i)));
		if(idx == ctrlNodesAll.size()){

			// if not included then add to newControlNodes list
			newCtrlNodes(cnt) = i;

			// adding new alpha to list
			for(dim =0; dim < (*mPtr).nDims; dim++){
				newAlpha(cnt, dim) = alpha(dim*N_c+i);
			}
			cnt++;

		// if its already included then add to accumulative sum
		}else{
			for(dim = 0; dim < (*mPtr).nDims; dim++){
				alphaTotal(idx,dim) += alpha(dim*N_c+i);
			}
		}
	}

	// resizing
	ctrlNodesAll.conservativeResize(ctrlNodesAll.size()+cnt);
	alphaTotal.conservativeResize(alphaTotal.rows()+cnt,(*mPtr).nDims);

	// adding
	ctrlNodesAll(Eigen::lastN(cnt)) = (*cPtr)(newCtrlNodes(Eigen::seqN(0,cnt)));
	alphaTotal(Eigen::lastN(cnt), Eigen::all) = newAlpha(Eigen::seqN(0,cnt), Eigen::all);

	// the deformation of the "internal" nodes
	delta(n.iNodesIdx,Eigen::all) += d;

	// deformation of the control nodes
	for(int dim = 0; dim < (*mPtr).nDims; dim++ ){
		delta(n.cNodesIdx,dim) += (defVec(Eigen::seqN(dim*N_c, N_c))).array();
	}
	std::cout << "Done\n";
	// set error of the prev lvl equal to current error for the next level
	errorPrevLvl = Eigen::ArrayXXd::Zero((*mPtr).N_m+(*mPtr).N_se+(*mPtr).N_ss, (*mPtr).nDims);
	if((*p).ptype){
		errorPrevLvl(n.iNodesIdx,Eigen::all) = errorPolarCylindrical;
	}else{
		errorPrevLvl(n.iNodesIdx,Eigen::all) = error;
	}

}

/* getAlphaVector
 *
 * sets up an array containing the interpolation coefficients of the selected control nodes of all levels of the multi-level algorithm
 * Is used in the update of the coordinates
 */

void greedy::getAlphaVector(){

	alphaGrdy.resize(alphaTotal.size());

	switch(alphaTotal.cols()){
		case 2:
			alphaGrdy << alphaTotal.col(0), alphaTotal.col(1);
			break;
		case 3:
			alphaGrdy << alphaTotal.col(0), alphaTotal.col(1), alphaTotal.col(2);
			break;
	}

}



/* setInitMaxErrorNodes
 *
 * Find the first control node based on the maximum prescribed deformation. A second one is found using the double edged selection
 *
 */

void greedy::setInitMaxErrorNodes(){

	Eigen::Array2i idxMax;
	Eigen::ArrayXd exctDispSquaredNorm;

	// Finding node with the maximum displacement
	exctDispSquaredNorm =(*exctDispPtr).rowwise().squaredNorm();
	exctDispSquaredNorm.maxCoeff(&idxMax(0));

	// for the multi level greedy the max error of the previous level is set, either based on max displacement or previous level error
	maxErrorPrevLvl = (*exctDispPtr).row(idxMax(0)).matrix().norm();


	// finding the second node using the double edged selection
	if((*p).doubleEdge){

		int N = (*mIdxPtr).size();

		// finding the error angles of the considered nodes
		errorAngle.resize(N, (*mPtr).nDims -1);

		for(int i = 0; i<N; i++){
			errorAngle(i,0) = atan2((*exctDispPtr)(i,1), (*exctDispPtr)(i,0));
			if((*mPtr).nDims ==3){
				errorAngle(i,1) = atan2(sqrt(pow((*exctDispPtr)(i,0),2) + pow((*exctDispPtr)(i,1),2)), (*exctDispPtr)(i,2));
			}
		}

		// determining the node with the max error and an error angle difference larger than 90 deg w.r.t. the first error node
		idxMax(1) = getDoubleEdgeError(idxMax(0), N, (*exctDispPtr));


		// if no node is found with the double edged selection, the node farthest away from the first node is selected as second node
		if (idxMax(1) == -1){
			Eigen::ArrayXXd movingCoords;
			movingCoords = (*mPtr).coords((*mIdxPtr),Eigen::all);

			Eigen::ArrayXd dist;
			dist = (movingCoords.rowwise() - (*mPtr).coords.row((*mIdxPtr)(idxMax(0)))).rowwise().norm();
			dist.maxCoeff(&idxMax(1));
		}
		maxErrorNodes.resize(2);
		maxErrorNodes << (*mIdxPtr)(idxMax);
	}
	// if not using double edged selection than only the node with the max error magnitude is selected
	else{
		maxErrorNodes.resize(1);
		maxErrorNodes << (*mIdxPtr)(idxMax(0));
	}




}




