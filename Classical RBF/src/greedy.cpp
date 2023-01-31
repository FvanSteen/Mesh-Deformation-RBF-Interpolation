#include "greedy.h"
#include "projection.h"
#include <iostream>
#include <Eigen/Dense>
#include <iterator>
#include <chrono>
#include <math.h>
greedy::greedy()
{}


void greedy::getError(Mesh& m, getNodeType& n, Eigen::ArrayXXd& d, double& maxError, Eigen::ArrayXi& maxErrorNodes, Eigen::ArrayXi& movingIndices, Eigen::ArrayXXd& exactDisp, Eigen::VectorXd& pnVec, projection* projPtr, bool multiLvl, int lvl){

	if(error.rows() == 0){
		error.resize(n.N_i,m.nDims);
	}

	// defining array for the error directions
	Eigen::ArrayXd errorAngle(n.N_i);


	if(multiLvl && lvl>0){

		getErrorMultiLvl(n,errorAngle,d,m,movingIndices,pnVec,projPtr);

	}else{
		getErrorSingleLvl(m,n,errorAngle,d,movingIndices,exactDisp,pnVec, projPtr);
	}


	// find index of largest error
	int idxMax;
	error.rowwise().squaredNorm().maxCoeff(&idxMax);

	// and the maximum error magnitude
	maxError = error.row(idxMax).matrix().norm();

	// find index of largest error where there is a 90 degree difference with the max magnitude direction
	int idxMaxAngle = getDoubleEdgeError(errorAngle, idxMax, n.N_i);

	// making array with the max error indices
	maxErrorNodes.resize(2);
	maxErrorNodes << (*n.iPtr)(idxMax), (*n.iPtr)(idxMaxAngle);

//	if(lvl==1 && n.N_mStd == 8){
//		std::cout << maxErrorNodes << std::endl;
//		m.coords(*n.iPtr,Eigen::all) += d;
//		m.writeMeshFile();
//		std::exit(0);
//	}



}

void greedy::getErrorSingleLvl(Mesh& m, getNodeType& n, Eigen::ArrayXd& errorAngle, Eigen::ArrayXXd& d, Eigen::ArrayXi& movingIndices, Eigen::ArrayXXd& exactDisp, Eigen::VectorXd& pnVec, projection* projPtr){

	int idx_m, idx_se, idx_ss, i;
	int edge;
	// for all of the boundary nodes the error will be determined
	for(i = 0; i< n.N_i; i++){

		// finding index of the node in consideration among the nodes with predescribed displacement
		idx_m = std::distance(std::begin(movingIndices), std::find(std::begin(movingIndices),std::end(movingIndices),(*n.iPtr)(i)));

		// finding index of the node in consideration among the sliding edge nodes
		idx_se = std::distance(std::begin(m.seNodes), std::find(std::begin(m.seNodes),std::end(m.seNodes),(*n.iPtr)(i)));

		// finding index of the node in consideration among the sliding edge nodes
		idx_ss = std::distance(std::begin(m.ssNodes), std::find(std::begin(m.ssNodes),std::end(m.ssNodes),(*n.iPtr)(i)));

		// if the node is part of the nodes with a prescribed displacement
		if(idx_m != movingIndices.size()){

//			std::cout <<"moving node: " <<  (*n.iPtr)(i) << std::endl;
			// the error is equal to the difference between the displacement and the prescribed (exact) displacement
			error.row(i) = d.row(i) - exactDisp.row(idx_m);
//			std::cout << d.row(i) << std::endl;
//			std::cout << exactDisp.row(idx_m) << std::endl;
//			std::cout << error.row(i) << std::endl;
		// if the node is part of the sliding edge nodes
		}else if(idx_se != m.N_se){
//			std::cout << "sliding edge node" << std::endl;
			if(m.nDims == 3){
				edge = 1;
			}else{
				edge = 0;
			}
			// projection is performed in which the projection itself is the error indication

			project(m, (*n.iPtr)(i), i, d, pnVec, projPtr, edge);
		}else if(idx_ss != m.N_ss){
//			std::cout << "sliding surf node" << std::endl;
			edge = 0;
			project(m, (*n.iPtr)(i), i, d, pnVec, projPtr, edge);
		}

		// if not the two above then the node is a static node with zero displacement. Therefore, its error is equal to the found displacement
		else{

			error.row(i) = d.row(i);
//			std::cout << (*n.iPtr)(i) << '\t' << error.row(i) << std::endl;
		}
		// finding the direction of the error.
		errorAngle(i) = atan2(error(i,1),error(i,0));
	}

//	std::cout << d << std::endl;
//	std::cout << error << std::endl;
//	std::cout << *n.iPtr << std::endl;



}

void greedy::getErrorMultiLvl(getNodeType& n, Eigen::ArrayXd& errorAngle,  Eigen::ArrayXXd& d, Mesh& m, Eigen::ArrayXi& movingIndices, Eigen::VectorXd& pnVec, projection* projPtr){
//	std::cout << "getting multi level error" << std::endl;
//	std::cout << displacement << std::endl;
//	std::cout << "\n" << displacement.rows() << '\t' << n.N_i << std::endl;

//	std::cout << "found delta: \n" << d << std::endl;

//	std::cout << displacement.row(111) - d.row(111) << '\n' << displacement.row(72) - d.row(72) << std::endl;


//	Eigen::ArrayXd errorAngle(n.N_i);
//	m.coords((*n.iPtr)(Eigen::seq(m.N_m,Eigen::last)), Eigen::all) += errorPrevLvl(Eigen::seq(m.N_m,Eigen::last),Eigen::all);
//	m.getMidPnts();
//	m.getVecs();
//
//
//	int idx_m, idx_se, i;
//
//	// for all of the boundary nodes the error will be determined
//	for(i = 0; i< n.N_i; i++){
//
//		// finding index of the node in consideration among the nodes with predescribed displacement
//		idx_m = std::distance(std::begin(movingIndices), std::find(std::begin(movingIndices),std::end(movingIndices),(*n.iPtr)(i)));
//
//		// finding index of the node in consideration among the sliding edge nodes
//		idx_se = std::distance(std::begin(m.seNodes), std::find(std::begin(m.seNodes),std::end(m.seNodes),(*n.iPtr)(i)));
//
//		// if the node is part of the nodes with a prescribed displacement
//		if(idx_m != movingIndices.size()){
//
//			// the error is equal to the difference between the displacement and the prescribed (exact) displacement
//			error.row(i) = d.row(i) - errorPrevLvl.row(i);
//
//
//		// if the node is part of the sliding edge nodes
//		}else if(idx_se != m.N_se){
//
//			// projection is performed in which the projection itself is the error indication
//			project(m, (*n.iPtr)(i),i, d, pnVec, projPtr);
//		}
//
//		// if not the two above then the node is a static node with zero displacement. Therefore, its error is equal to the found displacement
//		else{
//			error.row(i) = d.row(i);
//		}
//		// finding the direction of the error.
//		errorAngle(i) = atan2(error(i,1),error(i,0));
//	}
//	m.coords((*n.iPtr)(Eigen::seq(m.N_m,Eigen::last)), Eigen::all) -= errorPrevLvl(Eigen::seq(m.N_m,Eigen::last),Eigen::all);
	// differenece between actual displacement (displacement) and the found displacement (d)

//	m.coords(*n.iPtr,Eigen::all) += (d-error);
//	m.writeMeshFile();

	error = d - errorPrevLvl;
//	std::cout << d.row(67) << std::endl;
//	std::cout << errorPrevLvl.row(67) << std::endl;
//	std::cout << error.row(67) << std::endl;
//
//	m.writeMeshFile();

//	std::cout << "errors: \n\n" << std::endl;
//	std::cout << errorss << std::endl;

	// finding angle of the error directions
	for(int i = 0; i < n.N_i; i++){
		errorAngle(i) = atan2(error(i,1),error(i,0));
	}

	/*//identifying max error index
	int idxMax;
	error.rowwise().squaredNorm().maxCoeff(&idxMax);

	// value of max error
	e = error.row(idxMax).matrix().norm();
//	std::cout << errs.rowwise().norm() << std::endl;

//	std::cout << e << '\t' << idxMax << std::endl;

	// obtaining the double edge greedy maximum error

	int doubleEdgeMax = getDoubleEdgeError(errorAngle, idxMax ,  n.N_i);

	maxErrorNodes.resize(2);
	maxErrorNodes <<(*n.iPtr)(idxMax), (*n.iPtr)(doubleEdgeMax);*/
//	maxErrorNodes.resize(1);
//	maxErrorNodes << (*n.iPtr)(idxMax);
//	std::cout << maxErrorNodes << std::endl;



}

int greedy::getDoubleEdgeError(Eigen::ArrayXd& errorAngle, int idxMax, int N_i){
	// obtaining the error direction w.r.t. the direction of the max magnitude error
	errorAngle -= errorAngle(idxMax);



	// array containing the errors of the nodes that have an angle > 90 deg w.r.t. the max magnitude error direction
	Eigen::ArrayXd largeAngleError(N_i-1);
	// corresponding indices
	Eigen::ArrayXi largeAngleIdx(N_i-1);
	// Size of these arrays is unknown beforehand, but has a maximum of N_i-1 elements

	// counting the number of error with an relative angle > 90 def
	int cnt = 0;
	// reference angle of 90 deg.
	double refAngle = M_PI/2;
	// second reference angle of 360-90 = 270 deg.
	double refAngle2 = 3*refAngle;

	// for each boundary node
	for(int i=0; i<N_i; i++){
		// if the relative angle is between 90 and 270 degrees the squared norm is included in the array and its index is saved.
		if(abs(errorAngle(i)) >= refAngle && abs(errorAngle(i)) <= refAngle2){
			largeAngleError(cnt) = error.row(i).matrix().squaredNorm();
			largeAngleIdx(cnt) = i;
			cnt++;
		}
	}

	// resizing the largeAngleError array to only include the values included in the for-loop
	largeAngleError.conservativeResize(cnt);

	// finding the index where the errors is maximum
	int idxMaxLargeAngle;
	largeAngleError.maxCoeff(&idxMaxLargeAngle);

	// setting the integer equal to the corresponding index of all boundary nodes considered
	idxMaxLargeAngle = largeAngleIdx(idxMaxLargeAngle);

	return idxMaxLargeAngle;

}

void greedy::correction(Mesh& m, getNodeType& n, double& gamma){

	// the correction of the boundary nodes is equal to the error found previously.
	// for the control nodes the error is zero, so their position remains unchanged
	m.coords(*n.iPtr , Eigen::all) -= error;


	// The following part ensures the interpolation of the correction done at the boundary nodes into the internal volume
	// This is done by applying a local rbf interpolation of the displacement (correction) at the boundary nodes.
	// The correction radius is determined by the the multiplication of the maximum error with a factor gamma.
	// the factor gamma is set in the configuration file.

	// integer for storing the index where the error is largest
	int idxMax;

	// finding the node with the maximum error
	error.rowwise().norm().maxCoeff(&idxMax);

	// returning the largest error
	double maxError = error.row(idxMax).matrix().norm();

	// dist will store the distance between an internal node and the nearest boundary node
	double dist;

	// integer that stores the index of the nearest boundary node
	int idxNear;

	// looping through the internal nodes that were selected previously for applying the correction
//	for (int i = 0; i < m.intCorNodes.size(); i++){
	for (int i = 0; i < m.iNodes.size(); i++){

		// finding the nearest boundary node
//		getNearestNode(m, n, m.intCorNodes(i), idxNear, dist);
		getNearestNode(m, n, m.iNodes(i), idxNear, dist);

		// applying the interpolation
//		m.coords.row(m.intCorNodes(i)) -= error.row(idxNear)*rbfEval(dist,gamma*maxError);
		m.coords.row(m.iNodes(i)) -= error.row(idxNear)*rbfEval(dist,gamma*maxError);

		// keeping track of the progress
		std::cout << i << '\t' << m.iNodes.size() << std::endl;
	}
}

void greedy::getNearestNode(Mesh& m, getNodeType& n, int& node, int& idxMin, double& dist){

	// Array with the distance from the selected internal node to all the boundary nodes
	Eigen::ArrayXXd diff;
	diff = m.coords(*n.iPtr,Eigen::all);
	diff = diff.rowwise() - m.coords.row(node);


	// identifying the index to tho node that is nearest.
	diff.rowwise().squaredNorm().minCoeff(&idxMin);

	// dist from the internal node to the nearest boundary node
	dist = diff.row(idxMin).matrix().norm();
}

double greedy::rbfEval(double distance, double radius){
	double xi = distance/radius;	// distance scaled by support radius
	double f_xi;
	if(xi > 1){
		f_xi = 0;
	}else{
		f_xi = pow((1-(distance/radius)),4)*(4*(distance/radius)+1);
	}
	return f_xi;
}

void greedy::project(Mesh& m, int& node, int& idx, Eigen::ArrayXXd& disp,Eigen::VectorXd& pnVec, projection* projPtr, int edge){
	if(std::find(std::begin(m.staticNodes),std::end(m.staticNodes),node) != std::end(m.staticNodes)){
		std::cout << "static node: " << node << std::endl;
		std::cout << "displacement: \n" << disp.row(idx) << std::endl;
		std::cout << "PERIODIC VECTOR: \n" << pnVec << std::endl;

		error.row(idx) = disp.row(idx)*pnVec.transpose().array();
//		std::cout << error.row(idx) << std::endl;

	}else{

		Eigen::RowVectorXd err;
		Eigen::ArrayXXd dist;
		if(edge){
			dist = m.edgeMidPnts.rowwise() - (m.coords.row(node) + disp.row(idx));
		}else{
			dist = m.midPnts.rowwise() - (m.coords.row(node) + disp.row(idx));
		}


		projPtr->projectFun(m, err, dist,edge);

		error.row(idx) = -err;
	}
}


void greedy::setLevelParams(Mesh& m, getNodeType& n, int& lvl, int& lvlSize, Eigen::ArrayXXd& d, Eigen::VectorXd& alpha, Eigen::MatrixXd& Phi_imGreedy){
//	mNodesHist.conservativeResize(lvlSize,lvl+1);
//	mNodesHist.col(lvl) = *n.mPtr;

//	alphaHist.conservativeResize(lvlSize*m.nDims, lvl+1);
//	alphaHist.col(lvl) = alpha;



	if(lvl == 0){
		deltaInternal = Eigen::ArrayXXd::Zero(m.N_i, m.nDims);
		delta = Eigen::ArrayXXd::Zero(n.N_i, m.nDims);
	}

	delta += d;

	for(int dim = 0; dim < m.nDims; dim++){
		deltaInternal.col(dim) += (Phi_imGreedy*alpha(Eigen::seqN(dim*lvlSize,lvlSize))).array();
	}

	errorPrevLvl = -error;


}


