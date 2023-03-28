#include "greedy.h"

#include <iostream>
#include <Eigen/Dense>
#include <iterator>
#include <chrono>
#include <math.h>
#include "SPDS.h"
#include "nanoflann.hpp"
#include "SPDS.h"
greedy::greedy(Mesh& m, probParams& params,Eigen::ArrayXXd& disp, Eigen::ArrayXi& movingIndices, Eigen::VectorXd& alpha, Eigen::ArrayXXd& d, Eigen::VectorXd& perVec)
{


	mPtr = &m;
	exctDispPtr = &disp;
	mIdxPtr = &movingIndices;

	p = &params;
	pVec = &perVec;
	if((*p).multiLvl){
		alpha_step = &alphaGrdy;
		d_step = &delta;
		ctrlPtr = &ctrlNodesAll;
	}else{
		alpha_step = &alpha;
		d_step = &d;
	}

	setInitMaxErrorNodes(m);
	std::cout << "initial selected Nodes:\n" << maxErrorNodes << std::endl;


}


void greedy::getError(Mesh& m, getNodeType& n, Eigen::ArrayXXd& d, int lvl){

	error.resize(n.N_i,m.nDims);
	errorAngle.resize(n.N_i, m.nDims-1);
	// defining array for the error directions
//	Eigen::ArrayXXd errorAngle(n.N_i, m.nDims-1);



	if((*p).multiLvl && lvl>0){
		getErrorMultiLvl(n,d);
	}else{
		getErrorSingleLvl(m,n,d);
	}


	// find index of largest error
	int idxMax;
	error.rowwise().squaredNorm().maxCoeff(&idxMax);

	// and the maximum error magnitude
	maxError = error.row(idxMax).matrix().norm();


	// from here
	if((*p).doubleEdge){
		getErrorAngle(m.nDims, n.N_i);
		int idxMaxAngle = getDoubleEdgeError(idxMax, n.N_i, error);

		if(idxMaxAngle == -1){
			maxErrorNodes.resize(1);
			maxErrorNodes << (*n.iPtr)(idxMax);
		}else{
			maxErrorNodes.resize(2);
			maxErrorNodes << (*n.iPtr)(idxMax), (*n.iPtr)(idxMaxAngle);
		}
	}else{
		maxErrorNodes.resize(1);
		maxErrorNodes << (*n.iPtr)(idxMax);
	}
}

void greedy::getErrorMovingNodes(Eigen::ArrayXi* nodes, Eigen::ArrayXXd& d, size_t N){
	int idx_m;
	for(size_t i = 0; i < N; i++){
		idx_m = std::distance(std::begin(*mIdxPtr), std::find(std::begin(*mIdxPtr),std::end(*mIdxPtr),(*nodes)(i)));
		if(idx_m !=  (*mIdxPtr).size()){
			error.row(i) = d.row(i) - (*exctDispPtr).row(idx_m);
		}else{
			error.row(i) = d.row(i);
		}
	}
}




void greedy::getErrorAngle(size_t dims, size_t N){
	for(size_t i =0 ; i < N; i++){
		errorAngle(i,0) = atan2(error(i,1),error(i,0));
		if(dims == 3){
			errorAngle(i,1) = atan2(  sqrt(pow(error(i,0),2) + pow(error(i,1),2)), error(i,2));
		}
	}
}
void greedy::getErrorSingleLvl(Mesh& m, getNodeType& n, Eigen::ArrayXXd& d){


	size_t m_end, se_end, ss_end;
	m_end = m.N_m-n.N_c;
	se_end = m_end + m.N_se-n.N_se;

	getErrorMovingNodes(n.iPtr,  d,  m_end);
//	std::cout << "obtained error moving nodes\n";
	SPDS SPDSobj;
	if(m_end != size_t(n.N_i)){
		SPDSobj.projectEdge(m, n.iPtr, d, error, *pVec, m_end, se_end, 0);
	}
//	std::cout << "obtained error sliding edge nodes\n";
	if(se_end != size_t(n.N_i)){
		ss_end = se_end + m.N_ss-n.N_ss;
		SPDSobj.projectSurf(m, n.iPtr, d, error, *pVec, se_end, ss_end, 0);
	}

}

void greedy::getErrorMultiLvl(getNodeType& n, Eigen::ArrayXXd& d){
	error = d + errorPrevLvl(n.iNodesIdx, Eigen::all);
}

int greedy::getDoubleEdgeError(int idxMax, int N_i, Eigen::ArrayXXd& error){


	Eigen::ArrayXd maxErrorAngle = errorAngle.row(idxMax);
	// obtaining the error direction w.r.t. the direction of the max magnitude error
	errorAngle.rowwise() -= maxErrorAngle.transpose();

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
		if (errorAngle.cols() == 1){
			if(abs(errorAngle(i,0)) >= refAngle && abs(errorAngle(i,0)) <= refAngle2){
				largeAngleError(cnt) = error.row(i).matrix().squaredNorm();
				largeAngleIdx(cnt) = i;
				cnt++;
			}
		}
		else{
			if(abs(errorAngle(i,0)) >= refAngle && abs(errorAngle(i,0)) <= refAngle2 && abs(errorAngle(i,1)) >= refAngle && abs(errorAngle(i,1)) <= refAngle2){
				largeAngleError(cnt) = error.row(i).matrix().squaredNorm();
				largeAngleIdx(cnt) = i;
				cnt++;
			}
		}
//		std::cout << i << "\t" << errorAngle.row(i) << '\t' << error.row(i) << std::endl;
	}


	if (cnt == 0 && errorAngle.cols() == 2){
		for(int i=0; i<N_i; i++){
			if(abs(errorAngle(i,0)) >= refAngle && abs(errorAngle(i,0)) <= refAngle2){
				largeAngleError(cnt) = error.row(i).matrix().squaredNorm();
				largeAngleIdx(cnt) = i;
				cnt++;
			}
		}
	}


//	 in case there is no error with a relative angle larger than 90 deg.
//	if(cnt == 0){
//		refAngle -= M_PI/4;
//		refAngle2 += M_PI/4;
//		for(int i=0; i<N_i; i++){
//			// if the relative angle is between 90 and 270 degrees the squared norm is included in the array and its index is saved.
//			if (errorAngle.cols() == 1){
//				if(abs(errorAngle(i,0)) >= refAngle && abs(errorAngle(i,0)) <= refAngle2){
//					largeAngleError(cnt) = error.row(i).matrix().squaredNorm();
//					largeAngleIdx(cnt) = i;
//					cnt++;
//				}
//			}else{
//				if(abs(errorAngle(i,0)) >= refAngle && abs(errorAngle(i,0)) <= refAngle2 && abs(errorAngle(i,1)) >= refAngle && abs(errorAngle(i,1)) <= refAngle2){
//					largeAngleError(cnt) = error.row(i).matrix().squaredNorm();
//					largeAngleIdx(cnt) = i;
//					cnt++;
//				}
//			}
//
//	//		std::cout << i << "\t" << errorAngle.row(i) << '\t' << error.row(i) << std::endl;
//		}
//
//	}


//	std::cout << "here" << std::endl;
//	std::cout << largeAngleError(Eigen::seqN(0,cnt)) << std::endl;
//	std::cout << largeAngleIdx(Eigen::seqN(0,cnt)) << std::endl;

	int idxMaxLargeAngle;

	if(cnt!= 0){
		// resizing the largeAngleError array to only include the values included in the for-loop
		largeAngleError.conservativeResize(cnt);

		// finding the index where the errors is maximum
		largeAngleError.maxCoeff(&idxMaxLargeAngle);

		// setting the integer equal to the corresponding index of all boundary nodes considered
		idxMaxLargeAngle = largeAngleIdx(idxMaxLargeAngle);
	}else{
		std::cout << "no error found with a large angle" << std::endl;
		idxMaxLargeAngle = -1;
	}

//	some if statement to ensure that the same node cannot be added as doubleEdged error node

	return idxMaxLargeAngle;


}

void greedy::correction(Mesh& m, getNodeType& n, double& gamma, bool& multiLvl){

	Eigen::ArrayXXd* errorPtr;
	if(multiLvl){
		errorPtr = &errorPrevLvl;
	}else{
		errorPtr = &error;
	}
	m.coords(*n.iPtr , Eigen::all) -= *errorPtr;


	// integer for storing the index where the error is largest
	int idxMax;

	// finding the node with the maximum error
	error.rowwise().norm().maxCoeff(&idxMax);

	// returning the largest error
	double maxError = error.row(idxMax).matrix().norm();


	SPDS SPDS_obj;
	SPDS_obj.kdt_NNSearch(*n.iPtr, *n.iPtrGrdy,  m.coords, m.nDims, gamma, maxError, errorPtr);

}





void greedy::setLevelParams(Mesh& m, getNodeType& n, int lvl, Eigen::ArrayXXd& d, Eigen::VectorXd& alpha,Eigen::VectorXd& defVec, Eigen::ArrayXi* cPtr, int N_c){


	maxErrorPrevLvl = maxError;

	if(lvl == 0){
//		deltaInternal = Eigen::ArrayXXd::Zero(m.N_i, m.nDims);
		delta = Eigen::ArrayXXd::Zero(m.N_m+m.N_se+m.N_ss, m.nDims);
		alphaTotal.resize(0,0);
		ctrlNodesAll.resize(0);
//		alphaSum = Eigen::VectorXd::Zero(m.nDims*n.N_i);
	}


	Eigen::ArrayXi newCtrlNodes(N_c);
	Eigen::ArrayXXd newAlpha(N_c, m.nDims);

	int cnt = 0;
	int idx, dim;
	for(int i = 0; i < N_c; i++){
//		std::cout << i << '\t' << N_c << std::endl;
		idx = std::distance(std::begin(ctrlNodesAll), std::find(std::begin(ctrlNodesAll), std::end(ctrlNodesAll), (*cPtr)(i)));
		if(idx == ctrlNodesAll.size()){
//			std::cout << i << '\t' << (*n.mPtr)(i) << std::endl;
			newCtrlNodes(cnt) = i;

			for(dim =0; dim < m.nDims; dim++){
				newAlpha(cnt, dim) = alpha(dim*N_c+i);
			}
			cnt++;
		}else{
			for(dim = 0; dim < m.nDims; dim++){
				alphaTotal(idx,dim) += alpha(dim*N_c+i);
			}
		}
	}


	ctrlNodesAll.conservativeResize(ctrlNodesAll.size()+cnt);
	alphaTotal.conservativeResize(alphaTotal.rows()+cnt,m.nDims);

	ctrlNodesAll(Eigen::lastN(cnt)) = (*cPtr)(newCtrlNodes(Eigen::seqN(0,cnt)));
	alphaTotal(Eigen::lastN(cnt), Eigen::all) = newAlpha(Eigen::seqN(0,cnt), Eigen::all);


	delta(n.iNodesIdx,Eigen::all) += d;

	for(int dim = 0; dim < m.nDims; dim++ ){
		delta(n.cNodesIdx,dim) += (defVec(Eigen::seqN(dim*N_c, N_c))).array();
	}

	errorPrevLvl = Eigen::ArrayXXd::Zero(m.N_m+m.N_se+m.N_ss, m.nDims);
	errorPrevLvl(n.iNodesIdx,Eigen::all) = error;
//	std::cout << test(n.cNodesIdx) << std::endl;
}

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

void greedy::getAlphaIdx(Eigen::ArrayXi& mNodes, Eigen::ArrayXi* mNodesGrdy, int N, Eigen::ArrayXi& idxAlpha){
	int idx;
	idxAlpha.resize(N);

	for(int i = 0; i < N; i++){
		idx = std::distance(std::begin(mNodes), std::find(std::begin(mNodes), std::end(mNodes),(*mNodesGrdy)(i)));
		idxAlpha(i) = idx;
	}
}


void greedy::setInitMaxErrorNodes(Mesh& m){

	Eigen::Array2i idxMax;
	(*exctDispPtr).rowwise().norm().maxCoeff(&idxMax(0));
	std::cout << "max deformation: " << (*exctDispPtr).row(idxMax(0)).matrix().norm() << std::endl;
	maxErrorPrevLvl = (*exctDispPtr).row(idxMax(0)).matrix().norm();

	if((*p).doubleEdge){
		int N = (*mIdxPtr).size();
		errorAngle.resize(N, m.nDims -1);
		for(int i = 0; i<N; i++){
			errorAngle(i,0) = atan2((*exctDispPtr)(i,1), (*exctDispPtr)(i,0));
			if(m.nDims ==3){
				errorAngle(i,1) = atan2(sqrt(pow((*exctDispPtr)(i,0),2) + pow((*exctDispPtr)(i,1),2)), (*exctDispPtr)(i,2));
			}
		}



		idxMax(1) = getDoubleEdgeError(idxMax(0), N, (*exctDispPtr));

		if (idxMax(1) == -1){
			Eigen::ArrayXXd movingCoords;
			movingCoords = m.coords((*mIdxPtr),Eigen::all);

			Eigen::ArrayXd dist;
			dist = (movingCoords.rowwise() - m.coords.row((*mIdxPtr)(idxMax(0)))).rowwise().norm();
			dist.maxCoeff(&idxMax(1));
		}
		maxErrorNodes.resize(2);
		maxErrorNodes << (*mIdxPtr)(idxMax);
	}else{
		maxErrorNodes.resize(1);
		maxErrorNodes << (*mIdxPtr)(idxMax(0));
	}


}



