#include "greedy.h"
#include "projection.h"
#include <iostream>
#include <Eigen/Dense>
#include <iterator>
#include <chrono>
#include <math.h>
greedy::greedy(probParams& params, Eigen::VectorXd& alpha, Eigen::ArrayXXd& d)
{
	if(params.multiLvl){
		alpha_step = &alphaGrdy;
		d_step = &delta;
		ctrlPtr = &ctrlNodesAll;
	}else{
		alpha_step = &alpha;
		d_step = &d;
	}
}


void greedy::getError(Mesh& m, getNodeType& n, Eigen::ArrayXXd& d, double& maxError, Eigen::ArrayXi& maxErrorNodes, Eigen::ArrayXi& movingIndices, Eigen::ArrayXXd& exactDisp, Eigen::VectorXd& pnVec, projection& p, bool multiLvl, int lvl, bool doubleEdge){


	error.resize(n.N_i,m.nDims);

	// defining array for the error directions
	Eigen::ArrayXd errorAngle(n.N_i);



	if(multiLvl && lvl>0){
		getErrorMultiLvl(n,errorAngle,d,m,movingIndices,pnVec,p);

	}else{
		getErrorSingleLvl(m,n,errorAngle,d,movingIndices,exactDisp,pnVec, p);
	}

	// find index of largest error
	int idxMax;
	error.rowwise().squaredNorm().maxCoeff(&idxMax);

	// and the maximum error magnitude
	maxError = error.row(idxMax).matrix().norm();


	if(doubleEdge){
		int idxMaxAngle = getDoubleEdgeError(errorAngle, idxMax, n.N_i, error);
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

void greedy::getErrorSingleLvl(Mesh& m, getNodeType& n, Eigen::ArrayXd& errorAngle, Eigen::ArrayXXd& d, Eigen::ArrayXi& movingIndices, Eigen::ArrayXXd& exactDisp, Eigen::VectorXd& pnVec, projection& p){

	int idx_m, idx_se, idx_ss, i;
	int edge;
	// for all of the boundary nodes the error will be determined
	for(i = 0; i< n.N_i; i++){
//		std::cout << i << "/" << n.N_i << std::endl;
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
//			std::cout << "sliding edge node: " <<  (*n.iPtr)(i) << std::endl;
			if(m.nDims == 3){
				edge = 1;
			}else{
				edge = 0;
			}
			// projection is performed in which the projection itself is the error indication

			project(m, (*n.iPtr)(i), i, d, pnVec, p, edge);
		}else if(idx_ss != m.N_ss){
//			std::cout << "sliding surf node: " << (*n.iPtr)(i) <<  std::endl;
			edge = 0;
			project(m, (*n.iPtr)(i), i, d, pnVec, p, edge);

		}

		// if not the two above then the node is a static node with zero displacement. Therefore, its error is equal to the found displacement
		else{
//			std::cout << "static: " <<  (*n.iPtr)(i) << std::endl;
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

void greedy::getErrorMultiLvl(getNodeType& n, Eigen::ArrayXd& errorAngle,  Eigen::ArrayXXd& d, Mesh& m, Eigen::ArrayXi& movingIndices, Eigen::VectorXd& pnVec, projection& p){
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


	error = d + errorPrevLvl(n.iNodesIdx, Eigen::all);

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

int greedy::getDoubleEdgeError(Eigen::ArrayXd& errorAngle, int idxMax, int N_i, Eigen::ArrayXXd& error){

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


	// in case there is no error with a relative angle larger than 90 deg.
	if(cnt == 0){
		refAngle -= M_PI/4;
		refAngle2 += M_PI/4;
		for(int i=0; i<N_i; i++){
			if(abs(errorAngle(i)) >= refAngle && abs(errorAngle(i)) <= refAngle2){
				largeAngleError(cnt) = abs(error.row(i).matrix().squaredNorm());
				largeAngleIdx(cnt) = i;
				cnt++;
			}
		}
	}
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

	// the correction of the boundary nodes is equal to the error found previously.
	// for the control nodes the error is zero, so their position remains unchanged

	Eigen::ArrayXXd* errorPtr;
	if(multiLvl){
		errorPtr = &errorPrevLvl;
	}else{
		errorPtr = &error;
	}
	m.coords(*n.iPtr , Eigen::all) -= *errorPtr;

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
//		std::cout << i << '\t' << m.N_i << std::endl;
		// finding the nearest boundary node
//		getNearestNode(m, n, m.intCorNodes(i), idxNear, dist);
		getNearestNode(m, n, m.iNodes(i), idxNear, dist);

		// applying the interpolation
//		m.coords.row(m.intCorNodes(i)) -= error.row(idxNear)*rbfEval(dist,gamma*maxError);
		m.coords.row(m.iNodes(i)) -= (*errorPtr).row(idxNear)*rbfEval(dist,gamma*maxError);

		// keeping track of the progress

		if(i % 5000 == 0){
			std::cout << i << '\t' << m.iNodes.size() << std::endl;
		}
//		std::cout << i << '\t' << m.intCorNodes.size() << std::endl;
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

void greedy::project(Mesh& m, int& node, int& idx, Eigen::ArrayXXd& disp,Eigen::VectorXd& pnVec, projection& p, int edge){
	if(std::find(std::begin(m.periodicVerticesNodes),std::end(m.periodicVerticesNodes),node) != std::end(m.periodicVerticesNodes)){
//		std::cout << "static node: " << node << std::endl;
//		std::cout << "displacement: \n" << disp.row(idx) << std::endl;
//		std::cout << "PERIODIC VECTOR: \n" << pnVec << std::endl;

		error.row(idx) = disp.row(idx)*pnVec.transpose().array();
//		std::cout << error.row(idx) << std::endl;

	}else{
		int idxNode = std::distance(std::begin(m.periodicEdgeNodes), std::find(std::begin(m.periodicEdgeNodes),std::end(m.periodicEdgeNodes),node));
		if(idxNode != m.N_pe){
//			std::cout << idxNode << '\t' << m.N_pe << std::endl;
			error.row(idx)  = disp.row(idx)*m.n_ss.row(m.N_ss-m.N_pe+idxNode);
		}else{

			Eigen::RowVectorXd err;
			Eigen::ArrayXXd dist;
			if(edge){
				dist = m.edgeMidPnts.rowwise() - (m.coords.row(node) + disp.row(idx));
			}else{
				dist = m.midPnts.rowwise() - (m.coords.row(node) + disp.row(idx));
			}


			p.projectFun(m, err, dist,edge);

			error.row(idx) = -err;
		}
	}
}


void greedy::setLevelParams(Mesh& m, getNodeType& n, int lvl, Eigen::ArrayXXd& d, Eigen::VectorXd& alpha, double maxError, Eigen::VectorXd& defVec, Eigen::ArrayXi* cPtr, int N_c){

//	std::cout << n.N_i << std::endl;
//	std::cout << n.N_c << std::endl;
//	std::cout << "errorRows\t" << error.rows() ;
//

//	mNodesHist.conservativeResize(lvlSize,lvl+1);
//	mNodesHist.col(lvl) = *n.mPtr;

//	alphaHist.conservativeResize(lvlSize*m.nDims, lvl+1);
//	alphaHist.col(lvl) = alpha;
//	std::cout << *n.mPtr << std::endl;
//
//	std::cout << alpha << std::endl;
//	std::cout << *n.mPtr << std::endl;
//
//	Eigen::ArrayXi test(m.N_m+m.N_se+m.N_ss);
//	test << m.mNodes, m.seNodes,m.ssNodes;

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


void greedy::setInitMaxErrorNodes(Mesh& m, Eigen::ArrayXXd& coords, Eigen::ArrayXXd& disp, Eigen::ArrayXi& mIdx, Eigen::ArrayXi& maxErrorNodes, bool doubleEdge){

	Eigen::Array2i idxMax;
	disp.rowwise().norm().maxCoeff(&idxMax(0));
	std::cout << "max deformation: " << disp.row(idxMax(0)).matrix().norm() << std::endl;
	maxErrorPrevLvl = disp.row(idxMax(0)).matrix().norm();

	if(doubleEdge){
		int N = mIdx.size();
		Eigen::ArrayXd errorAngle(N);
		for(int i = 0; i<N; i++){
			errorAngle(i) = atan2(disp(i,1), disp(i,0));
		}



		idxMax(1) = getDoubleEdgeError(errorAngle, idxMax(0), N, disp);

		if (idxMax(1) == -1){
			Eigen::ArrayXXd movingCoords;
			movingCoords = coords(mIdx,Eigen::all);

			Eigen::ArrayXd dist;
			dist = (movingCoords.rowwise() - coords.row(mIdx(idxMax(0)))).rowwise().norm();
			dist.maxCoeff(&idxMax(1));
		}

		maxErrorNodes.resize(2);
		maxErrorNodes << mIdx(idxMax);
	}else{
		maxErrorNodes.resize(1);
		maxErrorNodes << mIdx(idxMax(0));
	}


}


//void greedy::setMaxErrorNodes(Mesh& m, Eigen::ArrayXi& maxErrorNodes){
//
//	Eigen::Array2i idxMax;
//	Eigen::ArrayXXd mError;
//	mError = error(Eigen::seqN(0,m.N_m), Eigen::all);
//	mError.rowwise().squaredNorm().maxCoeff(&idxMax(0));
//
//
//
//	Eigen::ArrayXd errorAngle(m.N_m);
//	for(int i= 0; i <m.N_m; i++){
//		errorAngle(i) = atan2(mError(i,1),mError(i,0));
//	}
//
//
//	idxMax(1) = getDoubleEdgeError(errorAngle, idxMax(0), m.N_m);
//
////	std::cout << idxMax << std::endl;
////
////	std::exit(0);
////
////	Eigen::ArrayXXd movingCoords;
////	movingCoords = m.coords(m.mNodes,Eigen::all);
////
////	Eigen::ArrayXd dist;
////	dist = (movingCoords.rowwise()-m.coords.row(m.mNodes(idxMax(0)))).rowwise().norm();
////	dist.maxCoeff(&idxMax(1));
////
////
////
////	maxErrorNodes.resize(2);
////	maxErrorNodes << m.mNodes(idxMax);
//
//}
