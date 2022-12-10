#include "greedy.h"
#include <iostream>
#include <Eigen/Dense>
#include <iterator>
#include <chrono>
greedy::greedy()  {
//
////	std::cout << ptr->mNodes << std::endl;
////	std::cout << ptr->iNodes <<std::endl;
////	std::cout << ptr->sNodes(2) << std::endl;
//	std::cout << "In the greedy class" << std::endl;
////	int idx = 2;
//
////	std::cout << n.sNodes << std::endl;
//
//	n.GreedyInit();
//
//	n.greedyNodes(rbf.m.intBdryNodes(0));
//
//
////	rbf.params.steps = 1;
//
////	std::cout <<rbf.params.steps << std::endl;
//	rbf.m.getVecs();
//	rbf.perform_rbf(n);
////	std::cout << rbf.m.coords << std::endl;
//	int idxMax = getError(n,rbf);
//	std::cout << n.iNodes(idxMax) << std::endl;
//	n.greedyNodes(n.iNodes(idxMax));
//
//
////	std::cout << n.iNodes << std::endl;
////	rbf.params.dataRed = false;
//
//	rbf.perform_rbf(n);
//	idxMax = getError(n,rbf);
//	std::cout << n.iNodes(idxMax) << std::endl;
//
////	std::cout << n.iNodes << std::endl;
//
////	rbf.params.dataRed = false;
//	n.greedyNodes(n.iNodes(idxMax));
//	rbf.perform_rbf(n);
//	idxMax = getError(n,rbf);
//	std::cout << n.iNodes(idxMax) << std::endl;
//
////	std::cout << n.iNodes << std::endl;
//
////	rbf.params.dataRed = false;
//	n.greedyNodes(n.iNodes(idxMax));
//	rbf.perform_rbf(n);
//	idxMax = getError(n,rbf);
//	std::cout << n.iNodes(idxMax) << std::endl;
//
////	rbf.params.dataRed = false;
//	n.greedyNodes(n.iNodes(idxMax));
//	rbf.perform_rbf(n);
//	idxMax = getError(n,rbf);
//	std::cout << n.iNodes(idxMax) << std::endl;
//
//
////	rbf.params.dataRed = false;
//	n.greedyNodes(n.iNodes(idxMax));
//	rbf.perform_rbf(n);
//	idxMax = getError(n,rbf);
//	std::cout << n.iNodes(idxMax) << std::endl;
//
//
////	rbf.params.dataRed = false;
//	n.greedyNodes(n.iNodes(idxMax));
//	rbf.perform_rbf(n);
//	idxMax = getError(n,rbf);
//	std::cout << n.iNodes(idxMax) << std::endl;
//
//
//
//	rbf.params.dataRed = false;
//	n.greedyNodes(n.iNodes(idxMax));
//	rbf.perform_rbf(n);
//
//	std::cout << n.mNodesStd << std::endl;
//	idxMax = getError(n,rbf);
//	std::cout << n.iNodes(idxMax) << std::endl;

//	std::cout << n.iNodes << std::endl;

//	std::cout << rbf.d.rows() << std::endl;
//	Eigen::ArrayXi freeNodesIdx;
//	freeNodesIdx <<

//	ptr->mNodes.resize(ptr->m.N_ib + ptr->m.N_es + 1);
//	ptr->mNodes << ptr->m.intBdryNodes, ptr->m.extStaticNodes, ptr->m.slidingEdgeNodes(idx);
//
//	ptr->sNodes.resize(1);
//	ptr->sNodes << ptr->m.slidingEdgeNodes(idx);
//
//	ptr->iNodes.resize(ptr->m.N_se-1);
//	ptr->iNodes << ptr->m.slidingEdgeNodes(Eigen::seqN(0,idx)), ptr->m.slidingEdgeNodes(Eigen::seqN(idx+1,ptr->m.N_se-idx-1));
//	std::cout << ptr->iNodes << std::endl;
//	ptr->N_m = ptr->mNodes.size();
//	ptr->N_s = ptr->sNodes.size();
//	ptr->N_i = ptr->iNodes.size();
//	std::cout << ptr->mNodes << std::endl;
//	std::cout << ptr->sNodes << std::endl;
//
//	std::cout << '\n' << ptr->iNodes << std::endl;

//	ptr->iNodes.conservativeResize(iNodes.size()+)
//	ptr->mNodes.conservativeResize

//	Eigen::ArrayXXd coordsInit;
//	coordsInit = ptr->m.coords(ptr->iNodes,Eigen::all);
//	ptr->m.getVecs();
//
//	ptr->perform_rbf();
//
//	Eigen::ArrayXXd diff;
//	diff = coordsInit-ptr->m.coords(ptr->iNodes,Eigen::all);
//	std::cout << diff << std::endl;
//	getError(diff,ptr);

}


void greedy::getError(getNodeType& n, Mesh& meshOb, Eigen::ArrayXXd& d, double& e, int& idxMax, std::string sMode, Eigen::ArrayXi& mIndex, Eigen::ArrayXXd& displacement){
//	std::cout << d << std::endl;
//	std::cout << std::endl;
//	std::cout << displacement << std::endl;

	error = Eigen::ArrayXXd::Zero(n.N_i,meshOb.nDims);

	std::cout << "obtaining error" << std::endl;
	int idx, i, dim;

	for(i = 0; i< n.N_i; i++){
		idx = std::distance(std::begin(mIndex), std::find(std::begin(mIndex),std::end(mIndex),(*n.iPtr)(i)));
		if(idx!=mIndex.size()){
 			for(dim = 0;dim<meshOb.nDims;dim++){
				error(i,dim) = d(i,dim)-displacement(idx,dim);
			}
		}
		else{
			for(dim = 0;dim<meshOb.nDims;dim++){
				error(i,dim) = d(i,dim);
			}
		}
	}



	error.rowwise().norm().maxCoeff(&idxMax);

	Eigen::VectorXd errorVec =  error.row(idxMax);
	e = errorVec.norm();
	idxMax = (*n.iPtr)(idxMax);

//	Eigen::VectorXd defVec;
//	int N = rbf.m.N_ib-n.N_ib;
//	defVec = Eigen::VectorXd::Zero(rbf.m.nDims*N);
////	std::cout << n.iNodes(Eigen::seq(rbf.m.N_i+rbf.m.N_p, rbf.m.N_i+rbf.m.N_p+rbf.m.N_ib-n.N_ib-1)) << std::endl;
//
//
//	Eigen::ArrayXi ibNodes = n.iNodes(Eigen::seq(rbf.m.N_i+rbf.m.N_p, rbf.m.N_i+rbf.m.N_p+rbf.m.N_ib-n.N_ib-1));
//	rbf.getDefVec(defVec, N, ibNodes);

//	Eigen::ArrayXd error;
//	error = Eigen::ArrayXd::Zero(n.N_i - meshOb.N_i-meshOb.N_p);
//	int i,j,k;
//	for(i=0; i<meshOb.N_ib-n.N_ib;i++){
//
//		for(int dim=0; dim<meshOb.nDims;dim++){
////			std::cout << rbf.d(rbf.m.N_i+rbf.m.N_p + i,dim)-defVec(N*dim+i) << std::endl;
//			error(i) += pow(d(meshOb.N_i+meshOb.N_p + i,dim)-exactDeformation((meshOb.N_ib-n.N_ib)*dim+i),2);
//
//		}
//		error(i) = sqrt(error(i));
////		std::cout << n.iNodes(rbf.m.N_i+rbf.m.N_p + i) << std::endl;
//
//	}
////	std::cout << "check" << std::endl;
//	for(j = 0; j < meshOb.N_es - n.N_es ; j++){
////		std::cout << i + j << std::endl;
//		for(int dim=0; dim<meshOb.nDims;dim++){
//			error(i+j) += pow(d(meshOb.N_i+meshOb.N_p + i + j,dim),2);
//		}
//		error(i+j) = sqrt(error(i+j));
//	}
////	std::cout << "check" << std::endl;
//
////	std::cout << rbf.m.n << std::endl;
////	std::cout< < meshOb.N_se << std::endl;
////	std::cout << meshOb.N_es << std::endl;
//
//
//
//	// todo include as function argument
////	std::string smode = "none";
//
//	int idx;
//	if(sMode == "ps"){
//		for(k=0; k< meshOb.N_se - n.N_s; k++){
////		std::cout << n.iNodes(rbf.m.N_i+rbf.m.N_p+i+j+k) << std::endl;
//			idx = std::distance(std::begin(meshOb.slidingEdgeNodes),std::find(std::begin(meshOb.slidingEdgeNodes),std::end(meshOb.slidingEdgeNodes),n.iNodes(meshOb.N_i+meshOb.N_p+i+j+k)));
//
////		std::cout << idx << std::endl;
////		std::cout << i+j+k << std::endl;
////		std::cout << rbf.m.n.row(k) << std::endl;
////		std::cout << rbf.d.row(rbf.m.N_i+rbf.m.N_p + i + j+k) << std::endl;
////		std::cout << std::abs(rbf.m.n.row(k).matrix().dot(rbf.d.row(rbf.m.N_i+rbf.m.N_p + i + j+k).matrix())) << std::endl;
//			error(i+j+k) = std::abs(meshOb.n.row(idx).matrix().dot(d.row(meshOb.N_i+meshOb.N_p + i + j + k).matrix()));
//		}
//	}
////	std::cout << n.iNodes << std::endl;
////	std::cout << d << std::endl;
////	std::cout << '\n';
//
//
//	if(sMode == "none"){
//		for(k=0; k< meshOb.N_se - n.N_s; k++){
//			// for loop through dimensions
////			std::cout << k << std::endl;
////			std::cout << meshOb.slidingEdgeNodes(k) << std::endl;
////			std::cout << d.row(meshOb.N_i+meshOb.N_p + i + j+k) << std::endl; // this is the error, since the displacement should actually by zero for the none sliding mode.
//			for(int dim=0; dim<meshOb.nDims;dim++){
//				error(i+j+k) += pow(d(meshOb.N_i+meshOb.N_p + i + j+ k,dim),2);
//			}
//			error(i+j+k) = sqrt(error(i+j+k));
//		}
//	}
//
////	std::cout << "check" << std::endl;
//
////	std::cout << "error Array: \n" << error << '\n' << std::endl;
////	std::exit(0);
////	std::exit(0);
//
//	error.maxCoeff(&idxMax);
//
////	std::cout << n.iNodes << std::endl;
//
////	std::cout << error(idxMax) << std::endl;
////	std::cout << meshOb.N_i+meshOb.N_p+idxMax << std::endl;
//	e = error(idxMax);
////	std::cout << rbf.d << std::endl;
////	std::cout << rbf.d << std::endl;
////	std::cout << defVec << std::endl;
////	Eigen::ArrayXXd error;
////
////	std::cout << ptr->m.n << std::endl;
//
////	error = ptr->m.n
//
//	idxMax = n.iNodes(meshOb.N_i+meshOb.N_p+idxMax);
//

}

void greedy::correction(Mesh& m, getNodeType& n){

	// doing the correction for the boundary nodes
	m.coords(*n.iPtr , Eigen::all) -= error;

	// obtaining the max error
	double maxError = getMaxError();

	// dist will store the distance between an internal node and the nearest boundary node
	double dist;
	// factor that determines the support radius
	double gamma = 2.5;

	// integer that stores the index of the nearest boundary node
	int idxNear;

	// looping through the internal nodes
	for (int i = 0; i < n.N_i_grdy; i++){

		// obtaining the index of the nearest node and the distance to it
		getNearestNode(m, n, (*n.iPtrGrdy)(i), idxNear, dist);

		// doing the correction
		m.coords.row((*n.iPtrGrdy)(i)) += error.row(idxNear)*rbfEval(dist,gamma*maxError);
	}
}

void greedy::getNearestNode(Mesh& m, getNodeType& n, int& node, int& idxMin, double& dist){

	// Array with the distance from the selected internal node to all the boundary nodes
	Eigen::ArrayXXd diff;
	diff = m.coords(*n.iPtr, Eigen::all).rowwise() - m.coords.row(node);

	// identifying the index to tho node that is nearest.
	diff.rowwise().norm().minCoeff(&idxMin);

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

double greedy::getMaxError(){
	// integer for storing the index where the error is largest
	int idxMax;

	// finding the node with the highest error
	error.rowwise().norm().maxCoeff(&idxMax);

	// returning the largest error
	return error.row(idxMax).matrix().norm();
}
