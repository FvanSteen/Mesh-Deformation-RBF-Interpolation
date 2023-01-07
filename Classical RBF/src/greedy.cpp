#include "greedy.h"
#include "projection.h"
#include <iostream>
#include <Eigen/Dense>
#include <iterator>
#include <chrono>
#include <math.h>
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


void greedy::getError(getNodeType& n, Mesh& meshOb, Eigen::ArrayXXd& d, double& e, Eigen::ArrayXi& maxErrorNodes, std::string sMode, Eigen::ArrayXi& mIndex, Eigen::ArrayXXd& displacement,Eigen::VectorXd& pnVec){
//	std::cout << d << std::endl;
//	std::cout << *n.iPtr << std::endl;
//	std::cout << std::endl;
//	std::cout << displacement << std::endl;

	error = Eigen::ArrayXXd::Zero(n.N_i,meshOb.nDims);
	Eigen::ArrayXd errorAngle(n.N_i);


	std::cout << "obtaining error" << std::endl;
	int idx_m, idx_se, i, dim, idxMax;
//	Eigen::ArrayXd dispRow;
	for(i = 0; i< n.N_i; i++){
		idx_m = std::distance(std::begin(mIndex), std::find(std::begin(mIndex),std::end(mIndex),(*n.iPtr)(i)));
		idx_se = std::distance(std::begin(meshOb.seNodes), std::find(std::begin(meshOb.seNodes),std::end(meshOb.seNodes),(*n.iPtr)(i)));
		if(idx_m!=mIndex.size()){

 			for(dim = 0;dim<meshOb.nDims;dim++){
				error(i,dim) = d(i,dim)-displacement(idx_m,dim);
			}

// 			std::cout << (*n.iPtr)(i) << " is an moving node" << std::endl;

		}else if(idx_se != meshOb.N_se){
//			std::cout << (*n.iPtr)(i) << " is an sliding edge node" << std::endl;
//			Eigen::ArrayXi& sNodes, Eigen::ArrayXXd& delta,Eigen::ArrayXXd& finalDef, Eigen::VectorXd& pVec

//			possibly not required due to the correction step!!
//			dispRow = d.row(i);
			project(meshOb, (*n.iPtr)(i),i, d, error,pnVec);


//			std::cout << (*n.iPtr)(i) << " is a sliding edge node, do something" << std::endl;
//			std::cout << "index in the seNodes array: " << idx_se << std::endl;
//			std::abs(meshOb.n.row(idx).matrix().dot(d.row(meshOb.N_i+meshOb.N_p + i + j + k).matrix()))
//			std::cout << meshOb.n.row(idx_se) << std::endl;
//			std::cout << d.row(i) << std::endl;
//			std::cout << meshOb.n.row(idx_se)*d.row(i) << std::endl;
//			error.row(i) = meshOb.n.row(idx_se)*d.row(i);

//			std::exit(0);


//			error.row(i) = d.row(i)*meshOb.n.row(idx_se);
//			std::cout << "\n\n" << "node: " << (*n.iPtr)(i) << "\n displacement\n"  << d.row(i) << "\n normal vec\n" << meshOb.n.row(idx_se) << "\nResulting error: \n" << error.row(i) << std::endl;
//			error.row(i) =  meshOb.n.row(idx_se).matrix().dot(d.row(i).matrix())*meshOb.n.row(idx_se);



		}
		else{
//			std::cout << (*n.iPtr)(i) << " should be a static node" << std::endl;
			for(dim = 0;dim<meshOb.nDims;dim++){
				error(i,dim) = d(i,dim);
			}
		}
		errorAngle(i) = atan2(error(i,1),error(i,0));
	}
//	std::cout << *n.iPtr << std::endl;
//	std::cout << "\n\n" <<  error << std::endl;


	error.rowwise().squaredNorm().maxCoeff(&idxMax);


	e = error.row(idxMax).matrix().norm();

	int doubleEdgeMax;
	getDoubleEdgeError(error, errorAngle, idxMax , n.N_i, doubleEdgeMax);

	maxErrorNodes.resize(2);
	maxErrorNodes << (*n.iPtr)(idxMax), (*n.iPtr)(doubleEdgeMax);
//	maxErrorNodes.resize(1);
//	maxErrorNodes << (*n.iPtr)(idxMax);
	std::cout << maxErrorNodes << std::endl;

//	std::exit(0);
//	std::cout << "Max error: " << e << std::endl;

//	meshOb.coords(*n.iPtr, Eigen::all) +=d;

//	std::cout <<  (error.col(1)/error.col(0)).atan()*180/M_PI << std::endl;








//	meshOb.writeMeshFile();


//	std::cout << error.rowwise().norm() << std::endl;
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

void greedy::getDoubleEdgeError(Eigen::ArrayXXd& error, Eigen::ArrayXd& errorAngle, int& idxMax, int& N_i, int& doubleEdgeMax){
//	std::cout << errorAngle*180/M_PI << std::endl;
//	std::cout << "error angle of max error: "<< errorAngle(idxMax)*180/M_PI << std::endl;

	// relative angle from max error vector
//	std::cout << "angle from max error vector: \n" << (errorAngle-errorAngle(idxMax))*180/M_PI << std::endl;

	errorAngle -= errorAngle(idxMax);


	Eigen::ArrayXd largeAngleError(N_i);
	Eigen::ArrayXi largeAngleIdx(N_i);


	int cnt = 0;
	double refAngle = M_PI/2;
	double refAngle2 = 3*refAngle;
	for(int i=0; i<N_i; i++){

		if(abs(errorAngle(i)) > refAngle && abs(errorAngle(i)) < refAngle2){
			largeAngleError(cnt) = error.row(i).matrix().squaredNorm();
			largeAngleIdx(cnt) = i;
			cnt++;
		}
	}
	largeAngleError.conservativeResize(cnt);
	largeAngleError.maxCoeff(&doubleEdgeMax);

	doubleEdgeMax = largeAngleIdx(doubleEdgeMax);
//	std::cout << doubleEdgeMax << '\t' << largeAngleIdx(doubleEdgeMax) << std::endl;

}

void greedy::correction(Mesh& m, getNodeType& n, double& gamma){


//	std::cout << "All the bdry points: " << std::endl;
//	std::cout << *n.iPtr << std::endl;
//	std::cout << "and their respective error: \n" << error << std::endl;
//	std::cout << m.intCorNodes << std::endl;
//	std::cout << m.intCorNodes.size() << std::endl;
//	std::exit(0);
	// doing the correction for the boundary nodes

	m.coords(*n.iPtr , Eigen::all) -= error;

//	m.writeMeshFile();
//	std::exit(0);



	// obtaining the max error
	double maxError = getMaxError();
//	std::cout << maxError << std::endl;
	// dist will store the distance between an internal node and the nearest boundary node
	double dist;
	// factor that determines the support radius
//	double gamma = 5;

	// integer that stores the index of the nearest boundary node
	int idxNear;

	// looping through the internal nodes
//	for (int i = 0; i < n.N_i_grdy; i++){
	for (int i = 0; i < m.intCorNodes.size(); i++){

		// obtaining the index of the nearest node and the distance to it
//		getNearestNode(m, n, (*n.iPtrGrdy)(i), idxNear, dist);
		getNearestNode(m, n, m.intCorNodes(i), idxNear, dist);
		std::cout << "node: " << m.intCorNodes(i) << "bdryNode "<< (*n.iPtr)(idxNear) << std::endl;
		std::cout << "dist: " << dist << std::endl;
//		std::cout << (*n.iPtrGrdy)(i) << '\t' << (*n.iPtr)(idxNear) << std::endl;
		// doing the correction
//		m.coords.row((*n.iPtrGrdy)(i)) -= error.row(idxNear)*rbfEval(dist,gamma*maxError);
		m.coords.row(m.intCorNodes(i)) -= error.row(idxNear)*rbfEval(dist,gamma*maxError);
//		std::cout << error.row(idxNear)*rbfEval(dist,gamma*maxError) << std::endl;
//		std::cout << i << '\t' << n.N_i_grdy << std::endl;
		std::cout << i << '\t' << m.intCorNodes.size() << std::endl;
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

void greedy::project(Mesh& m, int& node,int& idx, Eigen::ArrayXXd& disp, Eigen::ArrayXXd& error,Eigen::VectorXd& pnVec ){
//	std::cout << node << std::endl;

//	std::cout << "node in question: " << node << std::endl;
//	std::cout << "disp: \n" << disp.row(idx) << std::endl;


	if(std::find(std::begin(m.staticNodes),std::end(m.staticNodes),node) != std::end(m.staticNodes)){
		std::cout << "static node: " << node << std::endl;
		std::cout << "displacement: \n" << disp.row(idx) << std::endl;
		std::cout << "PERIODIC VECTOR: \n" << pnVec << std::endl;

		error.row(idx) = disp.row(idx)*pnVec.transpose().array();
		std::cout << error.row(idx) << std::endl;

//		error.row(idx) = disp.row(i)-disp.row(i)*pVec.transpose().array();
	}else{
		Eigen::RowVectorXd d;
		Eigen::ArrayXXd dist;
		int idxMin;
		int idxMinNew = -1;
		// distance to all midpoints
		dist = m.midPnts.rowwise()- (m.coords.row(node) + disp.row(idx));
		// find idx of
		dist.rowwise().norm().minCoeff(&idxMin);

//		std::cout << "initial position: " << m.coords.row(node) << std::endl;
//		std::cout << "new position: " << m.coords.row(node) + disp.row(idx) << std::endl;
//		std::cout << "dist to nearest midpoint: " << dist.row(idxMin) << std::endl;
//		std::cout << "Normal: " << m.midPntNormals.row(idxMin) << std::endl;

//		std::cout << "nearest midpoint: " << m.coords.row(node) + disp.row(idx) + dist.row(idxMin) << std::endl;
		error.row(idx) = -dist.row(idxMin).matrix().dot(m.midPntNormals.row(idxMin).matrix())*m.midPntNormals.row(idxMin);
//		std::cout << error.row(idx) << std::endl;
//		std::cout << "new position: " << m.coords.row(node) + disp.row(idx) - error.row(idx) << std::endl;

		while(idxMin != idxMinNew){
			dist = m.midPnts.rowwise() -(m.coords.row(node) + disp.row(idx) - error.row(idx));
			dist.rowwise().norm().minCoeff(&idxMinNew);



//			std::cout << idxMin << '\t' << idxMinNew << '\n';

			if(idxMin != idxMinNew){

				error.row(idx) -= dist.row(idxMinNew).matrix().dot(m.midPntNormals.row(idxMinNew).matrix())*m.midPntNormals.row(idxMinNew);


				idxMin = idxMinNew;
				idxMinNew = -1;
			}

		}

//		std::cout << error.row(idx) << std::endl;

	}


}


