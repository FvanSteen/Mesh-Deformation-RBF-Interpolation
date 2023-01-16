
#include "rbfGenFunc.h"
#include <iostream>
#include <Eigen/Dense>
#include <chrono>
#include <fstream>
#include <sstream>
rbfGenFunc::rbfGenFunc(Mesh& meshObject, struct probParams& probParamsObject)
//rbfGenFunc::rbfGenFunc(Mesh* meshPtr, Eigen::VectorXd& dVec, Eigen::RowVectorXd& rotPnt, Eigen::VectorXd& rotVec, const int& steps, const std::string& smode, const bool& curved, const std::string& pDir)
:m(meshObject), params(probParamsObject)
{

	std::cout << "Initialised the rbfGenFunc class" << std::endl;
	movingIndices.resize(m.ibIndices.size());
	exactDisp.resize(m.ibIndices.size(),m.nDims);
	readDisplacementFile();

	exactDisp = exactDisp/params.steps; // deformation per step is more usefull then the total deformation.
	getPeriodicParams();


}

void rbfGenFunc::getPhi(Eigen::MatrixXd& Phi, Eigen::ArrayXi& idxSet1, Eigen::ArrayXi& idxSet2){
//	double lambda = 0.05749995;
//	double lambda = 1;
//	double lambda = 0.045;
	Phi.resize(idxSet1.size(), idxSet2.size());
	double dist;
	for(int i=0; i<idxSet1.size();i++){
		for(int j=0; j<idxSet2.size();j++){
			if(m.nDims == 2){
				// Euclidian distance

//				dist = sqrt(pow(m.coords(idxSet1(i),0)-m.coords(idxSet2(j),0),2) + pow(m.coords(idxSet1(i),1)-m.coords(idxSet2(j),1),2));


				//todo following if statement is introduced to improve efficiency, check if anything else can be done
				if(m.pmode != "none"){
					dist = 0;
					for(int dim = 0; dim < m.nDims; dim++){
						if(pVec(dim)){
							//todo change function to allow for periodic length in stead of unit length
							dist += pow(m.lambda/M_PI*sin( (m.coords(idxSet1(i),dim)-m.coords(idxSet2(j),dim))*M_PI/m.lambda),2);
						}
						else{
							dist += pow(m.coords(idxSet1(i),dim)-m.coords(idxSet2(j),dim),2);
						}

					}
					dist = sqrt(dist);
				}
				else{
					//todo can the difference in coordinates by calculated for the whole row. then take the power of 2, sum and take square root??
					dist = sqrt(pow(m.coords(idxSet1(i),0)-m.coords(idxSet2(j),0),2) + pow(m.coords(idxSet1(i),1)-m.coords(idxSet2(j),1),2) );
				}
//				dist = sqrt(pow(m.coords(idxSet1(i),0)-m.coords(idxSet2(j),0),2) + pow(1/M_PI*sin( (m.coords(idxSet1(i),1)-m.coords(idxSet2(j),1))*M_PI/1),2));
//				std::cout << dist << std::endl;
//				std::cout << 'x' << '\t' << m.coords(idxSet1(i),0)-m.coords(idxSet2(j),0) << std::endl;
//				std::cout << 'y' << '\t' << m.coords(idxSet1(i),1)-m.coords(idxSet2(j),1) << std::endl;

			}
			//todo if statements can probably by removed if the calc is done with the previous for loop.
			else if(m.nDims == 3){
				dist = sqrt(pow(m.coords(idxSet1(i),0)-m.coords(idxSet2(j),0),2) + pow(m.coords(idxSet1(i),1)-m.coords(idxSet2(j),1),2) + pow(m.coords(idxSet1(i),2)-m.coords(idxSet2(j),2),2));
			}
			Phi(i,j) = pow((1-(dist/m.r)),4)*(4*(dist/m.r)+1);
		}
	}
//	std::exit(0);
}
//
void rbfGenFunc::getDefVec(Eigen::VectorXd& defVec, int& N, int& steps, Eigen::ArrayXi& movingNodes){
	//todo in this function the N_m from the getNodeType.cpp class should be taken instead of the one from the mesh.cpp



//	std::cout << '\n' << defVec.size() << '\t' << m.N_m << std::endl;


	int idx;
	//loop through the moving nodes
	for(int i = 0; i < N; i++){
		idx = std::distance(std::begin(movingIndices), std::find(std::begin(movingIndices), std::end(movingIndices),movingNodes(i)));
		if(idx!= movingIndices.size()){
//			std::cout << idx << std::endl;
//			std::cout << mIndex(idx) << std::endl;
//			std::cout << displacement.row(idx) << std::endl;
			for(int dim = 0; dim < m.nDims; dim++){
				defVec(dim*N+i) = exactDisp(idx,dim);
			}
		}
	}



//	std::cout << defVec << std::endl;
////	std::cout << "\n" << N << std::endl;
////	Eigen::MatrixXd intPnts(m.N_ib,m.nDims);
//	Eigen::MatrixXd intPnts(ibNodes.size(),m.nDims);
//	Eigen::MatrixXd rotDef;
//
////	intPnts = m.coords(m.intBdryNodes,Eigen::all);
//	intPnts = m.coords(ibNodes,Eigen::all);
//
//
//	if(m.nDims == 2){
//
//		rotDef = (rotMat*(intPnts.rowwise() - params.rotPnt).transpose()).transpose().rowwise() +params.rotPnt - intPnts;
//
//	}
//	else if(m.nDims == 3){
//		rotDef = Eigen::MatrixXd::Zero(m.N_ib,m.nDims);
//
//		if(params.rotVec[0] != 0){
//			rotDef+= (rotMatX*(intPnts.rowwise() - params.rotPnt).transpose()).transpose().rowwise() +params.rotPnt - intPnts;
//		}
//		if(params.rotVec[1] != 0){
//			rotDef+= (rotMatY*(intPnts.rowwise() - params.rotPnt).transpose()).transpose().rowwise() +params.rotPnt - intPnts;
//		}
//		if(params.rotVec[2] != 0){
//			rotDef+= (rotMatZ*(intPnts.rowwise() - params.rotPnt).transpose()).transpose().rowwise() +params.rotPnt - intPnts;
//		}
//	}
////	std::cout << "check" << std::endl;
////	std::cout << defVec << std::endl;
//
//	for(int dim = 0; dim < m.nDims; dim++){
////		defVec(Eigen::seqN(dim*N, m.N_ib)).array() += params.dVec(dim);
////		defVec(Eigen::seqN(dim*N, m.N_ib)) += rotDef.col(dim);
////		defVec(Eigen::seqN(dim*N, ibNodes.size())).array() += params.dVec(dim);
////		defVec(Eigen::seqN(dim*N, ibNodes.size())) += rotDef.col(dim);
//		defVec(dim*N + m.ibIndices).array() += params.dVec(dim);
//		defVec(dim*N + m.ibIndices) += rotDef.col(dim);
//	}
////	std::cout << "check" << std::endl;
////	std::cout << m.coords << std::endl;
//
//
////	std::exit(0);

}

void rbfGenFunc::readDisplacementFile(){
	std::cout << "Reading the displacement file :" << params.dispFile << std::endl;
	std::string line;							// string containing line obtained by getline() function
	std::ifstream file("C:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\Classical RBF\\defs\\" + params.dispFile);



	if(file.is_open()){
		int lineNo = 0;
		while(getline(file, line)){
			std::stringstream ss(line);
			if(m.nDims == 2){
				ss >> movingIndices(lineNo) >> exactDisp(lineNo,0) >> exactDisp(lineNo,1);
			}else if(m.nDims == 3){
				ss >> movingIndices(lineNo) >> exactDisp(lineNo,0) >> exactDisp(lineNo,1) >> exactDisp(lineNo,2);
			}
			lineNo++;
		}
	}else std::cout << "Unable to open the displacment file" << std::endl;

//	std::cout << displacement << std::endl;
//	std::cout << mIndex << std::endl;
}



//void rbfGenFunc::getNodeTypes(){
//	iNodes.resize(m.N_i+m.N_p);
//	iNodes << m.intNodes, m.periodicNodes;
//
//	if(m.smode== "none" || smode == "ps"){
//		mNodes.resize(m.N_ib + m.N_es + m.N_se);
//		mNodes << m.intBdryNodes, m.extStaticNodes, m.slidingEdgeNodes;
//
//		if(m.smode == "ps" && m.pmode == "moving"){
//			mNodesPro.resize(m.N_ib);
//			mNodesPro << m.intBdryNodes;
//
//			sNodes.resize(m.N_se + m.N_es);
//			sNodes << m.slidingEdgeNodes, m.extStaticNodes;
//
//		}else{
//			mNodesPro.resize(m.N_ib+m.N_es);
//			mNodesPro << m.intBdryNodes, m.extStaticNodes;
//
//			sNodes.resize(m.N_se); // sNodes is always the sliding nodes in 2D
//			sNodes << m.slidingEdgeNodes;
//		}
//		N_s = sNodes.size();
//		N_mPro = mNodesPro.size();
//
//	}else if(m.smode=="ds"){
//
//		if(m.pmode == "moving"){
//			mNodes.resize(m.N_ib);
//			mNodes << m.intBdryNodes;
//			sNodes.resize(m.N_se+m.N_es);
//			sNodes << m.slidingEdgeNodes, m.extStaticNodes;
//		}else{
//			mNodes.resize(m.N_ib+m.N_es);
//			mNodes << m.intBdryNodes, m.extStaticNodes;
//			sNodes.resize(m.N_se);
//			sNodes << m.slidingEdgeNodes;
//		}
//		if(curved){
//			// todo rename to make clearer
//			mNodes2.resize(m.N_ib+m.N_es+m.N_se);
//			mNodes2 << m.intBdryNodes,m.extStaticNodes, m.slidingEdgeNodes;
//			N_m2 = mNodes2.size();
//
//		}
//	}
//	N_i = iNodes.size();
//	N_s = sNodes.size();
//	N_m = mNodes.size();
//
//
//}


void rbfGenFunc::getPeriodicParams(){
	pVec.resize(m.nDims);
	pnVec.resize(m.nDims);
	if(m.nDims==2){
		if(m.pmode != "none"){
			if(params.pDir == "x"){
				pVec << 1,0;
				pnVec << 0,1;
			}
			else if(params.pDir == "y"){
				pVec << 0,1;
				pnVec << 1,0;
			}
		}
		else{
			pVec << 0,0;
		}

	// todo do the same for 3D
	}
	else if(m.nDims==3){

	}

//	std::cout << "vector in the periodic direction is: \n" << pVec << std::endl;

}

double rbfGenFunc::rbfEval(double distance){
//	double xi = distance/m.r;	// distance scaled by support radius

	double f_xi = pow((1-(distance/m.r)),4)*(4*(distance/m.r)+1);
	return f_xi;
}



void rbfGenFunc::getDefVecMultiGreedy(Eigen::VectorXd& defVec, getNodeType& n, Eigen::ArrayXXd& errors){
//	std::cout << "testing" << std::endl;
//	std::cout << *n.mPtr << std::endl;
//	std::cout << '\n' << *n.iPtr << std::endl;
	int idx;
	for(int i = 0; i < n.N_m;i++){
		std::cout << i << '\t' << (*n.mPtr)(i) << std::endl;
		idx = std::distance(std::begin(*n.iPtr), std::find(std::begin(*n.iPtr), std::end(*n.iPtr),(*n.mPtr)(i)));
//		std::cout << errors.row(idx) << std::endl;
		std::cout << "index: " << idx << std::endl;
		for(int j = 0; j < errors.cols();j++){
			defVec(n.N_m*j+i) = errors(idx,j);
		}

	}

//	std::cout << defVec << std::endl;




//	int idx;
//	//loop through the moving nodes
//	for(int i = 0; i < N; i++){
//		idx = std::distance(std::begin(mIndex), std::find(std::begin(mIndex), std::end(mIndex),movingNodes(i)));
//		if(idx!= mIndex.size()){
//
//			for(int dim = 0; dim < m.nDims; dim++){
//				defVec(dim*N+i) = displacement(idx,dim);
//			}
//		}
//	}
}
