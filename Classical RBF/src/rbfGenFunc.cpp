
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
				if(m.pmode != "none"){
					dist = 0;
					for(int dim = 0; dim < m.nDims; dim++){

						if(pVec(dim)){
							dist += pow(m.lambda/M_PI*sin( (m.coords(idxSet1(i),dim)-m.coords(idxSet2(j),dim))*M_PI/m.lambda),2);
						}else{
							dist += pow(m.coords(idxSet1(i),dim)-m.coords(idxSet2(j),dim),2);
						}

					dist = sqrt(dist);
					}
				}else{
//					std::cout << "check" << std::endl;
					dist = sqrt(pow(m.coords(idxSet1(i),0)-m.coords(idxSet2(j),0),2) + pow(m.coords(idxSet1(i),1)-m.coords(idxSet2(j),1),2) + pow(m.coords(idxSet1(i),2)-m.coords(idxSet2(j),2),2));
				}
			}
//			std::cout << m.coords(idxSet1(i),2)-m.coords(idxSet2(j),2) << std::endl;
//			std::cout << i << '\t' << j << '\t' << dist << std::endl;
			Phi(i,j) = pow((1-(dist/m.r)),4)*(4*(dist/m.r)+1);

		}
	}
//	std::exit(0);
}


void rbfGenFunc::getDefVec(Eigen::VectorXd& defVec, getNodeType& n, int lvl, Eigen::ArrayXXd& errorPrevLvl){

	int N;
	Eigen::ArrayXi* mPtr;
	if(params.smode == "ds" || (params.smode == "ps" && lvl > 0)){
		N = n.N_mStd;
		mPtr = n.mStdPtr;
	}else{
		N = n.N_m;
		mPtr = n.mPtr;
	}

	defVec = Eigen::VectorXd::Zero(N*m.nDims);
//
//	if(params.smode == "ds"){
//		defVec = Eigen::VectorXd::Zero((n.N_m+n.N_se)*m.nDims);
//	}else{
//		defVec = Eigen::VectorXd::Zero(n.N_m*m.nDims);
////		defVec = Eigen::VectorXd::Zero(n.N_mStd*m.nDims);
//	}



	if(params.multiLvl && lvl > 0){
		getDefVecMultiGreedy(defVec,n, errorPrevLvl,N, mPtr);
	}else{
		getDefVecStd(n, defVec);
	}

//	std::cout << '\n' << defVec.size() << '\t' << m.N_m << std::endl;



}

void rbfGenFunc::getDefVecStd(getNodeType& n, Eigen::VectorXd& defVec){
	std::cout << "determinging std defvec"  << std::endl;

//	std::cout << exactDisp << std::endl;
//	std::cout << *n.mPtr << std::endl;
//	std::cout << movingIndices << std::endl;

	int idx;
	//loop through the moving nodes
	for(int i = 0; i < n.N_m; i++){
		idx = std::distance(std::begin(movingIndices), std::find(std::begin(movingIndices), std::end(movingIndices),(*n.mPtr)(i)));
		if(idx!= movingIndices.size()){
//			std::cout << idx << std::endl;
//			std::cout << mIndex(idx) << std::endl;
//			std::cout << displacement.row(idx) << std::endl;
			for(int dim = 0; dim < m.nDims; dim++){
				defVec(dim*n.N_m+i) = exactDisp(idx,dim);
			}
		}
	}
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

	// todo define the pnVecs for 3D.
	}
	else if(m.nDims==3){
		if(m.pmode != "none"){
			if(params.pDir == "x"){
				pVec << 1,0,0;
			}else if(params.pDir == "y"){
				pVec << 0,1,0;
			}else if(params.pDir == "z"){
				pVec << 0,0,1;
			}
		}
	}

//	std::cout << "vector in the periodic direction is: \n" << pVec << std::endl;

}

double rbfGenFunc::rbfEval(double distance){
//	double xi = distance/m.r;	// distance scaled by support radius

	double f_xi = pow((1-(distance/m.r)),4)*(4*(distance/m.r)+1);
	return f_xi;
}



void rbfGenFunc::getDefVecMultiGreedy(Eigen::VectorXd& defVec, getNodeType& n, Eigen::ArrayXXd& errors, int N, Eigen::ArrayXi*& mPtr){
//	std::cout << "testing" << std::endl;
//	std::cout << *n.mPtr << std::endl;
//	std::cout << '\n' << *n.iPtr << std::endl;

	//todo check difference for pointers *n.mPtr/*n.mStdPtr for none/ps sliding and in the for loop
	int idx;
	for(int i = 0; i < N ;i++){
//		std::cout << i << '\t' << (*n.mPtr)(i) << std::endl;
		idx = std::distance(std::begin(*n.iPtr), std::find(std::begin(*n.iPtr), std::end(*n.iPtr),(*mPtr)(i)));
//		std::cout << errors.row(idx) << std::endl;
//		std::cout << "index: " << idx << std::endl;
		for(int j = 0; j < errors.cols();j++){
			defVec(N*j+i) = errors(idx,j);
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
