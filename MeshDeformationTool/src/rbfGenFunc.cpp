
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
	movingIndices.resize(m.N_nonzeroDisp);
	exactDisp.resize(m.N_nonzeroDisp,m.nDims);
	readDisplacementFile();

	dispIdx = &movingIndices;
	disp = &exactDisp;

	exactDisp = exactDisp/params.steps; // deformation per step is more usefull then the total deformation.

	PhiPtr = &Phis;
}


void rbfGenFunc::getPhis(getNodeType& n, int iter){

	if(params.dataRed && iter != 0){

		getPhisReduced(n);
	}else{
		getPhisFull(n);
	}

}

void rbfGenFunc::getPhisFull(getNodeType& n){

	// no sliding
	if(params.smode == "none"){
		getPhi(Phis.Phi_cc, n.mPtr, n.mPtr);
		getPhi(Phis.Phi_ic, n.iPtr, n.mPtr);

	}
	// pseudo sliding
	else if(params.smode == "ps"){

		getPhi(Phis.Phi_cc, n.mPtr, n.mPtr);
		getPhi(Phis.Phi_sc, n.sePtr, n.mPtr);

		getPhi(Phis.Phi_cs, n.mPtr, n.sePtr);
		getPhi(Phis.Phi_ss, n.sePtr, n.sePtr);


		Phis.Phi_bb.resize(n.N_m + n.N_se, n.N_m + n.N_se);
		Phis.Phi_bb << Phis.Phi_cc, Phis.Phi_cs, Phis.Phi_sc, Phis.Phi_ss;


		getPhi(Phis.Phi_ib, n.iPtr, n.bPtr);

	}
	// direct sliding
	else{
		if(m.nDims == 2){
			getPhi(Phis.Phi_cc, n.mPtr,n.mPtr);
			getPhi(Phis.Phi_cs, n.mPtr, n.sePtr);
			getPhi(Phis.Phi_sc, n.sePtr, n.mPtr);
			getPhi(Phis.Phi_ss, n.sePtr, n.sePtr);
//			getPhi(Phis.Phi_ic, n.iPtr, n.mPtr);
//			getPhi(Phis.Phi_is, n.iPtr, n.sePtr);

			Phis.Phi_bb.resize(n.N_b, n.N_b);
			Phis.Phi_bb << Phis.Phi_cc, Phis.Phi_cs,
							Phis.Phi_ss, Phis.Phi_sc;


			getPhi(Phis.Phi_ib, n.iPtr, n.bPtr);
//			Phis.Phi_ib.resize(n.N_i, n.N_b);
//			Phis.Phi_ib << Phis.Phi_ic, Phis.Phi_is;

		} else if (m.nDims == 3){
			//todo
			getPhi(Phis.Phi_cc ,n.mPtr,n.mPtr);
			getPhi(Phis.Phi_ce, n.mPtr, n.sePtr);
			getPhi(Phis.Phi_cs, n.mPtr, n.ssPtr);

			getPhi(Phis.Phi_ec, n.sePtr, n.mPtr);
			getPhi(Phis.Phi_ee, n.sePtr, n.sePtr);
			getPhi(Phis.Phi_es, n.sePtr, n.ssPtr);

			getPhi(Phis.Phi_sc, n.ssPtr, n.mPtr);
			getPhi(Phis.Phi_se, n.ssPtr, n.sePtr);
			getPhi(Phis.Phi_ss, n.ssPtr, n.ssPtr);

//			getPhi(Phis.Phi_ic, n.iPtr, n.mPtr);
//			getPhi(Phis.Phi_ie, n.iPtr, n.sePtr);
//			getPhi(Phis.Phi_is, n.iPtr, n.ssPtr);

			Phis.Phi_bb.resize(n.N_b,n.N_b);
			Phis.Phi_bb << Phis.Phi_cc, Phis.Phi_ce, Phis.Phi_cs,
						Phis.Phi_ec, Phis.Phi_ee, Phis.Phi_es,
						Phis.Phi_sc, Phis.Phi_se, Phis.Phi_ss;

//			Phis.Phi_ib.resize(n.N_i, n.N_b);
//			Phis.Phi_ib << Phis.Phi_ic, Phis.Phi_ie, Phis.Phi_is;


		}
	}
}

void rbfGenFunc::getPhisReduced(getNodeType& n){
	// 0 adds row, 1 adds column, 2 adds both
	// no sliding
	if(params.smode == "none"){

		getReducedPhi(Phis.Phi_ic, n);

		for(int i = 0; i < n.addedNodes.idx.size(); i++){
			if(n.addedNodes.type[i] == 0){
				getPhi(Phis.Phi_cc, n.mPtr, n.mPtr, n.addedNodes.idx[i], 2);
				getPhi(Phis.Phi_ic, n.iPtr, n.mPtr, n.addedNodes.idx[i], 1);
			}
		}
	}
	else{

		getReducedPhi(Phis.Phi_ib, n);

		for(int i = 0; i < n.addedNodes.idx.size(); i++){
			if(n.addedNodes.type[i] == 0){
				getPhi(Phis.Phi_cc, n.mPtr, n.mPtr, n.addedNodes.idx[i], 2);
				getPhi(Phis.Phi_sc, n.sePtr, n.mPtr, n.addedNodes.idx[i], 1);
				getPhi(Phis.Phi_cs, n.mPtr, n.sePtr, n.addedNodes.idx[i], 0);
				getPhi(Phis.Phi_ib, n.iPtr, n.bPtr, n.addedNodes.idx[i], 1);
			}else if(n.addedNodes.type[i] == 1){
				getPhi(Phis.Phi_ss, n.sePtr, n.sePtr, n.addedNodes.idx[i], 2);
				getPhi(Phis.Phi_sc, n.sePtr, n.mPtr, n.addedNodes.idx[i], 0);
				getPhi(Phis.Phi_cs, n.mPtr, n.sePtr, n.addedNodes.idx[i], 1);
				getPhi(Phis.Phi_ib, n.iPtr, n.bPtr, n.addedNodes.idx[i]+n.N_m, 1);
			}else{

			}


		}

		Phis.Phi_bb.resize(n.N_m + n.N_se, n.N_m + n.N_se);
		Phis.Phi_bb << Phis.Phi_cc, Phis.Phi_cs, Phis.Phi_sc, Phis.Phi_ss;
	}


}

void rbfGenFunc::getReducedPhi(Eigen::MatrixXd& Phi, getNodeType& n){

	for(int i = 0; i < n.addedNodes.idx_i.size(); i++){
		Phi.middleRows(n.addedNodes.idx_i[i],Phi.rows()-n.addedNodes.idx_i[i]-1) = Phi.bottomRows(Phi.rows()-n.addedNodes.idx_i[i]-1);
	}
	Phi.conservativeResize(n.N_i,Phi.cols());
}

void rbfGenFunc::getPhi(Eigen::MatrixXd& Phi, Eigen::ArrayXi* idxSet1, Eigen::ArrayXi* idxSet2, int newNode, int type){

		double distance;

		// adding a row
		if(type == 0){
			Phi.conservativeResize(Phi.rows()+1,Phi.cols());
			int lastRow = Phi.rows()-1;
			for(int j = 0; j < Phi.cols();j++){
				distance = getDistance((*idxSet1)[newNode],(*idxSet2)[j]);
				Phi(lastRow,j) = rbfEval(distance);
			}
		}
		// adding a column
		else if(type == 1){
			Phi.conservativeResize(Phi.rows(),Phi.cols()+1);
			if(newNode >= Phi.cols()-1){
				int col = Phi.cols()-1;
				for(int j = 0; j < Phi.rows();j++){
					distance = getDistance((*idxSet2)[newNode],(*idxSet1)[j]);
					Phi(j,col) = rbfEval(distance);
				}
			}else{
				for(int col = 1; col < Phi.cols()-(newNode); col++){
					Phi.col(Phi.cols()-col) = Phi.col(Phi.cols()-(col+1));
				}

				for(int j = 0; j < Phi.rows();j++){
					distance = getDistance((*idxSet2)[newNode],(*idxSet1)[j]);
					Phi(j,newNode) = rbfEval(distance);
				}
			}
		}
		// both row and column
		else{
			Phi.conservativeResize(Phi.rows()+1,Phi.cols()+1);
			int lastRow = Phi.rows()-1;
			for(int j = 0; j < Phi.cols();j++){
				distance = getDistance((*idxSet1)[newNode],(*idxSet2)[j]);
				Phi(lastRow,j) = rbfEval(distance);
			}
			Phi.col(lastRow) = Phi.row(lastRow);
		}
}

void rbfGenFunc::getPhi(Eigen::MatrixXd& Phi, Eigen::ArrayXi* idxSet1, Eigen::ArrayXi* idxSet2){
	Phi.resize((*idxSet1).size(), (*idxSet2).size());

	double distance;
	for(int i=0; i<(*idxSet1).size();i++){
		for(int j=0; j<(*idxSet2).size();j++){
			distance = getDistance((*idxSet1)[i],(*idxSet2)[j]);
			Phi(i,j) = rbfEval(distance);
		}
	}
}


double rbfGenFunc::getDistance(int node1, int node2){
	double dist;
	if(m.nDims == 2){
		// Euclidian distance

//				dist = sqrt(pow(m.coords(idxSet1(i),0)-m.coords(idxSet22(j),0),2) + pow(m.coords(idxSet1(i),1)-m.coords(idxSet2(j),1),2));


		//todo following if statement is introduced to improve efficiency, check if anything else can be done
		if(params.pmode != "none"){
			dist = 0;
			for(int dim = 0; dim < m.nDims; dim++){
				if(m.periodicVec(dim)){
					dist += pow(m.lambda/M_PI*sin( (m.coords(node1,dim)-m.coords(node2,dim))*M_PI/m.lambda),2);
				}
				else{
					dist += pow(m.coords(node1,dim)-m.coords(node2,dim),2);
				}
			}
			dist = sqrt(dist);
		}
		else{
			//todo can the difference in coordinates by calculated for the whole row. then take the power of 2, sum and take square root??
			dist = sqrt(pow(m.coords(node1,0)-m.coords(node2,0),2) + pow(m.coords(node1,1)-m.coords(node2,1),2) );
		}
//				dist = sqrt(pow(m.coords(idxSet1(i),0)-m.coords(idxSet2(j),0),2) + pow(1/M_PI*sin( (m.coords(idxSet1(i),1)-m.coords(idxSet2(j),1))*M_PI/1),2));
//				std::cout << dist << std::endl;
//				std::cout << 'x' << '\t' << m.coords(idxSet1(i),0)-m.coords(idxSet2(j),0) << std::endl;
//				std::cout << 'y' << '\t' << m.coords(idxSet1(i),1)-m.coords(idxSet2(j),1) << std::endl;

	}
	//todo if statements can probably by removed if the calc is done with the previous for loop.
	else if(m.nDims == 3){
		if(params.pmode != "none"){
			dist = 0;
//					std::cout << m.coords.row((*idxSet1)(i)) << std::endl;
//					std::cout << m.coords.row((*idxSet2)(j)) << std::endl;
			for(int dim = 0; dim < m.nDims; dim++){

				if(m.periodicVec(dim)){
					dist += pow(m.lambda/M_PI*sin( (m.coords(node1,dim)-m.coords(node2,dim))*M_PI/m.lambda),2);
//							std::cout << "here: " << pow(m.lambda/M_PI*sin( (m.coords((*idxSet1)(i),dim)-m.coords((*idxSet2)(j),dim))*M_PI/m.lambda),2) << std::endl;
				}else{
					dist += pow(m.coords(node1,dim)-m.coords(node2,dim),2);
//							std::cout << "local dist: " << pow(m.coords((*idxSet1)(i),dim)-m.coords((*idxSet2)(j),dim),2) << std::endl;
				}
//						std::cout << " squared distance: "<< dist << std::endl;


			}
			dist = sqrt(dist);
//					std::cout << "actual dist: " << dist << std::endl;
//					if(j>0){
//						std::exit(0);
//					}
		}else{
//					std::cout << "check" << std::endl;
			dist = sqrt(pow(m.coords(node1,0)-m.coords(node2,0),2) + pow(m.coords(node1,1)-m.coords(node2,1),2) + pow(m.coords(node1,2)-m.coords(node2,2),2));
		}
	}
	return dist;
}

//void rbfGenFunc::getDefVec(Eigen::VectorXd& defVec, getNodeType& n, int lvl, Eigen::ArrayXXd& errorPrevLvl){
//
//	int N;
//	Eigen::ArrayXi* mPtr;
//	if(params.smode == "ds" || (params.smode == "ps" && lvl > 0)){
//		N = n.N_mStd;
//		mPtr = n.mStdPtr;
//	}else{
//		N = n.N_m;
//		mPtr = n.mPtr;
//	}
//
//	defVec = Eigen::VectorXd::Zero(N*m.nDims);
//
//
//	if(lvl > 0){
//		getDefVecMultiGreedy(defVec,n, errorPrevLvl,N, mPtr);
//	}else{
//		setDefVec(defVec, n.N_m, n.mPtr);
//	}
//
////	std::cout << '\n' << defVec.size() << '\t' << m.N_m << std::endl;
//
//}

void rbfGenFunc::getDefVec(Eigen::VectorXd& defVec, int N_m, Eigen::ArrayXi* mPtr){
//	std::cout << "determinging std defvec"  << std::endl;
	defVec = Eigen::VectorXd::Zero(N_m*m.nDims);
	int idx;
	int size = (*mPtr).size();
	//loop through the control nodes
	for(int i = 0; i < size; i++){
		idx = std::distance(std::begin(*dispIdx), std::find(std::begin(*dispIdx), std::end(*dispIdx),(*mPtr)(i)));
		if(idx!= (*dispIdx).size()){

			for(int dim = 0; dim < m.nDims; dim++){
				defVec(dim*size+i) = (*disp)(idx,dim);
			}
		}
	}
}

void rbfGenFunc::getDefVec(Eigen::VectorXd& defVec_b, Eigen::VectorXd& defVec, getNodeType& n,Eigen::ArrayXXd& finalDef){
	defVec_b = Eigen::VectorXd::Zero(n.N_b*m.nDims);

	for(int dim = 0; dim< m.nDims; dim++){
		defVec_b(Eigen::seqN(dim*n.N_b,n.N_m)) = defVec(Eigen::seqN(dim*n.N_m,n.N_m));
		defVec_b(Eigen::seqN(dim*n.N_b+n.N_m,n.N_se+n.N_ss)) = finalDef.col(dim);
	}
}

void rbfGenFunc::readDisplacementFile(){
	std::cout << "Reading the displacement file :" << params.dispFile << std::endl;
	std::string line;							// string containing line obtained by getline() function
	std::ifstream file("C:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool\\defs\\" + params.dispFile);


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
	}else std::cout << "Unable to open the displacment file" << std::endl;
//	std::cout << displacement << std::endl;
//	std::cout << mIndex << std::endl;
}






double rbfGenFunc::rbfEval(double distance){
	double xi = distance/m.r;	// distance scaled by support radius
	double f_xi = 0;
	if(xi < 1){
		f_xi = pow((1-(distance/m.r)),4)*(4*(distance/m.r)+1);
	}
	return f_xi;
}


void rbfGenFunc::getDefVec(Eigen::VectorXd& defVec, getNodeType& n, Eigen::ArrayXXd& errors, int N_m){

	defVec.resize(m.nDims*N_m);
	for(int dim = 0; dim< m.nDims; dim++){
		defVec(Eigen::seqN(N_m*dim,N_m)) = -errors(n.cNodesIdx, dim);
	}
}

void rbfGenFunc::performRBF(Eigen::MatrixXd& Phi_cc, Eigen::MatrixXd& Phi_ic, Eigen::VectorXd& defVec, Eigen::ArrayXi* cNodes, Eigen::ArrayXi* iNodes, int& N){

	alpha.resize(N*m.nDims);

	if(params.dataRed){
		d.resize((*iNodes).size(),m.nDims);
	}
	for(int dim = 0; dim < m.nDims; dim++){
		alpha(Eigen::seqN(dim*N,N)) = Phi_cc.fullPivHouseholderQr().solve(defVec(Eigen::seqN(dim*N,N))); //fullPivHouseholderQr()
		if(params.dataRed){
			d.col(dim) = Phi_ic*alpha(Eigen::seqN(dim*N,N));
		}else{
			m.coords(*iNodes,dim) += (Phi_ic*alpha(Eigen::seqN(dim*N,N))).array();
			m.coords(*cNodes,dim) += defVec(Eigen::seqN(dim*N,N)).array();
		}
	}
}

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
			ptr = n.bPtr;
			N_m = n.N_b;
		}
	}


	getPhi(Phi_icGrdy,n.iPtrGrdy,ptr);


	m.coords(*n.iPtr, Eigen::all) += *d_step;



	for(int dim = 0; dim < m.nDims; dim++){
		if(params.multiLvl == false){
			m.coords(*ptr,dim) += (defVec(Eigen::seqN(dim*N_m,N_m))).array();
		}
		m.coords(*n.iPtrGrdy,dim) +=  (Phi_icGrdy*(*alpha_step)(Eigen::seqN(dim*N_m,N_m))).array();
	}
}
