
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

	getPeriodicParams();
	PhiPtr = &Phis;
}


//std getPhis
void rbfGenFunc::getPhis(Eigen::MatrixXd& Phi_cc, Eigen::MatrixXd& Phi_ic, Eigen::ArrayXi* cPtr, Eigen::ArrayXi* iPtr){
	getPhi(Phi_cc, cPtr, cPtr);
	getPhi(Phi_ic, iPtr, cPtr);
}

// pseudo getPhis
void rbfGenFunc::getPhis(Eigen::MatrixXd& Phi_cc, Eigen::MatrixXd& Phi_sc, Eigen::MatrixXd& Phi_bb, Eigen::MatrixXd& Phi_ib, Eigen::ArrayXi* cPtr, Eigen::ArrayXi* sPtr, Eigen::ArrayXi* bPtr, Eigen::ArrayXi* iPtr){
	getPhi(Phi_cc, cPtr, cPtr);

	getPhi(Phi_sc, sPtr, cPtr);
	getPhi(Phi_bb, bPtr, bPtr);
	getPhi(Phi_ib, iPtr, bPtr);
	//			getPhi(Phi_sm, *n.sePtr, *n.mPtr); 	FOR 2D
}

// direct getPhis 2D
void rbfGenFunc::getPhis(getNodeType& n){
	if(m.nDims == 2){
		getPhi(Phis.Phi_cc, n.cPtr,n.cPtr);
		getPhi(Phis.Phi_cs, n.cPtr, n.sPtr);
		getPhi(Phis.Phi_sc, n.sPtr, n.cPtr);
		getPhi(Phis.Phi_ss, n.sPtr, n.sPtr);
		getPhi(Phis.Phi_ic, n.iPtr, n.cPtr);
		getPhi(Phis.Phi_is, n.iPtr, n.sPtr);
	} else if (m.nDims == 3){
		getPhi(Phis.Phi_cc ,n.cPtr,n.cPtr);
		getPhi(Phis.Phi_ce, n.cPtr, n.sePtr);
		getPhi(Phis.Phi_cs, n.cPtr, n.ssPtr);

		getPhi(Phis.Phi_ec, n.sePtr, n.cPtr);
		getPhi(Phis.Phi_ee, n.sePtr, n.sePtr);
		getPhi(Phis.Phi_es, n.sePtr, n.ssPtr);

		getPhi(Phis.Phi_sc, n.ssPtr, n.cPtr);
		getPhi(Phis.Phi_se, n.ssPtr, n.sePtr);
		getPhi(Phis.Phi_ss, n.ssPtr, n.ssPtr);

		getPhi(Phis.Phi_ic, n.iPtr, n.cPtr);
		getPhi(Phis.Phi_ie, n.iPtr, n.sePtr);
		getPhi(Phis.Phi_is, n.iPtr, n.ssPtr);
	}


}

//void rbfGenFunc::getPhis(Eigen::MatrixXd& Phi_cc, Eigen::MatrixXd& Phi_ce, Eigen::MatrixXd& Phi_cs, Eigen::MatrixXd& Phi_ec, Eigen::MatrixXd& Phi_ee, Eigen::MatrixXd& Phi_es, Eigen::MatrixXd& Phi_sc, Eigen::MatrixXd& Phi_se, Eigen::MatrixXd& Phi_ss, Eigen::MatrixXd& Phi_ic, Eigen::MatrixXd& Phi_ie, Eigen::MatrixXd& Phi_is , getNodeType& n){
//	getPhi(Phi_cc ,n.cPtr,n.cPtr);
//	getPhi(Phi_ce, n.cPtr, n.sePtr);
//	getPhi(Phi_cs, n.cPtr, n.ssPtr);
//
//	getPhi(Phi_ec, n.sePtr, n.cPtr);
//	getPhi(Phi_ee, n.sePtr, n.sePtr);
//	getPhi(Phi_es, n.sePtr, n.ssPtr);
//
//	getPhi(Phi_sc, n.ssPtr, n.cPtr);
//	getPhi(Phi_se, n.ssPtr, n.sePtr);
//	getPhi(Phi_ss, n.ssPtr, n.ssPtr);
//
//	getPhi(Phi_ic, n.iPtr, n.cPtr);
//	getPhi(Phi_ie, n.iPtr, n.sePtr);
//	getPhi(Phi_is, n.iPtr, n.ssPtr);
//}



void rbfGenFunc::getPhi(Eigen::MatrixXd& Phi, Eigen::ArrayXi* idxSet1, Eigen::ArrayXi* idxSet2){
//	std::cout << "lambda: " << m.lambda << std::endl;
	Phi.resize((*idxSet1).size(), (*idxSet2).size());
	double dist;
	for(int i=0; i<(*idxSet1).size();i++){
		for(int j=0; j<(*idxSet2).size();j++){
			if(m.nDims == 2){
				// Euclidian distance

//				dist = sqrt(pow(m.coords(idxSet1(i),0)-m.coords(idxSet22(j),0),2) + pow(m.coords(idxSet1(i),1)-m.coords(idxSet2(j),1),2));


				//todo following if statement is introduced to improve efficiency, check if anything else can be done
				if(params.pmode != "none"){
					dist = 0;
					for(int dim = 0; dim < m.nDims; dim++){
						if(periodicVec(dim)){

							dist += pow(m.lambda/M_PI*sin( (m.coords((*idxSet1)(i),dim)-m.coords((*idxSet2)(j),dim))*M_PI/m.lambda),2);

						}
						else{
							dist += pow(m.coords((*idxSet1)(i),dim)-m.coords((*idxSet2)(j),dim),2);
						}

					}
					dist = sqrt(dist);
				}
				else{
					//todo can the difference in coordinates by calculated for the whole row. then take the power of 2, sum and take square root??
					dist = sqrt(pow(m.coords((*idxSet1)(i),0)-m.coords((*idxSet2)(j),0),2) + pow(m.coords((*idxSet1)(i),1)-m.coords((*idxSet2)(j),1),2) );
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

						if(periodicVec(dim)){
							dist += pow(m.lambda/M_PI*sin( (m.coords((*idxSet1)(i),dim)-m.coords((*idxSet2)(j),dim))*M_PI/m.lambda),2);
//							std::cout << "here: " << pow(m.lambda/M_PI*sin( (m.coords((*idxSet1)(i),dim)-m.coords((*idxSet2)(j),dim))*M_PI/m.lambda),2) << std::endl;
						}else{
							dist += pow(m.coords((*idxSet1)(i),dim)-m.coords((*idxSet2)(j),dim),2);
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
					dist = sqrt(pow(m.coords((*idxSet1)(i),0)-m.coords((*idxSet2)(j),0),2) + pow(m.coords((*idxSet1)(i),1)-m.coords((*idxSet2)(j),1),2) + pow(m.coords((*idxSet1)(i),2)-m.coords((*idxSet2)(j),2),2));
				}
			}
//			std::cout << m.coords(idxSet1(i),2)-m.coords(idxSet2(j),2) << std::endl;
//			std::cout << i << '\t' << j << '\t' << dist << std::endl;
			if(dist/m.r > 1){
				Phi(i,j) = 0;
			}else{
				Phi(i,j) = pow((1-(dist/m.r)),4)*(4*(dist/m.r)+1);
			}

		}
	}
//	std::exit(0);
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

void rbfGenFunc::getDefVec(Eigen::VectorXd& defVec, int N_c, Eigen::ArrayXi* cPtr){
//	std::cout << "determinging std defvec"  << std::endl;
	defVec = Eigen::VectorXd::Zero(N_c*m.nDims);
	int idx;
	int size = (*cPtr).size();
	//loop through the control nodes
	for(int i = 0; i < size; i++){
		idx = std::distance(std::begin(*dispIdx), std::find(std::begin(*dispIdx), std::end(*dispIdx),(*cPtr)(i)));
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
		defVec_b(Eigen::seqN(dim*n.N_b,n.N_c)) = defVec(Eigen::seqN(dim*n.N_c,n.N_c));
		defVec_b(Eigen::seqN(dim*n.N_b+n.N_c,n.N_s)) = finalDef.col(dim);
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




void rbfGenFunc::getPeriodicParams(){
	periodicVec.resize(m.nDims);
	periodicNormalVec1.resize(m.nDims);
	periodicNormalVec2.resize(m.nDims);
	if(m.nDims==2){
		if(params.pmode != "none"){
			if(params.pDir == "x"){
				periodicVec << 1,0;
				periodicNormalVec1 << 0,1;
			}
			else if(params.pDir == "y"){
				periodicVec << 0,1;
				periodicNormalVec1 << 1,0;
			}
		}
		else{
			periodicVec << 0,0;
		}

	// todo define the pnVecs for 3D.
	}
	else if(m.nDims==3){
		if(params.pmode != "none"){
			if(params.pDir == "x"){
				periodicVec << 1,0,0;
				periodicNormalVec1 << 0,1,0;
				periodicNormalVec2 << 0,0,1;
			}else if(params.pDir == "y"){
				periodicVec << 0,1,0;
				periodicNormalVec1 << 1,0,0;
				periodicNormalVec2 << 0,0,1;
			}else if(params.pDir == "z"){
				periodicVec << 0,0,1;
				periodicNormalVec1 << 1,0,0;
				periodicNormalVec2 << 0,1,0;
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



//void rbfGenFunc::getDefVecMultiGreedy(Eigen::VectorXd& defVec, getNodeType& n, Eigen::ArrayXXd& errors, int N, Eigen::ArrayXi* cPtr){
void rbfGenFunc::getDefVec(Eigen::VectorXd& defVec, getNodeType& n, Eigen::ArrayXXd& errors, int N_c){

	defVec.resize(m.nDims*N_c);
	for(int dim = 0; dim< m.nDims; dim++){
		defVec(Eigen::seqN(N_c*dim,N_c)) = -errors(n.cNodesIdx, dim);
	}
}

void rbfGenFunc::performRBF(Eigen::MatrixXd& Phi_cc, Eigen::MatrixXd& Phi_ic, Eigen::VectorXd& defVec, Eigen::ArrayXi* cNodes, Eigen::ArrayXi* iNodes, int& N){

	alpha.resize(N*m.nDims);

	if(params.dataRed){
		d.resize((*iNodes).size(),m.nDims);
	}

	for(int dim = 0; dim < m.nDims; dim++){
//		std::cout << "Solving for dimension: " << dim << std::endl;
//		auto start = std::chrono::high_resolution_clock::now();
		alpha(Eigen::seqN(dim*N,N)) = Phi_cc.fullPivHouseholderQr().solve(defVec(Eigen::seqN(dim*N,N))); //fullPivHouseholderQr()

//		auto stop = std::chrono::high_resolution_clock::now();
//		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);

//		std::cout << "Time to solve for alpha: \t"<<  duration.count()/1e6 << " seconds"<< std::endl;
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
		ptr = ctrlPtr;
		N_m = (*ctrlPtr).size();

//		m.coords(*n.iPtrGrdy,Eigen::all) += deltaInternal;
	}else{
		if(params.smode == "none"){

			ptr = n.cPtr;
			N_m = n.N_c;
		}else{
			ptr = n.bPtr;
			N_m = n.N_b;
		}
	}

//	std::cout << *alpha_step << std::endl;
	auto start = std::chrono::high_resolution_clock::now();


	getPhi(Phi_icGrdy,n.iPtrGrdy,ptr);


	m.coords(*n.iPtr, Eigen::all) += *d_step;



	for(int dim = 0; dim < m.nDims; dim++){
		if(params.multiLvl == false){
			m.coords(*ptr,dim) += (defVec(Eigen::seqN(dim*N_m,N_m))).array();
		}
		m.coords(*n.iPtrGrdy,dim) +=  (Phi_icGrdy*(*alpha_step)(Eigen::seqN(dim*N_m,N_m))).array();
	}

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);

	std::cout << "Time for updating internal nodes: \t"<<  duration.count()/1e6 << " seconds"<< std::endl;


//	auto stop = std::chrono::high_resolution_clock::now();
//	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
//	std::cout << "time: \t"<<  duration.count()/1e6 << " seconds"<< std::endl;

}
