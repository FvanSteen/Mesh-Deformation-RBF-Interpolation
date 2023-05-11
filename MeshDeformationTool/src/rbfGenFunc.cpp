
#include "rbfGenFunc.h"
#include <iostream>
#include <Eigen/Dense>
#include <chrono>
#include <fstream>
#include <sstream>
#include "CoordTransform.h"
rbfGenFunc::rbfGenFunc(Mesh& meshObject, struct probParams& probParamsObject)
//rbfGenFunc::rbfGenFunc(Mesh* meshPtr, Eigen::VectorXd& dVec, Eigen::RowVectorXd& rotPnt, Eigen::VectorXd& rotVec, const int& steps, const std::string& smode, const bool& curved, const std::string& pDir)
:m(meshObject), params(probParamsObject)
{

	std::cout << "Initialised the rbfGenFunc class" << std::endl;
	movingIndices.resize(m.N_nonzeroDisp);
	exactDisp.resize(m.N_nonzeroDisp,m.nDims);

	readDisplacementFile();
	exactDisp = exactDisp/params.steps; // deformation per step is more usefull then the total deformation.
	exactDispPtr = &exactDisp;
	dispIdx = &movingIndices;

	if(params.ptype){
		exactDisp_polar_cylindrical.resize(m.N_nonzeroDisp, m.nDims);
		CoordTransform transform;
		transform.vector_cart_to_polar_cylindrical(exactDisp, exactDisp_polar_cylindrical, movingIndices, m.coords);
		disp = &exactDisp_polar_cylindrical;

	}else{
		disp = &exactDisp;
	}

	PhiPtr = &Phis;
}


void rbfGenFunc::getPhis(getNodeType& n, int iter){

	if(params.dataRed && iter != 0){
		getPhisReduced(n);
	}else{
		getPhisFull(n);
	}

	if(params.smode == "ps"){

		Phis.Phi_mm.resize(n.N_m, n.N_m);
		Phis.Phi_mm = Phis.Phi_cc.block(0,0,n.N_m, n.N_m);

		Phis.Phi_em.resize(n.N_se, n.N_m);
		Phis.Phi_em = Phis.Phi_cc.block(n.N_m, 0, n.N_se, n.N_m);

		if(m.nDims == 3){
			Phis.Phi_meme.resize(n.N_m + n.N_se,n.N_m + n.N_se);
			Phis.Phi_meme = Phis.Phi_cc.block(0,0,n.N_m+n.N_se, n.N_m+n.N_se);

			Phis.Phi_sme.resize(n.N_ss, n.N_m+n.N_se);
			Phis.Phi_sme = Phis.Phi_cc.block(n.N_m+n.N_se, 0, n.N_ss, n.N_m+n.N_se);

			Phis.Phi_ic_reduced = Phis.Phi_ic.block(0,0,n.N_i- m.N_ss + n.N_ss, n.N_m+n.N_se);
		}

	}else if(params.smode == "ds"){

		Phis.Phi_mc.resize(n.N_m, n.N_c);
		Phis.Phi_mc = Phis.Phi_cc.block(0,0,n.N_m, n.N_c);

		Phis.Phi_ec.resize(n.N_se, n.N_c);
		Phis.Phi_ec = Phis.Phi_cc.block(n.N_m,0,n.N_se, n.N_c);

		Phis.Phi_sc.resize(n.N_ss, n.N_c);
		Phis.Phi_sc = Phis.Phi_cc.block(n.N_m+n.N_se,0,n.N_ss, n.N_c);
	}
}

void rbfGenFunc::getPhisFull(getNodeType& n){

	getPhi(Phis.Phi_cc, n.cPtr, n.cPtr);
	getPhi(Phis.Phi_ic, n.iPtr, n.cPtr);

}


void rbfGenFunc::getPhisReduced(getNodeType& n){
	// 0 adds row, 1 adds column, 2 adds both

	// removing rows from the Phi_ic matrix
	getReducedPhi(Phis.Phi_ic, n);

	// adjsuting sizes of the interpolation matrices
	adjustPhi(Phis.Phi_cc, n, 2);
	adjustPhi(Phis.Phi_ic, n, 1);

	// adding the new rows and/ or columns to the matrices
	for(int i = 0; i < n.addedNodes.idx.size(); i++){
		switch (n.addedNodes.type[i]){

		// in case of moving node
			case 0:
				getPhi(Phis.Phi_cc, n.cPtr, n.cPtr, n.addedNodes.idx[i], 2);
				getPhi(Phis.Phi_ic, n.iPtr, n.cPtr, n.addedNodes.idx[i], 1);
				break;

		// in case of edge node
			case 1:
				getPhi(Phis.Phi_cc, n.cPtr, n.cPtr, n.addedNodes.idx[i]+n.N_m, 2);
				getPhi(Phis.Phi_ic, n.iPtr, n.cPtr, n.addedNodes.idx[i]+n.N_m, 1);
				break;

		// in case of surface node
			case 2:
				getPhi(Phis.Phi_cc, n.cPtr, n.cPtr, n.addedNodes.idx[i]+n.N_m+n.N_se, 2);
				getPhi(Phis.Phi_ic, n.iPtr, n.cPtr, n.addedNodes.idx[i]+n.N_m+n.N_se, 1);
				break;
		}

	}
}

void rbfGenFunc::getReducedPhi(Eigen::MatrixXd& Phi, getNodeType& n){
	for(int i = 0; i < n.addedNodes.idx_i.size(); i++){
		Phi.middleRows(n.addedNodes.idx_i[i],Phi.rows()-n.addedNodes.idx_i[i]-1) = Phi.bottomRows(Phi.rows()-n.addedNodes.idx_i[i]-1);
	}
	Phi.conservativeResize(n.N_i,Phi.cols());
}

void rbfGenFunc::adjustPhi(Eigen::MatrixXd& Phi, getNodeType& n,  int type){

	int idx;
	switch(type){
		case 0:
			Phi.conservativeResize(Phi.rows()+n.addedNodes.idx.size(), Phi.cols());
			break;
		case 1:
			Phi.conservativeResize(Phi.rows(), Phi.cols()+n.addedNodes.idx.size());
			break;
		case 2:
			Phi.conservativeResize(Phi.rows()+n.addedNodes.idx.size(), Phi.cols()+n.addedNodes.idx.size());
			break;
	}


	for(int i = 0; i < n.addedNodes.idx.size(); i++){
		int idx_shift = 0;
		if(n.addedNodes.idx.size() > 1 && n.addedNodes.type[1] < n.addedNodes.type[0]){
			idx_shift = -1;
		}

		switch(n.addedNodes.type[i]){
		case 0:
			idx = n.addedNodes.idx[i];
			break;
		case 1:
			idx = n.addedNodes.idx[i] + n.N_m + idx_shift;
			break;
		case 2:
			idx = n.addedNodes.idx[i] + n.N_m + n.N_se + idx_shift;
			break;
		}

		if(idx < 0){
			idx = 0;
		}

		if( type > 0){
			int col = Phi.cols()-1;
			while(col > idx){
				Phi.col(col) = Phi.col(col-1);
				col--;
			}
		}
		if( type != 1){
			int row = Phi.rows()-1;
			while(row > idx){
				Phi.row(row) = Phi.row(row-1);
				row--;
			}
		}
	}

}


void rbfGenFunc::getPhi(Eigen::MatrixXd& Phi, Eigen::ArrayXi* idxSet1, Eigen::ArrayXi* idxSet2, int idx, int type){

	double distance;

	//depending on type add row, column or both
	switch (type) {

	// adding a row at row idx
		case 0:
			for(int col = 0; col < Phi.cols(); col++){
				distance = getDistance((*idxSet1)[idx],(*idxSet2)[col]);
				Phi(idx,col) = rbfEval(distance);
			}
			break;

	// adding column at column idx
		case 1:
			for(int row = 0; row < Phi.rows(); row++){
				distance = getDistance((*idxSet1)[row], (*idxSet2)[idx]);
				Phi(row,idx) = rbfEval(distance);
			}
			break;

	// adding both a row and column at idx
		case 2:
			for(int col = 0; col < Phi.cols(); col++){
				distance = getDistance((*idxSet1)[idx],(*idxSet2)[col]);
				Phi(idx,col) = rbfEval(distance);
			}

			Phi.col(idx) = Phi.row(idx);
			break;
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
	if(params.ptype){
		if(m.nDims == 2){

			dist = sqrt(pow(m.coords_polar_cylindrical(node1,0),2) + pow(m.coords_polar_cylindrical(node2,0),2) -2*m.coords_polar_cylindrical(node1,0)*m.coords_polar_cylindrical(node2,0)*cos(m.periodic_length/M_PI*sin( (m.coords_polar_cylindrical(node2,1)-m.coords_polar_cylindrical(node1,1))*M_PI/m.periodic_length) ));

		}else if(m.nDims == 3){
			dist = sqrt(pow(m.coords_polar_cylindrical(node1,0),2) + pow(m.coords_polar_cylindrical(node2,0),2) -2*m.coords_polar_cylindrical(node1,0)*m.coords_polar_cylindrical(node2,0)*cos(m.periodic_length/M_PI*sin( (m.coords_polar_cylindrical(node2,1)-m.coords_polar_cylindrical(node1,1))*M_PI/m.periodic_length)) + pow(m.coords_polar_cylindrical(node2,2) - m.coords_polar_cylindrical(node1,2),2) );
		}
	}else{
		if(m.nDims == 2){
			// Euclidian distance

	//				dist = sqrt(pow(m.coords(idxSet1(i),0)-m.coords(idxSet22(j),0),2) + pow(m.coords(idxSet1(i),1)-m.coords(idxSet2(j),1),2));


			//todo following if statement is introduced to improve efficiency, check if anything else can be done
			if(params.pmode != "none"){
				dist = 0;
				for(int dim = 0; dim < m.nDims; dim++){
					if(dim == params.pDir){
						dist += pow(m.periodic_length/M_PI*sin( (m.coords(node1,dim)-m.coords(node2,dim))*M_PI/m.periodic_length),2);
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

					if(dim == params.pDir){
						dist += pow(m.periodic_length/M_PI*sin( (m.coords(node1,dim)-m.coords(node2,dim))*M_PI/m.periodic_length),2);
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
	}
	return dist;
}


void rbfGenFunc::getDefVec(Eigen::VectorXd& defVec, int N_m, Eigen::ArrayXi* mPtr){

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

void rbfGenFunc::getDefVec(Eigen::VectorXd& defVec_b, Eigen::VectorXd& defVec, getNodeType& n,Eigen::ArrayXXd& finalDef, int N, int N_init){
	defVec_b = Eigen::VectorXd::Zero(N*m.nDims);

	for(int dim = 0; dim< m.nDims; dim++){
		defVec_b(Eigen::seqN(dim*N,N_init)) = defVec(Eigen::seqN(dim*N_init,N_init));
		defVec_b(Eigen::seqN(dim*N+N_init,N-N_init)) = finalDef.col(dim);
	}
}

void rbfGenFunc::readDisplacementFile(){
	std::cout << "Reading the displacement file :" << params.dispFile << std::endl;
	std::string line;							// string containing line obtained by getline() function
	// todo find file path automatically
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
	}else std::cout << "Unable to open the displacement file" << std::endl;

	//todo do the division by the number of steps here
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

void rbfGenFunc::performRBF(Eigen::MatrixXd& Phi_cc, Eigen::MatrixXd& Phi_ic, Eigen::VectorXd& defVec, Eigen::ArrayXi* cNodes, Eigen::ArrayXi* iNodes, int N){

	alpha.resize(N*m.nDims);

	if(params.dataRed){
		d.resize((*iNodes).size(),m.nDims);
	}
	for(int dim = 0; dim < m.nDims; dim++){

		alpha(Eigen::seqN(dim*N,N)) = Phi_cc.householderQr().solve(defVec(Eigen::seqN(dim*N,N))); //fullPivHouseholderQr()

		if(params.dataRed){
			d.col(dim) = Phi_ic*alpha(Eigen::seqN(dim*N,N));
		}else{
			(*m.ptrCoords)(*iNodes,dim) += (Phi_ic*alpha(Eigen::seqN(dim*N,N))).array();
			(*m.ptrCoords)(*cNodes,dim) += defVec(Eigen::seqN(dim*N,N)).array();
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
			ptr = n.cPtr;
			N_m = n.N_c;
		}
	}


	getPhi(Phi_icGrdy,n.iPtrGrdy,ptr);
	for(int dim = 0; dim < m.nDims; dim++){

		if(params.multiLvl == false){
			(*m.ptrCoords)(*ptr,dim) += (defVec(Eigen::seqN(dim*N_m,N_m))).array();
		}
		(*m.ptrCoords)(*n.iPtrGrdy,dim) +=  (Phi_icGrdy*(*alpha_step)(Eigen::seqN(dim*N_m,N_m))).array();
	}

	// todo call class more elegantly
	if(params.ptype){
		CoordTransform transform;
		transform.polar_cylindrical_to_cart(m.coords_polar_cylindrical, m.coords);
	}

	m.coords(*n.iPtr, Eigen::all) += *d_step;
}
