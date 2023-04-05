
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

	Phis.Phi_mm.resize(n.N_m, n.N_m);
	Phis.Phi_mm = Phis.Phi_cc.block(0,0,n.N_m, n.N_m);


	Phis.Phi_em.resize(n.N_se, n.N_m);
	Phis.Phi_em = Phis.Phi_cc.block(n.N_m, 0, n.N_se, n.N_m);

}

void rbfGenFunc::getPhisFull(getNodeType& n){

	// no sliding
	if(params.smode == "none"){
		getPhi(Phis.Phi_cc, n.cPtr, n.cPtr);
		getPhi(Phis.Phi_ic, n.iPtr, n.cPtr);

	}
	// pseudo sliding
	else if(params.smode == "ps"){
		// todo can probably be together with ds

		getPhi(Phis.Phi_cc, n.cPtr, n.cPtr);

//		Phis.Phi_mm.resize(n.N_m, n.N_m);
//		Phis.Phi_mm = Phis.Phi_cc.block(0,0,n.N_m, n.N_m);
//
//
//		Phis.Phi_em.resize(n.N_se, n.N_m);
//		Phis.Phi_em = Phis.Phi_cc.block(n.N_m, 0, n.N_se, n.N_m);
//		getPhi(Phis.Phi_mm, n.mPtr, n.mPtr);
//		getPhi(Phis.Phi_me, n.mPtr, n.sePtr);
//		getPhi(Phis.Phi_em, n.sePtr, n.mPtr);
//		getPhi(Phis.Phi_ee, n.sePtr, n.sePtr);


//		Phis.Phi_cc.resize(n.N_m + n.N_se, n.N_m + n.N_se);
//		Phis.Phi_cc << Phis.Phi_mm, Phis.Phi_me, Phis.Phi_em, Phis.Phi_ee;


		getPhi(Phis.Phi_ic, n.iPtr, n.cPtr);

	}
	// direct sliding todo can be added to previous ps
	else{
		//todo start here find Phi_cc and obtain smaller matrices from that
		// moving nodes

		getPhi(Phis.Phi_mm, n.mPtr,n.mPtr);

		// moving nodes and edge nodes
		getPhi(Phis.Phi_me, n.mPtr, n.sePtr);
		getPhi(Phis.Phi_em, n.sePtr, n.mPtr);


		// the real g's
		getPhi(Phis.Phi_mc, n.mPtr, n.cPtr);
		getPhi(Phis.Phi_ec, n.sePtr, n.cPtr);
		getPhi(Phis.Phi_sc, n.ssPtr, n.cPtr);

		// edge nodes
		getPhi(Phis.Phi_ee, n.sePtr, n.sePtr);

		// internal nodes and control nodes
		getPhi(Phis.Phi_ic, n.iPtr, n.cPtr);

		// moving nodes and surface nodes
		getPhi(Phis.Phi_ms, n.mPtr, n.ssPtr);
		getPhi(Phis.Phi_sm, n.ssPtr, n.mPtr);

		// edge nodes and surface nodes
		getPhi(Phis.Phi_es, n.sePtr, n.ssPtr);
		getPhi(Phis.Phi_se, n.ssPtr, n.sePtr);

		// surface nodes
		getPhi(Phis.Phi_ss, n.ssPtr, n.ssPtr);



		// control nodes
		Phis.Phi_cc.resize(n.N_c, n.N_c);
		if(params.curved){
			if(m.nDims == 2){
				Phis.Phi_cc << Phis.Phi_mm, Phis.Phi_me,
								Phis.Phi_em, Phis.Phi_ee;
			}else if(m.nDims == 3){
				Phis.Phi_cc << Phis.Phi_mm, Phis.Phi_me, Phis.Phi_ms,
							Phis.Phi_em, Phis.Phi_ee, Phis.Phi_es,
							Phis.Phi_sm, Phis.Phi_se, Phis.Phi_ss;
			}
		}
	}
}

void rbfGenFunc::adjustPhi(Eigen::MatrixXd& Phi, getNodeType& n,  int type){
//	std::cout << "NODES:\n" << n.addedNodes.idx << "\n\n\n" <<  std::endl;

//	std::cout << Phi << std::endl;
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
		switch(n.addedNodes.type[i]){
		case 0:
			idx = n.addedNodes.idx[i];
			break;
		case 1:
			idx = n.addedNodes.idx[i] + n.N_m;
			break;
		case 2:
			idx = n.addedNodes.idx[i] + n.N_m + n.N_se;
			break;
		}

		if( type > 0){
			int col = Phi.cols()-1;
		//		std::cout << "idx: " << idx << std::endl;
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
//		std::cout << "\n" << Phi << "\n" << std::endl;

	}

//	std::cout << Phi << std::endl;

}

void rbfGenFunc::getPhi2(Eigen::MatrixXd& Phi, Eigen::ArrayXi* idxSet1, Eigen::ArrayXi* idxSet2, int idx, int type){
//	std::cout << "\nGETPHI2\n" << std::endl;

//	Phi(idx,0) = 69;
//	Phi(0,idx) = 69;
//	std::cout << "\n" << Phi << std::endl;
//	std::cout << idx << std::endl;
	double distance;
	switch (type) {
		case 1:
			for(int row = 0; row < Phi.rows(); row++){
				distance = getDistance((*idxSet1)[row], (*idxSet2)[idx]);
				Phi(row,idx) = rbfEval(distance);
			}
			break;
		case 2:
			for(int col = 0; col < Phi.cols(); col++){
				distance = getDistance((*idxSet1)[idx],(*idxSet2)[col]);
				Phi(idx,col) = rbfEval(distance);
			}

			Phi.col(idx) = Phi.row(idx);
			break;
	}


//	std::cout << "\n" << Phi << "\n" << std::endl;
}

void rbfGenFunc::getPhisReduced(getNodeType& n){
	// 0 adds row, 1 adds column, 2 adds both
	// no sliding
//	std::cout << "newNodes:\n " << n.addedNodes.idx << std::endl;
	if(params.smode == "none"){

		getReducedPhi(Phis.Phi_ic, n);

		for(int i = 0; i < n.addedNodes.idx.size(); i++){
			if(n.addedNodes.type[i] == 0){
				getPhi(Phis.Phi_cc, n.cPtr, n.cPtr, n.addedNodes.idx[i], 2);
				getPhi(Phis.Phi_ic, n.iPtr, n.cPtr, n.addedNodes.idx[i], 1);
			}
		}
	}
	else{

		getReducedPhi(Phis.Phi_ic, n);

//		std::cout << "\n Initial Phi:\n" << Phis.Phi_cc << std::endl;
		Eigen::Array2i idx_order = {0,1};

		if(n.addedNodes.type[1] < n.addedNodes.type[0]){
			idx_order << 1,0;
		}

		adjustPhi(Phis.Phi_cc, n, 2);
		adjustPhi(Phis.Phi_ic, n, 1);

		for(int i = 0; i < n.addedNodes.idx.size(); i++){
			switch (n.addedNodes.type[i]){
				case 0:
					getPhi2(Phis.Phi_cc, n.cPtr, n.cPtr, n.addedNodes.idx[i], 2);
					getPhi2(Phis.Phi_ic, n.iPtr, n.cPtr, n.addedNodes.idx[i], 1);
					break;
				case 1:
					getPhi2(Phis.Phi_cc, n.cPtr, n.cPtr, n.addedNodes.idx[i]+n.N_m, 2);
					getPhi2(Phis.Phi_ic, n.iPtr, n.cPtr, n.addedNodes.idx[i]+n.N_m, 1);
					break;
				case 2:
					break;



			}



		}
/*
		for(int i = 0; i < n.addedNodes.idx.size(); i++){
//			std::cout << idx_order[i]<< std::endl;

			if(n.addedNodes.type[idx_order[i]] == 0){
//				std::cout << "moving\n";
				getPhi(Phis.Phi_cc, n.cPtr, n.cPtr, n.addedNodes.idx[idx_order[i]], 2);
//				std::cout << "\n" << Phis.Phi_cc << std::endl;
//				getPhi(Phis.Phi_mm, n.mPtr, n.mPtr, n.addedNodes.idx[i], 2);
//				getPhi(Phis.Phi_em, n.sePtr, n.mPtr, n.addedNodes.idx[i], 1);
//				getPhi(Phis.Phi_me, n.mPtr, n.sePtr, n.addedNodes.idx[i], 0);
//				getPhi(Phis.Phi_ms, n.mPtr, n.ssPtr, n.addedNodes.idx[i], 0);
//				getPhi(Phis.Phi_sm, n.ssPtr, n.mPtr, n.addedNodes.idx[i], 1);

				// required ones below
//				getPhi(Phis.Phi_mc, n.mPtr, n.cPtr, n.addedNodes.idx[i], 1);
//				getPhi(Phis.Phi_ec, n.sePtr, n.cPtr, n.addedNodes.idx[i], 1);
//				getPhi(Phis.Phi_sc, n.ssPtr, n.cPtr, n.addedNodes.idx[i], 1);
				getPhi(Phis.Phi_ic, n.iPtr, n.cPtr, n.addedNodes.idx[idx_order[i]], 1);

			}else if(n.addedNodes.type[idx_order[i]] == 1){
//				std::cout << "edge\n";

				getPhi(Phis.Phi_cc, n.cPtr, n.cPtr, n.addedNodes.idx[idx_order[i]]+n.N_m ,2);

//				std::cout << "\n" << Phis.Phi_cc << std::endl;
//				std::cout << Phis.Phi_mc <<std::endl;

//				std::cout << "\n" << Phis.Phi_mc <<std::endl;

//				getPhi(Phis.Phi_ee, n.sePtr, n.sePtr, n.addedNodes.idx[i], 2);
//				getPhi(Phis.Phi_em, n.sePtr, n.mPtr, n.addedNodes.idx[idx_order[i]], 0);
//				getPhi(Phis.Phi_me, n.mPtr, n.sePtr, n.addedNodes.idx[idx_order[i]],1);
//
//				getPhi(Phis.Phi_es, n.sePtr, n.ssPtr, n.addedNodes.idx[idx_order[i]], 0);
//				getPhi(Phis.Phi_se, n.ssPtr, n.sePtr, n.addedNodes.idx[idx_order[i]], 1);
//
//
//				// required ones below
//				getPhi(Phis.Phi_mc, n.mPtr, n.cPtr, n.addedNodes.idx[idx_order[i]] + n.N_m, 1);
//				getPhi(Phis.Phi_ec, n.sePtr, n.cPtr, n.addedNodes.idx[idx_order[i]] +n.N_m, 1);
//				getPhi(Phis.Phi_sc, n.ssPtr, n.cPtr, n.addedNodes.idx[idx_order[i]] +n.N_m, 1);
				getPhi(Phis.Phi_ic, n.iPtr, n.cPtr, n.addedNodes.idx[idx_order[i]]+n.N_m, 1);
			}else if(n.addedNodes.type[idx_order[i]] == 2){

				std::cout << "check the idx_shift\n";
				std::exit(0);

				getPhi(Phis.Phi_cc, n.cPtr, n.cPtr, n.addedNodes.idx[i]+n.N_m+n.N_se, 2);
//				std::cout << "surf\n";
//				std::cout << Phis.Phi_mc <<std::endl;
				getPhi(Phis.Phi_ms, n.mPtr, n.ssPtr, n.addedNodes.idx[i], 1);



//				std::cout << "\n" << Phis.Phi_mc <<std::endl;



				getPhi(Phis.Phi_sm, n.ssPtr, n.mPtr, n.addedNodes.idx[i], 0);

				getPhi(Phis.Phi_es, n.sePtr, n.ssPtr, n.addedNodes.idx[i], 1);

				getPhi(Phis.Phi_se, n.ssPtr, n.sePtr, n.addedNodes.idx[i], 0);
				getPhi(Phis.Phi_ss, n.ssPtr, n.ssPtr, n.addedNodes.idx[i],2);

				//required ones below
				getPhi(Phis.Phi_mc, n.mPtr, n.cPtr, n.addedNodes.idx[i] + n.N_m+ n.N_se, 1);
				getPhi(Phis.Phi_ec, n.sePtr, n.cPtr, n.addedNodes.idx[i] +n.N_m + n.N_se, 1);
				getPhi(Phis.Phi_sc, n.ssPtr, n.cPtr, n.addedNodes.idx[i] +n.N_m + n.N_se, 1);
				getPhi(Phis.Phi_ic, n.iPtr, n.cPtr, n.addedNodes.idx[i]+n.N_m+n.N_se, 1);
			}


//			std::cout << "\n" << Phis.Phi_mc << std::endl;
//			std::cout << "moving\n" << Phis.Phi_mm << "\n\nedge\n" << Phis.Phi_me << "\n\nsurf\n" << Phis.Phi_ms << std::endl;

		}

//		for(int i = 0; i < n.addedNodes.idx.size(); i++){
//			if(n.addedNodes.type[i] == 0){
//				getPhi(Phis.Phi_mc, n.mPtr, n.cPtr, n.addedNodes.idx[i], 0);
//			}else if(n.addedNodes.type[i] == 1){
//				getPhi(Phis.Phi_ec, n.sePtr, n.cPtr, n.addedNodes.idx[i], 0);
//			}else if(n.addedNodes.type[i] == 2){
//				getPhi(Phis.Phi_sc, n.ssPtr, n.cPtr, n.addedNodes.idx[i],0);
//
//			}
//		}
//		std::cout << "after adding potential row:\n" << Phis.Phi_mc << std::endl;

//		if(params.curved || params.smode == "ps"){
//			Phis.Phi_cc.resize(n.N_c, n.N_c);
//			if(m.nDims == 2){
//				Phis.Phi_cc << Phis.Phi_mm, Phis.Phi_me, Phis.Phi_em, Phis.Phi_ee;
//			}else if (m.nDims == 3){
//				Phis.Phi_cc << Phis.Phi_mm, Phis.Phi_me, Phis.Phi_ms,
//								Phis.Phi_em, Phis.Phi_ee, Phis.Phi_es,
//								Phis.Phi_sm, Phis.Phi_se, Phis.Phi_ss;
//			}
//		}


	}
//	if(Phis.Phi_mc.cols() == 14){
//		std::cout << Phis.Phi_mc << "\n\n" << Phis.Phi_me << "\n\n" << Phis.Phi_ms << std::endl;
//		std::exit(0);
//	}


//	Phis.Phi_mm.resize(n.N_m, n.N_m);
//	Phis.Phi_mm = Phis.Phi_cc.block(0,0,n.N_m, n.N_m);
//
//
//	Phis.Phi_em.resize(n.N_se, n.N_m);
//	Phis.Phi_em = Phis.Phi_cc.block(n.N_m, 0, n.N_se, n.N_m);
//	std::cout << "\n\n" << Phis.Phi_cc << std::endl;
//	std::cout << "\n\n" << Phis.Phi_em << std::endl;
//	std::cout << "\nresulting PHI_SS\n" << Phis.Phi_cc<< std::endl;

 */
	}
//	std::cout << "\nresulting PHI_ic\n" << Phis.Phi_ic.colwise().sum()<< std::endl;
//	getPhi(Phis.Phi_cc, n.cPtr, n.cPtr);
//	getPhi(Phis.Phi_mm, n.mPtr, n.mPtr);
//	getPhi(Phis.Phi_em, n.sePtr, n.mPtr);
//	getPhi(Phis.Phi_cc, n.cPtr, n.cPtr);
//	getPhi(Phis.Phi_ic, n.iPtr, n.cPtr);
//	std::cout << "\n\n" << Phis.Phi_em << std::endl;
//	std::cout << "\n\n" << Phis.Phi_cc << std::endl;
//	std::cout << "\nACTUAL PHI_cc\n" << Phis.Phi_ic.colwise().sum()<< std::endl;

//	if(Phis.Phi_cc.cols() == 12){
//		std::cout << "exit status reached\n";
//		std::exit(0);
//	}

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
//				std::cout << j << '\t' << (*idxSet1)[newNode] << '\t' << (*idxSet2)[j] << std::endl;
				distance = getDistance((*idxSet1)[newNode],(*idxSet2)[j]);
				Phi(lastRow,j) = rbfEval(distance);
			}

		}
		// adding a column
		else if(type == 1){
//			std::cout << "\n\n" << Phi << std::endl;
			Phi.conservativeResize(Phi.rows(),Phi.cols()+1);
			if(newNode>= Phi.cols()-1){
				int col = Phi.cols()-1;
				for(int j = 0; j < Phi.rows();j++){
					distance = getDistance((*idxSet2)[newNode],(*idxSet1)[j]);
					Phi(j,col) = rbfEval(distance);
				}
			}else{
				for(int col = 1; col < Phi.cols()-(newNode); col++){
					Phi.col(Phi.cols()-col) = Phi.col(Phi.cols()-(col+1));
//					std::cout << "\n\n"<<  Phi << std::endl;
				}

				for(int j = 0; j < Phi.rows();j++){
					distance = getDistance((*idxSet2)[newNode],(*idxSet1)[j]);
					Phi(j,newNode) = rbfEval(distance);
//					std::cout << "\n" << Phi(j,newNode) << "\n\n";
				}
			}
//			std::cout << "\n\n" << Phi << std::endl;
		}
		// both row and column
		else{


			Phi.conservativeResize(Phi.rows()+1,Phi.cols()+1);
			std::cout << "NEWNODE: " << newNode << std::endl;
//			int lastRow = Phi.rows()-1;
			std::cout << Phi << std::endl;
			int i = Phi.cols()-1;
			while(i > newNode){
//				std::cout << i << std::endl;
				Phi.col(i) = Phi.col(i-1);
//				std::cout << Phi << std::endl;
				i--;
			}
			i = Phi.rows()-1;
			while(i > newNode){
				Phi.row(i) = Phi.row(i-1);
//				std::cout << Phi << std::endl;
				i--;
			}
			std::cout << Phi << std::endl;

			for(int j = 0; j < Phi.cols();j++){
//				std::cout << (*idxSet1)[newNode] << '\t' << (*idxSet2)[j] << std::endl;
				distance = getDistance((*idxSet1)[newNode],(*idxSet2)[j]);
				Phi(newNode,j) = rbfEval(distance);
			}

			Phi.col(newNode) = Phi.row(newNode);
			std::cout << "AFTER ADJUSTING:\n"  << Phi << std::endl;
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
	defVec_b = Eigen::VectorXd::Zero(n.N_c*m.nDims);

	for(int dim = 0; dim< m.nDims; dim++){
		defVec_b(Eigen::seqN(dim*n.N_c,n.N_m)) = defVec(Eigen::seqN(dim*n.N_m,n.N_m));
		defVec_b(Eigen::seqN(dim*n.N_c+n.N_m,n.N_se+n.N_ss)) = finalDef.col(dim);
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
			ptr = n.cPtr;
			N_m = n.N_c;
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
