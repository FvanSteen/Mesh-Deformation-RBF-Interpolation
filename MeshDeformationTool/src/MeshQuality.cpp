#include "MeshQuality.h"
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <string>


MeshQuality::MeshQuality(probParams& p, Eigen::ArrayXXd& coords)
{

	nDim = coords.cols();
	if(p.generateQuality){

		fName = p.directory + "\\meshQuals\\" +  p.mesh_ifName.substr(0,p.mesh_ifName.size()-4) + "_qual.txt";

		alpha_ptr = &alpha_init;

		getInitialMeshQualParams(p, coords);



		bool existing = existTest(fName);

		if(!existing){

			getMeshQual();
			std::cout << "done\n";
			writeQualFile();
		}
		alpha_ptr = &alpha;
		fName = p.directory + "\\meshQuals\\" + p.mesh_ofName.substr(0,p.mesh_ofName.size()-4) + "_qual.txt";
	}
}

bool MeshQuality::existTest(std::string& fName){
	std::ifstream file(fName);
	return file.good();
}

void MeshQuality::getInitialMeshQualParams(probParams& p, Eigen::ArrayXXd& coords){

	getElemConnectivity(p);

	getQualParams(coords);
}

void MeshQuality::getDeformedMeshQual(Eigen::ArrayXXd& coords){
	getQualParams(coords);
	getMeshQual();
	writeQualFile();
}

void MeshQuality::getMeshQual(){
	Eigen::ArrayXd tau(nElem), f_size(nElem), f_skew(nElem);
	tau = (*alpha_ptr).rowwise().sum()/alpha_init.rowwise().sum();


	for(int i = 0; i < nElem; i++){
		f_size(i) = std::min(tau[i], 1/tau[i]);
		f_skew(i) = get_f_skew(elemType[i], i);
	}

	qual = sqrt(f_size)*f_skew;

	defQuals << qual.minCoeff(), qual.mean(), qual.maxCoeff();


}

double MeshQuality::get_f_skew(int type, int i){
	double f_skew;
	switch(type){

		case 5:
			f_skew = sqrt(3)*(*alpha_ptr)(i,0)/(lambda_11(i,0) + lambda_22(i,0) - lambda_12(i,0));
			break;

		case 9: {
			Eigen::ArrayXd ll_1(4), ll_2(4), alpha(4);
			ll_1 = lambda_11(i,Eigen::seqN(0,4));
			ll_2 = lambda_22(i,Eigen::seqN(0,4));
			alpha = (*alpha_ptr)(i,Eigen::seqN(0,4));

			if(quadClockWise){
				f_skew = -4/(sqrt(ll_1*ll_2)/alpha).sum();
			}else{
				f_skew = 4/(sqrt(ll_1*ll_2)/alpha).sum();
			}
		}break;

		case 10:
			f_skew = 3. * pow(((*alpha_ptr)(i,0) *sqrt(2.)), 2.0/3.0) / (1.5 * (lambda_11(i,0) + lambda_22(i,0) + lambda_33(i,0)) - (lambda_12(i,0) + lambda_23(i,0) + lambda_13(i,0)));
		break;

		case 12: {
			Eigen::ArrayXd ll_1(8), ll_2(8), ll_3(8), alpha(8);
			ll_1 = lambda_11(i,Eigen::seqN(0,8));
			ll_2 = lambda_22(i,Eigen::seqN(0,8));
			ll_3 = lambda_33(i,Eigen::seqN(0,8));

			alpha = (*alpha_ptr)(i,Eigen::seqN(0,8));

			f_skew = 8/ (sqrt(ll_1*ll_2*ll_3)/alpha).sum();
		}break;
	}

	return f_skew;
}


void MeshQuality::getQualParams(Eigen::ArrayXXd& coords){

	Eigen::MatrixXd A(nDim,nDim);
	Eigen::ArrayXXd tensor(nDim,nDim);
	int cols;
	for(int i = 0; i < nElem; i++){

		if( i == 0 || elemType[i] != elemType[i-1]){

			cols  = setElemTypeParams(elemType[i]);
			if(elemType[i] == 9 && !setQuadDir){
				getQuadRotation(i, coords);
			}
		}
		for(int ii = 0; ii < cols; ii++){
			A = coords(elems(i,k.col(ii)), Eigen::all);

			A = A.array().rowwise() -  coords(elems(i,ii), Eigen::all);

			(*alpha_ptr)(i,ii) = A.determinant();

			tensor = A*A.transpose();



			lambda_11(i,ii) = tensor(0,0);
			lambda_22(i,ii) = tensor(1,1);

			if(elemType[i] == 5){
				lambda_12(i,ii) = tensor(0,1);
			}

			if(elemType[i] == 10){
				lambda_12(i,ii) = tensor(0,1);
				lambda_23(i,ii) = tensor(1,2);
				lambda_13(i,ii) = tensor(0,2);
			}

			if(nDim == 3){
				lambda_33(i,ii) = tensor(2,2);
			}

		}
	}

}

int MeshQuality::setElemTypeParams(int type){
	int cols;
	if(nDim == 2){
		if(type == 5){
			k.resize(2,1);
			k.row(0) << 1;
			k.row(1) << 2;

		if((*alpha_ptr).cols() < 1){
			(*alpha_ptr).conservativeResize(nElem,1);
			lambda_11.conservativeResize(nElem,1);
			lambda_22.conservativeResize(nElem,1);
		}
		if(lambda_12.cols() < 1){
			lambda_12.conservativeResize(nElem,1);
		}

		cols = 1;
		}else if(type == 9){

			k.resize(2,4);
			k.row(0) << 1,2,3,0; // k+1
			k.row(1) << 3,0,1,2; // k+3

		if((*alpha_ptr).cols() < 4){
			(*alpha_ptr).conservativeResize(nElem,4);
			lambda_11.conservativeResize(nElem,4);
			lambda_22.conservativeResize(nElem,4);
		}
		cols = 4;
		}
	}else if(nDim == 3){
		if(type == 12){
			k.resize(3,8);
			k.row(0) << 1,2,3,0,7,4,5,6;  // k + 1
			k.row(1) << 3,0,1,2,5,6,7,4;  // k + 3
			k.row(2) << 4,5,6,7,0,1,2,3;  // k + 4

			if((*alpha_ptr).cols() != 8){
				(*alpha_ptr).conservativeResize(nElem,8);
				lambda_11.conservativeResize(nElem,8);
				lambda_22.conservativeResize(nElem,8);
				lambda_33.conservativeResize(nElem,8);
			}

			cols = 8;

		}else if(type == 10){
			k.resize(3, 1);
			k << 1,2,3;


			if((*alpha_ptr).cols() < 1){
				(*alpha_ptr).conservativeResize(nElem,1);
				lambda_11.conservativeResize(nElem,1);
				lambda_22.conservativeResize(nElem,1);
				lambda_33.conservativeResize(nElem,1);
			}

			if(lambda_12.cols() < 1){
				lambda_12.conservativeResize(nElem,1);
				lambda_23.conservativeResize(nElem,1);
				lambda_13.conservativeResize(nElem,1);
			}

			cols = 1;
		}

	}

	return cols;

}

void MeshQuality::writeQualFile(){
	std::ofstream outputF;		// Making an output stream class to operate on the output file
	outputF.precision(15);		// sets precision of the floats in the file

	outputF.open(fName);
	for(int i = 0; i < nElem; i++){
		outputF << qual[i] << std::endl;
	}

	outputF.close();
	std::cout << "Generated mesh quality file: " << fName << std::endl;

}

void MeshQuality::getElemConnectivity(probParams& p){
	std::string line;							// string containing line obtained by getline() function
	std::ifstream mFile(p.directory + "\\Meshes\\" + p.mesh_ifName); 	//opening file name stored in mFile object

	int startIdx;
	int lineNo = 0;
	bool elemData = false;
	if (mFile.is_open()){
		//Obtain line
		while (getline(mFile, line)){

			if (line.rfind("NELEM= ",0)==0){						// save nr of elements
				nElem = stoi(line.substr(7));
				startIdx = lineNo+1;
				elemData = true;

				elemType.resize(nElem);

			}
			else if(elemData && lineNo < startIdx + nElem){

				std::istringstream is(line);
				int data, lineElem = 0;

				while(is >> data){
					if(lineElem == 0){

						elemType(lineNo-startIdx) = data;

						if(data == 12 && elems.cols() < 8){
							elems.conservativeResize(nElem, 8);
						}else if((data == 9 || data == 10) && elems.cols() < 4){
							elems.conservativeResize(nElem, 4);
						}else if(data == 5 && elems.cols() < 3){
							elems.conservativeResize(nElem,3);
						}

					}else if(lineElem-1 < elems.cols()) {
						elems(lineNo-startIdx,lineElem-1) = data;
					}
					lineElem++;
				}
			}
			lineNo++;
		}
	}
}


void MeshQuality::getQuadRotation(int i, Eigen::ArrayXXd& coords){

	// the coordinates of the nodes of a single quadriliteral element
	Eigen::ArrayXXd elemCoords;
	elemCoords = coords(elems.row(i), Eigen::all);

	// defining the center of the element
	Eigen::ArrayXd midpnt;
	midpnt = elemCoords.colwise().sum()/elemCoords.rows();

	// angle of the first defined node w.r.t. the middle of the element
	double angle_0 = atan2(elemCoords(0,1)-midpnt(1), elemCoords(0,0)-midpnt(0));

	// angle of the second defined node w.r.t. the middle of the element
	double angle_1 = atan2(elemCoords(1,1)-midpnt(1), elemCoords(1,0)-midpnt(0));

	// setting the rotation in which the nodes are defined
	if(angle_1 > angle_0){
		quadClockWise = false;
	}else{
		quadClockWise = true;
	}

	// setting bool to prevent recalling of this function
	setQuadDir = true;
}
