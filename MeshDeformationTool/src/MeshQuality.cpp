#include "MeshQuality.h"
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <string>

/*
 * This class is able to determine the mesh quality of the mesh elements
 * This is done by using a set of algebraic mesh quality metrics as introduced by P.M. Knupp
 * See section 2.7 of the thesis report for more details
 */
MeshQuality::MeshQuality(probParams& p, Eigen::ArrayXXd& coords)
{
	// number of coordinates
	nDim = coords.cols();

	// in case a quality has to be generated
	if(p.generateQuality){
		// setting filename for the quality file of the input mesh
		fName = p.directory + "\\" +  p.mesh_ifName.substr(0,p.mesh_ifName.size()-4) + "_qual.txt";

		// set the pointer of the determinants of the Jacobian matrices
		alpha_ptr = &alpha_init;

		// Determine the mesh parameters of the initial mesh
		getInitialMeshQualParams(p, coords);

		// check if a mesh quality file already exists for the input mesh file
		bool existing = existTest(fName);

		// if it does not exist then create one
		if(!existing){
			// determine mesh quality
			getMeshQual();
			// write mesh quality file
			writeQualFile();
		}

		// set pointer to the alpha of the deformed mesh
		alpha_ptr = &alpha;
		// set mesh quality filename for the output mesh file
		fName = p.directory + "\\" + p.mesh_ofName.substr(0,p.mesh_ofName.size()-4) + "_qual.txt";
	}
}

/*
 * Check if mesh quality file exists in the directory
 */

bool MeshQuality::existTest(std::string& fName){
	std::ifstream file(fName);
	return file.good();
}


/*  getInitialMeshQualParams
 *
 * Establish the connectivity of the elements and finding the tensor elements and volumes of the undeformed mesh
 */

void MeshQuality::getInitialMeshQualParams(probParams& p, Eigen::ArrayXXd& coords){

	// establishing element connectivity
	getElemConnectivity(p);
	// Determining the volume and elements of the tensor elements of the mesh element using the Jacobian matrices
	getQualParams(coords);
}

/* getDEformedMeshQual
 *
 * Calling functions to obtain the mesh quality of the deformed mesh
 *
 */

void MeshQuality::getDeformedMeshQual(Eigen::ArrayXXd& coords, int i){
	getQualParams(coords);
	getMeshQual();
	writeQualFile();
}

/* getMeshQual
 * Computing the mesh quality using the size and skew metrics. See equations 2.43 through 2.49 of the manuscript
 */

void MeshQuality::getMeshQual(){
	Eigen::ArrayXd tau(nElem), f_size(nElem), f_skew(nElem);

	// ratio of the deformed and undeformed mesh element. See equation 2.43 of the manuscript
	tau = (*alpha_ptr).rowwise().sum()/alpha_init.rowwise().sum();

	// determine the size and skew quality metrics. See equations 2.44 and 2.45 through 2.48 of the manuscript.
	for(int i = 0; i < nElem; i++){
		f_size(i) = std::min(tau[i], 1/tau[i]);

		f_skew(i) = get_f_skew(elemType[i], i);
	}

	// computing quality
	qual = sqrt(f_size)*f_skew;

	// storing min, mean and max quality
	defQuals << qual.minCoeff(), qual.mean(), qual.maxCoeff();
}


/* get_f_skew
 *
 * Computing the skew quality metric. using equations 2.45 through 2.48 of the manuscript
 * A check has been done to see if the nodes of the elements are ordered by the clockwise or anticlockwise
 * In case the nodes are ordered anticlockwise, a minus has to be included to obtain a positive element volume
 *
 */

double MeshQuality::get_f_skew(int type, int i){
	double f_skew = 0.0;
	switch(type){

		// finding skew metric for triangular elements
		case 5:
			if(triClockWise){
				f_skew = -sqrt(3)*(*alpha_ptr)(i,0)/(lambda_11(i,0) + lambda_22(i,0) - lambda_12(i,0));
			}else{
				f_skew = sqrt(3)*(*alpha_ptr)(i,0)/(lambda_11(i,0) + lambda_22(i,0) - lambda_12(i,0));
			}
			break;

		// finding skew metric for quadriliteral elements
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

		// finding skew metric for tetrahedral elements
		case 10:
			f_skew = 3. * pow(((*alpha_ptr)(i,0) *sqrt(2.)), 2.0/3.0) / (1.5 * (lambda_11(i,0) + lambda_22(i,0) + lambda_33(i,0)) - (lambda_12(i,0) + lambda_23(i,0) + lambda_13(i,0)));
		break;


		// finding skew metric for hexahedral elements
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


/* getQualParams
 *
 * Determining the Jacobian matrices, the element volume and tensors
 *
 */

void MeshQuality::getQualParams(Eigen::ArrayXXd& coords){

	Eigen::MatrixXd A(nDim,nDim);
	Eigen::ArrayXXd tensor(nDim,nDim);
	int cols;
	for(int i = 0; i < nElem; i++){

		// if new type of 2D element is encountered then its ordering is determined
		if( i == 0 || elemType[i] != elemType[i-1]){
			cols  = setElemTypeParams(elemType[i]);
			if(elemType[i] == 9 && !setQuadDir){
				getRotation(i, coords, 9);
			}else if(elemType[i] == 5 && !setTriDir){
				getRotation(i, coords, 5);
			}
		}

		// finding Jacobian matrices A, element volume alpha and tensor elements lambda
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

/* setElemTypeParams
 *
 * Setting indices for determining the Jacobian matrices of the elements see section 2.7 of the manuscript for a detailed description
 *
 */

int MeshQuality::setElemTypeParams(int type){
	int cols = 0;
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

/* writeQualFile
 *
 * Writing of mesh quality file
 *
 */


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

/* getElemConnectivity
 *
 * This function established the element connectivity from the input mesh file
 */

void MeshQuality::getElemConnectivity(probParams& p){
	std::string line;											// string containing line obtained by getline() function
	std::ifstream mFile(p.directory + "\\" + p.mesh_ifName); 	//opening file name stored in mFile object

	int startIdx = -1;		// index at which the element data is stored in the mesh file
	int lineNo = 0;			// line number counter

	if (mFile.is_open()){
		//Obtain line
		while (getline(mFile, line)){
			// search for line containing the number of elements
			if (line.rfind("NELEM= ",0)==0){
				nElem = stoi(line.substr(7));				// save nr of elements
				startIdx = lineNo+1;						// starting index of the elements
				elemType.resize(nElem);						// resizing array containing the element types

			}

			// if line contains element data then store it
			else if(lineNo >= startIdx && lineNo < startIdx + nElem){

				std::istringstream is(line);
				int data, lineElem = 0;

				while(is >> data){
					if(lineElem == 0){
						// first row contains the element type
						elemType(lineNo-startIdx) = data;

						// resize the elems array according to the type of element
						if(data == 12 && elems.cols() < 8){
							elems.conservativeResize(nElem, 8);
						}else if((data == 9 || data == 10) && elems.cols() < 4){
							elems.conservativeResize(nElem, 4);
						}else if(data == 5 && elems.cols() < 3){
							elems.conservativeResize(nElem,3);
						}

					}
					// store the nodes of the element
					else if(lineElem-1 < elems.cols()) {
						elems(lineNo-startIdx,lineElem-1) = data;
					}
					lineElem++;
				}
			}
			lineNo++;
		}
	}
}


/* getRotation
 *
 * Function for determining whether an element is defined clockwise or anticlockwise
 */

void MeshQuality::getRotation(int i, Eigen::ArrayXXd& coords, int type){
	// the coordinates of the nodes of a single quadriliteral element
	Eigen::ArrayXXd elemCoords;

	// store coords of a single element
	if(type == 5)
		//triangular
		elemCoords = coords(elems(i,Eigen::seqN(0,3)), Eigen::all);
	else
		// quadriliteral
		elemCoords = coords(elems.row(i), Eigen::all);


	// sum over the edges of the element of (x_i*y_i+1 - x_i+1*y_i) (the sum equals twice the area of the polygon in question)
	// if the sum is positive then the nodes are order clockwise, else anticlockwise
	double sum = 0.0;
	for(int i = 0; i < elemCoords.rows()-1; i++){
		sum += elemCoords(i,0)*elemCoords(i+1,1) - elemCoords(i,1)*elemCoords(i+1,0);
	}
	sum += elemCoords(elemCoords.rows()-1,0)*elemCoords(0,1) - elemCoords(elemCoords.rows()-1,1)*elemCoords(0,0);


	// setting the rotation in which the nodes are defined
	if(type == 5){
		if(sum > 0)
			triClockWise = false;
		else
			triClockWise = true;

		// setting bool to prevent recalling of this function
		setTriDir = true;

	}else if(type == 9){
		if(sum > 0 )
			quadClockWise = false;
		else
			quadClockWise = true;

		// setting bool to prevent recalling of this function
		setQuadDir = true;

	}
}
