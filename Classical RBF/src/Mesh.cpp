#include "Mesh.h"

#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <string>
using Eigen::MatrixXd;
using namespace std;

Mesh::Mesh(std::string fileName, double supportRadius,int debugLvl)
:fName(fileName), nPnts(0), nDims(0),nElem(0),nBdryNodes(0),lvl(debugLvl), r(supportRadius)
{}

void Mesh::readMeshFile(vector<std::string> ibTags,vector<std::string> ebTags){
	if(lvl>=1){
		cout << "Reading mesh file: " << fName << endl;
	}

	int lineNo =0;	// keep track of the line number
	int nIntBdry = 0, nExtBdry = 0;		// ints for keeping track of number of internal/ external boundary points
	bool intBound = false, extBound = false;	// booleans to acknowledge if on int/ext boundary
	int cntInt = 0,cntExt = 0;					// counters for int/ext boundary nodes
	int point = 0;								// counter for number of points
	int pntsIdx;								// int that stores the line where "NPOIN= " is

	string line;	// string containing line obtained by getline() function
	ifstream mFile(fName); 	//opening file name stored in mFile object
	// Check if file is opened
	if (mFile.is_open()){
		//Obtain line
		while (getline(mFile, line)){

			// Following statements find information on various lines of the input file
			if (line.rfind("NDIME= ",0)==0){
				nDims = stoi(line.substr(7));
			}
			else if (line.rfind("NPOIN= ",0)==0){
				nPnts = stoi(line.substr(7));
				pntsIdx = lineNo;
				coords.resize(nPnts, nDims);
			}
			else if (line.rfind("NELEM= ",0)==0){
				nElem = stoi(line.substr(7));
			}

			// Checking whether provided boundary tags equals the amount in the mesh file
			else if (line.rfind("NMARK= ",0)==0 && lvl >= 1){
				try{
					if(stoi(line.substr(7)) == int(ibTags.size() + ebTags.size())){
						cout << "Number of boundary tags are matched" << endl;
					} else{
						throw(stoi(line.substr(7)));
					}
				}
				catch(int nTag){
					cout << "Number of tags provided (" << ibTags.size() + ebTags.size() << ") does not match number of tags found in file (" << nTag << ")." << endl;
					std::exit(0);
				}
			}

			// Finding tags of the boundaries
			else if (line.rfind("MARKER_TAG= ",0)==0){
				string tag =  line.substr(12);

				// Check if found tag is among provided internal boundary tags
				try{
					if(std::find(std::begin(ibTags), std::end(ibTags), tag) != std::end(ibTags)){
						if(lvl >=2){
							cout << "Finding nodes of internal boundary: " << tag << endl;
						}
						intBound = true;
						extBound = false;
					}
					// Check if found tag is among provided external boundary tags
					else if(std::find(std::begin(ebTags), std::end(ebTags), tag) != std::end(ebTags)){
						if(lvl >= 2){
							cout << "Finding nodes of external boundary: " << tag << endl;
						}
						intBound = false;
						extBound = true;
					} else throw(tag);
				}

				catch(string& tag){
					cout << "None of the provided tags match \"" << tag << "\" as found in mesh file" << endl;
					std::exit(0);
				}
			}

			/*
			 * Following piece of code should be adjusted so that all nodes are included in the array
			 * otherwise this will give problems in the 3d cases where the boundary not only consists of lines
			 * but also surfaces
			 */




			// Check whether line corresponds to internal boundary
			else if(intBound){
				istringstream is(line.substr(1)); // split line at '\t' and omitting the first number

				// check if line states the number of elements in the boundary
				if(line.rfind("MARKER_ELEMS= ",0)==0){
					nIntBdry += stoi(line.substr(13));	// updating number of int. boundary points
					intBdryNodes.resize(nIntBdry+1);	// resizing the array containing the int. boundary nodes

				}
				// check if line contains nodes corresponding to the boundary
				else if(isdigit(line[0])&& cntInt < nIntBdry-1){
					is >> intBdryNodes[cntInt];	//assigning first node to int. boundary array
					cntInt++; // update the int. node counter

				// in case the current boundary element is the final element of the boundary.
				// Then both nodes on that line should be included in the int. boundary array
				} else if(isdigit(line[0])&& cntInt < nIntBdry){
					is >> intBdryNodes[cntInt] >> intBdryNodes[cntInt+1];
					cntInt++;
					cout << cntInt << endl;
				}

			}

			else if(extBound){
				istringstream is(line.substr(1));
				if(line.rfind("MARKER_ELEMS= ",0)==0){
					nExtBdry += stoi(line.substr(13));
					extBdryNodes.resize(nExtBdry+1); //might need to use conservativeResize() if problematic
					// resize of the matrix containing the internal boundary nodes
				}
				else if(isdigit(line[0])&& cntExt < nExtBdry-1){
					is >> extBdryNodes[cntExt];
					cntExt++;

				} else if(isdigit(line[0])&& cntExt < nExtBdry){
					is >> extBdryNodes[cntExt] >> extBdryNodes[cntExt+1];
					cntExt++;
				}
			}
			if (lineNo > pntsIdx && point < nPnts){
				istringstream is(line);
				switch(nDims){
					case 2:
						is >> coords(point,0) >> coords(point,1);
						break;
					case 3:
						is >> coords(point,0) >> coords(point,1) >> coords(point,2);
						break;
				}
				point++;
			}
			lineNo++;
		}

		mFile.close();
	}
	else cout << "Not able to open input mesh file";

	nBdryNodes = nIntBdry+nExtBdry;

	 // temp statement for succeeding functions
	sort(begin(intBdryNodes), end(intBdryNodes));
	intBdryNodes = UniqueElems(intBdryNodes);
//	cout << intBdryNodes << endl;

	sort(begin(extBdryNodes), end(extBdryNodes));
	extBdryNodes = UniqueElems(extBdryNodes);

//	cout << extBdryNodes << endl;
//	cout << "ext \n" << extBdryNodes << endl;

//	cout << coords.size()/nDims << endl;
//	cout << nPnts << endl;
//	cout << coords << endl;

	bdryNodes.resize(intBdryNodes.size()+extBdryNodes.size());
	bdryNodes << intBdryNodes, extBdryNodes;

//	cout << bdryNodes << endl;

	sort(begin(bdryNodes), end(bdryNodes)); // I don;t want to sort this array, as i need to find idx again
//	cout << bdryNodes << endl;
	intNodes = obtainIntNodes();

//	cout << intNodes << endl;
//	cout << intNodes << endl;

	cout << nDims << '\t' << nPnts << '\t' << pntsIdx << '\t' << nElem << '\t' <<nBdryNodes << '\t'  <<endl;


	cout << "Mesh file read succesfully" << endl;


	cout << "Boundary nodes \n" << bdryNodes << "\nPoints \n" << coords << endl;

}

Eigen::VectorXi Mesh::UniqueElems(Eigen::ArrayXi arr){
	Eigen::VectorXi a(2*nBdryNodes);
	int cnt = 0;
	for(int x=0;x<arr.size()-1;x++){
		if(arr(x+1)!=arr(x)){
			a(cnt) = arr(x);
			cnt++;
		}
		if(x==arr.size()-2 && arr(x+1)==arr(x)){
			a(cnt) = arr(x);
		}
		else if(x==arr.size()-2 && arr(x+1)!=arr(x)){
			a(cnt) = arr(x+1);
		}

	}
	return a(Eigen::seq(0,cnt));
}

Eigen::ArrayXi Mesh::obtainIntNodes(){
	cout << "determining internal nodes " << endl;
	Eigen::ArrayXi arr(nPnts-intBdryNodes.size()-extBdryNodes.size());
//	cout << arr.size() << endl;


	int cnt = 0;
	int i=0, j=0;
	while(i<nPnts){// && j<(intBdryNodes.size()+extBdryNodes.size())){

		if(j<(intBdryNodes.size()+extBdryNodes.size())){

			if(i == bdryNodes(j)){

				i++; j++;
			}
			else if(i < bdryNodes(j)){
				arr(cnt) = i;
//				cout << arr(cnt) << endl;
				i++; cnt++;

			}
			else if(i > bdryNodes(j)){
				j++;

			}
		}
		else{
			arr(cnt) = i;
//			cout << arr(cnt) << endl;
			cnt++;
			i++;
		}

	}
	return arr;
}

Eigen::MatrixXd Mesh::interpMat(Eigen::ArrayXi idxSet1, Eigen::ArrayXi idxSet2){
	Eigen::MatrixXd Phi(idxSet1.size(), idxSet2.size());
	for(int i=0; i<idxSet1.size();i++){
		for(int j=0; j<idxSet2.size();j++){

			double dist = sqrt(pow(coords(idxSet1(i),0)-coords(idxSet2(j),0),2) + pow(coords(idxSet1(i),1)-coords(idxSet2(j),1),2));
			Phi(i,j) = rbfEval(dist);
		}
	}
	return Phi;
}



double Mesh::rbfEval(double distance){
	double xi = distance/r;	// distance scaled by support radius
	double f_xi = pow((1-xi),4)*(4*xi+1);
	return f_xi;
}

void Mesh::updateNodes(Eigen::VectorXd dxVec,Eigen::VectorXd dyVec, Eigen::VectorXd xDisp,Eigen::VectorXd yDisp){
	newCoords = coords;

	newCoords(bdryNodes,0) = coords(bdryNodes,0) + dxVec;
	newCoords(bdryNodes,1) = coords(bdryNodes,1) + dyVec;
	newCoords(intNodes,0) = coords(intNodes,0) + xDisp;
	newCoords(intNodes,1) = coords(intNodes,1) + yDisp;
//	return newCoords;

}


void Mesh::writeMeshFile(std::string ifName,std::string ofName){
	ofstream outputF;
	outputF.precision(15);		// sets precision of the floats in the file
	outputF.open(ofName, ios::out); // ios::out allows output to file

	ifstream inputF(ifName);
	string mystring;

	bool printFlag = false;
	int cnt = 0;

	while (getline(inputF, mystring)){
		if(printFlag && cnt < nPnts){
			outputF << newCoords(cnt,0)<< '\t' << newCoords(cnt,1) << '\t'<< cnt << endl;
			cnt++;

		} else outputF << mystring << endl;
		if (mystring.rfind("NPOIN= ",0)==0){
			printFlag = true;
		}
	}

	inputF.close();
	outputF.close();
}
