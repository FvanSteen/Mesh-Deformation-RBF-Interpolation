#include "Mesh.h"

#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <algorithm>
using Eigen::MatrixXd;
using namespace std;

Mesh::Mesh(std::string fileName, double supportRadius)
:fName(fileName), pntsIdx(0), nPnts(0), nDims(0),nElem(0),bcNodesIdx(0),nBdryNodes(0), r(supportRadius)
	{
}

void Mesh::findProblemChars(){
	int lineNo =0, bndyPnts=0;
	string line;
	ifstream meshFile(fName);

	if (meshFile.is_open()){
		while (getline(meshFile, line)){

			if (line.rfind("NDIME= ",0)==0){
				nDims = stoi(line.substr(7));
			}
			if (line.rfind("NPOIN= ",0)==0){
				nPnts = stoi(line.substr(7));
				pntsIdx = lineNo;
			}
			if (line.rfind("NELEM= ",0)==0){
				nElem = stoi(line.substr(7));
			}
			if (line.rfind("NMARK= ",0)==0){
//				bcNodesIdx = stoi(line.substr(7));
				bcNodesIdx = lineNo;
//				cout << lineNo << endl;
			}
			if (bcNodesIdx!= 0 && isdigit(line[0])){
				bndyPnts++;

			}
			lineNo++;
		}

		nBdryNodes = bndyPnts;

		meshFile.close();
	}
	else cout << "Not able to open file";


}

void Mesh::obtainCoords(){
	int lineNo =0;
	string line;
	ifstream meshFile(fName);
	int point = 0, bdryPnt =0;

	coords.resize(nPnts, nDims);
	bdryNodes.resize(2*nBdryNodes);
	if (meshFile.is_open()){
		while (getline(meshFile, line)){
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
			if (lineNo > bcNodesIdx && isdigit(line[0])){
				istringstream is( line.substr(1) );
				int t1, t2;
				switch(nDims){
					case 2:
						is >> t1 >> t2;
						bdryNodes[2*bdryPnt] = t1;
						bdryNodes[2*bdryPnt+1] = t2;
				}
				bdryPnt++;
			}

			lineNo++;
		}
		meshFile.close();
	}
	sort(begin(bdryNodes), end(bdryNodes));
	bdryNodes = UniqueElems();

	intNodes = obtainIntNodes();
}

Eigen::VectorXi Mesh::UniqueElems(){
	Eigen::VectorXi a(2*nBdryNodes);

	int cnt = 0;
	if(bdryNodes(0)==bdryNodes(1)){
		a(cnt) = bdryNodes(0);
		cnt++;
	}
	for(int x=1;x<bdryNodes.size();x++){
		if(bdryNodes(x)!=bdryNodes(x-1)){
			a(cnt) = bdryNodes(x);
			cnt++;
		}
	}

	return a(Eigen::seq(0,cnt-1));
}

Eigen::ArrayXi Mesh::obtainIntNodes(){
	Eigen::ArrayXi arr(nPnts-nBdryNodes);
	int cnt = 0;
	int i=0, j=0;
	while(i<nPnts && j<nBdryNodes){
		if(i < bdryNodes(j)){
			arr(cnt) = i;
			cnt++;
			i++;
		} else if(i > bdryNodes(j)){
			j++;
		} else if(i== bdryNodes(j)){
			i++;
			j++;
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

void Mesh::updateNodes(Eigen::VectorXd dVec, Eigen::VectorXd disp){
	newCoords = coords;

	newCoords(bdryNodes,0) = coords(bdryNodes,0) + dVec;
	newCoords(intNodes,0) = coords(intNodes,0) + disp;
//	return newCoords;

}

void Mesh::WriteMeshFile(std::string fName){

	ofstream meshFile;
	meshFile.precision(15);		// sets precision of the floats in the file
	meshFile.open(fName, ios::out); // ios::out allows output to file

	meshFile << "%" << endl;
	meshFile << "% Problem dimension" << endl;
	meshFile << "%" << endl;
	meshFile << "NDIME= " << nDims<< endl;
	meshFile << "%" << endl;
	meshFile << "% Node coordinates" << endl;
	meshFile << "%" << endl;
	meshFile << "NPOIN= " << nPnts << endl;
	for(int i=0; i<nPnts; i++){
		for(int j=0; j<nDims; j++){
			meshFile << newCoords(i,j) << "\t";
		}
		meshFile << i << "\n";
	}
	meshFile.close();
}

