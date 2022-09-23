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

Mesh::Mesh(std::string fileName, double supportRadius)
:fName(fileName), pntsIdx(0), nPnts(0), nDims(0),nElem(0),bcNodesIdx(0),nBdryNodes(0),intNodesIdx(0), r(supportRadius)
{}

void Mesh::findProblemChars(vector<std::string> ibTags,vector<std::string> ebTags){
	int lineNo =0, bndyPnts=0;
	string line;
	ifstream meshFile(fName);
	int ElemIdx;
	int nIntBdry = 0,nExtBdry;
	bool bdry = false, intBound = false, extBound = false;
	int cntInt = 0,cntExt = 0;
	int point = 0;

	if (meshFile.is_open()){
		while (getline(meshFile, line)){

			if (line.rfind("NDIME= ",0)==0){
				nDims = stoi(line.substr(7));
			}
			if (line.rfind("NPOIN= ",0)==0){
				nPnts = stoi(line.substr(7));
				pntsIdx = lineNo;
				coords.resize(nPnts, nDims);
			}
			if (line.rfind("NELEM= ",0)==0){
				nElem = stoi(line.substr(7));
				ElemIdx = lineNo;
			}
			if (line.rfind("NMARK= ",0)==0){
				int nTags = stoi(line.substr(7));
				if(nTags == int(ibTags.size() + ebTags.size())){
//					bcNodesIdx = lineNo;
					cout << "Number of tags are matched" << endl;
					bdry = true;
				} else cout << "Number of tags provided (" << ibTags.size() + ebTags.size() << ") does not match number of tags found in file (" << stoi(line.substr(7)) << ")." << endl;
			}


			if (bdry==true && isdigit(line[0])){

				bndyPnts++;

			}

			if (line.rfind("MARKER_TAG= ",0)==0){
				string tag =  line.substr(12);

				//check if it is an internal boundary
				if(std::find(std::begin(ibTags), std::end(ibTags), tag) != std::end(ibTags)){
					cout << "internal: " <<tag << endl;
					intBound = true;
					extBound = false;
				}
				// check whether its an external boundary
				else if(std::find(std::begin(ebTags), std::end(ebTags), tag) != std::end(ebTags)){
					cout << "external: " << tag << endl;
					intBound = false;
					extBound = true;
				}
				// if not both print out message
				else cout << "Tag not found in mesh file." << endl;

			}
			if(intBound){
				istringstream is(line.substr(1));
				if(line.rfind("MARKER_ELEMS= ",0)==0){
					nIntBdry += stoi(line.substr(13));
					intBdryNodes.resize(nIntBdry+1);
					// resize of the matrix containing the internal boundary nodes
				}
				else if(isdigit(line[0])&& cntInt < nIntBdry-1){
					is >> intBdryNodes[cntInt];
					cntInt++;


				} else if(isdigit(line[0])&& cntInt < nIntBdry){
					is >> intBdryNodes[cntInt] >> intBdryNodes[cntInt+1];
					cntInt++;
				}

			}

			if(extBound){
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

		nBdryNodes = bndyPnts;

		meshFile.close();
	}
	else cout << "Not able to open file";

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

//	cout << nDims << '\t' << nPnts << '\t' << pntsIdx << '\t' << nElem << '\t' << bcNodesIdx << '\t' <<nBdryNodes << '\t' <<intNodesIdx <<endl;
}

void Mesh::obtainCoords(){
	int lineNo =0;
	string line;
	ifstream meshFile(fName);
	int point = 0, bdryPnt =0;

	coords.resize(nPnts, nDims);
	bdryNodes.resize(2*nBdryNodes);
	intBdryNodes.resize(2*4);
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
			cout << "coords are assigned" << endl;
			if (lineNo > bcNodesIdx && isdigit(line[0])){
				istringstream is( line.substr(1) );
				int t1, t2;
				switch(nDims){
					case 2:
						is >> t1 >> t2;
						bdryNodes[2*bdryPnt] = t1;
						bdryNodes[2*bdryPnt+1] = t2;
						if (lineNo > intNodesIdx){
							intBdryNodes[2*(bdryPnt-20)] = t1;
							intBdryNodes[2*(bdryPnt-20)+1] = t2;
						}
				}
				bdryPnt++;
			}
			cout << "boundary points are assigned" << endl;

			lineNo++;
		}
		meshFile.close();
	}

	sort(begin(bdryNodes), end(bdryNodes));
	bdryNodes = UniqueElems(bdryNodes);


	sort(begin(intBdryNodes), end(intBdryNodes));
	intBdryNodes = UniqueElems(intBdryNodes);

	intNodes = obtainIntNodes();
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


//		else if(bdryNodes(j) != i && j<(intBdryNodes.size()+extBdryNodes.size())){
//			cout << "hier " << j << endl;
//			j++;
//		}
//		else if(bdryNodes(j) != i){
//			cout << "of hier " << endl;
//			arr(cnt) = i;
//			cnt++; i++;
//		}

//		else if (bdryNodes(j) == i){
//			i++;
//			cout << i << ' ' << j << endl;
//		}

//		if(i < bdryNodes(j)){
//			arr(cnt) = i;
//			cout<<i<<endl;
//			cnt++;
//			i++;
//		} else if(i > bdryNodes(j)){
//			j++;
//		} else if(i== bdryNodes(j)){
//			i++;
//			j++;
//		}

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


void Mesh::wmf(std::string ifName,std::string ofName){
	ofstream outputF;
	outputF.precision(15);		// sets precision of the floats in the file
	outputF.open(ofName, ios::out); // ios::out allows output to file


	ifstream inputF(ifName);
	string mystring;

	bool printFlag = false;
	int cnt = 0;
//	inputF >> mystring;
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
