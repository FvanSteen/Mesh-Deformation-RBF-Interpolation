
#include <iostream>
#include "Mesh.h"
#include <fstream>



using namespace std;


Mesh::Mesh(int nNodes, int mNodes, double hDomain, double wDomain)
:n(nNodes), m(mNodes), h(hDomain), w(wDomain)
{
	cout << "Members are initialized" << endl;
}



void Mesh::writeMeshFile(double x[],double y[]){
	ofstream meshFile;
	meshFile.precision(15);
	meshFile.open("mesh.su2", ios::out);

	meshFile << "%" << endl;
	meshFile << "% Node coordinates" << endl;
	meshFile << "%" << endl;
	meshFile << "NPOIN= " << n*m << endl;

	meshFile << x[0] << "\t" << y[0];
//	cout << x << " "<< y << endl;
	meshFile.close();

}

void Mesh::getMeshCoor(){
	int i, j, node = 0;
	double x[n*m], y[n*m];
	for(j=0; j<m; j++){
		for(i=0; i<n; i++){
//			cout << (double(i)/double(n-1))*w << "\t" << (double(j)/double(m-1))*h << "\t" << node << endl;
			x[node] = (double(i)/double(n-1))*w;

			y[node] = (double(j)/double(m-1))*h;
			node++;
		}
	}
	Mesh::writeMeshFile(x,y);
}
