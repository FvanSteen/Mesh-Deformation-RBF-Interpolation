
#ifndef MESH_H_
#define MESH_H_

class Mesh {
public:
	Mesh(int nNodes, int mNodes, double hDomain, double wDomain);
	void writeMeshFile(double x[],double y[]);
	void getMeshCoor();

private:
	int n, m;
	double h,w;


};
#endif /* MESH_H_ */
