#ifndef MESH_H_
#define MESH_H_
//#include "ReadConfigFile.h"
#include <string>
#include <Eigen/Dense>
#include <vector>
#include "probParams.h"
class Mesh
{
public:
//	Mesh(const std::string& inputFileName,const std::string& outputFileName, const std::vector<std::string>& Tags,const double& rFac,const int& debugLvl, const std::string& slidingMode, const std::string& periodicMode,const std::vector<std::string>& movingTags, const std::vector<std::string>& periodicTags);
	Mesh(probParams& params, const int& debugLvl);

	int nNodes, nDims;

	Eigen::ArrayXXd coords;

	Eigen::ArrayXXd surfMidPnts, edgeMidPnts;

	Eigen::ArrayXXd n, t;
	Eigen::ArrayXXd t_se, n1_se, n2_se, n_ss,t1_ss,t2_ss;


	Eigen::ArrayXXd surfMidPntNormals, edgeMidPntNormals1,edgeMidPntNormals2;
	Eigen::ArrayXi mNodes,seNodes,iNodes, ssNodes;

	Eigen::ArrayXi intCorNodes;

	Eigen::ArrayXi verticesNodes, periodicVerticesNodes;

	Eigen::ArrayXi periodicEdgeNodes;

	int N_pe;

//	struct Domains{
//		Eigen::ArrayXi d1;
//		Eigen::ArrayXi d2;
//		Eigen::ArrayXi d3;
//		Eigen::ArrayXi d4;
//	};
//
//	struct Domains subDoms;


	double r,lambda; // support radius
	int N_i, N_nonzeroDisp, N_se, N_ss, N_m, N_periodic_vertices;







	void readMeshFile(probParams& params);

	void getNodeTypes(probParams& params, int nMarker);
	void getIntNodes();
	void removeDuplicates(Eigen::ArrayXi& arr);

	void getEdgeConnectivity(std::string& pmode, Eigen::ArrayXXi& edgeConnectivity,Eigen::ArrayXi& seNodes, int size);
	void getSurfConnectivity();

	void getVecs();
	void getPerpVecs(Eigen::ArrayXXd& vecs, Eigen::ArrayXXd& p1, Eigen::ArrayXXd& p2);
	void getEdgeTan(Eigen::ArrayXXd& t,Eigen::ArrayXXi& edgeConnectivity, Eigen::ArrayXi& seNodes);
	void getSurfNormal();
	void getSurfNormalPeriodic();


	void getExtBdryEdgeSegments();
	void getMidPnts(probParams& params);

	void getSubDomains(Eigen::ArrayXi& subDomains, Eigen::ArrayXi& subDomLen, Eigen::ArrayXi& subDomBdry, Eigen::ArrayXi& subDomBdryLen);
	void getIntCorNodes(double& gamma, double& tol);


	void getCharPerLength(std::string& pDir);
	double getCharDomLength();

	void writeMeshFile(std::string& ifName, std::string& ofName);
	void findStringBounds(int& first, int& last, std::string& line);
//	void getInternalCorrectionNodes(Eigen::ArrayXi& subDomains,  Eigen::ArrayXi& subDomsLen, Eigen::ArrayXXd& bdryCoord);
private:
	Eigen::ArrayXXi edgeConnectivity, edgeConnectivityPeriodic, surfConnectivity, bdryNodesMat, extBdryEdgeSegments;
	Eigen::ArrayXi nrElemsBdry;
	std::vector<std::string> srtdTags;
	int lvl;
};

#endif /* MESH_H_ */
