#ifndef PROBPARAMS_H_
#define PROBPARAMS_H_
#include <Eigen/Dense>
#include <string>
#include <vector>
struct probParams{
	int steps;
	std::string smode;
	std::string pmode;
	int ptype;
	bool curved;
	int pDir;
	bool dataRed;
	double tol;
	std::string dispFile;
	double gamma;
	bool multiLvl;
	int lvlSize; // delete
	int lvlSizeInit; // delete
	std::string convHistFile;
	std::string mesh_ifName, mesh_ofName;
	std::vector<std::string> bdryTags, mTags, pTags;
	double rFac;
	bool doubleEdge;
	double tolCrit;
};


#endif /* PROBPARAMS_H_ */
