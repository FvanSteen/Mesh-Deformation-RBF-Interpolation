#ifndef PROBPARAMS_H_
#define PROBPARAMS_H_
#include <Eigen/Dense>
#include <string>
#include <vector>
struct probParams{
	int steps, ptype, pDir, lvlSize;

	bool curved, dataRed, multiLvl, doubleEdge, generateQuality;

	double tol, gamma, rFac,tolCrit;

	std::string smode, pmode, dispFile, mesh_ifName, mesh_ofName, mCrit, directory;

	std::vector<std::string> bdryTags, mTags, pTags, iTags;
};


#endif /* PROBPARAMS_H_ */
