#include "rbfps.h"

#include <iostream>
#include <chrono>
#include <Eigen/Dense>
#include "greedy.h"

//rbf_ps::rbf_ps(Mesh& meshObject, Eigen::VectorXd& dVec, Eigen::RowVectorXd& rotPnt, Eigen::VectorXd& rotVec, const int& steps, const std::string& smode, const bool& curved, const std::string& pDir)
rbf_ps::rbf_ps( struct probParams& probParamsObject, Mesh& meshObject, getNodeType& n)
:rbfGenFunc(meshObject, probParamsObject)
{
	if(params.dataRed){
		greedy g(m, params, exactDisp, movingIndices, alpha, d);
		perform_rbf(n,g);
	}else{
		perform_rbf(n);
	}

}

void rbf_ps::perform_rbf(getNodeType& n){
	std::cout<< "Performing RBF PS without data reduction" << std::endl;

	Eigen::VectorXd defVec, defVec_b;
	Eigen::ArrayXXd delta, finalDef;


	Eigen::VectorXd defVecS(n.N_m+n.N_se);
	Eigen::ArrayXXd deltaS, finalDefS;

	for(int i = 0; i < params.steps; i++){
		std::cout << "Deformation step: " << i+1 << " out of "<< params.steps << std::endl;

		if(params.curved || i==0){
			m.getMidPnts(params);
			m.getVecs(params);
		}

		getPhis(n, 0);

		if(i == 0){
			getDefVec(defVec, n.N_m, n.mPtr);
		}

		// find the finaldef of the edge nodes
		performRBF_PS(PhiPtr, defVec, delta, finalDef, n, n.N_se, n.N_m, 0, PhiPtr->Phi_mm, PhiPtr->Phi_em);

		// make a second deformation vector that also contains the displacement of the edges
		if(m.nDims  == 3){
			getDefVec(defVecS, defVec, n, finalDef, n.N_m+n.N_se, n.N_m);
			performRBF_PS(PhiPtr, defVecS, deltaS, finalDefS, n, n.N_ss, n.N_m+n.N_se, 1, PhiPtr->Phi_meme, PhiPtr->Phi_sme);
			getDefVec(defVec_b, defVecS, n, finalDefS, n.N_c, n.N_m+n.N_se);
		}else{
			getDefVec(defVec_b, defVec, n, finalDef, n.N_c, n.N_m); // check if right?
		}


		performRBF(PhiPtr->Phi_cc,PhiPtr->Phi_ic,defVec_b,n.cPtr,n.iPtr,n.N_c);

	}
}

void rbf_ps::perform_rbf(getNodeType& n, greedy& g){
	std::cout<< "Performing RBF PS with data reduction" << std::endl;

	Eigen::VectorXd defVec, defVec_b;
	Eigen::ArrayXXd delta, finalDef;


	int iter, lvl;
	bool iterating = true;


	for(int i = 0; i < params.steps; i++){
		std::cout << "Deformation step: " << i+1 << " out of "<< params.steps << std::endl;

		iter = 0;
		lvl = 0;

		if(params.curved || i==0){
			m.getMidPnts(params);
			m.getVecs(params);
		}

		while(iterating){

			n.addControlNodes(g.maxErrorNodes, params.smode, m);


			std::cout << "\n\n errornodes:\n" <<  g.maxErrorNodes << "\n" << std::endl;
			std::cout << "\n added node types:\n" << n.addedNodes.type << "\n" << std::endl;
			getPhis(n, iter);
//			if(Phis.Phi_cc.cols() >= 8){
//				std::cout << Phis.Phi_cc << std::endl;
//				std::exit(0);
//			}
			if(lvl > 0){
				getDefVec(defVec_b, n, g.errorPrevLvl, n.N_c);
			}else{
				getDefVec(defVec, n.N_m, n.mPtr);
			}

			if(lvl!=0){
				performRBF(Phis.Phi_cc, Phis.Phi_ic, defVec_b, n.cPtr, n.iPtr, n.N_c);
			}else{
				// finding displacement of the edges and moving nodes
				performRBF_PS(PhiPtr, defVec, delta, finalDef, n, n.N_se, n.N_m, 0, PhiPtr->Phi_mm, PhiPtr->Phi_em);

				getDefVec(defVec_b, defVec, n, finalDef, n.N_m+n.N_se, n.N_m);

//				std::cout << "\n" << defVec_b<< "\n" << std::endl;

				Eigen::MatrixXd p_cc = PhiPtr->Phi_cc.block(0,0,n.N_m+n.N_se,n.N_m+n.N_se);
				Eigen::MatrixXd p_ic = PhiPtr->Phi_ic.block(0,0,n.N_i- m.N_ss + n.N_ss, n.N_m+n.N_se);

				Eigen::ArrayXi iNodesReduced = (*n.iPtr)(Eigen::seqN(0,n.N_i- m.N_ss + n.N_ss));

				Eigen::ArrayXi* iPtrReduced = &iNodesReduced;
//				performRBF(p_cc,p_ic,defVec_b,n.cPtr,n.iPtr,n.N_c);

				performRBF(p_cc,p_ic,defVec_b,n.cPtr,iPtrReduced,n.N_m+n.N_se);

//				std::cout << "performed rbf" << std::endl;
//				std::cout << d << std::endl;
//				std::cout << "rows: " <<  d.rows() << std::endl;
			}

			g.getError(n,d, lvl);
			std::cout << "error: \t"<< g.maxError <<" at node: \t" << g.maxErrorNodes(0)<< std::endl;

			if(g.maxError < params.tol){
				iterating = false;
				if(params.multiLvl == false){
					g.maxErrorNodes.resize(0);
				}
			}

			if(params.multiLvl && (g.maxError/g.maxErrorPrevLvl < params.tolCrit || iterating == false)){

				g.setLevelParams(n,lvl, d, alpha, defVec_b, n.cPtr, n.N_c);

				std::cout << "LEVEL: " << lvl << " HAS BEEN DONE" << std::endl;
				lvl++;
				iter = -1;

				n.assignNodeTypesGrdy(m);

				if(iterating == false){
					g.getAlphaVector();
					g.setInitMaxErrorNodes();
				}

			}
			iter++;

		}
		std::cout << "OUTSIDE EDGE\n";



		int iter_surf = 0;
		Eigen::VectorXd defVecS(n.N_m + n.N_se);
		Eigen::ArrayXXd deltaS, finalDefS;
		if(m.nDims == 3){
			iterating = true;
			defVecS = defVec_b;

			performRBF_PS(PhiPtr, defVecS, deltaS, finalDefS, n, n.N_ss, n.N_m+n.N_se, 1, PhiPtr->Phi_meme, PhiPtr->Phi_sme);
			getDefVec(defVec_b, defVecS, n, finalDefS, n.N_c, n.N_m+n.N_se);
			performRBF(PhiPtr->Phi_cc,PhiPtr->Phi_ic,defVec_b,n.cPtr,n.iPtr,n.N_c);
//				std::cout << "DONE surface projection\n";


			g.getError(n,d,lvl);
		}

		while(iterating){


			n.addControlNodes(g.maxErrorNodes, params.smode, m);
			std::cout << "\n added node types:\n" << n.addedNodes.type << "\n" << std::endl;

			getPhis(n, iter);

			if(n.addedNodes.type.minCoeff() == 0){

				getDefVec(defVec, n.N_m, n.mPtr);
			}


			if(n.addedNodes.type.minCoeff() <= 1){

				performRBF_PS(PhiPtr, defVec, delta, finalDef, n, n.N_se, n.N_m, 0, PhiPtr->Phi_mm, PhiPtr->Phi_em);
				getDefVec(defVecS, defVec, n, finalDef, n.N_m+n.N_se, n.N_m);

			}
//			std::cout << "\n" << defVecS << std::endl;
//			getDefVec(defVec, n.N_m, n.mPtr);

			performRBF_PS(PhiPtr, defVecS, deltaS, finalDefS, n, n.N_ss, n.N_m+n.N_se, 1, PhiPtr->Phi_meme, PhiPtr->Phi_sme);
			getDefVec(defVec_b, defVecS, n, finalDefS, n.N_c, n.N_m+n.N_se);
			performRBF(PhiPtr->Phi_cc,PhiPtr->Phi_ic,defVec_b,n.cPtr,n.iPtr,n.N_c);

			g.getError(n,d,lvl);
//			std::cout << (g.error(Eigen::seqN(0,m.N_m-n.N_m+m.N_se-n.N_se), Eigen::all)).colwise().sum() << std::endl;
			std::cout << "error: \t"<< g.maxError <<" at node: \t" << g.maxErrorNodes(0)<< std::endl;


			if(g.maxError < params.tol){
				iterating = false;
				g.maxErrorNodes.resize(0);
//				std::cout << "error tol reached\n";
//
//				m.coords(*n.iPtr,Eigen::all) += (d);
//				for(int dim = 0 ; dim < m.nDims; dim++){
//					m.coords(*n.cPtr, dim) += (defVec_b(Eigen::seqN(dim*n.N_c,n.N_c))).array();
//				}

			}
			iter_surf++;
		}

//		std::cout << g.error <<std::endl;

//		std::cout << g.error.colwise().sum() << std::endl;
//		getDefVec(defVecS, defVec_b, n, finalDef,
		/*
		std::cout << "doing projection of surface nodes\n";

		std::cout << "done\n";



		Eigen::MatrixXd p_ic = PhiPtr->Phi_ic.block(m.N_m-n.N_m + m.N_se- n.N_se,0,m.N_ss - n.N_ss,n.N_m+n.N_se);

		Eigen::MatrixXd p_cc = PhiPtr->Phi_cc.block(0,0,n.N_m+n.N_se,n.N_m+n.N_se);
		Eigen::ArrayXi iNodesReduced = (*n.iPtr)(Eigen::seqN(m.N_m-n.N_m + m.N_se- n.N_se,m.N_ss-n.N_ss));

		Eigen::ArrayXi* iPtrReduced = &iNodesReduced;
//		performRBF(p_cc,p_ic,defVec_b,n.cPtr,iPtrReduced,n.N_c);

//		std::cout << d.rows() << '\t' << n.N_i << '\t' <<






		m.coords(*n.iPtr, Eigen::all) += d;
		for(int dim = 0 ; dim < m.nDims; dim++){
			m.coords(*n.cPtr, dim) += (defVec_b(Eigen::seqN(dim*n.N_c,n.N_c))).array();
		}
		m.writeMeshFile(params.mesh_ifName, params.mesh_ofName);
		std::exit(0);

		m.coords(*n.iPtr,Eigen::all) += (d-g.error);
		for(int dim = 0 ; dim < m.nDims; dim++){
			m.coords(*n.cPtr, dim) += (defVec_b(Eigen::seqN(dim*n.N_c,n.N_c))).array();
		}
		m.writeMeshFile(params.mesh_ifName, params.mesh_ofName);

		std::cout << "DONE\n";
		std::exit(0);
		*/
		updateNodes(n,defVec_b, g.d_step, g.alpha_step, g.ctrlPtr);
		g.correction(m,n, params.gamma, params.multiLvl);
		iterating = true;
//		m.writeMeshFile(params.mesh_ifName, params.mesh_ofName);
//		std::exit(0);
		std::cout << "number of control nodes: " << n.N_c << std::endl;

	}
}

void rbf_ps::performRBF_PS(PhiStruct* PhiPtr, Eigen::VectorXd& defVec,Eigen::ArrayXXd& delta, Eigen::ArrayXXd& finalDef, getNodeType& n, int N_i, int N_c, int shape_type, Eigen::MatrixXd& Phi_cc, Eigen::MatrixXd& Phi_ic){

	std::cout << "Performing rbf ps\n";

	delta.resize(N_i, m.nDims);
	finalDef.resize(N_i, m.nDims);

//	std::cout << N_c << '\t' << N_i << std::endl;
//	std::cout << delta.rows() << '\t' << Phi_ic.rows() << '\t' << Phi_ic.cols() << '\t' << Phi_cc.rows() << '\t' << Phi_cc.cols() << std::endl;
	for(int dim = 0; dim < m.nDims; dim++){
		delta(Eigen::seqN(0,N_i),dim) = (Phi_ic*(Phi_cc.fullPivLu().solve(defVec(Eigen::seqN(dim*N_c,N_c))))).array();
	}

//	p.project(m, n, delta, finalDef, m.periodicVec);
	if(shape_type == 0){
		p.projectEdge(m, n.sePtr, delta, finalDef, m.periodicVec, 0, N_i, 1);
	}else if(shape_type == 1){
		p.projectSurf(m, n.ssPtr, delta, finalDef, m.periodicVec, 0, N_i, 1);
	}

//	std::cout << "\nfinalDef:\n" << finalDef << "\n\n";

	//for 2D

//	delta.resize(n.N_se, m.nDims);
//	finalDef.resize(n.N_se,m.nDims);
//
//	for(int dim = 0; dim < m.nDims; dim++){
//		delta.col(dim) = (PhiPtr->Phi_em*(PhiPtr->Phi_mm.fullPivLu().solve(defVec(Eigen::seqN(dim*n.N_m,n.N_m))))).array();
////		delta.col(dim) = (PhiPtr->Phi_cc.block(n.N_m, 0, n.N_se, n.N_m)*(PhiPtr->Phi_cc.block(0,0,n.N_m,n.N_m).fullPivLu().solve(defVec(Eigen::seqN(dim*n.N_m,n.N_m))))).array();
//	}
////	p.project(m, n, delta, finalDef, m.periodicVec);
//	p.projectEdge(m, n.sePtr, delta, finalDef, m.periodicVec, 0, N_i, 1);

	// same for 3D fails as this projects on the other surfaces
	/*
	delta.resize(n.N_se+n.N_ss, m.nDims);
	finalDef.resize(n.N_se+n.N_ss,m.nDims);

	for(int dim = 0; dim < m.nDims; dim++){
		delta(Eigen::seqN(0,n.N_se),dim) = (PhiPtr->Phi_em*(PhiPtr->Phi_mm.fullPivLu().solve(defVec(Eigen::seqN(dim*n.N_m,n.N_m))))).array();
		delta(Eigen::seqN(n.N_se,n.N_ss),dim) = (PhiPtr->Phi_cc.block(n.N_m+n.N_se, 0, n.N_ss, n.N_m)*(PhiPtr->Phi_mm.fullPivLu().solve(defVec(Eigen::seqN(dim*n.N_m,n.N_m))))).array();
//		delta.col(dim) = (PhiPtr->Phi_cc.block(n.N_m, 0, n.N_se, n.N_m)*(PhiPtr->Phi_cc.block(0,0,n.N_m,n.N_m).fullPivLu().solve(defVec(Eigen::seqN(dim*n.N_m,n.N_m))))).array();
	}
*/

	/*
	delta.resize(n.N_se, m.nDims);
	finalDef.resize(n.N_se,m.nDims);

	for(int dim = 0; dim < m.nDims; dim++){
		delta(Eigen::seqN(0,n.N_se),dim) = (PhiPtr->Phi_em*(PhiPtr->Phi_mm.fullPivLu().solve(defVec(Eigen::seqN(dim*n.N_m,n.N_m))))).array();
//		delta(Eigen::seqN(n.N_se,n.N_ss),dim) = (PhiPtr->Phi_cc.block(n.N_m+n.N_se, 0, n.N_ss, n.N_m)*(PhiPtr->Phi_mm.fullPivLu().solve(defVec(Eigen::seqN(dim*n.N_m,n.N_m))))).array();
//		delta.col(dim) = (PhiPtr->Phi_cc.block(n.N_m, 0, n.N_se, n.N_m)*(PhiPtr->Phi_cc.block(0,0,n.N_m,n.N_m).fullPivLu().solve(defVec(Eigen::seqN(dim*n.N_m,n.N_m))))).array();
	}

	p.projectEdge(m, n.sePtr, delta, finalDef, m.periodicVec, 0, n.N_se, 1);
*/

	/*
	Eigen::VectorXd defVec2(m.nDims*(n.N_m+n.N_se));

	for(int dim = 0; dim < m.nDims; dim++){
		defVec2(Eigen::seqN(dim*(n.N_m+n.N_se), n.N_m)) = defVec(Eigen::seqN(dim*n.N_m,n.N_m));
		defVec2(Eigen::seqN(dim*(n.N_m+n.N_se)+n.N_m, n.N_se)) = finalDef.col(dim);
	}

	Eigen::ArrayXXd delta2;
	Eigen::ArrayXXd finalDef2;

	delta2.resize(n.N_ss, m.nDims);
	finalDef2.resize(n.N_ss,m.nDims);


	for(int dim = 0; dim < m.nDims; dim++){
			delta2(Eigen::seqN(0,n.N_ss),dim) = (PhiPtr->Phi_cc.block(n.N_m+n.N_se, 0, n.N_ss, n.N_m+n.N_se)*(PhiPtr->Phi_cc.block(0,0,n.N_m+n.N_se,n.N_m+n.N_se).fullPivLu().solve(defVec2(Eigen::seqN(dim*(n.N_m+n.N_se),n.N_m+n.N_se))))).array();
	//		delta(Eigen::seqN(n.N_se,n.N_ss),dim) = (PhiPtr->Phi_cc.block(n.N_m+n.N_se, 0, n.N_ss, n.N_m)*(PhiPtr->Phi_mm.fullPivLu().solve(defVec(Eigen::seqN(dim*n.N_m,n.N_m))))).array();
	//		delta.col(dim) = (PhiPtr->Phi_cc.block(n.N_m, 0, n.N_se, n.N_m)*(PhiPtr->Phi_cc.block(0,0,n.N_m,n.N_m).fullPivLu().solve(defVec(Eigen::seqN(dim*n.N_m,n.N_m))))).array();
	}


	p.projectSurf(m, n.ssPtr, delta2, finalDef2, m.periodicVec, 0, n.N_ss, 1);


	m.coords(*n.sePtr, Eigen::all) += finalDef;
	m.coords(*n.ssPtr, Eigen::all) += finalDef2;
	for(int dim =0;dim<m.nDims;dim++){
		m.coords(*n.mPtr, dim) += (defVec(Eigen::seqN(dim*n.N_m, n.N_m))).array();
	}
	m.writeMeshFile(params.mesh_ifName, params.mesh_ofName);
	std::cout << "done\n";
	std::exit(0);
*/
//
//	if(m.nDims == 3){
//		m.coords(*n.sePtr, Eigen::all) += finalDef;
//		for(int dim =0;dim<m.nDims;dim++){
//			m.coords(*n.mPtr, dim) += (defVec(Eigen::seqN(dim*n.N_m, n.N_m))).array();
//		}
//		m.writeMeshFile(params.mesh_ifName, params.mesh_ofName);
//		std::cout << "exit status reached\n";
//		std::exit(0);
//	}
//

//	for(int dim = 0; dim < m.nDims; dim++){
//		m.coords(*n.cPtr, dim) += (defVec(Eigen::seqN(dim*n.N_c, n.N_c))).array();
//	}






//	m.coords(*n.sPtr,Eigen::all) += finalDef;
//	m.writeMeshFile(params.mesh_ifName, params.mesh_ofName);
//	std::exit(0);



//	performRBF(PhiPtr->Phi_cc,PhiPtr->Phi_ic,defVec_b,n.cPtr,n.iPtr,n.N_c);
//	std::cout << "performed rbf" << std::endl;
//	std::cout << defVec_b << "\n" << std::endl;
//	std::cout << PhiPtr->Phi_cc << "\n\n" << PhiPtr->Phi_ic << std::endl;

}



