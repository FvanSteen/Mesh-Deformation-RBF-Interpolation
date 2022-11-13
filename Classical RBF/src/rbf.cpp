#include "Mesh.h"
#include "rbf.h"
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <string>
#include <chrono>
#include "TestingGround.h"

rbf::rbf(Mesh& meshOb, Eigen::VectorXd& displacementVector, Eigen::VectorXd& rVec, const int steps, Eigen::RowVectorXd& rotationPnt, const std::string& slidingMode, const std::string& periodicDirection, const bool& curved)
:m(meshOb), dVec(displacementVector), rotPnt(rotationPnt), steps(steps), smode(slidingMode), rotVec(rVec), curved(curved)
{
	// todo specifying the periodic direction(s)
	// todo smode is already in the mesh class -> replace all instances with m.smode


	getNodeTypes();

	getPeriodicParams(periodicDirection);

	getRotationalMat();



}

void rbf::RBFMain(){

	if(smode=="none"){
		RBF_standard();
	}
	else if(smode=="ds"){
		if(m.nDims ==3){
			RBF_DS_3D();
		}else if(m.nDims ==2){
			RBF_DS();
		}
	}
	else if(smode=="ps"){
		RBF_PS();
	}
}

void rbf::RBF_PS(){
	std::cout<< "Performing RBF PS" << std::endl;


	Eigen::MatrixXd Phi_mmPro, Phi_sm, Phi_mm, Phi_im; 	//In this case only the internal boundary and static ext bdry
	Eigen::VectorXd defVecPro, defVec;
	Eigen::ArrayXXd delta(N_s, m.nDims), finalDef(N_s,m.nDims);

	for(int i = 0; i < steps; i++){
		std::cout << "Deformation step: " << i+1 << std::endl;

		getPhi(Phi_mmPro, mNodesPro, mNodesPro);
		getPhi(Phi_sm, sNodes, mNodesPro);
		getPhi(Phi_mm, mNodes, mNodes);
		getPhi(Phi_im, iNodes, mNodes);

//		defVecPro = Eigen::VectorXd::Zero(N_mPro*m.nDims);
		getDefVec(defVecPro, N_mPro, N_mPro);

		if(curved || i==0){
			m.getMidPnts();
		}
		performRBF_PS(Phi_mmPro, Phi_sm, Phi_mm, Phi_im, defVecPro, delta, finalDef, defVec);
//		if(i==0){
//			m.writeMeshFile();
//			std::exit(0);
//		}
	}

	m.writeMeshFile();
}

void rbf::performRBF_PS(Eigen::MatrixXd& Phi_mmPro, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_mm, Eigen::MatrixXd& Phi_im, Eigen::VectorXd& defVecPro,Eigen::ArrayXXd& delta, Eigen::ArrayXXd& finalDef, Eigen::VectorXd& defVec){

	for(int dim = 0; dim < m.nDims; dim++){
		delta.col(dim) = (Phi_sm*(Phi_mmPro.fullPivLu().solve(defVecPro(Eigen::seqN(dim*N_mPro,N_mPro))))).array();
	}



	// todo make pVec an array instead of vector
//	m.getVecs();
//	for(int i=0; i<m.N_se; i++){
//		finalDef.row(i) = (delta.row(i)).matrix().transpose() - (delta.row(i).matrix()).dot((m.n.row(i)).matrix().transpose())*(m.n.row(i)).matrix().transpose();
////				std::cout << "finalDef \t" << finalDef.row(i) << std::endl;
////				std::cout << "delta \t" << delta.row(i) << std::endl;
//	}
	// updating the midpoints on the external boundary of the mesh
//	std::cout << delta << std::endl;
	project(delta,finalDef);

//	std::cout << finalDef << std::endl;

	//todo make an if statement for fixed vertices
	if(m.pmode == "moving"){
		for(int i=0; i < m.N_es ; i++){
			finalDef.row(m.N_se+i) = pVec.transpose().array()*delta.row(m.N_se+i);
//		std::cout << "finalDef \t" << finalDef.row(m.N_se+i) << std::endl;
//		std::cout << "delta \t" << delta.row(m.N_se+i) << std::endl;
//			m.coords.row(m.extStaticNodes(i)) +=  pVec.transpose().array()*delta.row(m.N_se+i);
		}
	}





//	m.coords(m.extStaticNodes,1) += delta(Eigen::seq(Eigen::last+1-4,Eigen::last), 1);
//	m.coords(sNodes,Eigen::all) = m.coords(sNodes,Eigen::all) + finalDef;
//	m.writeMeshFile();


	defVec = Eigen::VectorXd::Zero(N_m*m.nDims);
//	std::cout << N_m << '\t' << N_mPro <<  '\t' << m.N_se << std::endl;
	for(int dim = 0; dim< m.nDims; dim++){
		defVec(Eigen::seqN(dim*N_m,N_mPro)) = defVecPro(Eigen::seqN(dim*N_mPro,N_mPro));
//		defVec(Eigen::seqN(dim*N_m+N_mPro,m.N_se)) = finalDef.col(dim);
		defVec(Eigen::seqN(dim*N_m+N_mPro,N_s)) = finalDef.col(dim);
	}

//	std::cout << "finalDef: \n" << finalDef << "\n" << std::endl;
//	std::cout << "defVec: \n" << defVec << std::endl;

//	std::exit(0);
	performRBF(Phi_mm,Phi_im,defVec,mNodes,N_m);
}

void rbf::RBF_DS_3D(){
	std::cout << "Performing RBF DS in three dimensions" << std::endl;

	Eigen::MatrixXd Phi_mm, Phi_ms, Phi_sm, Phi_ss, Phi_im, Phi_is, Phi;
//	Eigen::ArrayXXd n(m.N_se, m.nDims), t(m.N_se, m.nDims);		// two column array containing normal vector components
	Eigen::VectorXd defVec, alpha(m.nDims*(N_m+m.N_se+m.N_ss));
//	Eigen::ArrayXXd t_se(m.N_se,m.nDims),n1_se(m.N_se,m.nDims), n2_se(m.N_se,m.nDims), n_ss(m.N_ss, m.nDims),t1_ss(m.N_ss, m.nDims),t2_ss(m.N_ss, m.nDims);
	Eigen::MatrixXd Phi_me, Phi_ee, Phi_es, Phi_em, Phi_se, Phi_ie;
	for (int i = 0; i < steps; i++){
		std::cout << "Deformation step: " << i+1 << std::endl;
		getPhi(Phi_mm, mNodes,mNodes);
		getPhi(Phi_me, mNodes, m.slidingEdgeNodes);
		getPhi(Phi_ms, mNodes, m.slidingSurfNodes);

		getPhi(Phi_em, m.slidingEdgeNodes, mNodes);
		getPhi(Phi_ee, m.slidingEdgeNodes, m.slidingEdgeNodes);
		getPhi(Phi_es, m.slidingEdgeNodes, m.slidingSurfNodes);

		getPhi(Phi_sm, m.slidingSurfNodes, mNodes);
		getPhi(Phi_se, m.slidingSurfNodes, m.slidingEdgeNodes);
		getPhi(Phi_ss, m.slidingSurfNodes, m.slidingSurfNodes);
		getPhi(Phi_im, m.intNodes, mNodes);
		getPhi(Phi_ie, m.intNodes, m.slidingEdgeNodes);
		getPhi(Phi_is, m.intNodes, m.slidingSurfNodes);




		std::cout << "Reached the point at which normal and tangential vectors should be found" << std::endl;
		m.getVecs();
//		m.getVecs3D(t_se, n1_se, n2_se, n_ss,t1_ss,t2_ss); // function that obtains all the normal and tan vectors


		defVec = Eigen::VectorXd::Zero((N_m+m.N_se+m.N_ss)*m.nDims);
		getDefVec(defVec,N_m,(N_m+m.N_se+m.N_ss));
//		std::cout << defVec << std::endl;

		getPhiDS_3D(Phi,Phi_mm, Phi_me, Phi_ms, Phi_em, Phi_ee, Phi_es, Phi_sm, Phi_se, Phi_ss);
//		std::cout << Phi << std::endl;
//		std::cout << "DONE" << std::endl;
//		std::exit(0);
		performRBF_DS_3D(Phi, Phi_im, Phi_ie, Phi_is, Phi_em, Phi_ee, Phi_es, Phi_sm, Phi_se, Phi_ss,defVec, alpha);
	}

	m.writeMeshFile();
}

void rbf::performRBF_DS_3D(Eigen::MatrixXd& Phi, Eigen::MatrixXd& Phi_im, Eigen::MatrixXd& Phi_ie, Eigen::MatrixXd& Phi_is, Eigen::MatrixXd& Phi_em, Eigen::MatrixXd& Phi_ee, Eigen::MatrixXd& Phi_es, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_se, Eigen::MatrixXd& Phi_ss, Eigen::VectorXd& defVec, Eigen::VectorXd& alpha){

	alpha = Phi.partialPivLu().solve(defVec);
	std::cout << "SOLUTION WAS FOUND" << std::endl;

	for (int dim = 0; dim < m.nDims; dim++){

		m.coords(m.intNodes, dim) += (Phi_im*alpha(Eigen::seqN(dim*(N_m+m.N_se+m.N_ss),N_m)) + Phi_ie*alpha(Eigen::seqN(dim*(N_m+m.N_se+m.N_ss)+N_m, m.N_se)) + Phi_is*alpha(Eigen::seqN(dim*(N_m+m.N_se+m.N_ss)+N_m+m.N_se, m.N_ss)) ).array();
		m.coords(m.slidingEdgeNodes, dim) += (Phi_em*alpha(Eigen::seqN(dim*(N_m+m.N_se+m.N_ss),N_m)) + Phi_ee*alpha(Eigen::seqN(dim*(N_m+m.N_se+m.N_ss)+N_m, m.N_se)) + Phi_es*alpha(Eigen::seqN(dim*(N_m+m.N_se+m.N_ss)+N_m+m.N_se, m.N_ss)) ).array();
		m.coords(m.slidingSurfNodes, dim) += (Phi_sm*alpha(Eigen::seqN(dim*(N_m+m.N_se+m.N_ss),N_m)) + Phi_se*alpha(Eigen::seqN(dim*(N_m+m.N_se+m.N_ss)+N_m, m.N_se)) + Phi_ss*alpha(Eigen::seqN(dim*(N_m+m.N_se+m.N_ss)+N_m+m.N_se, m.N_ss)) ).array();
		m.coords(m.intBdryNodes, dim) += (defVec(Eigen::seqN(dim*N_m,m.N_ib))).array();
		rotPnt(dim) += dVec(dim);
	}
}

void rbf::getPhiDS_3D(Eigen::MatrixXd& Phi, Eigen::MatrixXd& Phi_mm, Eigen::MatrixXd&  Phi_me, Eigen::MatrixXd&  Phi_ms, Eigen::MatrixXd& Phi_em, Eigen::MatrixXd& Phi_ee, Eigen::MatrixXd& Phi_es, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_se, Eigen::MatrixXd& Phi_ss){

	Phi = Eigen::MatrixXd::Zero(m.nDims*(N_m+m.N_se+m.N_ss),m.nDims*(N_m+m.N_se+m.N_ss));

	for(int dim = 0; dim< m.nDims; dim++){
		// blocks related to the known displacements
		Phi.block(dim*N_m, dim*(N_m+m.N_se+m.N_ss), N_m, N_m) = Phi_mm;
		Phi.block(dim*N_m, dim*(N_m+m.N_se+m.N_ss)+N_m, N_m, m.N_se) = Phi_me;
		Phi.block(dim*N_m, dim*(N_m+m.N_se+m.N_ss)+N_m+m.N_se, N_m, m.N_ss) = Phi_ms;

		//blocks realteð to the first zero normal displacement condition of the sliding edge nodes
		Phi.block(3*N_m, dim*(N_m+m.N_se+m.N_ss),					m.N_se, N_m ) = Phi_em.array().colwise() * m.n1_se.col(dim);
		Phi.block(3*N_m, dim*(N_m+m.N_se+m.N_ss) + N_m,				m.N_se, m.N_se ) = Phi_ee.array().colwise() * m.n1_se.col(dim);
		Phi.block(3*N_m, dim*(N_m+m.N_se+m.N_ss) + N_m + m.N_se,	m.N_se,	m.N_ss ) = Phi_es.array().colwise() * m.n1_se.col(dim);
		//blocks realteð to the second zero normal displacement condition of the sliding edge nodes
		Phi.block(3*N_m+m.N_se, dim*(N_m+m.N_se+m.N_ss),					m.N_se, N_m ) = Phi_em.array().colwise() * m.n2_se.col(dim);
		Phi.block(3*N_m+m.N_se, dim*(N_m+m.N_se+m.N_ss) + N_m,				m.N_se, m.N_se ) = Phi_ee.array().colwise() * m.n2_se.col(dim);
		Phi.block(3*N_m+m.N_se, dim*(N_m+m.N_se+m.N_ss) + N_m + m.N_se,		m.N_se,	m.N_ss ) = Phi_es.array().colwise() * m.n2_se.col(dim);
		// blocks related to the zero tangential contribution of the sliding edge nodes
		Phi.block(3*N_m+2*m.N_se, dim*(N_m+m.N_se+m.N_ss) + N_m,	m.N_se, m.N_se) = Eigen::MatrixXd(m.t_se.col(dim).matrix().asDiagonal());
		// blocks related to the zero normal displacement condition of the sliding surface nodes
		Phi.block(3*N_m+3*m.N_se, dim*(N_m+m.N_se+m.N_ss),		m.N_ss,N_m) = Phi_sm.array().colwise() * m.n_ss.col(dim);
		Phi.block(3*N_m+3*m.N_se, dim*(N_m+m.N_se+m.N_ss)+N_m,	m.N_ss,m.N_se) = Phi_se.array().colwise() * m.n_ss.col(dim);
		Phi.block(3*N_m+3*m.N_se, dim*(N_m+m.N_se+m.N_ss)+N_m+m.N_se,	m.N_ss,m.N_ss) = Phi_ss.array().colwise() * m.n_ss.col(dim);

		// blocks related to the zero tangential contribution of the sliding surface nodes.
		Phi.block(3*N_m+3*m.N_se+m.N_ss, dim*(N_m+m.N_se+m.N_ss)+N_m+m.N_se, m.N_ss,m.N_ss) = Eigen::MatrixXd(m.t1_ss.col(dim).matrix().asDiagonal());
		Phi.block(3*N_m+3*m.N_se+2*m.N_ss, dim*(N_m+m.N_se+m.N_ss)+N_m+m.N_se, m.N_ss,m.N_ss) = Eigen::MatrixXd(m.t2_ss.col(dim).matrix().asDiagonal());
	}


}


void rbf::RBF_DS(){
	std::cout << "Performing RBF DS " << std::endl;

	Eigen::MatrixXd Phi_mm, Phi_ms, Phi_sm, Phi_ss, Phi_im, Phi_is, Phi;
	Eigen::VectorXd defVec, alpha(m.nDims*(N_m+N_s));

	for (int i = 0; i < steps; i++){
		std::cout << "Deformation step: " << i+1 << std::endl;

		getPhi(Phi_mm, mNodes,mNodes);
		getPhi(Phi_ms, mNodes, sNodes);
		getPhi(Phi_sm, sNodes, mNodes);
		getPhi(Phi_ss, sNodes, sNodes);
		getPhi(Phi_im, iNodes, mNodes);
		getPhi(Phi_is, iNodes, sNodes);


		if(curved || i==0){
			// getVecs obtains average vector at the nodes
			m.getVecs();
			// getMidPnts obtains vectors at midpoint of boundary segments
			m.getMidPnts();


		}

//		defVec = Eigen::VectorXd::Zero((N_m+N_s)*m.nDims);
		getDefVec(defVec,N_m,(N_m+N_s));

		getPhiDS(Phi,Phi_mm,Phi_ms, Phi_sm, Phi_ss, m.n, m.t);

		performRBF_DS(Phi, Phi_im, Phi_is, Phi_sm, Phi_ss,defVec, alpha);
	}

	m.writeMeshFile();
}

void rbf::performRBF_DS(Eigen::MatrixXd& Phi, Eigen::MatrixXd& Phi_im, Eigen::MatrixXd& Phi_is, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_ss, Eigen::VectorXd& defVec, Eigen::VectorXd& alpha){

	alpha = Phi.fullPivLu().solve(defVec);

	if(curved){
		Eigen::ArrayXXd delta(N_s, m.nDims), finalDef(N_s,m.nDims);

		// find displacement
		for (int dim = 0; dim < m.nDims; dim++){
			delta.col(dim) = (Phi_sm*alpha(Eigen::seqN(dim*(N_m+N_s),N_m)) + Phi_ss*alpha(Eigen::seqN(dim*(N_m+N_s)+N_m, N_s))).array();

		}
//		m.coords(sNodes,Eigen::seq(0,1)) += delta;



		// calling project function to find the final deformation after the projection
		project(delta,finalDef);
//		m.coords(sNodes,Eigen::seq(0,1)) += finalDef;
//		m.writeMeshFile();
//

//		std::exit(0);
//		defVec = Eigen::VectorXd::Zero(N_m2*m.nDims);

		getDefVec(defVec,N_m2,N_m2);
//		std::cout << defVec << std::endl;

//		std::cout << defVec << std::endl;
		for(int dim = 0; dim< m.nDims; dim++){
//			defVec(Eigen::seqN(dim*(N_m2)+m.N_ib+m.N_es,N_s)) = finalDef.col(dim);
			defVec(Eigen::seqN(dim*(N_m2)+N_m,N_s)) = finalDef.col(dim);
			//todo difference in last two lines, for fixed/moving
		}




//		std::cout << finalDef << std::endl;
//		std::cout << defVec << std::endl;
//		std::exit(0);
		Eigen::MatrixXd Phi_mm2, Phi_im2;

		getPhi(Phi_mm2, mNodes2,mNodes2);
		getPhi(Phi_im2,iNodes,mNodes2);

		performRBF(Phi_mm2,Phi_im2,defVec,mNodes2,N_m2);

	}
	else{
		for (int dim = 0; dim < m.nDims; dim++){
			m.coords(iNodes, dim) += (Phi_im*alpha(Eigen::seqN(dim*(N_m+N_s),N_m)) + Phi_is*alpha(Eigen::seqN(dim*(N_m+N_s)+N_m, N_s))).array();
			m.coords(sNodes, dim) += (Phi_sm*alpha(Eigen::seqN(dim*(N_m+N_s),N_m)) + Phi_ss*alpha(Eigen::seqN(dim*(N_m+N_s)+N_m, N_s))).array();
			m.coords(m.intBdryNodes, dim) += (defVec(Eigen::seqN(dim*N_m,m.N_ib))).array();
			rotPnt(dim) += dVec(dim);
		}
	}
}



void rbf::getPhiDS(Eigen::MatrixXd& Phi,Eigen::MatrixXd& Phi_mm,Eigen::MatrixXd& Phi_ms, Eigen::MatrixXd& Phi_sm, Eigen::MatrixXd& Phi_ss, Eigen::ArrayXXd& n, Eigen::ArrayXXd& t){

	Phi = Eigen::MatrixXd::Zero(m.nDims*(N_m+N_s),m.nDims*(N_m+N_s));


	if(m.pmode == "moving"){
		n.conservativeResize(N_s,m.nDims);
		t.conservativeResize(N_s,m.nDims);
		for(int i = 0; i<m.N_es; i++){
			n.row(m.N_se+i) = pnVec;
			t.row(m.N_se+i) = pVec;
		}
	}


	for(int dim = 0; dim< m.nDims; dim++){
		// blocks related to the known displacements
		Phi.block(dim*N_m, dim*(N_m+N_s), N_m, N_m) = Phi_mm;
		Phi.block(dim*N_m, dim*(N_m+N_s)+N_m, N_m, N_s) = Phi_ms;

		// blocks related to the zero normal displacement condition
		Phi.block(2*N_m, dim*(N_m+N_s), N_s, N_m) = Phi_sm.array().colwise() * n.col(dim);
		Phi.block(2*N_m, dim*(N_m+N_s)+N_m, N_s, N_s) = Phi_ss.array().colwise() * n.col(dim);

		//blocks related to the zero tangential contribution condition
		Phi.block(2*N_m + N_s, dim*(N_m+N_s)+N_m, N_s, N_s) = Eigen::MatrixXd(t.col(dim).matrix().asDiagonal());
	}


}



void rbf::RBF_standard(){

	auto start = std::chrono::high_resolution_clock::now();
	Eigen::MatrixXd Phi_mm, Phi_im;
	Eigen::VectorXd defVec;

	for(int i=0; i<steps; i++){
		auto start2 = std::chrono::high_resolution_clock::now();
		std::cout << "Deformation step: " << i+1 << std::endl;

		std::cout << "Obtaining Phi_mm" << std::endl;
		auto starti = std::chrono::high_resolution_clock::now();
		getPhi(Phi_mm, mNodes,mNodes);
		auto stopi = std::chrono::high_resolution_clock::now();
		auto durationi = std::chrono::duration_cast<std::chrono::microseconds>(stopi-starti);
		std::cout << "Runtime duration obtaining Phi_mm: \t"<<  durationi.count()/1e6 << " seconds"<< std::endl;


		std::cout << "Obtaining Phi_im" << std::endl;
		starti = std::chrono::high_resolution_clock::now();
		getPhi(Phi_im, iNodes,mNodes);
		stopi = std::chrono::high_resolution_clock::now();
		durationi = std::chrono::duration_cast<std::chrono::microseconds>(stopi-starti);
		std::cout << "Runtime duration obtaining Phi_im: \t"<<  durationi.count()/1e6 << " seconds"<< std::endl;


		std::cout << "Obtaining deformation vector" << std::endl;
		starti = std::chrono::high_resolution_clock::now();

		getDefVec(defVec,N_m, N_m);
		stopi = std::chrono::high_resolution_clock::now();
		durationi = std::chrono::duration_cast<std::chrono::microseconds>(stopi-starti);
		std::cout << "Runtime duration obtaining defVec: \t"<<  durationi.count()/1e6 << " seconds"<< std::endl;

		std::cout << "Performing RBF" << std::endl;
		starti = std::chrono::high_resolution_clock::now();
		performRBF(Phi_mm, Phi_im, defVec,mNodes,N_m);
		stopi = std::chrono::high_resolution_clock::now();
		durationi = std::chrono::duration_cast<std::chrono::microseconds>(stopi-starti);
		std::cout << "Runtime duration obtaining solution and updating nodes: \t"<<  durationi.count()/1e6 << " seconds"<< std::endl;




		auto stop2 = std::chrono::high_resolution_clock::now();
		auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(stop2-start2);
		std::cout << "Runtime duration whole step: \t"<<  duration2.count()/1e6 << " seconds"<< std::endl;
//		TestingGround tg;
//		tg.testFunction();
//		std::exit(0);
	}
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);

	std::cout << "Runtime duration: \t"<<  duration.count()/1e6 << " seconds"<< std::endl;
	m.writeMeshFile();
}

void rbf::performRBF(Eigen::MatrixXd& Phi_mm, Eigen::MatrixXd& Phi_im, Eigen::VectorXd& defVec, Eigen::ArrayXi& movingNodes, int& N){
	for(int dim = 0; dim < m.nDims; dim++){
		std::cout << "Solving for dimension: " << dim << std::endl;
		m.coords(iNodes,dim) += (Phi_im*(Phi_mm.fullPivLu().solve(defVec(Eigen::seqN(dim*N,N))))).array();
		m.coords(movingNodes,dim) += defVec(Eigen::seqN(dim*N,N)).array();
		rotPnt(dim) += dVec(dim);
	}
}


void rbf::getPhi(Eigen::MatrixXd& Phi, Eigen::ArrayXi& idxSet1, Eigen::ArrayXi& idxSet2){
	Phi.resize(idxSet1.size(), idxSet2.size());
	double dist;
	for(int i=0; i<idxSet1.size();i++){
		for(int j=0; j<idxSet2.size();j++){
			if(m.nDims == 2){
				// Euclidian distance

//				dist = sqrt(pow(m.coords(idxSet1(i),0)-m.coords(idxSet2(j),0),2) + pow(m.coords(idxSet1(i),1)-m.coords(idxSet2(j),1),2));


				//todo following if statement is introduced to improve efficiency, check if anything else can be done
				if(m.pmode != "none"){
					dist = 0;
					for(int dim = 0; dim < m.nDims; dim++){
						if(pVec(dim)){
							//todo change function to allow for periodic length in stead of unit length
							dist += pow(1/M_PI*sin( (m.coords(idxSet1(i),dim)-m.coords(idxSet2(j),dim))*M_PI/1),2);
						}
						else{
							dist += pow(m.coords(idxSet1(i),dim)-m.coords(idxSet2(j),dim),2);
						}

					}
					dist = sqrt(dist);
				}
				else{
					//todo can the difference in coordinates by calculated for the whole row. then take the power of 2, sum and take square root??
					dist = sqrt(pow(m.coords(idxSet1(i),0)-m.coords(idxSet2(j),0),2) + pow(m.coords(idxSet1(i),1)-m.coords(idxSet2(j),1),2) );
				}
//				dist = sqrt(pow(m.coords(idxSet1(i),0)-m.coords(idxSet2(j),0),2) + pow(1/M_PI*sin( (m.coords(idxSet1(i),1)-m.coords(idxSet2(j),1))*M_PI/1),2));
//				std::cout << dist << std::endl;
//				std::cout << 'x' << '\t' << m.coords(idxSet1(i),0)-m.coords(idxSet2(j),0) << std::endl;
//				std::cout << 'y' << '\t' << m.coords(idxSet1(i),1)-m.coords(idxSet2(j),1) << std::endl;

			}
			//todo if statements can probably by removed if the calc is done with the previous for loop.
			else if(m.nDims == 3){
				dist = sqrt(pow(m.coords(idxSet1(i),0)-m.coords(idxSet2(j),0),2) + pow(m.coords(idxSet1(i),1)-m.coords(idxSet2(j),1),2) + pow(m.coords(idxSet1(i),2)-m.coords(idxSet2(j),2),2));
			}
			Phi(i,j) = pow((1-(dist/m.r)),4)*(4*(dist/m.r)+1);
		}
	}
//	std::exit(0);
}

void rbf::getDefVec(Eigen::VectorXd& defVec, int& N, int defVecLength){
	Eigen::MatrixXd intPnts(m.N_ib,m.nDims);
	Eigen::MatrixXd rotDef;
	defVec = Eigen::VectorXd::Zero(defVecLength*m.nDims);
	intPnts = m.coords(m.intBdryNodes,Eigen::all);

	if(m.nDims == 2){
		rotDef = (rotMat*(intPnts.rowwise() - rotPnt).transpose()).transpose().rowwise() +rotPnt - intPnts;
	}
	else if(m.nDims == 3){
		rotDef = Eigen::MatrixXd::Zero(m.N_ib,m.nDims);

		if(rotVec[0] != 0){
			rotDef+= (rotMatX*(intPnts.rowwise() - rotPnt).transpose()).transpose().rowwise() +rotPnt - intPnts;
		}
		if(rotVec[1] != 0){
			rotDef+= (rotMatY*(intPnts.rowwise() - rotPnt).transpose()).transpose().rowwise() +rotPnt - intPnts;
		}
		if(rotVec[2] != 0){
			rotDef+= (rotMatZ*(intPnts.rowwise() - rotPnt).transpose()).transpose().rowwise() +rotPnt - intPnts;
		}
	}

	for(int dim = 0; dim < m.nDims; dim++){
		defVec(Eigen::seqN(dim*N, m.N_ib)).array() += dVec(dim);
		defVec(Eigen::seqN(dim*N, m.N_ib)) += rotDef.col(dim);
	}

 }

void rbf::getRotationalMat(){
	if(m.nDims == 2){
		const double theta = rotVec[0]/steps*M_PI/180;
		rotMat << 	cos(theta), -sin(theta),
					sin(theta),	cos(theta);
	}
	if(m.nDims == 3){
		const double x_theta = rotVec[0]/steps*M_PI/180;
		const double y_theta = rotVec[1]/steps*M_PI/180;
		const double z_theta = rotVec[2]/steps*M_PI/180;
		rotMatX << 	1,	0,	0,
					0,	cos(x_theta),	-sin(x_theta),
					0,	sin(x_theta),	cos(x_theta);
		rotMatY << 	cos(y_theta),	0,	sin(y_theta),
					0,	1,	0,
					-sin(y_theta),	0,	cos(y_theta);
		rotMatZ <<	cos(z_theta),	-sin(z_theta),	 0,
					sin(z_theta),	cos(z_theta),	 0,
					0,	0,	1;

	}
}

void rbf::getNodeTypes(){
	iNodes.resize(m.N_i+m.N_p);
	iNodes << m.intNodes, m.periodicNodes;

	if(m.smode== "none" || smode == "ps"){
		mNodes.resize(m.N_ib + m.N_es + m.N_se);
		mNodes << m.intBdryNodes, m.extStaticNodes, m.slidingEdgeNodes;

		if(m.smode == "ps" && m.pmode == "moving"){
			mNodesPro.resize(m.N_ib);
			mNodesPro << m.intBdryNodes;

			sNodes.resize(m.N_se + m.N_es);
			sNodes << m.slidingEdgeNodes, m.extStaticNodes;

		}else{
			mNodesPro.resize(m.N_ib+m.N_es);
			mNodesPro << m.intBdryNodes, m.extStaticNodes;

			sNodes.resize(m.N_se); // sNodes is always the sliding nodes in 2D
			sNodes << m.slidingEdgeNodes;
		}
		N_s = sNodes.size();
		N_mPro = mNodesPro.size();

	}else if(m.smode=="ds"){

		if(m.pmode == "moving"){
			mNodes.resize(m.N_ib);
			mNodes << m.intBdryNodes;
			sNodes.resize(m.N_se+m.N_es);
			sNodes << m.slidingEdgeNodes, m.extStaticNodes;
		}else{
			mNodes.resize(m.N_ib+m.N_es);
			mNodes << m.intBdryNodes, m.extStaticNodes;
			sNodes.resize(m.N_se);
			sNodes << m.slidingEdgeNodes;
		}
		if(curved){
			// todo rename to make clearer
			mNodes2.resize(m.N_ib+m.N_es+m.N_se);
			mNodes2 << m.intBdryNodes,m.extStaticNodes, m.slidingEdgeNodes;
			N_m2 = mNodes2.size();

		}
	}
	N_i = iNodes.size();
	N_s = sNodes.size();
	N_m = mNodes.size();
}


void rbf::getPeriodicParams(const std::string& periodicDirection){
	pVec.resize(m.nDims);
	pnVec.resize(m.nDims);
	if(m.nDims==2){
		if(m.pmode != "none"){
			if(periodicDirection == "x"){
				pVec << 1,0;
				pnVec << 0,1;
			}
			else if(periodicDirection == "y"){
				pVec << 0,1;
				pnVec << 1,0;
			}
		}
		else{
			pVec << 0,0;
		}

	// todo do the same for 3D
	}
	else if(m.nDims==3){

	}

//	std::cout << "vector in the periodic direction is: \n" << pVec << std::endl;

}

double rbf::rbfEval(double distance){
//	double xi = distance/m.r;	// distance scaled by support radius

	double f_xi = pow((1-(distance/m.r)),4)*(4*(distance/m.r)+1);
	return f_xi;
}


void rbf::project(Eigen::ArrayXXd& delta,Eigen::ArrayXXd& finalDef){
	std::cout << "Doing Projection" << std::endl;

	Eigen::RowVectorXd d;
	Eigen::ArrayXd dist, projection;
	Eigen::ArrayXi index = Eigen::ArrayXi::LinSpaced(m.midPnts.rows(),0,m.midPnts.rows()-1);

	for(int i=0; i<N_s; i++){
		if(i<m.N_se){

			// distance to all midpoints
			dist = (m.midPnts.rowwise()-(m.coords.row(sNodes(i)) + delta.row(i))).rowwise().norm();

			// finding closest midpoint by sorting them in ascending order
			std::sort(index.begin(), index.end(),[&](const int& a, const int& b) {
				return (dist[a] < dist[b]);
				}
			);

			// vector with distance from closest midpoint to displaced sliding node.
			d = m.midPnts.row(index(0)) - (m.coords.row(sNodes(i)) + delta.row(i));

			// Projection required to bring node back on the external boundary
			projection = d.dot(m.midPntNormals.row(index(0)).matrix())*m.midPntNormals.row(index(0));

			// Final deformation is the initial displacement plus the projection to bring the node back on the boundary
			finalDef.row(i) = delta.row(i) + projection.transpose();
		}else{
//			std::cout << sNodes(i) << '\t' << i << std::endl;
//			std::cout << delta.row(i) << std::endl;
//			std::cout << pVec << std::endl;
//			std::cout << delta.row(i)*pVec.transpose().array() << std::endl;
			finalDef.row(i) = delta.row(i)*pVec.transpose().array();

		}
	}
	std::cout << "projection is done " << std::endl;
}


