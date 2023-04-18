#include "CoordTransform.h"

CoordTransform::CoordTransform()
{}


void CoordTransform::cart_to_polar_cylindrical(Eigen::ArrayXXd& in, Eigen::ArrayXXd& out){

	out.col(0) = in.leftCols(2).rowwise().norm();
	for(int row = 0; row < in.rows(); row++){
		out(row,1) = atan2(in(row,1), in(row,0));
	}
	if(in.cols() == 3){
		out.col(2) = in.col(2);
	}

}

void CoordTransform::polar_cylindrical_to_cart(Eigen::ArrayXXd& in, Eigen::ArrayXXd& out){

	out.col(0) = in.col(0)*cos(in.col(1));
	out.col(1) = in.col(0)*sin(in.col(1));
	if(in.cols() == 3){
		out.col(2) = in.col(2);
	}

//	std::cout << m.coords << std::endl;
}

void CoordTransform::vector_cart_to_polar_cylindrical(Eigen::ArrayXXd& in, Eigen::ArrayXXd& out, Eigen::ArrayXi& idx, Eigen::ArrayXXd& coords){

	Eigen::ArrayXXd initPos(in.rows(),in.cols());
	Eigen::ArrayXXd newPos(in.rows(), in.cols());

	initPos = coords(idx,Eigen::all);
	newPos = coords(idx, Eigen::all) + in;
	// in contains initial position in polar/spherical
	cart_to_polar_cylindrical(initPos, in);

	// out contains the new position in polar/speherical coordinates
	cart_to_polar_cylindrical(newPos, out);

	out -= in;

}

void CoordTransform::error_to_cart(Eigen::ArrayXXd& error, Mesh* m, getNodeType& n){
	Eigen::ArrayXXd updatedPosition(error.rows(), error.cols());
	updatedPosition = (*m).coords_polar_cylindrical((*n.iPtr)(Eigen::seqN(0,error.rows())), Eigen::all) + error;

	Eigen::ArrayXXd updatedPositionCartesian(error.rows(), error.cols());
	polar_cylindrical_to_cart(updatedPosition, updatedPositionCartesian);

	error = updatedPositionCartesian - (*m).coords((*n.iPtr)(Eigen::seqN(0,error.rows())), Eigen::all);

}





