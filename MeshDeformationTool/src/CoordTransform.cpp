#include "CoordTransform.h"

CoordTransform::CoordTransform()
{}


void CoordTransform::cart_to_polar_spherical(Eigen::ArrayXXd& in, Eigen::ArrayXXd& out){

	out.col(0) = in.rowwise().norm();
	for(int row = 0; row < in.rows(); row++){
		out(row,1) = atan2(in(row,1), in(row,0));
	}
	if(in.cols() == 3){
		out.col(2) = in.col(2);
	}

//	std::cout <<m.coords_polar_spherical << std::endl;
}

void CoordTransform::polar_spherical_to_cart(Eigen::ArrayXXd& in, Eigen::ArrayXXd& out){

	out.col(0) = in.col(0)*cos(in.col(1));
	out.col(1) = in.col(0)*sin(in.col(1));
	if(in.cols() == 3){
		out.col(2) = in.col(2);
	}

//	std::cout << m.coords << std::endl;
}

void CoordTransform::vector_cart_to_polar_spherical(Eigen::ArrayXXd& in, Eigen::ArrayXXd& out, Eigen::ArrayXi& idx, Eigen::ArrayXXd& coords){
	/*
	Eigen::ArrayXd theta(in.rows());
	for(int i = 0; i < in.rows(); i++){
		theta(i) = atan2(coords(idx(i),1),coords(idx(i),0));
	}
	std::cout << theta << std::endl;
	out.col(0) = in.col(0)*cos(theta) + in.col(1)*sin(theta);
	out.col(1) = in.col(1)*cos(theta) - in.col(0)*sin(theta);

	if(in.cols() == 3){
		out.col(2) = in.col(2);
	}
	*/
	Eigen::ArrayXXd initPos(in.rows(),in.cols());
	Eigen::ArrayXXd newPos(in.rows(), in.cols());

	initPos = coords(idx,Eigen::all);
	newPos = coords(idx, Eigen::all) + in;
	// in contains initial position in polar/spherical
	cart_to_polar_spherical(initPos, in);

	// out contains the new position in polar/speherical coordinates
	cart_to_polar_spherical(newPos, out);

	out -= in;

}





