/*
 * Levelset2d.h
 *
 *  Receive a point and dertermine if it is inside of detected grain
 *
 */
#ifndef LEVELSET2D_H_
#define LEVELSET2D_H_

#include "definitions.h"
class Levelset2d {

public:
	// Default constructor
	Levelset2d() {
		_xdim = 0; _ydim = 0;
	}
	// Constructor to save a grain LS values
	Levelset2d(const vector<double> & levelset, const unsigned int & xdim, const unsigned int & ydim):
		_levelset(levelset), _xdim(xdim), _ydim(ydim) {
		if (_xdim*_ydim != _levelset.size()) {
			cout << "ERROR: levelset size not consistent with dimensions" << endl;
		}
	}

	// Checks if there is penetration, if there is, finds the penetration amount and the normalized gradient
	bool isPenetration(const Vector2d & point, double & value, Vector2d & gradient) const {
		double x = point(0);
		double y = point(1);
		// Case 1: This point is outside the grid, return false
		if (x + 1 > (double)_xdim || y + 1 > (double)_ydim || x < 0 || y < 0){
			return false;
		}
		unsigned int xr = (unsigned int) round(x);
		unsigned int yr = (unsigned int)round(y);
		// Case 2: Phi value in the node closest to this point is larger than 1, return false
		if (getGridValue(xr, yr) > 1) {
			return false;
		}
		// Case 3: If inside the level set, do bi-linear interpolation
		unsigned int xf = (unsigned int) floor(x);
		unsigned int yf = (unsigned int) floor(y);
		unsigned int xc = (unsigned int) ceil(x);
		unsigned int yc = (unsigned int) ceil(y);
		double dx 		= x - xf;
		double dy 		= y - yf;
		double b1 		= getGridValue(xf, yf);
		double b2 		= getGridValue(xc, yf) - b1;
		double b3 		= getGridValue(xf, yc) - b1;
		double b4 		= -b2 - getGridValue(xf, yc) + getGridValue(xc, yc);
		value = b1 + b2*dx + b3*dy + b4*dx*dy;
		// Case 3.1: Phi value in this point (x, y) is actually larger than 0, return false
		if (value > 0) {
			return false;
		}
		// Case 3.2: Phi value in this point (x, y) is actually smaller than 0, return true and gradient
		gradient << b2 + b4*dy, b3 + b4*dx;
		gradient /= gradient.norm();
		return true;
	}



//	double getXdim() const {
//		return _xdim;
//	}
//	double getYdim() const {
//		return _ydim;
//	}
//	vector<double> getLevelset() const {
//		return _levelset;
//	}



private:

	double getGridValue(unsigned int & x, unsigned int & y) const {
		return _levelset[y*_xdim + x];
	}
	vector<double> _levelset;
	unsigned int _xdim;
	unsigned int _ydim;


};

#endif /* LEVELSET2D_H_ */


//void shearVertical(const double & gamma) {
//	for (size_t j = 0; j < _ydim; j++) {
//		for (size_t i = 0; i < _xdim; i++) {
////				_levelset[j*_xdim + i] += ((double)j - (double)_ydim/2.)*gamma;
//			_levelset[j*_xdim + i] += ((double)j - 45.)*gamma;
//		}
//	}
//}
//
//void shearHorizontal(const double & gamma) {
//	for (size_t j = 0; j < _ydim; j++) {
//		for (size_t i = 0; i < _xdim; i++) {
////				_levelset[j*_xdim + i] += ((double)i - (double)_xdim/2.)*gamma;
//			_levelset[j*_xdim + i] += ((double)i - 45.)*gamma;
//		}
//	}
//}
