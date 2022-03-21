/*
 * Wall2dPeriodicx.h
 *
 * Create a periodic, virtual wall
 *
 */

#ifndef WALL2DPERIODICX_H_
#define WALL2DPERIODICX_H_

#include "definitions.h"
#include "Grain2d.h"

class Wall2dPeriodicx {
public:
	// Create a default wall
	Wall2dPeriodicx() {_exists = false; _width = 0; _x = 0;}
	// Create a periodica wall start from x to x + width
	Wall2dPeriodicx(const double & x, const double & width):
		_x(x), _width(width) {
		_exists = true;
	}
	// Change width directly
	void changeWidth(double newWidth) {
		_width = newWidth;
	}
	// Change width by adding a length
	void moveWidth(double amt) {
		_width += amt;
	}
	// Switch to open or close the periodic wall
	const bool & exists() const {
		return _exists;
	}
	const double & getX() const {
		return _x;
	}
	const double & getWidth() const {
		return _width;
	}

private:
	bool _exists;
	double _x;
	double _width;
};


#endif /* WALL2DPERIODIC_H_ */
