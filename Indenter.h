/*
 * Indenter.h
 *
 *  Created on: Apr 30, 2015
 *      Author: Reid
 */

#ifndef INDENTER_H_
#define INDENTER_H_

#include "definitions.h"

class Indenter {
	
public:
	// constructor
	Indenter() {
		_angle = 0; _radius = 0; _kn = 0; _id = 0; _slope = 0; _tipHeight = 0;
	}
	Indenter(const Vector2d & tipPos, const double & radius, const double & angle, const double & kn, const Vector2d & velocity, const size_t & id): 
	_tipPos(tipPos), _radius(radius), _angle(angle), _kn(kn), _velocity(velocity), _id(id) {
		_tipHeight = _radius/tan(_angle);
		_slope = tan(_angle);
	}
	 
	// checks for penetration of a single point against the intender
	// penetration is negative, normal is the outward normal of the intender (to conform to convention of everything else)
	bool isPenetration(const Vector2d & pt, double & penetration, Vector2d & normal) const {
		// check against left wall: further right than left wall, further left than center, further up than base of tip
		if (pt(0) > _tipPos(0)-_radius && pt(0) < _tipPos(0) && pt(1) > _tipPos(1) + _tipHeight) {
			normal << -1, 0;
			penetration = _tipPos(0)-_radius - pt(0);
			return true;
		}
		// check against right wall: further left than right wall, further right than center, further up than base of tip (this is somewhat redundant but idgaf)
		else if (pt(0) < _tipPos(0)+_radius && pt(0) > _tipPos(0) && pt(1) > _tipPos(1) + _tipHeight ) {
			normal << 1,0;
			penetration = pt(0) - _tipPos(0)-_radius;
			return true;
		}
		// check left side of tip: over the line carved by the left side of the tip and inside the bounding box of the tip's left slope
		else if ( (_slope*pt(0) + pt(1) - _tipPos(1) - _slope*_tipPos(0))/sqrt(_slope*_slope+1) > 0 &&
					 pt(0) > _tipPos(0) - _radius && pt(0) < _tipPos(0) &&
					 pt(1) < _tipPos(1) + _tipHeight && pt(1) > _tipPos(1)
			) {
			normal << -cos(_angle), -sin(_angle);
			penetration = -(_slope*pt(0) + pt(1) - _tipPos(1) - _slope*_tipPos(0))/sqrt(_slope*_slope+1);
			
			return true;
		}
		// check right side of tip
		else if ( (-_slope*pt(0) + pt(1) - _tipPos(1) + _slope*_tipPos(0))/sqrt(_slope*_slope+1) > 0 &&
				 	 pt(0) < _tipPos(0) + _radius && pt(0) > _tipPos(0) &&
				 	 pt(1) < _tipPos(1) + _tipHeight && pt(1) > _tipPos(1)
			) {
			normal << cos(_angle), -sin(_angle);
			penetration = -(-_slope*pt(0) + pt(1) - _tipPos(1) + _slope*_tipPos(0))/sqrt(_slope*_slope+1);
			
			return true;
		}
		
		return false;
	}
	
	void moveTip(const double & amount) {
		_tipPos(1) += amount;
	}
	
	const Vector2d & getTipPos() const {
		return _tipPos;
	}
	const double & getRadius() const {
		return _radius;
	}
	const double & getAngle() const {
		return _angle;
	}
	const double & getKn() const {
		return _kn;
	}
	const Vector2d & getVelocity() const {
		return _velocity;
	}
	const size_t & getId() const {
		return _id;
	}
	
	
private:
	Vector2d _tipPos;
	double _radius;
	double _angle;	// angle of tip's arc with the vertical
	double _kn;
	Vector2d _velocity;
	size_t _id;
	double _tipHeight;
	double _slope;	// (positive) slope of tip
	
};



#endif /* INDENTER_H_ */
