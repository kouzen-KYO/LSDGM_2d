/*
 * WorldStates.h
 *
 * Save everything used in the simulations
 *
 */
#ifndef WORLDSTATES_H_
#define WORLDSTATES_H_

#include "definitions.h"

struct WorldState2d {
	WorldState2d() {}
	WorldState2d(const size_t & ngrains, const size_t & nwalls) {
		grainForces.resize(ngrains);
		grainMoments.resize(ngrains);
		wallForces.resize(nwalls);
		wallContacts.resize(nwalls);
	}
	WorldState2d(const Vector3d & svtemp, const vector<Vector2d> & gftemp, const vector<double> & gmtemp, const vector<Vector2d> & wtemp):
		stressVoigt(svtemp), grainForces(gftemp), grainMoments(gmtemp), wallForces(wtemp) {}

	Vector3d stressVoigt; // Macroscopic stress of assembly
	vector<Vector2d> grainForces; // Forces on grains
	vector<double> grainMoments; // Moments on grains
	vector<Vector2d> wallForces; // Forces on wall
	vector<size_t> wallContacts;

	void reset() {
		stressVoigt << 0,0,0;
		for (size_t i = 0; i < grainForces.size(); i++) {
			grainForces[i] << 0, 0;
			grainMoments[i] = 0;
		}
		for (size_t i = 0; i < wallForces.size(); i++) {
			wallForces[i] << 0, 0;
			wallContacts[i] = 0;
		}
	}

	void operator+=(const WorldState2d & w) {
		stressVoigt += w.stressVoigt;
		for (size_t i = 0; i < grainForces.size(); i++) {
			grainForces[i] += w.grainForces[i];
			grainMoments[i] += w.grainMoments[i];
		}
		for (size_t i = 0; i < wallForces.size(); i++) {
			wallForces[i] += w.wallForces[i];
			wallContacts[i] += w.wallContacts[i];
		}
	}
};

struct CData {
	// member variables
	vector<Vector2i> _cpairs;
	vector<size_t> _nodes; // nodes of _cpairs[i](0)
	vector<Vector2d> _forces;
	vector<Vector2d> _normals;
	vector<Vector2d> _clocs;

	void operator+=(const CData & c) {
		_cpairs.insert( _cpairs.end(),	c._cpairs.begin(),	c._cpairs.end());
		_nodes.insert(_nodes.end(),		c._nodes.begin(),		c._nodes.end());
		_forces.insert( _forces.end(),	c._forces.begin(),	c._forces.end());
		_normals.insert(_normals.end(),	c._normals.begin(),	c._normals.end());
		_clocs.insert(_clocs.end(),		c._clocs.begin(),		c._clocs.end());
	}
	size_t size() {
		return _cpairs.size();
	}
};



#endif /* WORLDSTATES_H_ */
