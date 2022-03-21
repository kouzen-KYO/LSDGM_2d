/*
 * World2dBiaxial.h
 *
 * Used to provide shear and compression effects
 *
 */

#ifndef World2dBiaxial_H_
#define World2dBiaxial_H_

#include "definitions.h"
#include "Wall2d.h"
#include "Wall2dPeriodicx.h"
#include "WorldStates.h"

class World2dBiaxial {

public:
	World2dBiaxial() {
		_gDamping = 0.; _dt = 0.; _ngrains = 0;
	}
	World2dBiaxial(const vector<Grain2d> & grains, vector<Wall2d> & walls, const double & gDamping, const double & dt):
		_grains(grains), _walls(walls), _gDamping(gDamping), _dt(dt) {
		_globalState = WorldState2d(_grains.size(), _walls.size() );
		_ngrains = _grains.size();
	}

	World2dBiaxial(const vector<Grain2d> & grains, vector<Wall2d> & walls, const Wall2dPeriodicx & wallPeriodicx, const double & gDamping, const double & dt):
		_grains(grains), _walls(walls), _wallPeriodicx(wallPeriodicx), _gDamping(gDamping), _dt(dt) {
		_globalState = WorldState2d(_grains.size(), _walls.size() );
		_ngrains = _grains.size();
	}

	// Compute the snapshot of the world (all grain forces/moments, all wall forces, total stress)
	void computeWorldState() {
		_globalState.reset();
		_ngrains = _grains.size();
		#pragma omp parallel default(none) shared(cout) // num_threads(7)
		{
			// Initialize local states which are reduced to the global state
			WorldState2d localState(_grains.size(), _walls.size() );
			localState.reset();
			Vector4d forcemoment;
			Vector2d force;
			double moment;
			Vector2d cmvec;
			Vector3d stress;
			size_t ncontacts;
			// Calculate forces/moments/stresses
			#pragma omp for schedule(dynamic, 100) //shared(cout)
			for (size_t i = 0; i < _grains.size(); i++) {
				// Interparticle
				for (size_t j = i + 1; j < _grains.size(); j++) {
					if (_grains[i].bcircleGrainIntersection(_grains[j])) {
						forcemoment = _grains[i].findInterparticleForceMoment(_grains[j], _dt);
						if (forcemoment(0) != 0.0 || forcemoment(1) != 0.0) {
							force << forcemoment(0), forcemoment(1);
							cmvec = _grains[j].getPosition() - _grains[i].getPosition();
							localState.grainForces[i] += force;
							localState.grainForces[j] -= force;
							localState.stressVoigt(0) += force(0)*cmvec(0);
							localState.stressVoigt(1) += force(1)*cmvec(1);
							localState.stressVoigt(2) += 0.5*(force(1)*cmvec(0) + force(0)*cmvec(1));
							localState.grainMoments[i] += forcemoment(2);
							localState.grainMoments[j] += forcemoment(3);
						}
					}
				}
				// Fixed walls
				_grains[i].changeWallContact(INT_MAX);
				for (size_t j = 0; j < _walls.size(); j++) {
					if (_walls[j].bcircleWallIntersection(_grains[i])) {
						if (_walls[j].findWallForceMoment(_grains[i], force, moment, ncontacts, stress, _dt)) {
							_grains[i].changeWallContact(j);
							localState.wallContacts[j] += ncontacts;
							// Compute things related to forces
							localState.grainForces[i] += force;
							localState.wallForces[j] -= force;
							// Compute things related to moments
							localState.grainMoments[i] += moment;
							localState.stressVoigt += stress;
						}
					}
				}
				// Periodic walls
				if (_wallPeriodicx.exists()) {
					// Move grain by width of wall and find contacts
					Grain2d grain(_grains[i]);
					grain.moveGrain(Vector2d(_wallPeriodicx.getWidth(), 0.) );
					for (size_t j = 0; j < _grains.size(); j++) {
						if (grain.bcircleGrainIntersection(_grains[j]) && i != j) {
							forcemoment = grain.findInterparticleForceMoment(_grains[j], _dt);
							if (forcemoment(0) != 0.0 || forcemoment(1) != 0.0) {
								force << forcemoment(0), forcemoment(1);
//								cout << "grains " << i << " and " << j << " with force " << force.transpose() << endl;
								cmvec = _grains[j].getPosition() - grain.getPosition();
								localState.grainForces[i] += force;
								localState.grainForces[j] -= force;
								localState.stressVoigt(0) += force(0)*cmvec(0);
								localState.stressVoigt(1) += force(1)*cmvec(1);
								localState.stressVoigt(2) += 0.5*(force(1)*cmvec(0) + force(0)*cmvec(1));
								localState.grainMoments[i] += forcemoment(2);
								localState.grainMoments[j] += forcemoment(3);
							}
						}
					}
				} // End if statement for periodic wall
			}
			// Reduce the localStates to the globalState
			#pragma omp critical
			{
				_globalState += localState;
			}
		}// Closes omp parallel
	}

	CData computeCstate() const {
		CData cDataRank;
		#pragma omp parallel default(none) shared(cout,cDataRank) //num_threads(1)
		{
			CData cDataThread;
			#pragma omp for schedule(dynamic, 40)
			for (size_t i = 0; i < _ngrains; i++) {
				// Loop over grains
				for (size_t j = i + 1; j < _ngrains; j++) {
					if (_grains[i].bcircleGrainIntersection(_grains[j])) {
						CData cDataContact = _grains[i].findContactData(_grains[j]);
						if (cDataContact._clocs.size() > 0) {
							cDataThread += cDataContact;
						}
					}
				} // Close grain subloop
				// Loop over walls
				for (size_t j = 0; j < _walls.size(); j++) {
					if (_walls[j].bcircleWallIntersection(_grains[i])) {
						CData cDataContact = _walls[j].findContactData(_grains[i]);
						if (cDataContact._clocs.size() > 0) {
							cDataThread += cDataContact;
						}
					}
				} // Close wall subloop
			} // End loop over grains
			#pragma omp critical
			{
				cDataRank += cDataThread;
			}
		} // Closes openmp parallel section
		return cDataRank;
	} // End computeCState method



	void takeTimestep() {
		#pragma omp for schedule(dynamic, 100)
		for (size_t grainid = 0; grainid < _grains.size(); grainid++) {
			_grains[grainid].takeTimestep(_globalState.grainForces[grainid], _globalState.grainMoments[grainid], _gDamping, _dt);

			if (_wallPeriodicx.exists()) {
				if (_grains[grainid].getPosition()(0) > _wallPeriodicx.getX() + _wallPeriodicx.getWidth() ) {
					_grains[grainid].moveGrain(Vector2d(-_wallPeriodicx.getWidth(),0.));
//					cout << "moving grainid " << grainid << " to left side" << endl;
				}
				if (_grains[grainid].getPosition()(0) < _wallPeriodicx.getX() ) {
					_grains[grainid].moveGrain(Vector2d(_wallPeriodicx.getWidth(),0.));
//					cout << "moving grainid " << grainid << " to right side" << endl;
				}
			}
		}
	}

	void takeTimestepShear(double v) {
		#pragma omp for schedule(dynamic, 40)
		for (size_t grainid = 0; grainid < _grains.size(); grainid++) {

			if (_grains[grainid].getWallContact() == INT_MAX) {
				_grains[grainid].takeTimestep(_globalState.grainForces[grainid], _globalState.grainMoments[grainid], _gDamping, _dt);
			}
			else if (_grains[grainid].getWallContact() == 0) { // bottom wall
				Vector2d vel(-v,0);
				// _globalState.grainForces[grainid]
				_grains[grainid].changeVel(vel);
				_grains[grainid].changeOmega(0);
				_grains[grainid].takeTimestep(Vector2d(0,0), 0, 0, _dt);
			}
			else { // top wall
				Vector2d vel(v,0);
				_grains[grainid].changeVel(vel);
				_grains[grainid].changeOmega(0);
				_grains[grainid].takeTimestep(Vector2d(0,0), 0, 0, _dt);
			}

			if (_wallPeriodicx.exists()) {
				if (_grains[grainid].getPosition()(0) > _wallPeriodicx.getX() + _wallPeriodicx.getWidth() ) {
					_grains[grainid].moveGrain(Vector2d(-_wallPeriodicx.getWidth(),0.));
//					cout << "moving grainid " << grainid << " to left side" << endl;
				}
				if (_grains[grainid].getPosition()(0) < _wallPeriodicx.getX() ) {
					_grains[grainid].moveGrain(Vector2d(_wallPeriodicx.getWidth(),0.));
//					cout << "moving grainid " << grainid << " to right side" << endl;
				}
			}
		}
	}

	void addGravityForce(const double & grav) {
		for (size_t grainid = 0; grainid < _grains.size(); grainid++) {
			_globalState.grainForces[grainid](1) -= grav*_grains[grainid].getMass();
		}
	}

	void moveProportional(const double & x, const double & y, const double & w, const double & h) {
		#pragma omp for schedule(dynamic, 100)
		for (size_t grainid = 0; grainid < _grains.size(); grainid++) {
			double dx = x*_grains[grainid].getPosition()(0)/w;
			double dy = y*_grains[grainid].getPosition()(1)/h;
			_grains[grainid].moveGrain(Vector2d(dx,dy));
		}
	}

	const WorldState2d & getWorldState() const {
		return _globalState;
	}
	const vector<Grain2d> & getGrains() const {
		return _grains;
	}
	const vector<Wall2d> & getWalls() const {
		return _walls;
	}
	const Wall2dPeriodicx & getWallPeriodic() const {
		return _wallPeriodicx;
	}

	// non-const methods
	vector<Grain2d> & getGrainsNC() {
		return _grains;
	}
	vector<Wall2d> & getWallsNC() {
		return _walls;
	}
	Wall2dPeriodicx & getWallPeriodicNC() {
		return _wallPeriodicx;
	}

private:

	vector<Grain2d> _grains;
	vector<Wall2d> _walls;
	Wall2dPeriodicx _wallPeriodicx;
//	Indenter _ind;
	size_t _ngrains;
	double _gDamping;
	double _dt;
	WorldState2d _globalState;
};



#endif /* World2dBiaxial_H_ */
