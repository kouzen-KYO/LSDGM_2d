/*
 * Grain2d.h
 *
 *  2D grain class including basic properties and some functional methods such as contact calculation
 *
 */
#ifndef GRAIN2D_H_
#define GRAIN2D_H_

#include "definitions.h"
#include "WorldStates.h"
//#include "Wall2dFixed.h"
//#include "Indenter.h"

struct UtkContactInfo {//Contact information container

	UtkContactInfo(const size_t & ngrains) {
//		_ncontacts = 0;
		_ngrains = ngrains;
		_sForceVec.resize(ngrains);
		_nForceVec.resize(ngrains);
		_grainsVec.resize(ngrains);
		_cLocVec.resize(ngrains);
	}

	void clear() {
		for (size_t g = 0; g < _ngrains; g++) {
			_sForceVec[g].clear();
			_nForceVec[g].clear();
			_grainsVec[g].clear();
			_cLocVec[g].clear();
		}
		_sForceVec.reserve(7);
		_nForceVec.reserve(7);
		_grainsVec.reserve(7);
		_cLocVec.reserve(7);
	}

	size_t _ngrains;		// Number of grains in assembly
//	size_t _ncontacts;	// Number of contacts in assembly
	vector<vector<Vector2d> > _sForceVec;	// Normal force
	vector<vector<Vector2d> > _nForceVec;	// Shear force
	vector<vector<size_t> > _grainsVec;		// Grain indices of contact
	vector<vector<Vector2d> > _cLocVec;		// Contact location
};

class Grain2d {
public:
	// Constructors
	Grain2d() {
		_radius = 0; _rsq = 0; _kn = 0; _ks = 0; _mass = 0; _momentInertia = 0; _mu = 0; _omega = 0; _theta = 0; _id = 0; _density = 1; _ncontacts = 0; _wallContact = INT_MAX;
		_corr = 1;
	}

	Grain2d(const double & mass, const Vector2d & position, const Vector2d & velocity,
			  const double & momentInertia, const double & theta, const double & omega,
			  const Vector2d & cmLset, const vector<Vector2d> & pointList,
			  const double & radius, const shared_ptr<Levelset2d> lset, const double & kn, const double & ks, const double & mu, const size_t & id):
			  _mass(mass), _position(position), _velocity(velocity),  _momentInertia(momentInertia),
			  _theta(theta), _omega(omega), _cmLset(cmLset),
			  _pointList(pointList), _radius(radius), _lset(lset), _kn(kn), _ks(ks), _mu(mu), _id(id) {
		_rsq = _radius*_radius;
		_ncontacts = 0;
		_density = 1;
		_wallContact = INT_MAX;
		_nodeShears.resize(_pointList.size());
		_nodeContact.resize(_pointList.size());
		_corr = 1.;
		for (size_t i = 0; i < _pointList.size(); i++) {
			_nodeShears[i] = 0;
			_nodeContact[i] = 0;
		}
		// assemble kd tree
//		vector<Vector2d> tempPtList;
//		_pointChildren.resize(_pointList.size());
//		createKD(_pointList, tempPtList, _pointChildren, 0);
//		_pointList = tempPtList;
	}

//	size_t createKD(const vector<Vector2d> & ptList, vector<Vector2d> & newPtList, vector<Vector3i> & children, const int & depth) {
//		int axis = depth%2;
//		vector<std::pair<double,size_t> > temp(ptList.size());
//		for (size_t i = 0; i < temp.size(); i++) {
//			temp[i] = std::pair<double,size_t>(ptList[i](axis), i);
//		}
//		std::nth_element(temp.begin(), temp.begin() + temp.size()/2, temp.end());
//		size_t ptListIdx = temp[temp.size()/2].second;
//		newPtList.push_back( ptList[ ptListIdx ] );
//		size_t curIdx = newPtList.size()-1;
//		children[curIdx](2) = axis;
//		// left/bottom child
//		if (temp.size()/2 > 0) {
//			vector<Vector2d> vLeft(temp.size()/2);
//			for (size_t i = 0; i < temp.size()/2; i++) {
//				vLeft[i] = ptList[ temp[i].second];
//			}
//			children[curIdx](0) = createKD(vLeft, newPtList, children, depth+1);
//		}
//		else {
//			children[curIdx](0) = -1;
//		}
//		// right/top child
//		if (temp.size()/2+1 < temp.size()) {
//			vector<Vector2d> vRight( temp.size() - temp.size()/2 - 1 );
//			for (size_t i = temp.size()/2+1; i < temp.size(); i++) {
//				vRight[i-(temp.size()/2+1)] = ptList[ temp[i].second];
//			}
//			children[curIdx](1) = createKD(vRight, newPtList, children, depth+1);
//		}
//		else {
//			children[curIdx](1) = -1;
//		}
//		return curIdx;
//	}

	// Checks if bounding circles intersect between this and other
	bool bcircleGrainIntersection(const Grain2d & other) const {
		// First-order contact check based on the bounding radii
		double rsum = other.getRadius() + _radius;
		return (other.getPosition() - _position).squaredNorm() < rsum*rsum;
	}


	// Finds interparticle force and moment based on the points of this and the level set of other
	Vector4d findInterparticleForceMoment(const Grain2d & other, const double & dt) {
		// Declare temporary variables
		Vector4d		forceMoment(0.0, 0.0, 0.0, 0.0); // [fx, fy, moment on this, moment on other]
		Vector2d		ptOtherCM; 		// Point wrt the center of mass of other in real space
		Vector2d		ptOtherLset; 	// Point in the reference config of other's level set
		double			penetration;	// Penetration amount (is negative by convention of the level set)
		Vector2d		normal; 		// Surface normal
		const double cos2 = cos(other.getTheta());
		const double sin2 = sin(other.getTheta());
		Vector2d		ptThisCM; 		// Point wrt the center of mass of this in real space
		Vector2d		df; 			// Force increment from a single point
		Vector2d		tangent; 		// Surface tangent
		double			sdot; 			// Relative velocity of a point of contact
//		double		ds;
		Vector2d		Fs; 			// vector of frictional/shear force
		Vector2d        relativeVel;
		double 			Fsmag;
//		double cres = .6;
		double GamaN = -2*sqrt(_kn*_mass*other.getMass()/(_mass+other.getMass()))*log(_corr)/sqrt(M_PI*M_PI + log(_corr)*log(_corr));
//		double GamaN = 0;

		for (size_t ptidx = 0; ptidx < _pointList.size(); ptidx++) {
			ptOtherCM = _pointList[ptidx] - other.getPosition();
			// Check point against bounding radius of other
			if ( ptOtherCM.squaredNorm() < other.getRsq() ) {
				ptThisCM = _pointList[ptidx] - _position;
				ptOtherLset(0) =  ptOtherCM(0)*cos2 + ptOtherCM(1)*sin2;
				ptOtherLset(1) = -ptOtherCM(0)*sin2 + ptOtherCM(1)*cos2;
				ptOtherLset += other.getCmLset();
				//	If there is penetration, finds forces and moments
				if ( other.getLset()->isPenetration(ptOtherLset, penetration, normal) ) {
					// rotate the normal from the reference config of 2's level set to real space
					normal << normal(0)*cos2 - normal(1)*sin2, normal(0)*sin2 + normal(1)*cos2;
					// Find the tangent
					tangent << -normal(1), normal(0);
					// Update force and moments: normal force contribution
					relativeVel << other.getVelocity()(0) - other.getOmega()*ptOtherCM(1) - (_velocity(0) - _omega*ptThisCM(1)),
										other.getVelocity()(1) + other.getOmega()*ptOtherCM(0) - (_velocity(1) + _omega*ptThisCM(0));
					Vector2d relativeVelNormal = normal.dot(relativeVel)*normal;
					df = penetration*normal*_kn - GamaN*relativeVelNormal; // normal.dot(relativeVel)*normal
					forceMoment.block<2,1>(0,0) -= df;
					forceMoment(3) -= df(0)*ptOtherCM(1) - df(1)*ptOtherCM(0);
					forceMoment(2) += df(0)*ptThisCM(1) - df(1)*ptThisCM(0);
					// force/moment calculations based on friction
					//sdot = tangent.dot(_velocity - other.getVelocity()) - (_omega*ptThisCM.norm() + other.getOmega()*ptOtherCM.norm());
					sdot = tangent.dot(-relativeVel);
//					ds = sdot*dt;
					_nodeContact[ptidx] = other.getId();
					_nodeShears[ptidx] += _ks*sdot*dt;
					if (_nodeShears[ptidx] > 0) {
						Fsmag = min(_nodeShears[ptidx], df.norm()*_mu );
					}
					else {
						Fsmag = max(_nodeShears[ptidx], -df.norm()*_mu );
					}
					_nodeShears[ptidx] = Fsmag; // Cuts the length of the spring if _ks*ds > Fsmag (slip)

					Fs = -tangent*Fsmag;
					// Damp Fs with COR
					Fs -= GamaN*(relativeVel-relativeVelNormal);
					forceMoment.block<2,1>(0,0) += Fs;
					forceMoment(2) -= Fs(0)*ptThisCM(1) - Fs(1)*ptThisCM(0);
					forceMoment(3) += Fs(0)*ptOtherCM(1) - Fs(1)*ptOtherCM(0);
				}
				else if (_nodeContact[ptidx] == other.getId() ){ // If this surface point is inside other grain, clear it
					_nodeContact[ptidx] = 0;
					_nodeShears[ptidx] = 0;
	//				_nodeNormals[ptidx] << 0,0;
				}
			}
		} // End loop over points
		return forceMoment;
	}

//	// note: doesn't work.  need to take into account theta when performing left/right/up/down checks.  Also, it's slower than checking all of the points.
//	Vector4d findInterparticleForceMoment2(const Grain2d & other, const double & dt) {
//		// declare temporary variables
//		Vector4d		forceMoment(0.0, 0.0, 0.0, 0.0); // [fx, fy, moment on this, moment on other]
//		Vector2d		ptOtherCM; 		// point wrt the center of mass of other in real space
//		Vector2d		ptOtherLset; 	// point in the reference config of other's level set
//		double		penetration;	// penetration amount (is negative by convention of the level set)
//		Vector2d		normal; 			// surface normal
//		const double cos2 = cos(other.getTheta());
//		const double sin2 = sin(other.getTheta());
//		Vector2d		ptThisCM; 		// point wrt the center of mass of this in real space
//		Vector2d		df; 				// force increment from a single point
//		Vector2d		tangent; 		// surface tangent
//		double		sdot; 			// relative velocity of a point of contact
//		Vector2d		Fs; 				// vector of frictional/shear force
//		Vector2d    relativeVel;
//		double 		Fsmag;
//		Vector3i		child;
//		vector<size_t> stack( int(ceil(log2(_pointList.size()))) );
//		stack[0] = 0;
//		int stackidx = 0;
//		bool isDone = false;
//		while (!isDone) {
//			size_t ptidx = stack[stackidx];
//			// look at point on top of the stack
//			ptOtherCM = _pointList[ptidx] - other.getPosition();
//			cout << _pointList[ptidx].transpose() << endl;
//			if ( ptOtherCM.squaredNorm()<other.getRsq() ) {
//				ptThisCM = _pointList[ptidx] - _position;
//				ptOtherLset(0) =  ptOtherCM(0)*cos2 + ptOtherCM(1)*sin2;
//				ptOtherLset(1) = -ptOtherCM(0)*sin2 + ptOtherCM(1)*cos2;
//				ptOtherLset += other.getCmLset();
//				if ( other.getLset()->isPenetration(ptOtherLset, penetration, normal) ) {
////					cout << _pointList[ptidx].transpose() << endl;
//					// rotate the normal from the reference config of 2's level set to real space
//					normal << normal(0)*cos2 - normal(1)*sin2, normal(0)*sin2 + normal(1)*cos2;
//					// find the tangent
//					tangent << -normal(1), normal(0);
//					// update force and moments: normal force contribution
//					relativeVel << other.getVelocity()(0) - other.getOmega()*ptOtherCM(1) - (_velocity(0) - _omega*ptThisCM(1)),
//										other.getVelocity()(1) + other.getOmega()*ptOtherCM(0) - (_velocity(1) + _omega*ptThisCM(0));
//
//					df = penetration*normal*_kn; // - GamaN*normal.dot(relativeVel)*normal;
//					forceMoment.block<2,1>(0,0) -= df;
//					forceMoment(3) -= df(0)*ptOtherCM(1) - df(1)*ptOtherCM(0);
//					forceMoment(2) += df(0)*ptThisCM(1) - df(1)*ptThisCM(0);
//					// force/moment calculations based on friction
//					sdot = tangent.dot(-relativeVel);
//					_nodeContact[ptidx] = other.getId();
//					_nodeShears[ptidx] += _ks*sdot*dt;
//					if (_nodeShears[ptidx] > 0) {
//						Fsmag = min(_nodeShears[ptidx], df.norm()*_mu );
//					}
//					else {
//						Fsmag = max(_nodeShears[ptidx], -df.norm()*_mu );
//					}
//					_nodeShears[ptidx] = Fsmag; // cuts the length of the spring if _ks*ds > Fsmag (slip)
//
//					Fs = -tangent*Fsmag;
//					forceMoment.block<2,1>(0,0) += Fs;
//					forceMoment(2) -= Fs(0)*ptThisCM(1) - Fs(1)*ptThisCM(0);
//					forceMoment(3) += Fs(0)*ptOtherCM(1) - Fs(1)*ptOtherCM(0);
//				}
//				else if (_nodeContact[ptidx] == other.getId() ){
//					_nodeContact[ptidx] = 0;
//					_nodeShears[ptidx] = 0;
//				}
//			} // end check for this current point
//
//
//			// remove point from stack as it has been checked
//			stackidx--;
//			child = _pointChildren[ptidx];
//			// check both children
//			if (child(0)>=0 && ptOtherCM(child(2))>-other.getRadius()) {
//				stackidx++;
//				stack[stackidx] = _pointChildren[ptidx](0);
//			}
//			if (child(1)>=0 && ptOtherCM(child(2))<other.getRadius()) {
//				stackidx++;
//				stack[stackidx] = _pointChildren[ptidx](1);
//			}
//			if (stackidx < 0) {
//				isDone = true;
//			}
//		} // end while loop
//
//		return forceMoment;
//	}

	CData findContactData(const Grain2d & other) const {
		CData cData;			// Output
		Vector2d force;			// Total force
		Vector2d ptThisCM; 		// Point wrt the center of mass of *this in real space
		Vector2d ptOtherCM; 	// Point wrt the center of mass of other in real space
		Vector2d ptOtherLset; 	// Point in the reference config of other's level set
		double	penetration;	// Penetration amount
		Vector2d normal; 		// Surface normal pointing out of other in the reference config (initially) and then pointing out of *this in real space (after)
		Vector2d tangent;
//		double Fsmag = 0;
		const double cos2 = cos(other.getTheta());
		const double sin2 = sin(other.getTheta());
		double Fsmag;

		// Iterate through all of the points/nodes of *this and check for contact for each one
		for (size_t ptidx = 0; ptidx < _pointList.size(); ptidx++) {
			// Move the surface point to the reference space of other grain
//			ptThisCM = _pointList[ptidx] - _position;
			ptOtherCM = _pointList[ptidx] - other.getPosition();
			ptOtherLset(0) =  ptOtherCM(0)*cos2 + ptOtherCM(1)*sin2;
			ptOtherLset(1) = -ptOtherCM(0)*sin2 + ptOtherCM(1)*cos2;
			ptOtherLset += other.getCmLset();
			// Calculate phi value this surface point and decide if continue or not
			if ( other.getLset()->isPenetration(ptOtherLset, penetration, normal) ) {
				normal << normal(0)*cos2 - normal(1)*sin2, normal(0)*sin2 + normal(1)*cos2; // outward normal of other
				tangent << -normal(1), normal(0);
				force = penetration*normal*_kn; // Penetration is negative, along with outward normal of other, this is the force of *this on other
				if (_nodeShears[ptidx] > 0) {
					Fsmag = min(_nodeShears[ptidx], force.norm()*_mu );
				}
				else {
					Fsmag = max(_nodeShears[ptidx], -force.norm()*_mu );
				}
				force += Fsmag*tangent;
				cData._cpairs.push_back(Vector2i(_id, other.getId()));
				cData._forces.push_back(force);
				cData._normals.push_back(normal);
				cData._clocs.push_back(_pointList[ptidx]);
				cData._nodes.push_back(ptidx);
			}
		}
		return cData;
	}

	double computeKineticEnergy() const {
		double ke = 0;
		ke += .5*_mass*_velocity.squaredNorm();
		ke += .5*_momentInertia*_omega*_omega;
		return ke;
	}


	void takeTimestep(const Vector2d & force, const double & moment, const double & gDamping, const double & dt) {
		_velocity = 1/(1+gDamping*dt/2)*( (1-gDamping*dt/2)*_velocity + dt*force/_mass   );
		_omega = 1/(1+gDamping*dt/2)*( (1-gDamping*dt/2)*_omega + dt*moment/_momentInertia);
		double cosd = cos(_omega*dt);
		double sind = sin(_omega*dt);
		// Must update the points
		for (size_t ptid = 0; ptid < _pointList.size(); ptid++) {
			_pointList[ptid] << (_pointList[ptid](0)-_position(0))*cosd - (_pointList[ptid](1)-_position(1))*sind,
									  (_pointList[ptid](0)-_position(0))*sind + (_pointList[ptid](1)-_position(1))*cosd;
			_pointList[ptid] += _position + _velocity*dt;
		}
		_position = _position + dt*_velocity;
		_theta = _theta + dt*_omega;
	}

	void moveGrain(const Vector2d & amount) {
		_position = _position + amount;
		for (size_t ptid = 0; ptid < _pointList.size(); ptid++) {
			_pointList[ptid] += amount;
		}
	}

	// Change methods
	void changeMu(const double & newmu) {
		_mu = newmu;
	}

	void changePos(const Vector2d & pos) {
		Vector2d disp = pos - _position;
		for (size_t ptid = 0; ptid < _pointList.size(); ptid++) {
			_pointList[ptid] += disp;
		}
		_position = pos;
	}

	void changeVel(const Vector2d & vel) {
		_velocity = vel;
	}

	void changeRot(const double & rot) {
		double dtheta = rot - _theta;
		_theta = rot;
		double cosd = cos(dtheta);
		double sind = sin(dtheta);
		for (size_t ptid = 0; ptid < _pointList.size(); ptid++) {
			_pointList[ptid] << (_pointList[ptid](0)-_position(0))*cosd - (_pointList[ptid](1)-_position(1))*sind,
									  (_pointList[ptid](0)-_position(0))*sind + (_pointList[ptid](1)-_position(1))*cosd;
			_pointList[ptid] += _position;
		}
	}

	void changeOmega(const double & om) {
		_omega = om;
	}

	void rotateGrain(const double & rot) {
		double dtheta = rot;
		_theta += rot;
		double cosd = cos(dtheta);
		double sind = sin(dtheta);
		for (size_t ptid = 0; ptid < _pointList.size(); ptid++) {
			_pointList[ptid] << (_pointList[ptid](0)-_position(0))*cosd - (_pointList[ptid](1)-_position(1))*sind,
									  (_pointList[ptid](0)-_position(0))*sind + (_pointList[ptid](1)-_position(1))*cosd;
			_pointList[ptid] += _position;
		}
	}

	void changeKn(const double & kn) {
		_kn = kn;
	}
	void changeKs(const double & ks) {
		_ks = ks;
	}
	void changeDensity(const double & density) {
		_mass *= density/_density;
		_momentInertia *= density/_density;
		_density = density;
	}
	void changeId(const size_t & id) {
		_id = id;
	}
	void changeWallContact(size_t w) {
		_wallContact = w;
	}
//	void changeShearHist(const size_t nodeid, )


	// Helper methods
	const double & getMass() const {
		return _mass;
	}
	const Vector2d & getPosition() const {
		return _position;
	}
	const Vector2d & getVelocity() const {
		return _velocity;
	}
	const double & getTheta() const {
		return _theta;
	}
	const double & getOmega() const {
		return _omega;
	}
	const Vector2d & getCmLset() const {
		return _cmLset;
	}
	const double & getRadius() const {
		return _radius;
	}
	const double & getRsq() const {
		return _rsq;
	}
	const shared_ptr<Levelset2d> & getLset() const {
		return _lset;
	}
	const vector<Vector2d> & getPointList() const {
		return _pointList;
	}
	const vector<Vector3i> & getPointChildren() const {
		return _pointChildren;
	}
	const double & getKn() const {
		return _kn;
	}
	const double & getKs() const {
		return _ks;
	}
	const double & getMu() const {
		return _mu;
	}
	const size_t & getId() const {
		return _id;
	}
	const size_t & getWallContact() const {
		return _wallContact;
	}
	const vector<double> & getNodeShear() const {
			return _nodeShears;
	}

	// Non const
	vector<double> & getNodeShearsNonConst() {
		return _nodeShears;
	}
	vector<size_t> & getNodeContactNonConst() {
		return _nodeContact;
	}

private:

	double 		_mass;
	Vector2d 	_position; // Location of center of mass in real space
	Vector2d 	_velocity;
	double 		_momentInertia;
	double 		_theta;
	double 		_omega;
	Vector2d		_cmLset; // Center of mass wrt the level set reference configuration
	vector<Vector2d> _pointList; // List of points comprising the grain in real space (translated and rotated)
	vector<Vector3i> _pointChildren; // List of children of points (kd tree implementation)
	double 		_radius;
	double		_rsq; 		// squared radius
	shared_ptr<Levelset2d> 	_lset; // Smart POINTER to a level set
	double		_kn;		// Normal stiffness
	double		_ks;		// Shear stiffness
	double		_mu;		// Interparticle friction
	double		_density;
	double		_corr; 		// Coefficient of resitution

	size_t _id;

	size_t _wallContact; // Index of wall the grain is touching (int_max = not touching a wall)

	size_t _ncontacts; // Number of contacts of the grain (walls, other grains, everything)

	vector<double>		_nodeShears;	// Shears at each node
	vector<size_t>		_nodeContact;  // Index of grain the node is contacting
//	vector<Vector2d>	_nodeNormals;
};

#endif /* Grain2D_H_ */
