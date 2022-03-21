/*
 * main.cpp
 *
 *  Created on: May 30, 2014
 *      Author: Reid
 */
 
#include "definitions.h"
#include "Levelset2d.h"
#include "Grain2d.h"
#include "readInputFile.h"
#include "World2dBiaxial.h"
 
int main() {
	omp_set_dynamic(0);     // Explicitly disable dynamic threading
	omp_set_num_threads(6);
	// file name 
	string filename("caicos.dat");
	vector<Grain2d> grains = generateGrainsFromFile(filename);
	vector<Grain2d> temp = grains;
	for (size_t i = 0; i < 63; i++) {
		grains.insert(grains.end(), temp.begin(), temp.end());
	}
	size_t ngrains = grains.size();
	cout <<"Number of grains: " << ngrains << endl;
	// grain properties
	double mu = 0.0;
	double kn = 30000;
	double density = 1;
	// assembly properties
	double height = 645.3*2;
	double width = 645.3 + 70; 
	// world properties
	double gDamping = .2; // .4 or .2
	double dt = 0.2*(2*sqrt(.5*200/kn));
	Vector2d trans(545.07320305-634.37729014-5.5, 784.81266662-656.22997566+5.5);
	Vector2d trans2(0, height/2);
	
	vector<Vector2d> pos = readPositionFile("positions_caicos.dat",ngrains/2);
//	vector<Vector2d> temp2 = pos;
	pos.insert(pos.end(), pos.begin(), pos.end());
	vector<double> rot = readRotationFile("rotations_caicos.dat",ngrains/2);
//	vector<double> temp3 = rot;
	rot.insert(rot.end(), rot.begin(), rot.end());
	
	/* set things up */
	for (size_t i = 0; i < ngrains; i++) {
		grains[i].changePos(pos[i]-trans);
		if (i >= ngrains/2) {
			grains[i].changePos( grains[i].getPosition() + trans2 );
		}
		grains[i].changeRot(rot[i]);
		grains[i].changeMu(mu);
		grains[i].changeKn(kn);
		grains[i].changeDensity(density);
		grains[i].changeId(i);
	}
	
	// create walls
	vector<Wall2d> walls; // 0bot 1top
	walls.resize(2);
	// bottom wall
	walls[0] = Wall2d(Vector2d(0,0), Vector2d(0,1), kn, 0.9*kn, 0.0, INT_MAX-3);
	// top wall
	walls[1] = Wall2d(Vector2d(0,height), Vector2d(0,-1), kn, 0.9*kn, 0.0, INT_MAX-2);
	Wall2dPeriodicx wallp(0, width);
	
	string dir("test");
	mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if (chdir(dir.c_str())) {
		cout << "Invalid folder" << endl;
	}
	string cinfo("cinfo");
	mkdir(cinfo.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
//	FILE * props;
	FILE * positions;
	FILE * rotations;
	FILE * wallPos;
	FILE * stress;
//	props 	 = fopen("props.dat","w");
	positions = fopen("positions_caicos_2.dat","w");
	rotations = fopen("rotations_caicos_2.dat","w");
	wallPos   = fopen("wallPos_caicos_2.dat","w");
	stress    = fopen("stress_caicos_2.dat","w");
	
	vector<Grain2d> grainList;
	double duration;
	double start = omp_get_wtime();
	
	int tiso = 40000;
	int tshear = 1;
	double pressure = 40;
	double v = .137/2.;
	
	// time integration for isotropic compression
//	size_t iter = 0;
	cout <<"Starting iso stage" << endl;
	World2dBiaxial world(grains, walls, wallp, gDamping, dt);
	
	for (int t = 0; t < tiso; t++) {
		world.computeWorldState();
		// non output timestep
		if (t % 1000 != 0) {
		}
		// output timestep
		else {
			duration = omp_get_wtime() - start;
			printf( "iso timestep %d of %d (%4.2f%% complete, %4.2f minutes, ~%4.2f remaining)\n", t+1, tiso, 100*double(t+1)/double(tiso), duration/60., duration/60.*tiso/(t+1) - duration/60.);
			fflush(stdout);
			double volume = world.getWalls()[1].getPosition()(1)*world.getWallPeriodic().getWidth();
			cout << "stress = " << world.getWorldState().stressVoigt.transpose()/volume << endl;
			cout << "volume = " << volume << endl;
		}
		WorldState2d ws = world.getWorldState();
		double inForce = pressure*world.getWallPeriodic().getWidth();
		double volume = world.getWalls()[1].getPosition()(1)*world.getWallPeriodic().getWidth();
		double outForce = -world.getWorldState().stressVoigt(1)/volume*world.getWallPeriodic().getWidth();
		if ( .99*fabs(world.getWorldState().stressVoigt(1)) > fabs(world.getWorldState().stressVoigt(0))) {
			world.getWallPeriodicNC().moveWidth(-0.001);
		}
		if ( 1.01*fabs(world.getWorldState().stressVoigt(1)) < fabs(world.getWorldState().stressVoigt(0))) {
			world.getWallPeriodicNC().moveWidth(0.001);
		}
		world.getWallsNC()[1].takeForceTimestep(inForce-outForce, ws.wallContacts[1]);
		world.takeTimestep();
	}
	
	// shear
	cout <<"Starting shear stage" << endl;
	start = omp_get_wtime();
	size_t station = 0;
	for (int t = 0; t < tshear; t++) {
		world.computeWorldState();
		// non output timestep
		if (t % 1000 != 0) {
		}
		// output timestep
		else {
			CData cData = world.computeCstate();
			string fname = cinfo + string("/") + to_string(station) + string(".dat");
			FILE * cinfofile;
			cinfofile = fopen(fname.c_str(), "w");
			for (size_t j = 0; j < cData.size(); j++) {
				fprintf(cinfofile, "%d %d %.3f %.3f %.4f %.4f %.3f %.3f %lu\n", cData._cpairs[j](0), cData._cpairs[j](1), cData._forces[j](0), cData._forces[j](1),
																					 cData._normals[j](0), cData._normals[j](1), cData._clocs[j](0), cData._clocs[j](1), cData._nodes[j]);
			}
			station++;
			duration = omp_get_wtime() - start;
			printf( "iso timestep %d of %d (%4.2f%% complete, %4.2f minutes, ~%4.2f remaining)\n", t+1, tshear, 100*double(t+1)/double(tshear), duration/60., duration/60.*tshear/(t+1) - duration/60.);
			fflush(stdout);
			double volume = world.getWalls()[1].getPosition()(1)*world.getWallPeriodic().getWidth();
			cout << "stress = " << world.getWorldState().stressVoigt.transpose()/volume << endl;
			cout << "volume = " << volume << endl;
			// print positions and velocities to file
			grainList = world.getGrains();
		//	double ke = 0;
			for (size_t i = 0; i < world.getGrains().size(); i++) {
				fprintf(positions, "%.3f %.3f\n", grainList[i].getPosition()(0), grainList[i].getPosition()(1));
				fprintf(rotations, "%.3f\n", grainList[i].getTheta());
			}
			fprintf(wallPos,"%.3f %.3f\n", world.getWalls()[1].getPosition()(1), world.getWallPeriodic().getWidth() );
		}
		// shear
		world.takeTimestepShear(v);
		// apply pressure
		WorldState2d ws = world.getWorldState();
		double inForce = pressure*world.getWallPeriodic().getWidth();
		double volume = world.getWalls()[1].getPosition()(1)*world.getWallPeriodic().getWidth();
		double outForce = -world.getWorldState().stressVoigt(1)/volume*world.getWallPeriodic().getWidth();
		world.getWallsNC()[1].takeForceTimestep(inForce-outForce, ws.wallContacts[1]);
	}
	
	fclose(positions);
	fclose(rotations);
	fclose(wallPos);
	fclose(stress);
	
	duration = omp_get_wtime() - start;
	printf("Time taken: %.2fs (%.2f minutes)\n", duration, duration/60.);
	cout << endl << "program complete!\n";
	return 0;
}
