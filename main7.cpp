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
	// file name
	string filename("caicosBiaxial.dat");
	vector<Grain2d> grains = generateGrainsFromFile(filename);
	size_t ngrains = grains.size();
	// grain properties
	double mu = 0.0;
	double kn = 30000;
	double density = 1;
	// assembly properties
	double height = 645.3;
	double width = 645.3; 
	// world properties
	double gDamping = .2; // .4 or .2
	double dt = 0.2*(2*sqrt(.5*200/kn));
	int tiso = 5000;
	int tshear = 200001;
	double pressure1 = 500;
	double pressure2 = 40;
	
	Vector2d trans(545.07320305-634.37729014-5.5, 784.81266662-656.22997566+5.5);
	//double topWall = 784.81266662;
	//double rightWall = 545.07320305;
	//double height = 656.22997566; 
	//double width = 634.37729014; 
//	int tshear = 100001;
	
	vector<Vector2d> pos = readPositionFile("positions_caicos.dat",ngrains);
	vector<double> rot = readRotationFile("rotations_caicos.dat",ngrains);
	
	/* set things up */
	for (size_t i = 0; i < grains.size(); i++) {
		grains[i].changePos(pos[i]-trans);
		grains[i].changeRot(rot[i]);
		grains[i].changeMu(mu);
		grains[i].changeKn(kn);
		grains[i].changeDensity(density);
	}
	
	// create walls
	vector<Wall2d> walls; // 0bot 1left 2right
	walls.resize(4);
	// bottom wall
	walls[0] = Wall2d(Vector2d(0,0), Vector2d(0,1), kn, 0.9*kn, 0.5, INT_MAX-3);
	// top wall
	walls[1] = Wall2d(Vector2d(0,height), Vector2d(0,-1), kn, 0.9*kn, 0.5, INT_MAX-2);
	// left wall
	walls[2] = Wall2d(Vector2d(0,0), Vector2d(1,0), kn, 0.9*kn, 0.5, INT_MAX-1);
	// right wall
	walls[3] = Wall2d(Vector2d(width,0), Vector2d(-1,0), kn, 0.9*kn, 0.5, INT_MAX-0);
	
	
	
	FILE * positions;
	FILE * rotations;
	FILE * wallPos;
	FILE * stress;
	positions = fopen ("positions_caicos_3.dat","w");
	rotations = fopen ("rotations_caicos_3.dat","w");
	wallPos   = fopen("wallPos_caicos_3.dat","w");
	stress    = fopen("stress_caicos_3.dat","w");
	
	vector<Grain2d> grainList;
	double duration;
	double start = omp_get_wtime();
//	UtkContactInfo info(world.getGrains().size());
	
	// time integration for isotropic compression
//	size_t iter = 0;
	cout <<"Starting iso stage" << endl;
	World2dBiaxial world(grains, walls, gDamping, dt);
	double pressure;
	for (size_t i = 0; i < 2; i++) {
		if (i == 1 || i == 3) {
			pressure = pressure2;
		}
		else {
			pressure = pressure1;
		}
		for (int t = 0; t < tiso; t++) {
			// non output timestep
			world.computeWorldState();
			if (t % 1000 != 0) {
				
			}
			// output timestep
			else {
				duration = omp_get_wtime() - start;
				printf( "iso timestep %d of %d (%4.2f%% complete, %4.2f minutes, ~%4.2f remaining)\n", t+1, tiso, 100*double(t+1)/double(tiso), duration/60., duration/60.*tiso/(t+1) - duration/60.);
				fflush(stdout);
				double volume = world.getWalls()[1].getPosition()(1)*world.getWalls()[3].getPosition()(0);
				cout << "stress = " << world.getWorldState().stressVoigt.transpose()/volume << endl;
				cout << "volume = " << volume << endl;
			}
			WorldState2d ws = world.getWorldState();
			double inForce = pressure*world.getWalls()[3].getPosition()(0);
			double outForce = ws.wallForces[1](1);
			world.getWallsNC()[1].takeForceTimestep(inForce-outForce, ws.wallContacts[1]);
			inForce = pressure*world.getWalls()[1].getPosition()(1);
			outForce = ws.wallForces[3](0);
			world.getWallsNC()[3].takeForceTimestep(inForce-outForce, ws.wallContacts[3]);
			world.takeTimestep();
		}
	}
	
	mu = 0.5;
	for (size_t i = 0; i < grains.size(); i++) {
		grains[i].changeMu(mu);
	}
	// shear
	for (int t = 0; t < tshear; t++) {
		// non output timestep
		world.computeWorldState();
		if (t % 1000 != 0) {
			
		}
		// output timestep
		else {
			duration = omp_get_wtime() - start;
			printf( "iso timestep %d of %d (%4.2f%% complete, %4.2f minutes, ~%4.2f remaining)\n", t+1, tshear, 100*double(t+1)/double(tshear), duration/60., duration/60.*tshear/(t+1) - duration/60.);
			fflush(stdout);
			double volume = world.getWalls()[1].getPosition()(1)*world.getWalls()[3].getPosition()(0);
			cout << "stress = " << world.getWorldState().stressVoigt.transpose()/volume << endl;
			cout << "volume = " << volume << endl;
		}
		
		// shear
		world.getWallsNC()[2].rotateWall(-30.*M_PI/180./tshear);
		world.getWallsNC()[3].rotateWall(-30.*M_PI/180./tshear);
		// top pressure
		WorldState2d ws = world.getWorldState();
		double inForce = pressure*world.getWalls()[3].getPosition()(0);
		double outForce = ws.wallForces[1](1);
		world.getWallsNC()[1].takeForceTimestep(inForce-outForce, ws.wallContacts[1]);
//		inForce = pressure*world.getWalls()[1].getPosition()(1);
//		outForce = ws.wallForces[3](0);
//		world.getWallsNC()[3].takeForceTimestep(inForce-outForce, ws.wallContacts[3]);
		world.takeTimestep();
	}
	

	// print positions and velocities to file
	grainList = world.getGrains();
//	double ke = 0;
	for (size_t i = 0; i < world.getGrains().size(); i++) {
		fprintf(positions, "%.4f %.4f\n", grainList[i].getPosition()(0), grainList[i].getPosition()(1));
		fprintf(rotations, "%.4f\n", grainList[i].getTheta());
	}
	fprintf(wallPos,"%.4f %.4f\n", world.getWalls()[3].getPosition()(0), world.getWalls()[1].getPosition()(1));
	
	fclose(positions);
	fclose(rotations);
	fclose(wallPos);
	fclose(stress);
	
	duration = omp_get_wtime() - start;
	printf("Time taken: %.2fs (%.2f minutes)\n", duration, duration/60.);
	cout << endl << "program complete!";
	return 0;
}
