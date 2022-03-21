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
 
int main(int argc, char*argv[]) {
//	omp_set_dynamic(0);     // Explicitly disable dynamic threading
//	omp_set_num_threads(6);
	// file name 
	string filename;
	if (argc == 3) {
		filename = string(argv[1]) + ".dat";
	}
	else {
		filename = "disc.dat";
	}
	vector<Grain2d> grains = generateGrainsFromFile(filename);
	vector<Grain2d> temp = grains;
	for (size_t i = 0; i < 779; i++) {
		grains.insert(grains.end(), temp.begin(), temp.end());
	}
	size_t ngrains = grains.size();
	cout <<"Number of grains: " << ngrains << endl;
	// grain properties
	double mu;
	if (argc == 3) {
		mu = atof(argv[2]);
	}
	else {
		mu = 0.3;
	}
	cout << "Type of grain: " << argv[1] << endl;
	cout << "mu: " << mu << endl;
	double kn = 30000;
	double density = 1;
	// assembly properties
	vector<double> wallp = readRotationFile("walls_19500.dat",2);
	double widthInit = wallp[0]; 
	double heightInit = wallp[1];
	
	// world properties
	double gDamping = .7; // .4 or .2
	double dt = 0.2*(2*sqrt(.5*200/kn));
	
	vector<Vector2d> pos = readPositionFile("positions_19500.dat",ngrains);
	vector<double> rot = readRotationFile("rotations_19500.dat",ngrains);
	
	/* set things up */
	for (size_t i = 0; i < ngrains; i++) {
		grains[i].changePos(pos[i]);
		grains[i].changeRot(rot[i]);
		grains[i].changeMu(mu);
		grains[i].changeKn(kn);
		grains[i].changeKs(.9*kn);
		grains[i].changeDensity(density);
		grains[i].changeId(i);
	}
	
	// create walls
	vector<Wall2d> walls; // 0bot 1top
	walls.resize(4);
	// bottom wall
	walls[0] = Wall2d(Vector2d(0,0), Vector2d(0,1), kn, 0.9*kn, 0.0, INT_MAX-3);
	// top wall
	walls[1] = Wall2d(Vector2d(0,heightInit), Vector2d(0,-1), kn, 0.9*kn, 0.0, INT_MAX-2);
	// left wall
	walls[2] = Wall2d(Vector2d(0,0), Vector2d(1,0), kn, 0.9*kn, 0.0, INT_MAX-1);
	// right wall
	walls[3] = Wall2d(Vector2d(widthInit,0), Vector2d(-1,0), kn, 0.9*kn, 0.0, INT_MAX-0);
	
	string dir = string(argv[1]) + "_mu=" + string(argv[2]);
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
	positions = fopen("positions_iso.dat","w");
	rotations = fopen("rotations_iso.dat","w");
	wallPos   = fopen("wallPos_iso.dat","w");
	stress    = fopen("stress_iso.dat","w");
	
	vector<Grain2d> grainList;
	double duration;
	double start = omp_get_wtime();
	
	int tiso = 200001;
	double pressure = 40;
	double dx = -0.1;
	
	// time integration for isotropic compression
//	size_t iter = 0;
	cout <<"Starting iso stage" << endl;
	World2dBiaxial world(grains, walls, gDamping, dt);
	double width, height, volume;
	for (int t = 0; t < tiso; t++) {
		world.computeWorldState();
		width = (world.getWalls()[3].getPosition()(0)-world.getWalls()[2].getPosition()(0));
		height = (world.getWalls()[1].getPosition()(1)-world.getWalls()[0].getPosition()(1));
		volume = width*height;
							 
		// non output timestep
		if (t % 5000 != 0) {
		}
		// output timestep
		else {
			duration = omp_get_wtime() - start;
			printf( "iso timestep %d of %d (%4.2f%% complete, %4.2f minutes, ~%4.2f remaining)\n", t+1, tiso, 100*double(t+1)/double(tiso), duration/60., duration/60.*tiso/(t+1) - duration/60.);
			fflush(stdout);
			cout << "stress = " << world.getWorldState().stressVoigt.transpose()/volume << endl;
			cout << "volume = " << volume << endl;
		}
		if ( (fabs(world.getWorldState().stressVoigt(1)) + fabs(world.getWorldState().stressVoigt(0)))/volume < 60) {
			// move walls
//			world.getWallsNC()[0].moveWall(Vector2d(0, 0.1));
			world.getWallsNC()[1].moveWall(Vector2d(0,dx));
//			world.getWallsNC()[2].moveWall(Vector2d( 0.1,0));
			world.getWallsNC()[3].moveWall(Vector2d(dx,0));
			// move grains accordingly
			world.moveProportional(dx,dx,width,height);
		}
		else {
			WorldState2d ws = world.getWorldState();
			double inForce = pressure*width;
			double outForce = -world.getWorldState().stressVoigt(1)/volume*width;
			world.getWallsNC()[1].takeForceTimestep(inForce-outForce, ws.wallContacts[1]);
			inForce = pressure*height;
			outForce = -world.getWorldState().stressVoigt(0)/volume*height;
			world.getWallsNC()[3].takeForceTimestep(inForce-outForce, ws.wallContacts[3]);
		}
		world.takeTimestep();
	}
	
	
	
	// output at end
	CData cData = world.computeCstate();
	string fname = cinfo + string("/") + to_string(0) + string(".dat");
	FILE * cinfofile;
	cinfofile = fopen(fname.c_str(), "w");
	for (size_t j = 0; j < cData.size(); j++) {
		fprintf(cinfofile, "%d %d %.3f %.3f %.4f %.4f %.3f %.3f %lu\n", cData._cpairs[j](0), cData._cpairs[j](1), cData._forces[j](0), cData._forces[j](1),
																			 cData._normals[j](0), cData._normals[j](1), cData._clocs[j](0), cData._clocs[j](1), cData._nodes[j]);
	}
	fflush(stdout);
	// print positions and velocities to file
	grainList = world.getGrains();
//	double ke = 0;
	for (size_t i = 0; i < world.getGrains().size(); i++) {
		fprintf(positions, "%.3f %.3f\n", grainList[i].getPosition()(0), grainList[i].getPosition()(1));
		fprintf(rotations, "%.3f\n", grainList[i].getTheta());
	}
	fprintf(wallPos,"%.3f %.3f\n", world.getWalls()[1].getPosition()(1), world.getWalls()[3].getPosition()(0) );
	Vector3d stressVoigt = world.getWorldState().stressVoigt;
	fprintf(stress,"%.4f %.4f %.4f", stressVoigt(0)/volume, stressVoigt(1)/volume, stressVoigt(2)/volume);
	
	fclose(positions);
	fclose(rotations);
	fclose(wallPos);
	fclose(stress);
	
	duration = omp_get_wtime() - start;
	printf("Time taken: %.2fs (%.2f minutes)\n", duration, duration/60.);
	cout << endl << "program complete!\n";
	return 0;
}
