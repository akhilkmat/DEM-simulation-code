#include <iostream>
#include <stdlib.h> 
#include <vector>
#include <cmath>
#include<algorithm>
#include <fstream>
#include <chrono>
#include <omp.h>

#include "Functions.h"


using namespace std;
using namespace std::chrono;


int main()
{
	float  binl = 0.01, ins_reg_low_x = 0, ins_reg_up_x = 0.05, ins_reg_low_y = 0, ins_reg_up_y = 0.05, ins_reg_low_z = 0, ins_reg_up_z = 1, dom_low_x = 0, dom_up_x = 0.05, dom_low_y = 0, dom_up_y = 0.05, dom_low_z = 0, dom_up_z = 1, rad_p1 = 0.005;
	float dlx = dom_low_x, dux = dom_up_x, dly = dom_low_y, duy = dom_up_y, dlz = dom_low_z, duz = dom_up_z;
	int no_particles, ghost_particles = 0, trialcounterlimit = 20, guessbincount = 100, no_timesteps = 5000000, dumpevery=1000,tcounter = 1, exitparticles = 0, ins_choice = 1, boundx = 3, boundy = 3, boundz = 1;
	double  Kn, Cn, Kt, Ct, G, pi = 2 * acos(0), timestep = pow(10, -6);
	double Ym = 5 * pow(10, 8), Y, e = 0.6, poi = 0.2, m, R, density = 2500, gx =0 , gy = 0.0, gz =-1 , g = 9.81;

	Functions F;

	R = rad_p1 / 2;
	m = (4.0 / 3.0) * pi * pow(rad_p1, 3) * density;
	Y = pow((1 - pow(poi, 2)) * 2 / Ym, -1);
	Kn = (4.0 / 3.0) * Y * sqrt(R);
	Cn = 2 * sqrt(5.0 / 6.0) * (log(e) / sqrt(log(e) * log(e) + pi * pi)) * sqrt(2 * Y * sqrt(R) * m / 2);
	G = Ym / (4 * (2 - poi) * (1 + poi));
	Kt = 8 * G * sqrt(R);
	Ct = 2 * sqrt(5.0 / 6.0) * (log(e) / sqrt(log(e) * log(e) + pi * pi)) * sqrt(8 * G * sqrt(R) * m / 2);
	std::cout << m << "  " << Kn << "  " << Cn << "  " << endl;
	std::cout << " enter no of particles to be inserted" << endl;

	cin >> no_particles;

	auto start = high_resolution_clock::now(); // to measure the time span of simulation

	ofstream oFile;
	if (boundx == 3)  // checking for periodic conditions
	{
		dom_low_x = dom_low_x - binl;
		dom_up_x = dom_up_x + binl;
	}
	if (boundy == 3)
	{
		dom_low_y = dom_low_y - binl;
		dom_up_y = dom_up_y + binl;
	}
	if (boundz == 3)
	{
		dom_low_z = dom_low_z - binl;
		dom_up_z = dom_up_z + binl;
	}



	oFile.open(("dumpih.atom")); // opoening the dump file
	if (!oFile)
	{
		cerr << "Unable to open file dumpih.atom";
		exit(1);
	}
	oFile.precision(7);



	vector<vector<double>> P(4 * no_particles, vector<double>(23, 0)); // main particle array/vector (P) intialization.
	// P: flag type xpos ypos zpos vx vy vz fx fy fz wx wy wz Tx Ty Tz cxx cyy czz cxy cxz cyz




	F.insert_particles(ins_choice, ins_reg_low_x, ins_reg_up_x, ins_reg_low_y, ins_reg_up_y, ins_reg_low_z, ins_reg_up_z, rad_p1, no_particles, trialcounterlimit, P);

	cout << " check";
	int binx = ceil((dom_up_x - dom_low_x) / binl); // calculations for total no of bins in x y z directions
	int biny = ceil((dom_up_y - dom_low_y) / binl);
	int binz = ceil((dom_up_z - dom_low_z) / binl);
	int  bincount;
	double cx, cy, cz, Vx, Vy, Vz, Vtx, Vty, Vtz, Vnx, Vny, Vnz, Vtwx, Vtwy, Vtwz, Vdot_t, sx, sy, sz, s, Vttx, Vtty, Vttz, Vtt;
	double d, cd, Fn, nx, ny, nz, Vrn, fx, fy, fz, ft, Ft, mu = 0.2, temp, temp1, temp2, Fx, Fy, Fz;

	bincount = pow(binl, 3) / ((4.0 / 3.0) * pi * pow(rad_p1, 3)); // how many particles can be contained in a cube of size binl


	bincount = min(int(15 * bincount), guessbincount);
	cout << "bincout: " << bincount << endl;

	vector<vector<vector<vector<float>>>> cbin(binx, vector<vector<vector<float>>> // contact bin array
		(biny, vector<vector<float>>
			(binz, vector<float>
				(bincount, 0))));

	vector<vector<vector<double>>> PP(no_particles, vector<vector<double>> // particle particle array ,to store the history of tangentiual contact
		(4 * no_particles, vector<double>
			(4, 0)));


	while (tcounter <= no_timesteps)
	{
		ghost_particles = 0;
		if (tcounter > 1000000)
		{
			gx = sin(0.3);
			gz = -cos(0.3);
		}

        #pragma omp parallel for  private(temp, temp1, temp2,cx, cy, cz, Vx, Vy, Vz, Vtx, Vty, Vtz, Vnx, Vny, Vnz, Vtwx, Vtwy, Vtwz, Vdot_t, sx, sy, sz, s, Vttx, Vtty, Vttz, Vtt,d, cd, Fn, nx, ny, nz, Vrn, fx, fy, fz, ft, Ft)
		for (int i = 0; i < no_particles; i++)
		{

			if (P[i][0] == -1)
			{
				continue;
			}
			P[i][0] = 0;

			P[i][2] = P[i][2] + P[i][5] * timestep + 0.5 * (P[i][8] / m) * timestep * timestep;  // position updation velocity verlet scheme.
			P[i][3] = P[i][3] + P[i][6] * timestep + 0.5 * (P[i][9] / m) * timestep * timestep;
			P[i][4] = P[i][4] + P[i][7] * timestep + 0.5 * (P[i][10] / m) * timestep * timestep;

			P[i][5] = P[i][5] + 0.5 * timestep * P[i][8] / m; // half velocity updation VV scheme.
			P[i][6] = P[i][6] + 0.5 * timestep * P[i][9] / m;
			P[i][7] = P[i][7] + 0.5 * timestep * P[i][10] / m;

			P[i][11] = P[i][11] + 0.5 * timestep * P[i][14] / ((2 * m * rad_p1 * rad_p1) / 5);  //  half omega updation, VV scheme.
			P[i][12] = P[i][12] + 0.5 * timestep * P[i][15] / ((2 * m * rad_p1 * rad_p1) / 5);
			P[i][13] = P[i][13] + 0.5 * timestep * P[i][16] / ((2 * m * rad_p1 * rad_p1) / 5);

			switch (boundx)
			{
			case(2):
				if (P[i][2] >= dux)
				{
					P[i][0] = -1;
					exitparticles++;
				}
				else  if (P[i][2] < dlx)
				{
					P[i][0] = -1;
					exitparticles++;
				}break;
			case(3):
				if (P[i][2] >= dux)
				{
					P[i][2] = (P[i][2] - dux) + dlx;
				}
				else  if (P[i][2] < dlx)
				{
					P[i][2] = (P[i][2] - dlx) + dux;
				}break;
			default:
				break;

			}


			switch (boundy)
			{
			case(2):
				if (P[i][3] >= duy)
				{
					P[i][0] = -1;
					exitparticles++;
				}
				else  if (P[i][3] < dly)
				{
					P[i][0] = -1;
					exitparticles++;
				}break;
			case(3):
				if (P[i][3] >= duy)
				{
					P[i][3] = (P[i][3] - duy) + dly;
				}
				else  if (P[i][3] < dly)
				{
					P[i][3] = (P[i][3] - dly) + duy;
				}break;
			default:
				break;

			}


			switch (boundz)
			{
			case(2):
				if (P[i][4] >= duz)
				{
					P[i][0] = -1;
					exitparticles++;
				}
				else  if (P[i][4] < dlz)
				{
					P[i][0] = -1;
					exitparticles++;
				}break;
			case(3):
				if (P[i][4] >= duz)
				{
					P[i][4] = (P[i][4] - duz) + dlz;
				}
				else  if (P[i][4] < dlz)
				{
					P[i][4] = (P[i][4] - dlz) + duz;
				}break;
			default:
				break;

			}


			P[i][8] = m * gx * g;  // intializing forces with gravity in the current timestep
			P[i][9] = m * gy * g;
			P[i][10] = m * gz * g;

			P[i][14] = 0;  // intializing torque = 0
			P[i][15] = 0;
			P[i][16] = 0;



			int xb = floor((P[i][2] - dom_low_x) / binl);
			int yb = floor((P[i][3] - dom_low_y) / binl);
			int zb = floor((P[i][4] - dom_low_z) / binl);


			switch (boundx)
			{
			case(1):
			{

				if (xb == 0) // near to the lower x domain
				{
					d = rad_p1 - (P[i][2] - dom_low_x);
					if (d > 0)
					{
						Vx = P[i][5];
						Fn = Kn * sqrt(2) * pow(d, 1.5) + Cn * sqrt(2 * sqrt(2)) * pow(d, 0.25) * Vx;
						P[i][8] += Fn;

						nx = (P[i][2] - dom_low_x) / abs(P[i][2] - dom_low_x);
						ny = 0;
						nz = 0;

						Vtx = 0; // subtracting normal componentrs to get tangential component of center of mass velocity.
						Vty = P[i][6];
						Vtz = P[i][7];

						Vtwx = 0; // Vtwx= Vtw2x-Vt1wx ; micro sliding due to relative rotation after taking v = w x(cross) r of each particle  to get relative tangential velocity.
						Vtwy = -rad_p1 * nx * P[i][13];
						Vtwz = rad_p1 * nx * P[i][12];


						Vttx = Vtx + Vtwx; // total micro sliding components in the tangential directions
						Vtty = Vty + Vtwy;
						Vttz = Vtz + Vtwz;

						Vtt = sqrt(pow(Vttx, 2) + pow(Vtty, 2) + pow(Vttz, 2)); //magnitude of total tangential velocity



						sx = Vttx * timestep; // tangential overlap 
						sy = Vtty * timestep;
						sz = Vttz * timestep;

						if (PP[i][i + no_particles][0] == tcounter - 1) // if wall is present, periodic particle allocation ids can be used for other purposes.
						{											// here, paticle to wall id is taken as particle to the same periodic particle in that axis.
							PP[i][i + no_particles][1] = PP[i][i + no_particles][1] + sx;
							PP[i][i + no_particles][2] = PP[i][i + no_particles][2] + sy;
							PP[i][i + no_particles][3] = PP[i][i + no_particles][3] + sz;
							PP[i][i + no_particles][0] = tcounter;
						}
						else
						{
							PP[i][i + no_particles][1] = sx;
							PP[i][i + no_particles][2] = sy;
							PP[i][i + no_particles][3] = sz;
							PP[i][i + no_particles][0] = tcounter;
						}

						s = sqrt(pow(PP[i][i + no_particles][1], 2) + pow(PP[i][i + no_particles][2], 2) + pow(PP[i][i + no_particles][3], 2)); // magnitude of total tangential overlap
						if (s != 0)
						{

							ft = -(Kt * pow(d, 0.5) * s + Ct * pow(d, 0.25) * Vtt); //tangential force

							if ((abs(ft) >= mu * abs(Fn)) && Vtt > 0)
							{
								Ft = -mu * abs(Fn);

								PP[i][i + no_particles][1] = 0;  // slip
								PP[i][i + no_particles][2] = 0;
								PP[i][i + no_particles][3] = 0;


								P[i][8] += Ft * Vttx / Vtt;    // tangential force components // cannot usse sx etc.. since its overlap for current timestep
								P[i][9] += Ft * Vtty / Vtt;
								P[i][10] += Ft * Vttz / Vtt;


								temp = -rad_p1 * (ny * Ft * Vttz / Vtt - nz * Ft * Vtty / Vtt);   // torque 
								temp1 = rad_p1 * (nx * Ft * Vttz / Vtt - nz * Ft * Vttx / Vtt);
								temp2 = -rad_p1 * (nx * Ft * Vtty / Vtt - ny * Ft * Vttx / Vtt);

								P[i][14] += temp;   // torque updation
								P[i][15] += temp1;
								P[i][16] += temp2;


							}
							else
							{
								Ft = -min(abs(ft), mu * abs(Fn));   //  tangnetial force with limit

								fx = -(Kt * pow(d, 0.5) * PP[i][i + no_particles][1] + Ct * pow(d, 0.25) * Vttx);
								fy = -(Kt * pow(d, 0.5) * PP[i][i + no_particles][2] + Ct * pow(d, 0.25) * Vtty);
								fz = -(Kt * pow(d, 0.5) * PP[i][i + no_particles][3] + Ct * pow(d, 0.25) * Vttz);


								P[i][8] += Ft * fx / ft;    // tangential force components // cannot usse sx etc.. since its overlap for current timestep
								P[i][9] += Ft * fy / ft;
								P[i][10] += Ft * fz / ft;


								temp = -rad_p1 * (ny * Ft * fz / ft - nz * Ft * fy / ft);   // torque updation
								temp1 = rad_p1 * (nx * Ft * fz / ft - nz * Ft * fx / ft);
								temp2 = -rad_p1 * (nx * Ft * fy / ft - ny * Ft * fx / ft);

								P[i][14] += temp;
								P[i][15] += temp1;
								P[i][16] += temp2;
							}

						
						}


					}
				}
				else if (xb == binx - 1) // near to the upper x domain
				{
					d = rad_p1 - (dom_up_x - P[i][2]);
					if (d > 0)
					{
						Vx = P[i][5];
						Fn = -Kn * sqrt(2) * pow(d, 1.5) + Cn * sqrt(2 * sqrt(2)) * pow(d, 0.25) * Vx;
						P[i][8] += Fn;


						nx = (P[i][2] - dom_up_x) / abs(P[i][2] - dom_up_x);


						Vtx = 0; // subtracting normal componentrs to get tangential component of center of mass velocity.
						Vty = P[i][6];
						Vtz = P[i][7];

						Vtwx = 0; // Vtwx= Vtw2x-Vt1wx ; micro sliding due to relative rotation after taking v = w x(cross) r of each particle  to get relative tangential velocity.
						Vtwy = -rad_p1 * nx * P[i][13];
						Vtwz = rad_p1 * nx * P[i][12];


						Vttx = Vtx + Vtwx; // total micro sliding components in the tangential directions
						Vtty = Vty + Vtwy;
						Vttz = Vtz + Vtwz;

						Vtt = sqrt(pow(Vttx, 2) + pow(Vtty, 2) + pow(Vttz, 2)); //magnitude of total tangential velocity

						sx = Vttx * timestep; // tangential overlap 
						sy = Vtty * timestep;
						sz = Vttz * timestep;

						if (PP[i][i + no_particles][0] == tcounter - 1) // if wall is present, periopdic particle allocation ids can be used for other purposes.
						{											// here, paticle to wall id is taken as particle to the same periodic particle in that axis.
							PP[i][i + no_particles][1] = PP[i][i + no_particles][1] + sx;
							PP[i][i + no_particles][2] = PP[i][i + no_particles][2] + sy;
							PP[i][i + no_particles][3] = PP[i][i + no_particles][3] + sz;
							PP[i][i + no_particles][0] = tcounter;
						}
						else
						{
							PP[i][i + no_particles][1] = sx;
							PP[i][i + no_particles][2] = sy;
							PP[i][i + no_particles][3] = sz;
							PP[i][i + no_particles][0] = tcounter;
						}

						s = sqrt(pow(PP[i][i + no_particles][1], 2) + pow(PP[i][i + no_particles][2], 2) + pow(PP[i][i + no_particles][3], 2)); // magnitude of total tangential overlap

						if (s != 0)
						{
							ft = -(Kt * pow(d, 0.5) * s + Ct * pow(d, 0.25) * Vtt); //tangential force

							if ((abs(ft) >= mu * abs(Fn)) && Vtt > 0)
							{
								Ft = -mu * abs(Fn);

								PP[i][i + no_particles][1] = 0; // slip
								PP[i][i + no_particles][2] = 0;
								PP[i][i + no_particles][3] = 0;


								P[i][8] += Ft * Vttx / Vtt;    // tangential force components // cannot usse sx etc.. since its overlap for current timestep
								P[i][9] += Ft * Vtty / Vtt;
								P[i][10] += Ft * Vttz / Vtt;


								temp = -rad_p1 * (ny * Ft * Vttz / Vtt - nz * Ft * Vtty / Vtt);   // torque 
								temp1 = rad_p1 * (nx * Ft * Vttz / Vtt - nz * Ft * Vttx / Vtt);
								temp2 = -rad_p1 * (nx * Ft * Vtty / Vtt - ny * Ft * Vttx / Vtt);

								P[i][14] += temp;   // torque updation
								P[i][15] += temp1;
								P[i][16] += temp2;


							}
							else
							{
								Ft = -min(abs(ft), mu * abs(Fn));   //  tangnetial force with limit

								fx = -(Kt * pow(d, 0.5) * PP[i][i + no_particles][1] + Ct * pow(d, 0.25) * Vttx);
								fy = -(Kt * pow(d, 0.5) * PP[i][i + no_particles][2] + Ct * pow(d, 0.25) * Vtty);
								fz = -(Kt * pow(d, 0.5) * PP[i][i + no_particles][3] + Ct * pow(d, 0.25) * Vttz);


								P[i][8] += Ft * fx / ft;    // tangential force components // cannot usse sx etc.. since its overlap for current timestep
								P[i][9] += Ft * fy / ft;
								P[i][10] += Ft * fz / ft;


								temp = -rad_p1 * (ny * Ft * fz / ft - nz * Ft * fy / ft);   // torque updation
								temp1 = rad_p1 * (nx * Ft * fz / ft - nz * Ft * fx / ft);
								temp2 = -rad_p1 * (nx * Ft * fy / ft - ny * Ft * fx / ft);

								P[i][14] += temp;
								P[i][15] += temp1;
								P[i][16] += temp2;
							}
							
						}
					}
				}

			}break;

			case(3):
			{
				if (xb == 1) // xb = 0 is ghost cell
				{
					P[no_particles + i][0] = -1; // ghost particles only impart forces and we dont calculate forces on them so this is kept as -1. the impart forces will be calculated as neighbour cell loop is executed.
					P[no_particles + i][1] = 2; // just a type asignment for ghost particles
					P[no_particles + i][2] = P[i][2] + (dux - dlx);
					P[no_particles + i][3] = P[i][3];
					P[no_particles + i][4] = P[i][4];
					P[no_particles + i][5] = P[i][5];
					P[no_particles + i][6] = P[i][6];
					P[no_particles + i][7] = P[i][7];
					P[no_particles + i][11] = P[i][11];
					P[no_particles + i][12] = P[i][12];
					P[no_particles + i][13] = P[i][13];
					
					

						if (tcounter != cbin[binx - 1][yb][zb][1]) // to implement conditions when a bin is populated 1st time within a particular timestep iteration.
						{
							#pragma omp critical
							{
								cbin[binx - 1][yb][zb][0] = 0;
								cbin[binx - 1][yb][zb][1] = tcounter;
							}
						}
						
						#pragma omp critical
						{
							cbin[binx - 1][yb][zb][0]++;
							cbin[binx - 1][yb][zb][cbin[binx - 1][yb][zb][0] + 1] = no_particles + i;
							ghost_particles++;
						}
					    
				}
				else if (xb == binx - 2)
				{
					P[no_particles + i][0] = -1;
					P[no_particles + i][1] = 2;
					P[no_particles + i][2] = P[i][2] - (dux - dlx);
					P[no_particles + i][3] = P[i][3];
					P[no_particles + i][4] = P[i][4];
					P[no_particles + i][5] = P[i][5];
					P[no_particles + i][6] = P[i][6];
					P[no_particles + i][7] = P[i][7];
					P[no_particles + i][11] = P[i][11];
					P[no_particles + i][12] = P[i][12];
					P[no_particles + i][13] = P[i][13];
					
					
						if (tcounter != cbin[0][yb][zb][1]) // to implement conditions when a bin is populated 1st time within a particular timestep iteration.
						{
							#pragma omp critical
							{	
								cbin[0][yb][zb][0] = 0;
								cbin[0][yb][zb][1] = tcounter;
							}
						}
					
						#pragma omp critical
						{
							cbin[0][yb][zb][0]++;
							cbin[0][yb][zb][cbin[0][yb][zb][0] + 1] = no_particles + i;
							ghost_particles++;
						}
				        
				}

			}break;
			default:
				break;
			}


			switch (boundy)
			{
			case(1):
			{

				if (yb == 0)
				{
					d = rad_p1 - (P[i][3] - dom_low_y);
					if (d > 0)
					{
						Vy = P[i][6];
						Fn = Kn * sqrt(2) * pow(d, 1.5) + Cn * sqrt(2 * sqrt(2)) * pow(d, 0.25) * Vy;
						P[i][9] += Fn;

						ny = (P[i][3] - dom_low_y) / abs(P[i][3] - dom_low_y);
						nx = 0;
						nz = 0;


						Vtx = P[i][5]; // subtracting normal componentrs to get tangential component of center of mass velocity.
						Vty = 0;
						Vtz = P[i][7];

						Vtwx = rad_p1 * ny * P[i][13]; // Vtx= Vt1x-Vt2x ; micro sliding due to rotation after taking v = w x(cross) r of each particle and subtracting to get relative tangential velocity.
						Vtwy = 0;
						Vtwz = -rad_p1 * P[i][11] * ny;

						Vttx = Vtx + Vtwx; // total micro sliding components in the tangential directions
						Vtty = Vty + Vtwy;
						Vttz = Vtz + Vtwz;

						Vtt = sqrt(pow(Vttx, 2) + pow(Vtty, 2) + pow(Vttz, 2)); //magnitude of total tangential velocity

						sx = Vttx * timestep; // tangential overlap 
						sy = Vtty * timestep;
						sz = Vttz * timestep;

						if (PP[i][i + 2 * no_particles][0] == tcounter - 1) // if wall is present, periopdic particle allocation ids can be used for other purposes.
						{											// here, paticle to wall id is taken as particle to the same periodic particle in that axis.
							PP[i][i + 2 * no_particles][1] = PP[i][i + 2 * no_particles][1] + sx;
							PP[i][i + 2 * no_particles][2] = PP[i][i + 2 * no_particles][2] + sy;
							PP[i][i + 2 * no_particles][3] = PP[i][i + 2 * no_particles][3] + sz;
							PP[i][i + 2 * no_particles][0] = tcounter;
						}
						else
						{
							PP[i][i + 2 * no_particles][1] = sx;
							PP[i][i + 2 * no_particles][2] = sy;
							PP[i][i + 2 * no_particles][3] = sz;
							PP[i][i + 2 * no_particles][0] = tcounter;
						}

						s = sqrt(pow(PP[i][i + 2 * no_particles][1], 2) + pow(PP[i][i + 2 * no_particles][2], 2) + pow(PP[i][i + 2 * no_particles][3], 2)); // magnitude of total tangential overlap

						if (s != 0)
						{

							ft = -(Kt * pow(d, 0.5) * s + Ct * pow(d, 0.25) * Vtt); //tangential force

							if ((abs(ft) >= mu * abs(Fn)) && Vtt > 0)
							{
								Ft = -mu * abs(Fn);

								PP[i][i + 2 * no_particles][1] = 0; // slip
								PP[i][i + 2 * no_particles][2] = 0;
								PP[i][i + 2 * no_particles][3] = 0;

								P[i][8] += Ft * Vttx / Vtt;    // tangential force components // cannot usse sx etc.. since its overlap for current timestep
								P[i][9] += Ft * Vtty / Vtt;
								P[i][10] += Ft * Vttz / Vtt;


								temp = -rad_p1 * (ny * Ft * Vttz / Vtt - nz * Ft * Vtty / Vtt);   // torque 
								temp1 = rad_p1 * (nx * Ft * Vttz / Vtt - nz * Ft * Vttx / Vtt);
								temp2 = -rad_p1 * (nx * Ft * Vtty / Vtt - ny * Ft * Vttx / Vtt);

								P[i][14] += temp;   // torque updation
								P[i][15] += temp1;
								P[i][16] += temp2;


							}
							else
							{
								Ft = -min(abs(ft), mu * abs(Fn));   //  tangnetial force with limit

								fx = -(Kt * pow(d, 0.5) * PP[i][i + 2 * no_particles][1] + Ct * pow(d, 0.25) * Vttx);
								fy = -(Kt * pow(d, 0.5) * PP[i][i + 2 * no_particles][2] + Ct * pow(d, 0.25) * Vtty);
								fz = -(Kt * pow(d, 0.5) * PP[i][i + 2 * no_particles][3] + Ct * pow(d, 0.25) * Vttz);


								P[i][8] += Ft * fx / ft;    // tangential force components // cannot usse sx etc.. since its overlap for current timestep
								P[i][9] += Ft * fy / ft;
								P[i][10] += Ft * fz / ft;


								temp = -rad_p1 * (ny * Ft * fz / ft - nz * Ft * fy / ft);   // torque updation
								temp1 = rad_p1 * (nx * Ft * fz / ft - nz * Ft * fx / ft);
								temp2 = -rad_p1 * (nx * Ft * fy / ft - ny * Ft * fx / ft);

								P[i][14] += temp;
								P[i][15] += temp1;
								P[i][16] += temp2;
							}

							
						}

					}
				}
				else if (yb == biny - 1)
				{
					d = rad_p1 - (dom_up_y - P[i][3]);
					if (d > 0)
					{
						Vy = P[i][6];
						Fn = -Kn * sqrt(2) * pow(d, 1.5) + Cn * sqrt(2 * sqrt(2)) * pow(d, 0.25) * Vy;
						P[i][9] += Fn;

						ny = (P[i][3] - dom_up_y) / abs(P[i][3] - dom_up_y);


						Vtx = P[i][5]; // subtracting normal componentrs to get tangential component of center of mass velocity.
						Vty = 0;
						Vtz = P[i][7];

						Vtwx = rad_p1 * ny * P[i][13]; // Vtx= Vt1x-Vt2x ; micro sliding due to rotation after taking v = w x(cross) r of each particle and subtracting to get relative tangential velocity.
						Vtwy = 0;
						Vtwz = -rad_p1 * P[i][11] * ny;



						Vttx = Vtx + Vtwx; // total micro sliding components in the tangential directions
						Vtty = Vty + Vtwy;
						Vttz = Vtz + Vtwz;

						Vtt = sqrt(pow(Vttx, 2) + pow(Vtty, 2) + pow(Vttz, 2)); //magnitude of total tangential velocity

						sx = Vttx * timestep; // tangential overlap 
						sy = Vtty * timestep;
						sz = Vttz * timestep;

						if (PP[i][i + 2 * no_particles][0] == tcounter - 1) // if wall is present, periopdic particle allocation ids can be used for other purposes.
						{											// here, paticle to wall id is taken as particle to the same periodic particle in that axis.
							PP[i][i + 2 * no_particles][1] = PP[i][i + 2 * no_particles][1] + sx;
							PP[i][i + 2 * no_particles][2] = PP[i][i + 2 * no_particles][2] + sy;
							PP[i][i + 2 * no_particles][3] = PP[i][i + 2 * no_particles][3] + sz;
							PP[i][i + 2 * no_particles][0] = tcounter;
						}
						else
						{
							PP[i][i + 2 * no_particles][1] = sx;
							PP[i][i + 2 * no_particles][2] = sy;
							PP[i][i + 2 * no_particles][3] = sz;
							PP[i][i + 2 * no_particles][0] = tcounter;
						}

						s = sqrt(pow(PP[i][i + 2 * no_particles][1], 2) + pow(PP[i][i + 2 * no_particles][2], 2) + pow(PP[i][i + 2 * no_particles][3], 2)); // magnitude of total tangential overlap

						if (s != 0)
						{
							ft = -(Kt * pow(d, 0.5) * s + Ct * pow(d, 0.25) * Vtt); //tangential force

							if ((abs(ft) >= mu * abs(Fn)) && Vtt > 0)
							{
								Ft = -mu * abs(Fn);

								PP[i][i + 2 * no_particles][1] = 0; //slip
								PP[i][i + 2 * no_particles][2] = 0;
								PP[i][i + 2 * no_particles][3] = 0;


								P[i][8] += Ft * Vttx / Vtt;    // tangential force components // cannot usse sx etc.. since its overlap for current timestep
								P[i][9] += Ft * Vtty / Vtt;
								P[i][10] += Ft * Vttz / Vtt;


								temp = -rad_p1 * (ny * Ft * Vttz / Vtt - nz * Ft * Vtty / Vtt);   // torque 
								temp1 = rad_p1 * (nx * Ft * Vttz / Vtt - nz * Ft * Vttx / Vtt);
								temp2 = -rad_p1 * (nx * Ft * Vtty / Vtt - ny * Ft * Vttx / Vtt);

								P[i][14] += temp;   // torque updation
								P[i][15] += temp1;
								P[i][16] += temp2;


							}
							else
							{
								Ft = -min(abs(ft), mu * abs(Fn));   //  tangnetial force with limit

								fx = -(Kt * pow(d, 0.5) * PP[i][i + 2 * no_particles][1] + Ct * pow(d, 0.25) * Vttx);
								fy = -(Kt * pow(d, 0.5) * PP[i][i + 2 * no_particles][2] + Ct * pow(d, 0.25) * Vtty);
								fz = -(Kt * pow(d, 0.5) * PP[i][i + 2 * no_particles][3] + Ct * pow(d, 0.25) * Vttz);


								P[i][8] += Ft * fx / ft;    // tangential force components // cannot usse sx etc.. since its overlap for current timestep
								P[i][9] += Ft * fy / ft;
								P[i][10] += Ft * fz / ft;


								temp = -rad_p1 * (ny * Ft * fz / ft - nz * Ft * fy / ft);   // torque updation
								temp1 = rad_p1 * (nx * Ft * fz / ft - nz * Ft * fx / ft);
								temp2 = -rad_p1 * (nx * Ft * fy / ft - ny * Ft * fx / ft);

								P[i][14] += temp;
								P[i][15] += temp1;
								P[i][16] += temp2;
							}
							
						}


					}
				}

			}break;
			case(3):
			{
				if (yb == 1)
				{
					P[2 * no_particles + i][0] = -1;
					P[2 * no_particles + i][1] = 2;
					P[2 * no_particles + i][2] = P[i][2];
					P[2 * no_particles + i][3] = P[i][3] + (duy - dly);
					P[2 * no_particles + i][4] = P[i][4];
					P[2 * no_particles + i][5] = P[i][5];
					P[2 * no_particles + i][6] = P[i][6];
					P[2 * no_particles + i][7] = P[i][7];
					P[2 * no_particles + i][11] = P[i][11];
					P[2 * no_particles + i][12] = P[i][12];
					P[2 * no_particles + i][13] = P[i][13];
					
					
						if (tcounter != cbin[xb][biny - 1][zb][1]) // to implement conditions when a bin is populated 1st time within a particular timestep iteration.
						{	
							#pragma omp critical
							{	
								cbin[xb][biny - 1][zb][0] = 0;
								cbin[xb][biny - 1][zb][1] = tcounter;
							}
						}
						
						#pragma omp critical
						{
							cbin[xb][biny - 1][zb][0]++;
							cbin[xb][biny - 1][zb][cbin[xb][biny - 1][zb][0] + 1] = 2 * no_particles + i;
							ghost_particles++;
						}
					    
				}
				else if (yb == biny - 2)
				{
					P[2 * no_particles + i][0] = -1;
					P[2 * no_particles + i][1] = 2;
					P[2 * no_particles + i][2] = P[i][2];
					P[2 * no_particles + i][3] = P[i][3] - (duy - dly);
					P[2 * no_particles + i][4] = P[i][4];
					P[2 * no_particles + i][5] = P[i][5];
					P[2 * no_particles + i][6] = P[i][6];
					P[2 * no_particles + i][7] = P[i][7];
					P[2 * no_particles + i][11] = P[i][11];
					P[2 * no_particles + i][12] = P[i][12];
					P[2 * no_particles + i][13] = P[i][13];
					
					
						if (tcounter != cbin[xb][0][zb][1]) // to implement conditions when a bin is populated 1st time within a particular timestep iteration.
						{	
							#pragma omp critical
							{
								cbin[xb][0][zb][0] = 0;
								cbin[xb][0][zb][1] = tcounter;
							}
						}
						
						#pragma omp critical
						{
							cbin[xb][0][zb][0]++;
							cbin[xb][0][zb][cbin[xb][0][zb][0] + 1] = 2 * no_particles + i;
							ghost_particles++;
						}
				}

			}break;

			default:
				break;
			}

			switch (boundz)
			{
			case(1):
			{

				if (zb == 0)
				{
					d = rad_p1 - (P[i][4] - dom_low_z);
					if (d > 0)
					{
						Vz = P[i][7];
						Fn = Kn * sqrt(2) * pow(d, 1.5) + Cn * sqrt(2 * sqrt(2)) * pow(d, 0.25) * Vz;
						P[i][10] += Fn;

						nz = (P[i][4] - dom_low_z) / abs(P[i][4] - dom_low_z);
						nx = 0;
						ny = 0;


						P[i][5] = 0; //no velocity in x  direction at z ground surface bumpy inclined plane
						P[i][6] = 0;

						Vtx = P[i][5]; // subtracting normal componentrs to get tangential component of center of mass velocity.
						Vty = P[i][6];
						Vtz = 0;

						Vtwx = -rad_p1 * P[i][12] * nz; // Vtx= Vt1x-Vt2x ; micro sliding due to rotation after taking v = w x(cross) r of each particle and subtracting to get relative tangential velocity.
						Vtwy = rad_p1 * P[i][11] * nz;
						Vtwz = 0;

						

						Vttx = Vtx + Vtwx; // total micro sliding components in the tangential directions
						Vtty = Vty + Vtwy;
						Vttz = Vtz + Vtwz;

						Vtt = sqrt(pow(Vttx, 2) + pow(Vtty, 2) + pow(Vttz, 2)); //magnitude of total tangential velocity

						sx = Vttx * timestep; // tangential overlap 
						sy = Vtty * timestep;
						sz = Vttz * timestep;

						if (PP[i][i + 3 * no_particles][0] == tcounter - 1) // if wall is present, periopdic particle allocation ids can be used for other purposes.
						{											// here, paticle to wall id is taken as particle to the same periodic particle in that axis.
							PP[i][i + 3 * no_particles][1] = PP[i][i + 3 * no_particles][1] + sx;
							PP[i][i + 3 * no_particles][2] = PP[i][i + 3 * no_particles][2] + sy;
							PP[i][i + 3 * no_particles][3] = PP[i][i + 3 * no_particles][3] + sz;
							PP[i][i + 3 * no_particles][0] = tcounter;
						}
						else
						{
							PP[i][i + 3 * no_particles][1] = sx;
							PP[i][i + 3 * no_particles][2] = sy;
							PP[i][i + 3 * no_particles][3] = sz;
							PP[i][i + 3 * no_particles][0] = tcounter;
						}

						s = sqrt(pow(PP[i][i + 3 * no_particles][1], 2) + pow(PP[i][i + 3 * no_particles][2], 2) + pow(PP[i][i + 3 * no_particles][3], 2)); // magnitude of total tangential overlap

						if (s != 0)
						{
							ft = -(Kt * pow(d, 0.5) * s + Ct * pow(d, 0.25) * Vtt); //tangential force

							if ((abs(ft) >= mu * abs(Fn)) && Vtt > 0) 
							{
								Ft = -mu * abs(Fn);

								PP[i][i + 3 * no_particles][1] = 0; // slip
								PP[i][i + 3 * no_particles][2] = 0;
								PP[i][i + 3 * no_particles][3] = 0;


								P[i][8] += Ft * Vttx / Vtt;    // tangential force components // cannot usse sx etc.. since its overlap for current timestep
								P[i][9] += Ft * Vtty / Vtt;
								P[i][10] += Ft * Vttz / Vtt;


								temp = -rad_p1 * (ny * Ft * Vttz / Vtt - nz * Ft * Vtty / Vtt);   // torque 
								temp1 = rad_p1 * (nx * Ft * Vttz / Vtt - nz * Ft * Vttx / Vtt);
								temp2 = -rad_p1 * (nx * Ft * Vtty / Vtt - ny * Ft * Vttx / Vtt);

								P[i][14] += temp;   // torque updation
								P[i][15] += temp1;
								P[i][16] += temp2;


							}
							else
							{
								Ft = -min(abs(ft), mu * abs(Fn));   //  tangnetial force with limit

								fx = -(Kt * pow(d, 0.5) * PP[i][i + 3 * no_particles][1] + Ct * pow(d, 0.25) * Vttx);
								fy = -(Kt * pow(d, 0.5) * PP[i][i + 3 * no_particles][2] + Ct * pow(d, 0.25) * Vtty);
								fz = -(Kt * pow(d, 0.5) * PP[i][i + 3 * no_particles][3] + Ct * pow(d, 0.25) * Vttz);


								P[i][8] += Ft * fx / ft;    // tangential force components // cannot usse sx etc.. since its overlap for current timestep
								P[i][9] += Ft * fy / ft;
								P[i][10] += Ft * fz / ft;


								temp = -rad_p1 * (ny * Ft * fz / ft - nz * Ft * fy / ft);   // torque updation
								temp1 = rad_p1 * (nx * Ft * fz / ft - nz * Ft * fx / ft);
								temp2 = -rad_p1 * (nx * Ft * fy / ft - ny * Ft * fx / ft);

								P[i][14] += temp;
								P[i][15] += temp1;
								P[i][16] += temp2;
							}
							
						}
			
					}
				}
				else if (zb == binz - 1)
				{
					d = rad_p1 - (dom_up_z - P[i][4]);
					if (d > 0)
					{
						Vz = P[i][7];
						Fn = -Kn * sqrt(2) * pow(d, 1.5) + Cn * sqrt(2 * sqrt(2)) * pow(d, 0.25) * Vz;
						P[i][10] += Fn;

						nz = (P[i][4] - dom_up_z) / abs(P[i][4] - dom_up_z);


						Vtx = P[i][5]; // subtracting normal componentrs to get tangential component of center of mass velocity.
						Vty = P[i][6];
						Vtz = 0;

						Vtwx = -rad_p1 * P[i][12] * nz; // Vtx= Vt1x-Vt2x ; micro sliding due to rotation after taking v = w x(cross) r of each particle and subtracting to get relative tangential velocity.
						Vtwy = rad_p1 * P[i][11] * nz;
						Vtwz = 0;

						Vttx = Vtx + Vtwx; // total micro sliding components in the tangential directions
						Vtty = Vty + Vtwy;
						Vttz = Vtz + Vtwz;

						Vtt = sqrt(pow(Vttx, 2) + pow(Vtty, 2) + pow(Vttz, 2)); //magnitude of total tangential velocity

						sx = Vttx * timestep; // tangential overlap 
						sy = Vtty * timestep;
						sz = Vttz * timestep;

						if (PP[i][i + 3 * no_particles][0] == tcounter - 1) // if wall is present, periopdic particle allocation ids can be used for other purposes.
						{											// here, paticle to wall id is taken as particle to the same periodic particle in that axis.
							PP[i][i + 3 * no_particles][1] = PP[i][i + 3 * no_particles][1] + sx;
							PP[i][i + 3 * no_particles][2] = PP[i][i + 3 * no_particles][2] + sy;
							PP[i][i + 3 * no_particles][3] = PP[i][i + 3 * no_particles][3] + sz;
							PP[i][i + 3 * no_particles][0] = tcounter;
						}
						else
						{
							PP[i][i + 3 * no_particles][1] = sx;
							PP[i][i + 3 * no_particles][2] = sy;
							PP[i][i + 3 * no_particles][3] = sz;
							PP[i][i + 3 * no_particles][0] = tcounter;
						}

						s = sqrt(pow(PP[i][i + 3 * no_particles][1], 2) + pow(PP[i][i + 3 * no_particles][2], 2) + pow(PP[i][i + 3 * no_particles][3], 2)); // magnitude of total tangential overlap

						if (s != 0)
						{
							ft = -(Kt * pow(d, 0.5) * s + Ct * pow(d, 0.25) * Vtt); //tangential force

							if ((abs(ft) >= mu * abs(Fn)) && Vtt > 0)
							{
								Ft = -mu * abs(Fn);

								PP[i][i + 3 * no_particles][1] = 0; //slip
								PP[i][i + 3 * no_particles][2] = 0;
								PP[i][i + 3 * no_particles][3] = 0;


								P[i][8] += Ft * Vttx / Vtt;    // tangential force components // cannot usse sx etc.. since its overlap for current timestep
								P[i][9] += Ft * Vtty / Vtt;
								P[i][10] += Ft * Vttz / Vtt;


								temp = -rad_p1 * (ny * Ft * Vttz / Vtt - nz * Ft * Vtty / Vtt);   // torque 
								temp1 = rad_p1 * (nx * Ft * Vttz / Vtt - nz * Ft * Vttx / Vtt);
								temp2 = -rad_p1 * (nx * Ft * Vtty / Vtt - ny * Ft * Vttx / Vtt);

								P[i][14] += temp;   // torque updation
								P[i][15] += temp1;
								P[i][16] += temp2;


							}
							else
							{
								Ft = -min(abs(ft), mu * abs(Fn));   //  tangnetial force with limit

								fx = -(Kt * pow(d, 0.5) * PP[i][i + 3 * no_particles][1] + Ct * pow(d, 0.25) * Vttx);
								fy = -(Kt * pow(d, 0.5) * PP[i][i + 3 * no_particles][2] + Ct * pow(d, 0.25) * Vtty);
								fz = -(Kt * pow(d, 0.5) * PP[i][i + 3 * no_particles][3] + Ct * pow(d, 0.25) * Vttz);


								P[i][8] += Ft * fx / ft;    // tangential force components // cannot usse sx etc.. since its overlap for current timestep
								P[i][9] += Ft * fy / ft;
								P[i][10] += Ft * fz / ft;


								temp = -rad_p1 * (ny * Ft * fz / ft - nz * Ft * fy / ft);   // torque updation
								temp1 = rad_p1 * (nx * Ft * fz / ft - nz * Ft * fx / ft);
								temp2 = -rad_p1 * (nx * Ft * fy / ft - ny * Ft * fx / ft);

								P[i][14] += temp;
								P[i][15] += temp1;
								P[i][16] += temp2;
							}


						}
					}
				}

			}break;
			case(3):
			{
				if (zb == 1)
				{
					P[3 * no_particles + i][0] = -1;
					P[3 * no_particles + i][1] = 2;
					P[3 * no_particles + i][2] = P[i][2];
					P[3 * no_particles + i][3] = P[i][3];
					P[3 * no_particles + i][4] = P[i][4] + (duz - dlz);
					P[3 * no_particles + i][5] = P[i][5];
					P[3 * no_particles + i][6] = P[i][6];
					P[3 * no_particles + i][7] = P[i][7];
					P[3 * no_particles + i][11] = P[i][11];
					P[3 * no_particles + i][12] = P[i][12];
					P[3 * no_particles + i][13] = P[i][13];
					
					
						if (tcounter != cbin[xb][yb][binz - 1][1]) // to implement conditions when a bin is populated 1st time within a particular timestep iteration.
						{
							#pragma omp critical
							{
								cbin[xb][yb][binz - 1][0] = 0;
								cbin[xb][yb][binz - 1][1] = tcounter;
							}
						}
                       
						#pragma omp critical
						{
							cbin[xb][yb][binz - 1][0]++;
							cbin[xb][yb][binz - 1][cbin[xb][yb][binz - 1][0] + 1] = 3 * no_particles + i;
							ghost_particles++;
						}
						
				}
				else if (zb == binz - 2)
				{
					P[3 * no_particles + i][0] = -1;
					P[3 * no_particles + i][1] = 2;
					P[3 * no_particles + i][2] = P[i][2];
					P[3 * no_particles + i][3] = P[i][3];
					P[3 * no_particles + i][4] = P[i][4] - (duz - dlz);
					P[3 * no_particles + i][5] = P[i][5];
					P[3 * no_particles + i][6] = P[i][6];
					P[3 * no_particles + i][7] = P[i][7];
					P[3 * no_particles + i][11] = P[i][11];
					P[3 * no_particles + i][12] = P[i][12];
					P[3 * no_particles + i][13] = P[i][13];
					
					
						if (tcounter != cbin[xb][yb][0][1]) // to implement conditions when a bin is populated 1st time within a particular timestep iteration.
						{	
							#pragma omp critical
							{	
								cbin[xb][yb][0][0] = 0;
								cbin[xb][yb][0][1] = tcounter;
							}
						}
						#pragma omp critical
						{
							cbin[xb][yb][0][0]++;
							cbin[xb][yb][0][cbin[xb][yb][0][0] + 1] = 3 * no_particles + i;
							ghost_particles++;
						}
				}

			}break;

			default:
				break;
			}

			if (P[i][0] == -1) // to check if the particle have gone out in the current/on going iteration.
			{
				continue;
			}
			
			
				if (tcounter != cbin[xb][yb][zb][1]) // to implement conditions when a bin is populated 1st time within a particular timestep iteration.
				{
					#pragma omp critical
					{
						cbin[xb][yb][zb][0] = 0;
						cbin[xb][yb][zb][1] = tcounter;
					}
				}
				#pragma omp critical			
				{	 
					 cbin[xb][yb][zb][0]++;     // filling particles in the contact detection array - particle count	
					 cbin[xb][yb][zb][cbin[xb][yb][zb][0] + 1] = i; // 
				}

		}
	
		
       #pragma omp barrier
       #pragma omp parallel for  private(Fx,Fy,Fz,temp, temp1, temp2,cx, cy, cz, Vx, Vy, Vz, Vtx, Vty, Vtz, Vnx, Vny, Vnz, Vtwx, Vtwy, Vtwz, Vdot_t, sx, sy, sz, s, Vttx, Vtty, Vttz, Vtt,d, cd, Fn, nx, ny, nz, Vrn, fx, fy, fz, ft, Ft)
		for (int q = 0; q < no_particles; q++)
		{
			if (P[q][0] == 1 || P[q][0] == -1) // to check if the particle is already considered or out of domain.
			{
				continue;
			}

			P[q][0] = 1;  // assign already considered.

			int i = floor((P[q][2] - dom_low_x) / binl);
			 int j = floor((P[q][3] - dom_low_y) / binl);
			 int k = floor((P[q][4] - dom_low_z) / binl);

			 if (tcounter % dumpevery == 0)
			 {
				 P[q][17] = 0;  // collissional stress tensor initialising;
				 P[q][18] = 0;
				 P[q][19] = 0;
				 P[q][20] = 0;
				 P[q][21] = 0;
				 P[q][22] = 0;

			 }

			 



				for (int ii = i - 1; ii <= i + 1; ii++)
				{
					if (ii < 0 || ii >= binx) // to check whether ii represents a cell outside the domain
					{
						continue;
					}
					for (int jj = j - 1; jj <= j + 1; jj++)
					{
						if (jj < 0 || jj >= biny)
						{
							continue;
						}
						for (int kk = k - 1; kk <= k + 1; kk++)
						{
							if (kk < 0 || kk >= binz)
							{
								continue;
							}
							if (tcounter == cbin[ii][jj][kk][1]) // checking whether the adjacent cell is newly populated
							{


								if (ii == i && jj == j && kk == k) // checking whether the loop is at the primary mid cell
								{


										for (int p2 = 2; p2 <= cbin[i][j][k][0] + 1; p2++) // going through all the particles in the primary mid cell 1 by 1
										{

											 int cp2 = cbin[i][j][k][p2];  // id of the contact particle (cp2)

											if (cp2 == q) // to check if the 2nd particle is same as 1st particle
											{
												continue;
											}

											cx = P[q][2] - P[cp2][2];
											cy = P[q][3] - P[cp2][3];
											cz = P[q][4] - P[cp2][4];
											cd = sqrt(pow(cx, 2) + pow(cy, 2) + pow(cz, 2));  // distance between 2 particles


											d = 2 * rad_p1 - cd; // overlap

											if (d > 0)
											{
												nx = cx / cd;   // directioion unit vector
												ny = cy / cd;
												nz = cz / cd;

												Vx = P[q][5] - P[cp2][5];  // relative velocities
												Vy = P[q][6] - P[cp2][6];
												Vz = P[q][7] - P[cp2][7];

												Vrn = nx * Vx + ny * Vy + nz * Vz; // relative velocity in the normal direction

												Fn = Kn * pow(d, 1.5) + Cn * pow(d, 0.25) * Vrn;   // normal force
												fx = nx * Fn;
												fy = ny * Fn;
												fz = nz * Fn;

												P[q][8] += fx;
												P[q][9] += fy;
												P[q][10] += fz;

												if (tcounter % dumpevery == 0)
												{
													Fx = fx;  // to later take the total sum of normal and tangential forces.
													Fy = fy;
													Fz = fz;
												}
												Vtx = Vx - Vrn * nx; // subtracting normal componentrs to get tangential component of center of mass velocity.
												Vty = Vy - Vrn * ny;
												Vtz = Vz - Vrn * nz;

												Vtwx = -rad_p1 * ((P[cp2][12] * nz - ny * P[cp2][13]) + (P[q][12] * nz - ny * P[q][13])); // Vtx= Vt1x-Vt2x ; micro sliding due to rotation after taking v = w x(cross) r of each particle and subtracting to get relative tangential velocity.
												Vtwy = rad_p1 * ((P[cp2][11] * nz - nx * P[cp2][13]) + (P[q][11] * nz - nx * P[q][13]));
												Vtwz = -rad_p1 * ((P[cp2][12] * ny - nx * P[cp2][12]) + (P[q][11] * ny - nx * P[q][12]));



												Vttx = Vtx + Vtwx; // total micro sliding components in the tangential directions
												Vtty = Vty + Vtwy;
												Vttz = Vtz + Vtwz;

												Vtt = sqrt(pow(Vttx, 2) + pow(Vtty, 2) + pow(Vttz, 2)); //magnitude of total tangential velocity

												sx = Vttx * timestep; // tangential overlap 
												sy = Vtty * timestep;
												sz = Vttz * timestep;

												if (PP[q][cp2][0] == tcounter - 1)  // to check ongoing contact
												{

													PP[q][cp2][1] = PP[q][cp2][1] + sx;
													PP[q][cp2][2] = PP[q][cp2][2] + sy;
													PP[q][cp2][3] = PP[q][cp2][3] + sz;
													PP[q][cp2][0] = tcounter;

												
												}
												else                    // if new contact
												{
													PP[q][cp2][1] = sx;
													PP[q][cp2][2] = sy;
													PP[q][cp2][3] = sz;
													PP[q][cp2][0] = tcounter;

												
												}

												s = sqrt(pow(PP[q][cp2][1], 2) + pow(PP[q][cp2][2], 2) + pow(PP[q][cp2][3], 2)); // magnitude of total tangential overlap

												if (s != 0)
												{

													ft = -(Kt * pow(d, 0.5) * s + Ct * pow(d, 0.25) * Vtt); //tangential force

													if ((abs(ft) >= mu * abs(Fn)) && Vtt > 0)
													{
														Ft = -mu * abs(Fn);

														PP[q][cp2][1] = 0; // slip occurs, so tangential displacement = 0.
														PP[q][cp2][2] = 0;
														PP[q][cp2][3] = 0;



														P[q][8] += Ft * Vttx / Vtt;    // tangential force components // cannot usse sx etc.. since its overlap for current timestep
														P[q][9] += Ft * Vtty / Vtt;
														P[q][10] += Ft * Vttz / Vtt;

														if (tcounter % dumpevery == 0)
														{

															Fx += Ft * Vttx / Vtt; // summing both normal force and tangential force for a particular interaction.
															Fy += Ft * Vtty / Vtt;
															Fz += Ft * Vttz / Vtt;

															P[q][17] += 0.5 * Fx * cx;  // collissional stress tensor;
															P[q][18] += 0.5 * Fy * cy;
															P[q][19] += 0.5 * Fz * cz;
															P[q][20] += 0.5 * Fx * cy;
															P[q][21] += 0.5 * Fx * cz;
															P[q][22] += 0.5 * Fy * cz;

														}


														temp = -rad_p1 * (ny * Ft * Vttz / Vtt - nz * Ft * Vtty / Vtt);   // torque updation
														temp1 = rad_p1 * (nx * Ft * Vttz / Vtt - nz * Ft * Vttx / Vtt);
														temp2 = -rad_p1 * (nx * Ft * Vtty / Vtt - ny * Ft * Vttx / Vtt);



														P[q][14] += temp;
														P[q][15] += temp1;
														P[q][16] += temp2;



													}
													else
													{
														Ft = -min(abs(ft), mu * abs(Fn));   //  tangnetial force with limit

														fx = -(Kt * pow(d, 0.5) * PP[q][cp2][1] + Ct * pow(d, 0.25) * Vttx);
														fy = -(Kt * pow(d, 0.5) * PP[q][cp2][2] + Ct * pow(d, 0.25) * Vtty);
														fz = -(Kt * pow(d, 0.5) * PP[q][cp2][3] + Ct * pow(d, 0.25) * Vttz);


														P[q][8] += Ft * fx / ft;    // tangential force components // cannot usse sx etc.. since its overlap for current timestep
														P[q][9] += Ft * fy / ft;
														P[q][10] += Ft * fz / ft;

														if (tcounter % dumpevery == 0)
														{

															Fx += Ft * fx / ft; // summing both normal force and tangential force for a particular interaction.
															Fy += Ft * fy / ft;
															Fz += Ft * fz / ft;

															P[q][17] += 0.5 * Fx * cx;  // collissional stress tensor;
															P[q][18] += 0.5 * Fy * cy;
															P[q][19] += 0.5 * Fz * cz;
															P[q][20] += 0.5 * Fx * cy;
															P[q][21] += 0.5 * Fx * cz;
															P[q][22] += 0.5 * Fy * cz;

														}

														temp = -rad_p1 * (ny * Ft * fz / ft - nz * Ft * fy / ft);   // torque updation
														temp1 = rad_p1 * (nx * Ft * fz / ft - nz * Ft * fx / ft);
														temp2 = -rad_p1 * (nx * Ft * fy / ft - ny * Ft * fx / ft);



														P[q][14] += temp;
														P[q][15] += temp1;
														P[q][16] += temp2;


													

													}

												
												}
											}
										}
									
								}
								else
								{

									for (int p2 = 2; p2 <= cbin[ii][jj][kk][0] + 1; p2++)
									{

										int cp2 = cbin[ii][jj][kk][p2];
										
										cx = P[q][2] - P[cp2][2];
										cy = P[q][3] - P[cp2][3];
										cz = P[q][4] - P[cp2][4];
										cd = sqrt(pow(cx, 2) + pow(cy, 2) + pow(cz, 2));


										d = 2 * rad_p1 - cd;

										if (d > 0)
										{
											nx = cx / cd;
											ny = cy / cd;
											nz = cz / cd;

											Vx = P[q][5] - P[cp2][5];
											Vy = P[q][6] - P[cp2][6];
											Vz = P[q][7] - P[cp2][7];

											Vrn = nx * Vx + ny * Vy + nz * Vz;

											Fn = Kn * pow(d, 1.5) + Cn * pow(d, 0.25) * Vrn;

											fx = nx * Fn;
											fy = ny * Fn;
											fz = nz * Fn;

											if (tcounter % dumpevery == 0)
											{
												Fx = fx;  // to later take the total sum of normal and tangential forces.
												Fy = fy;
												Fz = fz;
											}

											P[q][8] += fx;
											P[q][9] += fy;
											P[q][10] += fz;

											Vtx = Vx - Vrn * nx; // subtracting normal components to get tangential component of center of mass velocity.
											Vty = Vy - Vrn * ny;
											Vtz = Vz - Vrn * nz;

											Vtwx = -rad_p1 * ((P[cp2][12] * nz - ny * P[cp2][13]) + (P[q][12] * nz - ny * P[q][13])); // Vtwx= Vtw2x-Vt1wx ; micro sliding due to relative rotation after taking v = w x(cross) r of each particle  to get relative tangential velocity.
											Vtwy = rad_p1 * ((P[cp2][11] * nz - nx * P[cp2][13]) + (P[q][11] * nz - nx * P[q][13]));
											Vtwz = -rad_p1 * ((P[cp2][12] * ny - nx * P[cp2][12]) + (P[q][11] * ny - nx * P[q][12]));


											Vttx = Vtx + Vtwx; // total micro sliding components in the tangential directions
											Vtty = Vty + Vtwy;
											Vttz = Vtz + Vtwz;

											Vtt = sqrt(pow(Vttx, 2) + pow(Vtty, 2) + pow(Vttz, 2)); //magnitude of total tangential velocity

											sx = Vttx * timestep; // tangential overlap 
											sy = Vtty * timestep;
											sz = Vttz * timestep;

											if (PP[q][cp2][0] == tcounter - 1)
											{
												PP[q][cp2][1] = PP[q][cp2][1] + sx;
												PP[q][cp2][2] = PP[q][cp2][2] + sy;
												PP[q][cp2][3] = PP[q][cp2][3] + sz;
												PP[q][cp2][0] = tcounter;

												if (cp2 >= no_particles)   // for periodic conditions
												{
													temp = cp2 % no_particles;
													PP[q][temp][1] = PP[q][temp][1] + sx;
													PP[q][temp][2] = PP[q][temp][2] + sy;
													PP[q][temp][3] = PP[q][temp][3] + sz;
													PP[q][temp][0] = tcounter;

												}

											}
											else
											{
												PP[q][cp2][1] = sx;
												PP[q][cp2][2] = sy;
												PP[q][cp2][3] = sz;
												PP[q][cp2][0] = tcounter;

												if (cp2 >= no_particles) // periodic conditions
												{
													temp = cp2 % no_particles;
													PP[q][temp][1] = sx;
													PP[q][temp][2] = sy;
													PP[q][temp][3] = sz;
													PP[q][temp][0] = tcounter;

												}

											}

											s = sqrt(pow(PP[q][cp2][1], 2) + pow(PP[q][cp2][2], 2) + pow(PP[q][cp2][3], 2)); // magnitude of total tangential overlap


											if (s != 0)
											{
												ft = -(Kt * pow(d, 0.5) * s + Ct * pow(d, 0.25) * Vtt); //tangential force

												if ((abs(ft) >= mu * abs(Fn)) && Vtt > 0)
												{
													Ft = -mu * abs(Fn);

													PP[q][cp2][1] = 0; // slip occurs, so tangential displacement = 0.
													PP[q][cp2][2] = 0;
													PP[q][cp2][3] = 0;

													if (cp2 >= no_particles)   // for periodic conditions 
													{
														temp = cp2 % no_particles;
														PP[q][temp][1] = 0;
														PP[q][temp][2] = 0;
														PP[q][temp][3] = 0;


													}


													P[q][8] += Ft * Vttx / Vtt;    // tangential force components // cannot usse sx etc.. since its overlap for current timestep
													P[q][9] += Ft * Vtty / Vtt;
													P[q][10] += Ft * Vttz / Vtt;

													if (tcounter % dumpevery == 0)
													{

														Fx += Ft * Vttx / Vtt; // summing both normal force and tangential force for a particular interaction.
														Fy += Ft * Vtty / Vtt;
														Fz += Ft * Vttz / Vtt;

														P[q][17] += 0.5 * Fx * cx;  // collissional stress tensor;
														P[q][18] += 0.5 * Fy * cy;
														P[q][19] += 0.5 * Fz * cz;
														P[q][20] += 0.5 * Fx * cy;
														P[q][21] += 0.5 * Fx * cz;
														P[q][22] += 0.5 * Fy * cz;

													}


													temp = -rad_p1 * (ny * Ft * Vttz / Vtt - nz * Ft * Vtty / Vtt);   // torque 
													temp1 = rad_p1 * (nx * Ft * Vttz / Vtt - nz * Ft * Vttx / Vtt);
													temp2 = -rad_p1 * (nx * Ft * Vtty / Vtt - ny * Ft * Vttx / Vtt);

													P[q][14] += temp;   // torque updation
													P[q][15] += temp1;
													P[q][16] += temp2;


												}
												else
												{
													Ft = -min(abs(ft), mu * abs(Fn));   //  tangnetial force with limit

													fx = -(Kt * pow(d, 0.5) * PP[q][cp2][1] + Ct * pow(d, 0.25) * Vttx);
													fy = -(Kt * pow(d, 0.5) * PP[q][cp2][2] + Ct * pow(d, 0.25) * Vtty);
													fz = -(Kt * pow(d, 0.5) * PP[q][cp2][3] + Ct * pow(d, 0.25) * Vttz);


													P[q][8] += Ft * fx / ft;    // tangential force components // cannot usse sx etc.. since its overlap for current timestep
													P[q][9] += Ft * fy / ft;
													P[q][10] += Ft * fz / ft;

													if (tcounter % dumpevery == 0)
													{

														Fx += Ft * fx / ft; // summing both normal force and tangential force for a particular interaction.
														Fy += Ft * fy / ft;
														Fz += Ft * fz / ft;

														P[q][17] += 0.5 * Fx * cx;  // collissional stress tensor;
														P[q][18] += 0.5 * Fy * cy;
														P[q][19] += 0.5 * Fz * cz;
														P[q][20] += 0.5 * Fx * cy;
														P[q][21] += 0.5 * Fx * cz;
														P[q][22] += 0.5 * Fy * cz;

													}


													temp = -rad_p1 * (ny * Ft * fz / ft - nz * Ft * fy / ft);   // torque updation
													temp1 = rad_p1 * (nx * Ft * fz / ft - nz * Ft * fx / ft);
													temp2 = -rad_p1 * (nx * Ft * fy / ft - ny * Ft * fx / ft);

													P[q][14] += temp;
													P[q][15] += temp1;
													P[q][16] += temp2;
												}


											}


										}
									}


								}
							}
						}
					}

				}

				P[q][5] = P[q][5] + 0.5 * timestep * P[q][8] / m;  // other half velocity updation, VV scheme.
				P[q][6] = P[q][6] + 0.5 * timestep * P[q][9] / m;
				P[q][7] = P[q][7] + 0.5 * timestep * P[q][10] / m;

				P[q][11] = P[q][11] + 0.5 * timestep * P[q][14] / ((2 * m * rad_p1 * rad_p1) / 5);  // other half omega updation, VV scheme.
				P[q][12] = P[q][12] + 0.5 * timestep * P[q][15] / ((2 * m * rad_p1 * rad_p1) / 5);
				P[q][13] = P[q][13] + 0.5 * timestep * P[q][16] / ((2 * m * rad_p1 * rad_p1) / 5);



			
		}
		
		if (tcounter % dumpevery == 0)
		{
			#pragma omp barrier
			
			cout << "timestep: ";
			cout << tcounter << "\n";
			oFile << "ITEM: TIMESTEP" << endl;
			oFile << tcounter << "\n";

			oFile << "ITEM: NUMBER OF ATOMS\n";
			oFile << no_particles - exitparticles << "\n";
			oFile << "ITEM: BOX BOUNDS ff ff ff\n";
			oFile << dlx << "  " << dux << "\n";
			oFile << dly << "  " << duy << "\n";
			oFile << dlz << "  " << duz << "\n";
			oFile << "ITEM: ATOMS id type x y z vx vy vz fx fy fz wx wy wz Tx Ty Tz cxx cyy czz cxy cxz cyz radius\n";


			for (int i = 0; i < no_particles; i++)
			{
				if (P[i][0] == -1)
				{
					continue;
				}
				oFile << i + 1 << "  ";
				for (int j = 1; j < 23; j++)
				{
					//if (i == 16 && j == 1)     // checking a particular particle
					//{
					//	oFile << 2 << "  ";
					//	continue;
					//}
					oFile << P[i][j] << "  ";

				}


				if (i != no_particles - 1)
				{
					oFile << rad_p1 << "\n";
				}
				else
				{
					oFile << rad_p1 << endl;  // to clear buffer ie.. to print after each timestep, 'endl' must be used otherwise it gets printed totaly at end of program.
				}
			}

		}
		tcounter++;
	}
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);
	cout << floor(duration.count() / 60.0) << " min and " << duration.count() % 60 << " sec" << endl;

	return(0);
}