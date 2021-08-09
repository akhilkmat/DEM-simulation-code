
#include <iostream>
#include <fstream>
#include <cmath>
#include <string> 
#include <vector>
#include <limits>
#include <cstdlib>
#include <cstdio>
#include <sstream>

using namespace std;

int main()
{

	float id = 0, x = 0, y = 0, z = 0, ux = 0, uy = 0, uz = 0, us = 0, vx = 0, vy = 0, vz = 0, m = 0, s = 0, s_old = 0, fx = 0, fy = 0, fz = 0, wx = 0, wy = 0, wz = 0,Tx,Ty,Tz, cxx = 0, cyy = 0, czz = 0, cxy = 0, cxz = 0, cyz = 0, r = 0, count = 0;
	//float xb = 0.02, yb = 0.001, zb = 0.001, xl = -0.01, xu = 0.01, yl = 0.0, yu = 0.03, zl = 0.0, zu = 0.03, tstart = 80000, tstop = 100000,Vc=0,rho=2500;
	float xb = 0.05, yb = 0.05, zb = 0.01, xl = 0, xu = 0.05, yl = 0, yu = 0.05, zl = 0.0, zu = 1, tstart = 1000000, tstop = 5000000, Vc = 0, rho = 2500;
	int flag = 0,mark=0, bin = 0, i = 0, j = 0, k = 0, lc = 0, binold = 0, type = 0;
	double pi = 2 * acos(0.0);
	ifstream inFile;
	ofstream oFile;
	int tbx = 0, tby = 0, tbz = 0, tbin = 0;

	/*cout << endl<<"bin size in x ddirection(SI units) : ";
	cin >> xb;
	cout << endl << "bin size in y ddirection : ";
	cin >> yb;
	cout << endl << "bin size in z ddirection : ";
	cin >> zb;
	cout << endl << "enter lower limit in x direction : ";
	cin >> xl;
	cout << endl << "enter upper limit in x direction : ";
	cin >> xu;
	cout << endl << "enter lower limit in y direction : ";
	cin >> yl;
	cout << endl << "enter upper limit in y direction : ";
	cin >> yu;
	cout << endl << "enter lower limit in z direction : ";
	cin >> zl;
	cout << endl << "enter upper limit in z direction : ";
	cin >> zu;
	cout << " enter the time step after which averaging for granular temperature should be done : ";
	cin >> tstart;
	cout << " enter the time step at which averaging should be stopped. (includes the timestep too) : ";
	cin >> tstop;
	cout << endl << "enter the density of the particle : ";
	cin >> rho;
	*/
	tbx = ceil((xu - xl) / xb);
	tby = ceil((yu - yl) / yb);
	tbz = ceil((zu - zl) / zb);

	tbin = int(tbx * tby * tbz);
	cout << tbin << "  " << tbx << "  " << tby << "  " << tbz;

	Vc = xb * yb * zb;
	inFile.open("dump.atom");
	if (!inFile)
	{
		cerr << "Unable to open file dump.atom";
		exit(1);
	}

	//oFile.open("codeinf_%tstart_.dat");
	oFile.open(("long 0.5cof 0.4cor " + std::to_string(tstart) + ".dat")); // write file name can be changed here eg: "0.4cor_0.4cof_cdata"
	if (!oFile)
	{
		cerr << "Unable to open file codeinf.dat";
		exit(1);
	}
	vector<vector<vector<float>>> Vbin(tbin, vector<vector<float>> // for storing data associated with vx, vy, vz, speed, particle density , radius, and stresses respectively.
		(7, vector<float>
			(200, 0)));  // adjust (x,0) no according to what would be the approximate no of particles per bin.





	float** inf2 = new float* [tbin];
	for (i = 0; i < tbin; i++)
	{
		inf2[i] = new float[19];
	}
	for (i = 0; i < tbin; i++)
	{
		for (j = 0; j < 19; j++)
		{
			inf2[i][j] = 0;
		}
	}




	string dummyLine;
	float timestep = 0, No_of_atoms = 0;
	flag = 1;


	while (flag)
	{

		for (i = 0; i < tbin; i++)
		{
			for (j = 0; j < 7; j++)
			{
				for (k = 0; k < 100; k++)
				{
					Vbin[i][j][k] = 0;
				}
			}
		}



		getline(inFile, dummyLine);
		if (inFile.eof())
			flag = 0;
		inFile >> timestep;
		getline(inFile, dummyLine);
		getline(inFile, dummyLine);
		inFile >> No_of_atoms;
		getline(inFile, dummyLine);

		for (int i_2 = 1; i_2 <= 5; i_2++)
			getline(inFile, dummyLine);

		if (timestep > tstop)
		{
			flag = 0;
		}
		else if (timestep < tstart)
		{
			for (int i_2 = 0; i_2 < No_of_atoms; i_2++)
			{
				inFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			}
		}
		else if (flag)
		{
			count = count + 1; // no of timesteps processed.
			for (int i_2 = 0; i_2 < No_of_atoms; i_2++)
			{
				inFile >> mark >> type >> x >> y >> z >> vx >> vy >> vz >> fx >> fy >> fz >>wx>>wy>>wz>>Tx>>Ty>>Tz>> cxx >> cyy >> czz >> cxy >> cxz >> cyz >> r; // reading values separated by spaces from input file,edit these quantities according to your input dump file.
				if (x < xl || y < yl || z < zl)
				{
					continue;
				}
				bin = ceil((x - xl) / xb) + floor((y - yl) / yb) * tbx + floor((z - zl) / zb) * tbx * tby; // binning of particles.
				//cout << timestep << "  " << id << "  " << bin << endl;
				s = sqrt(pow(vx, 2) + pow(vy, 2) + pow(vz, 2)); //speed
				m = rho * (4 / 3) * pi * pow(r, 3); // mass

				Vbin[bin - 1][0][0] += vx; // bin cumilative velocity of different particles in that bin.
				Vbin[bin - 1][1][0] += vy;
				Vbin[bin - 1][2][0] += vz;
				Vbin[bin - 1][3][0] += s;

				Vbin[bin - 1][3][1] += 1; // count of particles in that bin, same for vx,vy,vz bins too and hence left blank for middle indexes 0,1 and 2.

				Vbin[bin - 1][0][2] = Vbin[bin - 1][0][0] / Vbin[bin - 1][3][1]; // average speed of particles in that bin
				Vbin[bin - 1][1][2] = Vbin[bin - 1][1][0] / Vbin[bin - 1][3][1];
				Vbin[bin - 1][2][2] = Vbin[bin - 1][2][0] / Vbin[bin - 1][3][1];
				Vbin[bin - 1][3][2] = Vbin[bin - 1][3][0] / Vbin[bin - 1][3][1];

				lc = 4 + Vbin[bin - 1][3][1]; // means ,4 + no of particles, making a variable for indexing while storing speed info of each particle.

			   /*if (lc >= 99)
				{
					cout << endl<<lc-4;
					cout << endl<<timestep << "  " << id << "  " << bin << endl;
				}*/

				Vbin[bin - 1][0][lc] = vx; // storing speed of each particles in separate locations in inf array row.
				Vbin[bin - 1][1][lc] = vy;
				Vbin[bin - 1][2][lc] = vz;
				Vbin[bin - 1][3][lc] = s;
				if (type == 1)
				{
					Vbin[bin - 1][4][lc] = 2500;   // defining density for the particle type.
				}
				Vbin[bin - 1][5][lc] = r; // define radius for each particle stored.

				//cxx += m * pow(vx, 2); // only contact stresses remains after these operations (for LIGGGHTS only).
				//cyy += m * pow(vy, 2);
				//czz += m * pow(vz, 2);
				//cxy += m * vx * vy;
				//cxz += m * vx * vz;
				//cyz += m * vy * vz;

				Vbin[bin - 1][3][4] += (cxx + cyy + czz); // bin cumilative contact pressure of all particles in the bin , only in the speed bin.

				Vbin[bin - 1][6][0] += cxx;
				Vbin[bin - 1][6][1] += cyy;
				Vbin[bin - 1][6][2] += czz;
				Vbin[bin - 1][6][3] += cxy;
				Vbin[bin - 1][6][4] += cxz;
				Vbin[bin - 1][6][5] += cyz;

			}
			getline(inFile, dummyLine);//to clear the end of a timestep's read data.

			for (i = 0; i < tbin; i++)
			{
				if (Vbin[i][3][1] != 0) // checking bin is not empty
				{
					for (j = 5; j <= (Vbin[i][3][1] + 4); j++)
					{
						ux = Vbin[i][0][j] - Vbin[i][0][2]; //delta v
						uy = Vbin[i][1][j] - Vbin[i][1][2];
						uz = Vbin[i][2][j] - Vbin[i][2][2];
						us = Vbin[i][3][j] - Vbin[i][3][2];

						Vbin[i][0][3] += pow(ux, 2); // bin cumilative delta v squared.
						Vbin[i][1][3] += pow(uy, 2);
						Vbin[i][2][3] += pow(uz, 2);
						Vbin[i][3][3] += pow(us, 2);

						m = Vbin[i][4][j] * (4 / 3) * pi * pow(Vbin[i][5][j], 3); // mass of that particle.

						Vbin[i][6][0] -= m * ux * ux; // including streaming stresses
						Vbin[i][6][1] -= m * uy * uy;
						Vbin[i][6][2] -= m * uz * uz;
						Vbin[i][6][3] -= m * ux * uy;
						Vbin[i][6][4] -= m * ux * uz;
						Vbin[i][6][5] -= m * uy * uz;

					}

					inf2[i][0] = ((count - 1) * inf2[i][0] + (Vbin[i][0][3] / Vbin[i][3][1])) / count; //  average granular temperatures x,y,z and total
					inf2[i][1] = ((count - 1) * inf2[i][1] + (Vbin[i][1][3] / Vbin[i][3][1])) / count;
					inf2[i][2] = ((count - 1) * inf2[i][2] + (Vbin[i][2][3] / Vbin[i][3][1])) / count;
					inf2[i][3] = ((count - 1) * inf2[i][3] + (Vbin[i][3][3] / Vbin[i][3][1])) / count;

				}
				else
				{
					inf2[i][0] = ((count - 1) * inf2[i][0] + 0) / count; // this step is needed since if no particles are present dividing by particle count will give nan values.
					inf2[i][1] = ((count - 1) * inf2[i][1] + 0) / count;  // - average granular temperatures x,y,z and total
					inf2[i][2] = ((count - 1) * inf2[i][2] + 0) / count;
					inf2[i][3] = ((count - 1) * inf2[i][3] + 0) / count;
				}

				inf2[i][4] = ((count - 1) * inf2[i][4] + (Vbin[i][3][1] * 4 * pi * pow(r, 3) / (3 * Vc))) / count; // avg volume fraction

				inf2[i][5] = inf2[i][4] * rho; // avg density

				inf2[i][6] = ((count - 1) * inf2[i][6] + Vbin[i][0][2]) / count; // avg velocities
				inf2[i][7] = ((count - 1) * inf2[i][7] + Vbin[i][1][2]) / count;
				inf2[i][8] = ((count - 1) * inf2[i][8] + Vbin[i][2][2]) / count;
				inf2[i][9] = ((count - 1) * inf2[i][9] + Vbin[i][3][2]) / count;

				inf2[i][10] = ((count - 1) * inf2[i][10] + -1 * (Vbin[i][3][4] / (3 * Vc))) / count; // avg contact pressure

				inf2[i][11] = ((count - 1) * inf2[i][11] + 1 * (Vbin[i][6][0]) / Vc) / count;  // xx total stress
				inf2[i][12] = ((count - 1) * inf2[i][12] + 1 * (Vbin[i][6][1]) / Vc) / count;  // yy total stress
				inf2[i][13] = ((count - 1) * inf2[i][13] + 1 * (Vbin[i][6][2]) / Vc) / count;  // zz total stress
				inf2[i][14] = ((count - 1) * inf2[i][14] + 1 * (Vbin[i][6][3]) / Vc) / count;  // xy total stress
				inf2[i][15] = ((count - 1) * inf2[i][15] + 1 * (Vbin[i][6][4]) / Vc) / count;  // xz total stress
				inf2[i][16] = ((count - 1) * inf2[i][16] + 1 * (Vbin[i][6][5]) / Vc) / count;  // yz total stress

				inf2[i][17] = -(inf2[i][11] + inf2[i][12] + inf2[i][13]) / 3; // avg total pressure

				inf2[i][18] = inf2[i][17] - inf2[i][10]; // avg streaming pressure

			}

		}

	}

	oFile << "TITLE = LIGGHTS" << endl;
	oFile << "VARIABLES = Bin X Y Z GTx GTy GTz GT VolF Density Vx Vy Vz V Sxx Syy Szz Sxy Sxz Syz CP SP TP" << endl;
	oFile << "ZONE T = \"" << tstart << "\" I = " << tbx << " J = " << tby << " K = " << tbz << " F = POINT" << endl;

	for (i = 1; i <= tbin; i++)
	{
		oFile << i << "  ";
		if ((i % (tby * tbx)) == 0)
		{
			oFile << xl + (tbx - 0.5) * xb << "  " << yl + (tby - 0.5) * yb << "  " << zl + ((i / (tby * tbx)) - 0.5) * zb << "  " << inf2[i - 1][0] << "  " << inf2[i - 1][1] << "  " << inf2[i - 1][2] << " " << inf2[i - 1][3] << " " << inf2[i - 1][4] << " " << inf2[i - 1][5] << " " << inf2[i - 1][6] << " " << inf2[i - 1][7] << " " << inf2[i - 1][8] << " " << inf2[i - 1][9] << " " << inf2[i - 1][11] << " " << inf2[i - 1][12] << " " << inf2[i - 1][13] << " " << inf2[i - 1][14] << " " << inf2[i - 1][15] << " " << inf2[i - 1][16] << " " << inf2[i - 1][10] << " " << inf2[i - 1][18] << " " << inf2[i - 1][17] << " " << endl;
		}
		else if (((i % (tby * tbx)) % tbx) == 0)
		{
			oFile << xl + (tbx - 0.5) * xb << "  " << yl + (((i % (tby * tbx)) / tbx) - 0.5) * yb << "  " << zl + (floor(i / (tby * tbx)) + 0.5) * zb << "  " << inf2[i - 1][0] << "  " << inf2[i - 1][1] << "  " << inf2[i - 1][2] << " " << inf2[i - 1][3] << " " << inf2[i - 1][4] << " " << inf2[i - 1][5] << " " << inf2[i - 1][6] << " " << inf2[i - 1][7] << " " << inf2[i - 1][8] << " " << inf2[i - 1][9] << " " << inf2[i - 1][11] << " " << inf2[i - 1][12] << " " << inf2[i - 1][13] << " " << inf2[i - 1][14] << " " << inf2[i - 1][15] << " " << inf2[i - 1][16] << " " << inf2[i - 1][10] << " " << inf2[i - 1][18] << " " << inf2[i - 1][17] << " " << endl;
		}
		else
		{
			oFile << xl + (((i % (tby * tbx)) % tbx) - 0.5) * xb << "  " << yl + (floor((i % (tby * tbx)) / tbx) + 0.5) * yb << "  " << zl + (floor(i / (tby * tbx)) + 0.5) * zb << "  " << inf2[i - 1][0] << "  " << inf2[i - 1][1] << "  " << inf2[i - 1][2] << " " << inf2[i - 1][3] << " " << inf2[i - 1][4] << " " << inf2[i - 1][5] << " " << inf2[i - 1][6] << " " << inf2[i - 1][7] << " " << inf2[i - 1][8] << " " << inf2[i - 1][9] << " " << inf2[i - 1][11] << " " << inf2[i - 1][12] << " " << inf2[i - 1][13] << " " << inf2[i - 1][14] << " " << inf2[i - 1][15] << " " << inf2[i - 1][16] << " " << inf2[i - 1][10] << " " << inf2[i - 1][18] << " " << inf2[i - 1][17] << " " << endl;
		}



	}
	inFile.close();
	oFile.close();


	for (i = 0; i < tbin; i++)
	{
		delete[] inf2[i];
	}
	delete[] inf2;


	exit(0);



}