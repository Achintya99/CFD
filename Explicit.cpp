

#include <iostream>
#include<vector>
#include<fstream>
#include<math.h>
#include<time.h>
#include<iomanip>

using namespace std;

int main()
{
	clock_t ct;
	
	int nx, ny;
	double dx, dy, dt;
	dt = 0.0005333;
	double ts = 0.777 / dt;

	nx = ny = 7;
	dx = 1 / (double(nx) - 1);
	dy = 1 / (double(ny) - 1);
	
	double alpha = 1.0;
	double lambda = alpha * dt / pow(dx, 2.0);
	cout << lambda << endl;
	double w = 1.9;


	vector<vector<double>> Tg(ny, vector<double>(nx)), Tn(ny, vector<double>(nx)), To(ny, vector<double>(nx));

	//initializing 
	for (int i = 0;i < ny;i++)
	{
		for (int j = 0;j < nx;j++)
		{
			To[i][j] = Tn[i][j] = Tg[i][j] = 0.0;
		}
	}


	//boundary values
	//left
	for (int i = 0;i < ny;i++)
	{
		Tg[i][0] = Tn[i][0] = 300.0;
	}
	//top
	for (int j = 0;j < nx;j++)
	{
		Tg[ny - 1][j] = Tn[ny - 1][j] = 420.0;
	}
	//right
	for (int i = 0;i < ny;i++)
	{
		Tg[i][nx - 1] = Tn[i][nx - 1] = 400.0;
	}
	//bottom
	for (int j = 0;j < nx;j++)
	{
		Tg[0][j] = Tn[0][j] = 325.0;
	}
	//corners
	Tn[0][0] = Tg[0][0] = double(625.0 / 2);
	Tn[ny - 1][0] = Tg[ny - 1][0] = double(720 / 2);
	Tn[ny - 1][nx - 1] = Tg[ny - 1][nx - 1] = 410.0;
	Tn[0][nx - 1] = Tg[0][nx - 1] = double(725.0 / 2);

	double errm = 0.0, err;
	fstream tfile;
	tfile.open("C:/Users/ACHINTYA/Desktop/Mini_Project/Explicit_time.dat", ios::out);
	ct = clock();
	int g = 0;
	//time loop here
	for (int t = 0;t < ts;t++)
	{
		for (int i = 0;i < ny;i++)
		{
			for (int j = 0;j < nx;j++)
			{
				To[i][j] = Tn[i][j];
			}
		}
		//interior nodes GS SOR
		do
		{
			err = 0.0;
			errm = 0.000001;
			for (int i = 1;i < ny - 1;i++)
			{
				for (int j = 1;j < nx - 1;j++)
				{
					Tn[i][j] = (To[i][j] * (1 - (4 * lambda))) + (lambda * (To[i][j + 1] + To[i + 1][j] + To[i - 1][j] + To[i][j - 1]));
					Tn[i][j] = (Tn[i][j] * w) + ((1 - w) * Tg[i][j]);
					err = abs(Tn[i][j] - Tg[i][j]);
					if (err > errm)
						errm = err;
					Tg[i][j] = Tn[i][j];
				}
			}

			g++;
		} while (errm >= 0.00001);


		tfile << t * dt << " " << Tn[(ny - 1) / 2][(nx - 1) / 2]<<endl;
		//end of time loop
	}

	tfile.close();
	fstream ofile;
	ofile.open("C:/Users/ACHINTYA/Desktop/Mini_Project/TrialE.csv", ios::out);

	for (int i = 0;i < ny;i++)
	{
		for (int j = 0;j < nx;j++)
		{
			ofile << Tn[i][j] << ",";
		}
		ofile << endl;
	}
	ct = clock() - ct;
	cout << setprecision(5)<<Tn[(ny - 1) / 2][(nx - 1) / 2] << endl;

	cout << "Time is " << double(ct / CLOCKS_PER_SEC) << " and iterations= " << g << endl;
	ofile.close();
	//putting in a dat file as well
	fstream op;
	op.open("C:/Users/ACHINTYA/Desktop/Mini_Project/TrialE.dat", ios::out);
	for (int j = 0;j < nx;j++)
	{
		for (int i = 0;i < ny;i++)
		{
			double x = j * dx;
			double y = i * dy;
			op << x << " " << y << " " << Tn[i][j];
			op << endl;
		}
		op << endl;
	}
	op.close();
}




