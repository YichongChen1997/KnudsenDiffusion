#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <time.h>
#include <cstdlib>

using namespace std;
const float PI = 3.141592653589793238463;

double fRand(double fMin, double fMax);

// DATA: https://ww2.chemistry.gatech.edu/~lw26/structure/small_molecules/index.html

int main()
{
    string line;
    double x, y, z, xLo, xHi, yLo, yHi, zLo, zHi;
    double rCut, thickness, H, mG, sigmaG, mW, rho;
    int atoms, typ, i, j, k, m, n, id, mol = 0, molId = 0, nWallCount = 0, co;
    int nAtomsTot = 0;

    // reference scales (real to SI unit)
    double refLength = 1e-10;
    double refTime = 1e-15;
    double refVelocity = refLength / refTime;
    double kB = 1.38064852e-23;
    double zRigidRange = 3;   // rigid bound of the wall
    double zDynamicRange = 1; // dynamic range of the wall

    //***********************************
    int Nlines = 7; // find out the No. of lines of parameters    "wc -l < Parameters.dat"
    ifstream Parameters("Parameters.dat", ios::in);
    for (int n = 1; n < Nlines + 1; n++)
    {
        if (n == 1)
        {
            Parameters >> xLo >> xHi;
            getline(Parameters, line);
        }
        if (n == 2)
        {
            Parameters >> yLo >> yHi;
            getline(Parameters, line);
        }
        if (n == 3)
        {
            Parameters >> H >> rCut >> thickness;
            getline(Parameters, line);
        }
        if (n == 4)
        {
            Parameters >> rho;
            getline(Parameters, line);
        }
        if (n == 5)
        {
            Parameters >> sigmaG;
            getline(Parameters, line);
        }
        if (n == 6)
        {
            Parameters >> mW;
            getline(Parameters, line);
        }
        if (n == 7)
        {
            Parameters >> mG;
            getline(Parameters, line);
        }
    }

    srand(time(NULL));

    double Lx = xHi - xLo;
    double Ly = yHi - yLo;
    zLo = -thickness;
    zHi = H + thickness;

    double Lsigma = 3.4; // distance between different sheet
    int nSheet = thickness / Lsigma + 1;
    cout << "No. of graphene sheet = " << nSheet << endl;

    double d = 1.42; // bond length Angstroms
    double theta = 120 * PI / 180.0;
    double dx = d * sin(0.5 * theta);
    double dy = d * cos(0.5 * theta);

    cout << "initial Lx = " << Lx << endl;
    cout << "initial Ly = " << Ly << endl;

    int nX = int(Lx / (2 * dx));
    Lx = double(nX) * 2 * dx;

    int nY = int(Ly / (2 * dy + 2 * d));
    Ly = double(nY) * (2 * dy + 2 * d);

    cout << "nX = " << nX << endl;
    cout << "nY = " << nY << endl;

    cout << "new Lx = " << Lx << endl;
    cout << "new Ly = " << Ly << endl;

    xHi = Lx;
    yHi = Ly;

    int nWall = nSheet * 2 * nX * nY * 4;
    nAtomsTot += nWall;

    // create gas molecules
    double volume = (xHi - xLo) * (yHi - yLo) * H * refLength * refLength * refLength;
    int nGasPr = rho * volume;

    double offsetX = 3;
    double offsetY = 3;
    double offsetZ = 5;
    int nPx, nPy, nPz, nP;

    vector<double> xG(nGasPr);
    vector<double> yG(nGasPr);
    vector<double> zG(nGasPr);

    // Define molecular box
    double rhoG = nGasPr / (Lx * Ly * H);

    int nBx = int(pow((rhoG * Lx * Lx * Lx), 1.0 / 3));
    int nBy = int(pow((rhoG * Ly * Ly * Ly), 1.0 / 3));
    int nBz = int(pow((rhoG * H * H * H), 1.0 / 3));

    cout << "------------------------------------" << endl;
    cout << "nBx = " << nBx << endl;
    cout << "nBy = " << nBy << endl;
    cout << "nBz = " << nBz << endl;

    nPx = int(((Lx - offsetX) * nBx / Lx) + 1);
    nPy = int(((Ly - offsetY) * nBy / Ly) + 1);
    nPz = int(((H  - offsetZ) * nBz / H) + 1);
    
    nP = 0;
    bool search = true;
    while (search)
    {
        nP = nPx * nPy * nPz;
        if (nP > nGasPr)
        {
            search = false;
        }
        else
        {
            nPx++;
            nPz++;
            nPy++;
        }
    }

    double uCX = (Lx - offsetX) / double(nPx);
    double uCY = (Ly - offsetY) / double(nPy);
    double uCZ = (H  - 2 * offsetZ) / double(nPz);

    cout << "nPx = " << nPx << endl
         << "nPy = " << nPy << endl
         << "nPz = " << nPz << endl;

    cout << "uCX = " << uCX << endl
         << "uCY = " << uCY << endl
         << "uCZ = " << uCZ << endl;

    co = 0;

    for (j = 0; j < nPy; j++)
    {
        for (k = 0; k < nPz; k++)
        {
            for (i = 0; i < nPx; i++)
            {
                if (co < nGasPr)
                {
                    xG[co] = xLo + i * uCX + offsetX;
                    yG[co] = yLo + j * uCY + offsetY;
                    zG[co] = 0.0 + k * uCZ + offsetZ;
                    co++;
                }
            }
        }
    }

    nAtomsTot += co;

    cout << "------------------------------------" << endl
         << "nGas  mols  to insert = " << nGasPr << endl
         << "nWall atoms to insert = " << nWall << endl
         << "Total atoms to insert = " << nAtomsTot << endl
         << "------------------------------------" << endl;

    ofstream writer("data.dat");

    writer << "LAMMPS data file" << endl;
    writer << endl;
    writer << nAtomsTot << '\t' << "atoms" << endl;
    writer << endl;
    writer << "3	atom types" << endl;
    writer << endl;
    writer << xLo << '\t' << xHi << "	xlo xhi" << endl;
    writer << yLo << '\t' << yHi << "	ylo yhi" << endl;
    writer << zLo << '\t' << zHi << "	zlo zhi" << endl;

    writer << endl;
    writer << "Masses" << endl;
    writer << endl;
    writer << "1 " << mG << endl; // gas
    writer << "2 " << mW << endl; // wall bottom
    writer << "3 " << mW << endl; // wall top 
    writer << endl;
    writer << "Atoms" << endl;
    writer << endl;

    typ = 1;
    for (i = 0; i < nGasPr; i++)
    {
        molId++;
        mol++;
        writer << mol << '\t' << molId << '\t' << typ << '\t' << 0.0 << '\t' << xG[i] << '\t' << yG[i] << '\t' << zG[i] << endl;
    }

    cout << "Number of GasPr atoms placed = " << mol << endl;
    double refNdenstiy = mol / (Lx * Ly * H * pow(refLength, 3));

    int nCountWallAtoms = 0;

    typ = 2;
    for (n = 0; n < nSheet; n++)
    {
        for (i = 0; i < nX; i++)
        {
            for (j = 0; j < nY; j++)
            {
                // if (n == nSheet -1)
                // {
                //     typ = 4;
                // }
                mol++;
                // First lattice point (0, 0, 0)
                x = 0.0 + 2 * dx * i;
                y = 0.0 + (2 * dy + 2 * d) * j;
                z = 0.0 - n * Lsigma;

                x += (n % 2) * dx;
                y += (n % 2) * 0.5 * dy;
                if (x > Lx)
                {
                    x -= Lx;
                }
                if (y > Ly)
                {
                    y -= Ly;
                }
                writer << mol << '\t' << mol << '\t' << typ << '\t' << 0.0 << '\t' << x << '\t' << y << '\t' << z << endl;

                mol++;
                // Second lattice point (dx, dy, 0)
                x = dx + 2 * dx * i;
                y = dy + (2 * dy + 2 * d) * j;
                z = 0.0 - n * Lsigma;

                x += (n % 2) * dx;
                y += (n % 2) * 0.5 * dy;
                if (x > Lx)
                {
                    x -= Lx;
                }
                if (y > Ly)
                {
                    y -= Ly;
                }
                writer << mol << '\t' << mol << '\t' << typ << '\t' << 0.0 << '\t' << x << '\t' << y << '\t' << z << endl;

                mol++;
                // Third lattice point (dx, dy + d, 0)
                x = dx + 2 * dx * i;
                y = dy + d + (2 * dy + 2 * d) * j;
                z = 0.0 - n * Lsigma;

                x += (n % 2) * dx;
                y += (n % 2) * 0.5 * dy;
                if (x > Lx)
                {
                    x -= Lx;
                }
                if (y > Ly)
                {
                    y -= Ly;
                }
                writer << mol << '\t' << mol << '\t' << typ << '\t' << 0.0 << '\t' << x << '\t' << y << '\t' << z << endl;

                mol++;
                // Fourth lattice point (0, 2dy + d, 0)
                x = 0.0 + 2 * dx * i;
                y = 2 * dy + d + (2 * dy + 2 * d) * j;
                z = 0.0 - n * Lsigma;

                x += (n % 2) * dx;
                y += (n % 2) * 0.5 * dy;
                if (x > Lx)
                {
                    x -= Lx;
                }
                if (y > Ly)
                {
                    y -= Ly;
                }
                writer << mol << '\t' << mol << '\t' << typ << '\t' << 0.0 << '\t' << x << '\t' << y << '\t' << z << endl;

                nCountWallAtoms += 4;
            }
        }
    }

    typ = 3;
    for (n = 0; n < nSheet; n++)
    {
        for (i = 0; i < nX; i++)
        {
            for (j = 0; j < nY; j++)
            {
                // if (n == nSheet - 1)
                // {
                //     typ = 4;
                // }
                mol++;
                x = 0.0 + 2 * dx * i;
                y = 0.0 + (2 * dy + 2 * d) * j;
                z = H + n * Lsigma;

                x += (n % 2) * dx;
                y += (n % 2) * 0.5 * dy;
                if (x > Lx)
                {
                    x -= Lx;
                }
                if (y > Ly)
                {
                    y -= Ly;
                }
                writer << mol << '\t' << mol << '\t' << typ << '\t' << 0.0 << '\t' << x << '\t' << y << '\t' << z << endl;

                mol++;
                x = dx + 2 * dx * i;
                y = dy + (2 * dy + 2 * d) * j;
                z = H + n * Lsigma;

                x += (n % 2) * dx;
                y += (n % 2) * 0.5 * dy;
                if (x > Lx)
                {
                    x -= Lx;
                }
                if (y > Ly)
                {
                    y -= Ly;
                }
                writer << mol << '\t' << mol << '\t' << typ << '\t' << 0.0 << '\t' << x << '\t' << y << '\t' << z << endl;

                mol++;
                x = dx + 2 * dx * i;
                y = dy + d + (2 * dy + 2 * d) * j;
                z = H + n * Lsigma;

                x += (n % 2) * dx;
                y += (n % 2) * 0.5 * dy;
                if (x > Lx)
                {
                    x -= Lx;
                }
                if (y > Ly)
                {
                    y -= Ly;
                }
                writer << mol << '\t' << mol << '\t' << typ << '\t' << 0.0 << '\t' << x << '\t' << y << '\t' << z << endl;

                mol++;
                x = 0.0 + 2 * dx * i;
                y = 2 * dy + d + (2 * dy + 2 * d) * j;
                z = H + n * Lsigma;

                x += (n % 2) * dx;
                y += (n % 2) * 0.5 * dy;
                if (x > Lx)
                {
                    x -= Lx;
                }
                if (y > Ly)
                {
                    y -= Ly;
                }
                writer << mol << '\t' << mol << '\t' << typ << '\t' << 0.0 << '\t' << x << '\t' << y << '\t' << z << endl;

                nCountWallAtoms += 4;
            }
        }
    }
    cout << "Number of total atoms placed = " << mol << endl
         << "------------------------------------" << endl;

    writer << endl;
    writer.close();

    double redNdensity = refNdenstiy * (PI * pow(sigmaG, 3) * pow(refLength, 3) / 6);
    double confRatio = H / sigmaG;
    double lambda = 1 / (sqrt(2) * PI * pow(sigmaG, 2) * pow(refLength, 2) * refNdenstiy);
    double Kn = lambda / (H * refLength);

    cout << "Reference number density = " << refNdenstiy << endl;
    cout << "Reduced number density   = " << redNdensity << endl;
    cout << "Confinemen ratio R       = " << confRatio << endl;
    cout << "Mean free path           = " << lambda << endl;
    cout << "Knudsen number (Kn)      = " << Kn << endl;
    cout << "------------------------------------" << endl;

    return 0;
}

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}
