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

    double rCut, thickness, H, xLf, xL, xDelta, xLv, mGas, sigmaG, mW, rho;

    int atoms, typ, i, j, k, m, n, id, mol = 0, molId = 0, nWallCount = 0, co;
    int nAtomsTot = 0;
    int nSheet;

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
        if (n == 1)  {Parameters >> xLo >> xDelta >> xLf >> xL >> xLv;    getline(Parameters, line);}
        if (n == 2)  {Parameters >> yLo >> yHi;                           getline(Parameters, line);}
        if (n == 3)  {Parameters >> zLo >> zHi;                           getline(Parameters, line);}
        if (n == 4)  {Parameters >> H >> rCut;                            getline(Parameters, line);}
        if (n == 5)  {Parameters >> rho;                                  getline(Parameters, line);}
        if (n == 6)  {Parameters >> sigmaG >> mGas;                       getline(Parameters, line);}
        if (n == 7)  {Parameters >> mW;                                   getline(Parameters, line);}
    }

    srand(time(NULL));

    xHi = xLo + xDelta + xLf + xL + xLv;
    double Ly = yHi - yLo;

    zHi += rCut;
    zLo += rCut;
    double Lz = zHi - zLo;
    cout  << "------------------------------------" << endl;
    cout << "Adjusted zLo = " << zLo << endl;
    cout << "Adjusted zHi = " << zHi << endl; 

    /* 
    Parameters related to the setup of the wall 
    */
    double LZsigma = 3.4;   // distance between different sheet
    double epsilon = 1.8;   // distance for preventing gas leakage

    thickness = rCut;
    nSheet = ceil(thickness / LZsigma);

    double d = 1.42;        // bond length Angstroms
    double theta = 120 * PI / 180.0;
    double dx = d * sin(0.5 * theta);
    double dy = d * cos(0.5 * theta);

    cout  << "------------------------------------" << endl;
    cout << "initial xL = " << xL << endl;
    cout << "initial Ly = " << Ly << endl;

    int nX = int(xL / (2 * dx));
    xL = double(nX) * 2 * dx;

    int nY = int(Ly / (2 * dy + 2 * d));
    Ly = double(nY) * (2 * dy + 2 * d);

    cout << "nX = " << nX << endl;
    cout << "nY = " << nY << endl;

    cout << "new xL = " << xL << endl;
    cout << "new Ly = " << Ly << endl;

    // create gas molecules in gas reservoir
    double volume = xLf * (yHi - yLo) * (zHi - zLo) * refLength * refLength * refLength;

    // account for gas molecules within the channel
    double channelVolume = xL * (yHi - yLo) * H * refLength * refLength * refLength;

    int nGasPr = rho * volume + rho * channelVolume / 2;
    
    cout  << "------------------------------------" << endl;
    cout  << "Number of GasPr atoms to insert = " << nGasPr << endl;

    double offsetX = 3;
    double offsetY = 3;
    double offsetZ = 3;
    int nPx, nPy, nPz, nP;

    vector<double> xG(nGasPr);
    vector<double> yG(nGasPr);
    vector<double> zG(nGasPr);

    // Define molecular box
    double rhoG = nGasPr / (xLf * (yHi - yLo) * (zHi - zLo));

    // number of moleucles in each direction
    int nBx = int(pow((rhoG * xLf * xLf * xLf), 1.0 / 3));
    int nBy = int(pow((rhoG * (yHi - yLo)* (yHi - yLo)* (yHi - yLo)), 1.0 / 3));
    int nBz = int(pow((rhoG * (zHi - zLo)* (zHi - zLo)* (zHi - zLo)), 1.0 / 3));

    nPx = int(((xLf - offsetX) * nBx / xLf) + 1);
    nPy = int((((yHi - yLo) - offsetY) * nBy / (yHi - yLo)) + 1);
    nPz = int((((zHi - zLo) - offsetZ) * nBz / (zHi - zLo)) + 1);
    
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

    double uCX = (xLf - offsetX) / double(nPx);
    double uCY = ((yHi - yLo) - offsetY) / double(nPy);
    double uCZ = ((zHi - zLo) - offsetZ) / double(nPz);

    cout << "------------------------------------" << endl
         << "Initial number of atoms in each direction: " << endl 
         << "nPx = " << nPx << endl
         << "nPy = " << nPy << endl
         << "nPz = " << nPz << endl
         << "------------------------------------" << endl;

    cout << "------------------------------------" << endl
         << "Initial space interval in each direction: " << endl
         << "uCX = " << uCX << endl
         << "uCY = " << uCY << endl
         << "uCZ = " << uCZ << endl
         << "------------------------------------" << endl;

    co = 0;

    for (j = 0; j < nPy; j++)
    {
        for (k = 0; k < nPz; k++)
        {
            for (i = 0; i < nPx; i++)
            {
                if (co < nGasPr)
                {
                    xG[co] = xLo + i * uCX + offsetX + xDelta;
                    yG[co] = yLo + j * uCY + offsetY;
                    zG[co] = zLo + k * uCZ + offsetZ;
                    co++;
                }
            }
        }
    }
    
    cout << "Number of Gas molecules added = " << co << endl
         << "------------------------------------" << endl;
         
    nAtomsTot += co;

    int nWallAtoms = 0;

    double channelShift = xLf + xDelta;

    cout << "The value of channelShift = " << channelShift << endl;

    // Define z and x cutoff ranges to exclude inner atoms, keeping only near-surface atoms.
    double zMinInterface = - thickness;  // near lower surface
    double zMaxInterface = H   + thickness;  // near upper surface
    // double xMinInterface = channelShift + thickness;       // near left side wall
    // double xMaxInterface = channelShift + xL - thickness;  // near right side wall
    double xMinInterface = channelShift - epsilon;       // near left side wall
    double xMaxInterface = channelShift + xL + epsilon;  // near right side wall

    typ = 2;
    for (n = 0; n < nSheet; n++)
    {
        for (i = 0; i < nX; i++)
        {
            for (j = 0; j < nY; j++)
            {
                // First lattice point (0, 0, 0)
                x = 0.0 + 2 * dx * i + channelShift;
                y = 0.0 + (2 * dy + 2 * d) * j;
                z = 0.0 - n * LZsigma;

                x += (n % 2) * dx;
                y += (n % 2) * 0.5 * dy;
                if (x > xL + channelShift)
                {
                    x -= xL;
                }
                if (y > Ly)
                {
                    y -= Ly;
                }
                if (!(z < zMinInterface && x > xMinInterface && x < xMaxInterface)) {
                    nWallAtoms ++;
                    }

                // Second lattice point (dx, dy, 0)
                x = dx + 2 * dx * i + channelShift;
                y = dy + (2 * dy + 2 * d) * j;
                z = 0.0 - n * LZsigma;

                x += (n % 2) * dx;
                y += (n % 2) * 0.5 * dy;
                if (x > xL + channelShift)
                {
                    x -= xL;
                }
                if (y > Ly)
                {
                    y -= Ly;
                }
                if (!(z < zMinInterface && x > xMinInterface && x < xMaxInterface)) {
                    nWallAtoms ++;
                    }

                // Third lattice point (dx, dy + d, 0)
                x = dx + 2 * dx * i + channelShift;
                y = dy + d + (2 * dy + 2 * d) * j;
                z = 0.0 - n * LZsigma;

                x += (n % 2) * dx;
                y += (n % 2) * 0.5 * dy;
                if (x > xL + channelShift)
                {
                    x -= xL;
                }
                if (y > Ly)
                {
                    y -= Ly;
                }
                if (!(z < zMinInterface && x > xMinInterface && x < xMaxInterface)) {
                    nWallAtoms ++;
                    }

                // Fourth lattice point (0, 2dy + d, 0)
                x = 0.0 + 2 * dx * i + channelShift;
                y = 2 * dy + d + (2 * dy + 2 * d) * j;
                z = 0.0 - n * LZsigma;

                x += (n % 2) * dx;
                y += (n % 2) * 0.5 * dy;
                if (x > xL + channelShift)
                {
                    x -= xL;
                }
                if (y > Ly)
                {
                    y -= Ly;
                }
                if (!(z < zMinInterface && x > xMinInterface && x < xMaxInterface)) {
                    nWallAtoms ++;
                    }
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
                // First lattice point
                x = 0.0 + 2 * dx * i + channelShift;
                y = 0.0 + (2 * dy + 2 * d) * j;
                z = H + n * LZsigma;

                x += (n % 2) * dx;
                y += (n % 2) * 0.5 * dy;
                if (x > xL + channelShift)
                {
                    x -= xL;
                }
                if (y > Ly)
                {
                    y -= Ly;
                }
                if (!(z > zMaxInterface && x > xMinInterface && x < xMaxInterface)) {
                    nWallAtoms ++;
                    }

                // Second lattice point
                x = dx + 2 * dx * i + channelShift;
                y = dy + (2 * dy + 2 * d) * j;
                z = H + n * LZsigma;

                x += (n % 2) * dx;
                y += (n % 2) * 0.5 * dy;
                if (x > xL + channelShift)
                {
                    x -= xL;
                }
                if (y > Ly)
                {
                    y -= Ly;
                }
                if (!(z > zMaxInterface && x > xMinInterface && x < xMaxInterface)) {
                    nWallAtoms ++;
                    }


                // Third lattice point
                x = dx + 2 * dx * i + channelShift;
                y = dy + d + (2 * dy + 2 * d) * j;
                z = H + n * LZsigma;

                x += (n % 2) * dx;
                y += (n % 2) * 0.5 * dy;
                if (x > xL + channelShift)
                {
                    x -= xL;
                }
                if (y > Ly)
                {
                    y -= Ly;
                }
                if (!(z > zMaxInterface && x > xMinInterface && x < xMaxInterface)) {
                    nWallAtoms ++;
                    }


                // Fourth lattice point
                x = 0.0 + 2 * dx * i + channelShift;
                y = 2 * dy + d + (2 * dy + 2 * d) * j;
                z = H + n * LZsigma;

                x += (n % 2) * dx;
                y += (n % 2) * 0.5 * dy;
                if (x > xL + channelShift)
                {
                    x -= xL;
                }
                if (y > Ly)
                {
                    y -= Ly;
                }
                if (!(z > zMaxInterface && x > xMinInterface && x < xMaxInterface)) {
                    nWallAtoms ++;
                    }
            }
        }
    }
    

    nAtomsTot += nWallAtoms;

    cout << "------------------------------------" << endl
         << "nGas  mols  to insert = " << nGasPr << endl
         << "nWall atoms to insert = " << nWallAtoms << endl
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
    writer << "1 " << mGas << endl;  // gas
    writer << "2 " << mW << endl; // Bottom Graphene
    writer << "3 " << mW << endl; // Top    Graphene
    writer << endl;
    writer << "Atoms" << endl;
    writer << endl;

    typ = 1;
    for (i = 0; i < nGasPr; i++)
    {
        mol++;
        molId++;
        writer << mol << '\t' << molId << '\t' << typ << '\t' << 0.0 << '\t' << xG[i] << '\t' << yG[i] << '\t' << zG[i] << endl;
    }

    cout << "Number of GasPr atoms placed = " << mol << endl;
    double refNdenstiy = mol / (xLf * (yHi - yLo) * (zHi - zLo) * pow(refLength, 3));


    typ = 2;  // bottom walls
    for (n = 0; n < nSheet; n++)
    {
        for (i = 0; i < nX; i++)
        {
            for (j = 0; j < nY; j++)
            {
                // First lattice point (0, 0, 0)
                x = 0.0 + 2 * dx * i + channelShift;
                y = 0.0 + (2 * dy + 2 * d) * j;
                z = 0.0 - n * LZsigma;

                x += (n % 2) * dx;
                y += (n % 2) * 0.5 * dy;
                if (x > xL + channelShift)
                {
                    x -= xL;
                }
                if (y > Ly)
                {
                    y -= Ly;
                }
                if (!(z < zMinInterface && x > xMinInterface && x < xMaxInterface)) {
                    mol++;
                    molId++;
                    writer << mol << '\t' << molId << '\t' << typ << '\t' << 0.0 << '\t' << x << '\t' << y << '\t' << z << endl;
                    }

                // Second lattice point (dx, dy, 0)
                x = dx + 2 * dx * i + channelShift;
                y = dy + (2 * dy + 2 * d) * j;
                z = 0.0 - n * LZsigma;

                x += (n % 2) * dx;
                y += (n % 2) * 0.5 * dy;
                if (x > xL + channelShift)
                {
                    x -= xL;
                }
                if (y > Ly)
                {
                    y -= Ly;
                }
                if (!(z < zMinInterface && x > xMinInterface && x < xMaxInterface)) {
                    mol++;
                    molId++;
                    writer << mol << '\t' << molId << '\t' << typ << '\t' << 0.0 << '\t' << x << '\t' << y << '\t' << z << endl;
                    }


                // Third lattice point (dx, dy + d, 0)
                x = dx + 2 * dx * i + channelShift;
                y = dy + d + (2 * dy + 2 * d) * j;
                z = 0.0 - n * LZsigma;

                x += (n % 2) * dx;
                y += (n % 2) * 0.5 * dy;
                if (x > xL + channelShift)
                {
                    x -= xL;
                }
                if (y > Ly)
                {
                    y -= Ly;
                }
                if (!(z < zMinInterface && x > xMinInterface && x < xMaxInterface)) {
                    mol++;
                    molId++;
                    writer << mol << '\t' << molId << '\t' << typ << '\t' << 0.0 << '\t' << x << '\t' << y << '\t' << z << endl;
                    }

                // Fourth lattice point (0, 2dy + d, 0)
                x = 0.0 + 2 * dx * i + channelShift;
                y = 2 * dy + d + (2 * dy + 2 * d) * j;
                z = 0.0 - n * LZsigma;

                x += (n % 2) * dx;
                y += (n % 2) * 0.5 * dy;
                if (x > xL + channelShift)
                {
                    x -= xL;
                }
                if (y > Ly)
                {
                    y -= Ly;
                }
                if (!(z < zMinInterface && x > xMinInterface && x < xMaxInterface)) {
                    mol++;
                    molId++;
                    writer << mol << '\t' << molId << '\t' << typ << '\t' << 0.0 << '\t' << x << '\t' << y << '\t' << z << endl;
                    }
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
                // First lattice point
                x = 0.0 + 2 * dx * i + channelShift;
                y = 0.0 + (2 * dy + 2 * d) * j;
                z = H + n * LZsigma;

                x += (n % 2) * dx;
                y += (n % 2) * 0.5 * dy;
                if (x > xL + channelShift)
                {
                    x -= xL;
                }
                if (y > Ly)
                {
                    y -= Ly;
                }
                if (!(z > zMaxInterface && x > xMinInterface && x < xMaxInterface)) {
                    mol++;
                    molId++;
                    writer << mol << '\t' << molId << '\t' << typ << '\t' << 0.0 << '\t' << x << '\t' << y << '\t' << z << endl;
                    }


                // Second lattice point
                x = dx + 2 * dx * i + channelShift;
                y = dy + (2 * dy + 2 * d) * j;
                z = H + n * LZsigma;

                x += (n % 2) * dx;
                y += (n % 2) * 0.5 * dy;
                if (x > xL + channelShift)
                {
                    x -= xL;
                }
                if (y > Ly)
                {
                    y -= Ly;
                }
                if (!(z > zMaxInterface && x > xMinInterface && x < xMaxInterface)) {
                    mol++;
                    molId++;
                    writer << mol << '\t' << molId << '\t' << typ << '\t' << 0.0 << '\t' << x << '\t' << y << '\t' << z << endl;
                    }


                // Third lattice point
                x = dx + 2 * dx * i + channelShift;
                y = dy + d + (2 * dy + 2 * d) * j;
                z = H + n * LZsigma;

                x += (n % 2) * dx;
                y += (n % 2) * 0.5 * dy;
                if (x > xL + channelShift)
                {
                    x -= xL;
                }
                if (y > Ly)
                {
                    y -= Ly;
                }
                if (!(z > zMaxInterface && x > xMinInterface && x < xMaxInterface)) {
                    mol++;
                    molId++;
                    writer << mol << '\t' << molId << '\t' << typ << '\t' << 0.0 << '\t' << x << '\t' << y << '\t' << z << endl;
                    }


                // Fourth lattice point
                x = 0.0 + 2 * dx * i + channelShift;
                y = 2 * dy + d + (2 * dy + 2 * d) * j;
                z = H + n * LZsigma;

                x += (n % 2) * dx;
                y += (n % 2) * 0.5 * dy;
                if (x > xL + channelShift)
                {
                    x -= xL;
                }
                if (y > Ly)
                {
                    y -= Ly;
                }
                if (!(z > zMaxInterface && x > xMinInterface && x < xMaxInterface)) {
                    mol++;
                    molId++;
                    writer << mol << '\t' << molId << '\t' << typ << '\t' << 0.0 << '\t' << x << '\t' << y << '\t' << z << endl;
                    }
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