#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <time.h>
#include <cstdlib>

using namespace std;
const float PI = 3.1415;

//DATA: https://ww2.chemistry.gatech.edu/~lw26/structure/small_molecules/index.html

int main()
{
    string line;
    double mGas, xLo, xHi, yLo, yHi, zLo, zHi, rho;
    double sigmaG;
    int typ, i, j, k, id, mol = 0, molId = 0, co;
    int nAtomsTot = 0;

    // reference scales (real to SI unit)
    double refLength = 1e-10;
    double refForce = 6.9477E-11;
    double refTime = 1e-15;
    double refVelocity = refLength / refTime;
    double kB = 1.38064852e-23;

    int Nlines = 6;           // find out the No. of lines of parameters    "wc -l < Parameters.dat"
    ifstream Parameters("Parameters.dat", ios::in);
    for (int n = 1; n < Nlines + 1; n++)
    {
        if (n == 1)  {Parameters >> xLo >> xHi;    getline(Parameters, line);}
        if (n == 2)  {Parameters >> yLo >> yHi;    getline(Parameters, line);}
        if (n == 3)  {Parameters >> zLo >> zHi;    getline(Parameters, line);}
        if (n == 4)  {Parameters >> rho;           getline(Parameters, line);}
        if (n == 5)  {Parameters >> sigmaG;        getline(Parameters, line);}
        if (n == 6)  {Parameters >> mGas;          getline(Parameters, line);}
    }

    srand(time(NULL));

    double Lx = xHi - xLo;
    double Ly = yHi - yLo;
    double Lz = zHi - zLo;

    double H = zHi - zLo;

    // create gas molecules
    double volume = (xHi - xLo) * (yHi - yLo) * H * refLength * refLength * refLength;
    int nGasPr = rho * volume;
    
    cout  << "------------------------------------" << endl;
    cout  << "Number of GasPr atoms to insert = " << nGasPr << endl
          << "------------------------------------" << endl;

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
    
    cout << "------------------------------------" << endl
         << "Number of Gas molecules added = " << co << endl;
         
    nAtomsTot += co;


    ofstream writer("data.dat");

    writer << "LAMMPS data file" << endl;
    writer << endl;
    writer << nAtomsTot << '\t' << "atoms" << endl;
    writer << endl;
    writer << "1	atom types" << endl;
    writer << endl;
    writer << xLo << '\t' << xHi << "	xlo xhi" << endl;
    writer << yLo << '\t' << yHi << "	ylo yhi" << endl;
    writer << zLo << '\t' << zHi << "	zlo zhi" << endl;

    writer << endl;
    writer << "Masses" << endl;
    writer << endl;
    writer << "1  "<<  mGas << endl;
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

    cout << "Number of total atoms placed = " << mol << endl;
    double refNdenstiy = mol / (Lx * Ly * H * pow(refLength, 3));

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
    cout << "------------------------------------"  << endl;

    return 0;
}