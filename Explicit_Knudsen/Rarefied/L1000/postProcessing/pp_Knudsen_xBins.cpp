#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>

using namespace std;

int main()
{
    string line;
    int id, i, j, k, n, m, t, nAtoms, currentTimeStep, typ;
    double x, y, z, vx, vy, vz, tau1, tau2, tau3, tau4, tau5, tau6, PE, KE, fx, fy, fz;

    double xLo, xHi, yLo, yHi, zLo, zHi;
    double Lx, Ly, Lz, volume;

    double kB = 1.38064852e-23;
    double refMass = 1.66054e-27;
    double refLength = 1e-10;
    double refTime = 1e-15;
    double refVelocity = refLength / refTime; // 1 Angstrom / 1 fs
    double refPressure = 101325;
    double refVolume = refLength * refLength * refLength;

    //***************************
    // Specify the parameters for the post-processing
    int nTimeSteps;
    double deltaT, tSkip, skipTimeStep;
    double rCut, H;
    double xDelta, xLf, xL, xLv;
    double Tg, Tw;
    double mG, mW, mGas, mWall;

    int Nlines = 6;        // find out the No. of lines of parameters  "wc -l < Parameters.dat"
    ifstream Specification("Specification.dat", ios::in);
    if (!Specification.is_open())
    {
        cerr << "Error opening Specification.dat file" << endl;
        return 1;
    }
    for (int n = 1; n < Nlines + 1; n++)
    {
        if (n == 1)  {Specification >> nTimeSteps;                                 getline(Specification, line);}
        if (n == 2)  {Specification >> deltaT >> tSkip >> skipTimeStep;            getline(Specification, line);}
        if (n == 3)  {Specification >> xLo >> xDelta >> xLf >> xL;                 getline(Specification, line);}
        if (n == 4)  {Specification >> rCut >> H;                                  getline(Specification, line);}
        if (n == 5)  {Specification >> Tg >> Tw;                                   getline(Specification, line);}
        if (n == 6)  {Specification >> mG >> mW;                                   getline(Specification, line);}
    }
    mGas  = mG * refMass;
    mWall = mW * refMass;
    //***************************

    double binWidth = 0.5;   // coarse
    double binWidth2 = 0.25; // fine

    int nBins = 0, nBins2 = 0;
    //**************************

    ifstream data("dump_meas_gas.lammpstrj", ios::in);

    vector<double> molFieldTemp;
    vector<double> velFieldTemp;
    vector<double> keFieldTemp;
    vector<double> molField;
    vector<double> keField;
    vector<double> velField;
    vector<double> stressField;

    vector<double> molFieldTemp2;
    vector<double> velFieldTemp2;
    vector<double> keFieldTemp2;
    vector<double> molField2;
    vector<double> keField2;
    vector<double> velField2;
    vector<double> stressField2;

    int count = 0;
    double nAvSteps = 0.0;

    int bx, bx2;
    int dN = floor(nTimeSteps / 100);

    for (t = 0; t < nTimeSteps; t++)
    {
        for (n = 1; n < 10; n++)
        {
            if (n == 2)
            {
                data >> currentTimeStep;
                count++;
                if (count >= dN)
                {
                    cout << "currentTimeStep = " << currentTimeStep
                         << "; t = " << t << " [ "
                         << 100 * float(t + 1) / float(nTimeSteps)
                         << "% ]" << endl;

                    count = 0;
                }
            }
            if (n == 4)
            {
                data >> nAtoms;
                //                 cout << "nAtoms = " << nAtoms << endl;
            }
            if (n == 6)
            {
                data >> xLo >> xHi;
            }

            if (n == 7)
            {
                data >> yLo >> yHi;
            }

            if (n == 8)
            {
                data >> zLo >> zHi;
            }
            getline(data, line);
        }

        if (t == 0)
        {
            // I'm trying to make profiles symmetric
            Lx = xHi - xLo;               // average the entire domain
            Ly = yHi - yLo;
            Lz = H;                       // channel height
            nBins = ceil(Lx / binWidth);

            int halfnBins = ceil(nBins / 2);
            nBins = 2 * halfnBins;
            binWidth = Lx / double(nBins);

            cout << "binwidth update = " << binWidth << endl;
            cout << "nBins = " << nBins << endl;

            molFieldTemp.resize(nBins, 0.0);
            molField.resize(nBins, 0.0);
            stressField.resize(nBins, 0.0);

            nBins2 = ceil(Lx / binWidth2);
            int halfnBins2 = ceil(nBins2 / 2);
            nBins2 = 2 * halfnBins2;
            binWidth2 = Lx / double(nBins2);

            cout << "binwidth update = " << binWidth2 << endl;
            cout << "nBins = " << nBins2 << endl;

            molFieldTemp2.resize(nBins2, 0.0);
            molField2.resize(nBins2, 0.0);
            stressField2.resize(nBins2, 0.0);
        }

        // store the lagrangian properties first because we need to calculate velocity first
        vector<double> xR(nAtoms, 0.0);
        // vector<double> yR(nAtoms, 0.0);
        vector<double> zR(nAtoms, 0.0);
        vector<double> vxR(nAtoms, 0.0);
        vector<double> vyR(nAtoms, 0.0);
        vector<double> vzR(nAtoms, 0.0);
        vector<double> tau1R(nAtoms, 0.0);
        vector<double> tau2R(nAtoms, 0.0);
        vector<double> tau3R(nAtoms, 0.0);

        for (n = 0; n < nAtoms; n++)
        {
            data >> id >> typ >> x >> y >> z >> vx >> vy >> vz >> tau1 >> tau2 >> tau3 >> KE >> PE >> fx >> fy >> fz;
            // data >> id >> typ >> x >> y >> z >> vx >> vy >> vz;

            xR[n] = x;
            // yR[n]=y;
            zR[n] = z;
            vxR[n] = vx;
            vyR[n] = vy;
            vzR[n] = vz;
            tau1R[n] = tau1;
            tau2R[n] = tau2;
            tau3R[n] = tau3;
        }

        // now calculate other properties
        for (n = 0; n < nAtoms; n++)
        {
            if (t >= skipTimeStep)
            {
                x = xR[n];
                // y = yR[n];
                z = zR[n];
                vx = vxR[n];
                vy = vyR[n];
                vz = vzR[n];
                tau1 = tau1R[n];
                tau2 = tau2R[n];
                tau3 = tau3R[n];

                if ((x >= xLo) && (x <= xHi) && (z >= 0) && (z <= Lz))
                {
                        bx = floor((x - xLo) / binWidth);

                        if (bx < 0)
                        {
                            bx = 0;
                        }
                        if (bx >= nBins)
                        {
                            bx = nBins - 1;
                        }

                        molField[bx] += 1.0;
                        stressField[bx] += (tau1 + tau2 + tau3);

                        // fine measurements
                        bx2 = floor((x - xLo) / binWidth2);

                        if (bx2 < 0)
                        {
                            bx2 = 0;
                        }
                        if (bx2 >= nBins2)
                        {
                            bx2 = nBins2 - 1;
                        }

                        molField2[bx2] += 1.0;
                        stressField2[bx2] += (tau1 + tau2 + tau3);
                }
            }
        }

        if (t >= skipTimeStep)
        {
            nAvSteps += 1.0;
        }

        //****
        getline(data, line);
    }

    // bin measurements
    cout << "bin measurements" << endl;

    {
        ofstream binFile1("xbins_channel_coarse.txt", ios::out);

        double bin, rho = 0.0, vel = 0.0, temp = 0.0, binVol, binVol2, pressure;

        for (i = 0; i < nBins; i++)
        {
            bin = xLo + binWidth * 0.5 + binWidth * i; // coordinate of bin

            binVol = binWidth * Ly * Lz * refVolume; // volume of each bin

            rho = molField[i] * mGas / (binVol * nAvSteps); // average density of each bin over time

            pressure = -(stressField[i] / (3.0 * binVol * nAvSteps / refVolume)) * refPressure / 1e6;

            binFile1 << bin << '\t'
                     << rho << '\t'
                     << molField[i] / (binVol * nAvSteps) << '\t' // number density
                     << pressure << '\t'
                     << endl;
        }
    }

    {
        ofstream binFile2("xbins_channel_fine.txt", ios::out);

        double bin, rho = 0.0, vel = 0.0, temp = 0.0, binVol, binVol2, pressure;

        for (i = 0; i < nBins2; i++)
        {
            bin = xLo + binWidth2 * 0.5 + binWidth2 * i;

            binVol2 = binWidth2 * Ly * Lz * refVolume;

            rho = molField2[i] * mGas / (binVol2 * nAvSteps);

            pressure = -(stressField2[i] / (3.0 * binVol2 * nAvSteps / refVolume)) * refPressure / 1e6;

            binFile2 << bin << '\t'
                     << rho << '\t'
                     << molField2[i] / (binVol2 * nAvSteps) << '\t' // number density
                     << pressure << '\t'
                     << endl;
        }
    }

    return 0;
}