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
    double zbinWidth = 0.5;   // coarse, across the channel height
    double zbinWidth2 = 0.25; // fine

    double xbinWidth = 50; // coarse, across the channel length
    double xbinWidth2 = 25;

    int nzBins = 0, nzBins2 = 0;
    int nxBins = 0, nxBins2 = 0;
    //**************************

    ifstream data("dump_meas_gas.lammpstrj", ios::in);

    vector<vector<double> > molFieldTemp;
    vector<vector<double> > velFieldTemp;
    vector<vector<double> > keFieldTemp;
    vector<vector<double> > molField;
    vector<vector<double> > keField;
    vector<vector<double> > velField;
    vector<vector<double> > stressField;

    vector<vector<double> > molFieldTemp2;
    vector<vector<double> > velFieldTemp2;
    vector<vector<double> > keFieldTemp2;
    vector<vector<double> > molField2;
    vector<vector<double> > keField2;
    vector<vector<double> > velField2;
    vector<vector<double> > stressField2;

    int count = 0;
    double nAvSteps = 0.0;

    int bz, bz2;
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
            Lx = xL; // average the channel region
            Ly = yHi - yLo;
            Lz = zHi - zLo;

            nxBins = ceil(Lx / xbinWidth);
            nzBins = ceil(Lz / zbinWidth);

            int halfnzBins = ceil(nzBins / 2);
            nzBins = 2 * halfnzBins;
            zbinWidth = Lz / double(nzBins);

            int halfnxBins = ceil(nxBins / 2);
            nxBins = 2 * halfnxBins;
            xbinWidth = Lx / double(nxBins);

            cout << "x binwidth update = " << xbinWidth << endl;
            cout << "z binwidth update = " << zbinWidth << endl;
            cout << "nxBins = " << nxBins << endl;
            cout << "nzBins = " << nzBins << endl;

            molFieldTemp.resize(nxBins, vector<double>(nzBins, 0.0));
            velFieldTemp.resize(nxBins, vector<double>(nzBins, 0.0));
            keFieldTemp.resize(nxBins, vector<double>(nzBins, 0.0));

            molField.resize(nxBins, vector<double>(nzBins, 0.0));
            keField.resize(nxBins, vector<double>(nzBins, 0.0));
            velField.resize(nxBins, vector<double>(nzBins, 0.0));
            stressField.resize(nxBins, vector<double>(nzBins, 0.0));

            nxBins2 = ceil(Lx / xbinWidth2);
            nzBins2 = ceil(Lz / zbinWidth2);

            int halfnxBins2 = ceil(nxBins2 / 2);
            nxBins2 = 2 * halfnxBins2;
            xbinWidth2 = Lx / double(nxBins2);

            int halfnzBins2 = ceil(nzBins2 / 2);
            nzBins2 = 2 * halfnzBins2;
            zbinWidth2 = Lz / double(nzBins2);

            cout << "x binwidth2 update = " << xbinWidth2 << endl;
            cout << "z binwidth2 update = " << zbinWidth2 << endl;
            cout << "nxBins2 = " << nxBins2 << endl;
            cout << "nzBins2 = " << nzBins2 << endl;

            molFieldTemp2.resize(nxBins2, vector<double>(nzBins2, 0.0));
            velFieldTemp2.resize(nxBins2, vector<double>(nzBins2, 0.0));
            keFieldTemp2.resize(nxBins2, vector<double>(nzBins2, 0.0));

            molField2.resize(nxBins2, vector<double>(nzBins2, 0.0));
            keField2.resize(nxBins2, vector<double>(nzBins2, 0.0));
            velField2.resize(nxBins2, vector<double>(nzBins2, 0.0));
            stressField2.resize(nxBins2, vector<double>(nzBins2, 0.0));
        }

        // store the lagrangian properties first because we need to calculate velocity first
        vector<double> xR(nAtoms, 0.0);
        vector<double> yR(nAtoms, 0.0);
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
            yR[n] = y;
            zR[n] = z;
            vxR[n] = vx;
            vyR[n] = vy;
            vzR[n] = vz;
            tau1R[n] = tau1;
            tau2R[n] = tau2;
            tau3R[n] = tau3;
        }

        // reset fields before use
        for (i = 0; i < nxBins; i++)
        {
            for (j = 0; j < nzBins; j++)
            {
                molFieldTemp[i][j] = 0.0;
                velFieldTemp[i][j] = 0.0;
            }
        }

        for (i = 0; i < nxBins2; i++)
        {
            for (j = 0; j < nzBins2; j++)
            {
                molFieldTemp2[i][j] = 0.0;
                velFieldTemp2[i][j] = 0.0;
            }
        }

        // first calculate mean velocity in bin
        for (n = 0; n < nAtoms; n++)
        {
            x = xR[n];
            z = zR[n];
            vx = vxR[n];

            if (t >= skipTimeStep)
            {
                if ((x >= xLo + xDelta + xLf) && (x <= xLo + xDelta + xLf + xL)) // within the channel
                {
                    if ((z >= zLo) && (z <= zHi))
                    {
                        bx = floor((x - xLo - xDelta - xLf) / xbinWidth);
                        bz = floor((z - zLo) / zbinWidth);

                        if (bx < 0)
                        {
                            bx = 0;
                        }
                        if (bx >= nxBins)
                        {
                            bx = nxBins - 1;
                        }

                        if (bz < 0)
                        {
                            bz = 0;
                        }
                        if (bz >= nzBins)
                        {
                            bz = nzBins - 1;
                        }

                        molFieldTemp[bx][bz] += 1.0;
                        velFieldTemp[bx][bz] += vx;

                        // fine measurements
                        bx2 = floor((x - xLo - xDelta - xLf) / xbinWidth2);
                        bz2 = floor((z - zLo) / zbinWidth2);

                        if (bx2 < 0)
                        {
                            bx2 = 0;
                        }
                        if (bx2 >= nxBins2)
                        {
                            bx2 = nxBins2 - 1;
                        }

                        if (bz2 < 0)
                        {
                            bz2 = 0;
                        }
                        if (bz2 >= nzBins2)
                        {
                            bz2 = nzBins2 - 1;
                        }

                        molFieldTemp2[bx2][bz2] += 1.0;
                        velFieldTemp2[bx2][bz2] += vx;
                    }
                }
            }
        }

        for (i = 0; i < nxBins; i++)
        {
            for (j = 0; j < nzBins; j++)
            {
                if (molFieldTemp[i][j] > 0)
                {
                    velFieldTemp[i][j] /= molFieldTemp[i][j]; // mean velocity in bin
                }
            }
        }

        for (i = 0; i < nxBins2; i++)
        {
            for (j = 0; j < nzBins2; j++)
            {
                if (molFieldTemp2[i][j] > 0)
                {
                    velFieldTemp2[i][j] /= molFieldTemp2[i][j];
                }
            }
        }

        // now calculate other properties
        for (n = 0; n < nAtoms; n++)
        {
            if (t >= skipTimeStep)
            {
                x = xR[n];
                y = yR[n];
                z = zR[n];
                vx = vxR[n];
                vy = vyR[n];
                vz = vzR[n];
                tau1 = tau1R[n];
                tau2 = tau2R[n];
                tau3 = tau3R[n];

                if ((z >= zLo) && (z <= zHi))
                {
                    if ((x >= xLo + xDelta + xLf) && (x <= xLo + xDelta + xLf + xL)) // within the channel
                    {
                        bx = floor((x - xLo - xDelta - xLf) / xbinWidth);
                        bz = floor((z - zLo) / zbinWidth);

                        if (bx < 0)
                        {
                            bx = 0;
                        }
                        if (bx >= nxBins)
                        {
                            bx = nxBins - 1;
                        }

                        if (bz < 0)
                        {
                            bz = 0;
                        }
                        if (bz >= nzBins)
                        {
                            bz = nzBins - 1;
                        }

                        molField[bx][bz] += 1.0;
                        keField[bx][bz] += mGas * (vx * vx + vy * vy + vz * vz - velFieldTemp[bx][bz] * velFieldTemp[bx][bz]);

                        velField[bx][bz] += vx;
                        stressField[bx][bz] += (tau1 + tau2 + tau3);

                        // fine measurements
                        bx2 = floor((x - xLo - xDelta - xLf) / xbinWidth2);
                        bz2 = floor((z - zLo) / zbinWidth2);

                        if (bx2 < 0)
                        {
                            bx2 = 0;
                        }
                        if (bx2 >= nxBins2)
                        {
                            bx2 = nxBins2 - 1;
                        }

                        if (bz2 < 0)
                        {
                            bz2 = 0;
                        }
                        if (bz2 >= nzBins2)
                        {
                            bz2 = nzBins2 - 1;
                        }

                        molField2[bx2][bz2] += 1.0;
                        keField2[bx2][bz2] += mGas * (vx * vx + vy * vy + vz * vz - velFieldTemp2[bx2][bz2] * velFieldTemp2[bx2][bz2]);

                        velField2[bx2][bz2] += vx;
                        stressField2[bx2][bz2] += (tau1 + tau2 + tau3);
                    }
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
        ofstream binFile1("doubleBins_channel_coarse.txt", ios::out);

        vector<vector<double> > bin;
        vector<vector<double> > rho;
        vector<vector<double> > vel;
        vector<vector<double> > temp;
        vector<vector<double> > pressure;

        bin.resize(nxBins, vector<double>(nzBins, 0.0));
        rho.resize(nxBins, vector<double>(nzBins, 0.0));
        vel.resize(nxBins, vector<double>(nzBins, 0.0));
        temp.resize(nxBins, vector<double>(nzBins, 0.0));
        pressure.resize(nxBins, vector<double>(nzBins, 0.0));

        double binVol;

        for (i = 0; i < nxBins; i++)
        {
            binVol = xbinWidth * zbinWidth * Ly * refVolume;
            for (j = 0; j < nzBins; j++)
            {
                bin[i][j] = zLo + zbinWidth * 0.5 + zbinWidth * j;       // coordinate of bin
                rho[i][j] = molField[i][j] * mGas / (binVol * nAvSteps); // average density of each bin over time

                vel[i][j] = 0.0;
                temp[i][j] = 0.0;

                if (molField[i][j] > 0.0)
                {
                    vel[i][j] = velField[i][j] * refVelocity / molField[i][j];                            // average velocity of each bin over time
                    temp[i][j] = keField[i][j] * refVelocity * refVelocity / (3.0 * kB * molField[i][j]); // average temperature
                }

                pressure[i][j] = -(stressField[i][j] / (3.0 * binVol * nAvSteps / refVolume)) * refPressure / 1e6;

                binFile1 << bin[i][j] << '\t'
                         << rho[i][j] << '\t'
                         << molField[i][j] / (binVol * nAvSteps) << '\t' // number density
                         << vel[i][j] << '\t'
                         << temp[i][j] << '\t'
                         << pressure[i][j] << '\t'
                         << endl;
            }
        }
    }

    {
        ofstream binFile2("doubleBins_channel_fine.txt", ios::out);

        vector<vector<double> > bin2;
        vector<vector<double> > rho2;
        vector<vector<double> > vel2;
        vector<vector<double> > temp2;
        vector<vector<double> > pressure2;

        bin2.resize(nxBins2, vector<double>(nzBins2, 0.0));
        rho2.resize(nxBins2, vector<double>(nzBins2, 0.0));
        vel2.resize(nxBins2, vector<double>(nzBins2, 0.0));
        temp2.resize(nxBins2, vector<double>(nzBins2, 0.0));
        pressure2.resize(nxBins2, vector<double>(nzBins2, 0.0));

        double binVol2;

        for (i = 0; i < nxBins2; i++)
        {
            binVol2 = xbinWidth2 * zbinWidth2 * Ly * refVolume;

            for (j = 0; j < nzBins2; j++)
            {
                bin2[i][j] = zLo + zbinWidth2 * 0.5 + zbinWidth2 * j;       // coordinate of bin
                rho2[i][j] = molField2[i][j] * mGas / (binVol2 * nAvSteps); // average density of each bin over time

                vel2[i][j] = 0.0;
                temp2[i][j] = 0.0;

                if (molField2[i][j] > 0.0)
                {
                    vel2[i][j] = velField2[i][j] * refVelocity / molField2[i][j];                            // average velocity of each bin over time
                    temp2[i][j] = keField2[i][j] * refVelocity * refVelocity / (3.0 * kB * molField2[i][j]); // average temperature
                }

                pressure2[i][j] = -(stressField2[i][j] / (3.0 * binVol2 * nAvSteps / refVolume)) * refPressure / 1e6;

                binFile2 << bin2[i][j] << '\t'
                         << rho2[i][j] << '\t'
                         << molField2[i][j] / (binVol2 * nAvSteps) << '\t' // number density
                         << vel2[i][j] << '\t'
                         << temp2[i][j] << '\t'
                         << pressure2[i][j] << '\t'
                         << endl;
            }
        }
    }
    return 0;
}