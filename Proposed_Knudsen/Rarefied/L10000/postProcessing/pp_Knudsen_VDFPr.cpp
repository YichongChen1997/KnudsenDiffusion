#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <cstdlib>

using namespace std;

vector<double> linspace(double min, double max, size_t N);

int main()
{
    string line;
    int id, n, i, t, nAtoms, currentTimeStep, typ, bo;

    double x, y, z, vx, vy, vz, tau1, tau2, tau3, tau4, tau5, tau6, PE, KE, fx, fy, fz;
    double xLo, xHi, yLo, yHi, zLo, zHi, Lz;

    double kB = 1.38064852e-23;
    double refMass = 1.66054e-27;
    double refLength = 1e-10;
    double refTime = 1e-15;
    double refVelocity = refLength / refTime; // 1 Angstrom / 1 fs
    double refPressure = 101325;

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

    double vM = sqrt(2 * kB * Tw / mGas); // most probable speed
    double binWidth = 10;                 // binwidth of velocity
    int maxVout = vM * 4;                 // max output velocity
    int nBins = ceil(maxVout / binWidth);
    cout << "The number of bins " << nBins << endl;
    // **************************

    ifstream data("dump_meas_gas.lammpstrj", ios::in);

    // **************************

    int nWallCollisions = 0;    // collect total number of collided molecules with wall
    int nStartedCollisions = 0; // collect number of collisions that have started
    vector<int> nCollisions;    // per molecule collect how many times molecule collided

    vector<bool> startFromMiddle;
    vector<bool> crossedStart;

    // keep track of molecule positions, velocities and time
    // only when they are in the middle of the domain (stop when molecules leave)
    vector<double> tM;       // time
    vector<double> xM;       // position x
    vector<double> yM;       // position y
    vector<double> zM;       // position z
    vector<double> vxM;      // velocity x
    vector<double> vyM;      // velocity y
    vector<double> vzM;      // velocity z
    vector<int> leftFromTop; // 1 = top; 0 = bottom

    // the following is collision data ("start" = when the molecule enters rCut and
    // "end" = when the molecule has left rCut). Every row is a single atomic collision.
    // Note that one molecule can hit a wall multiple times. so it can be in this list many times

    vector<double> tEnter;
    vector<double> xEnter;
    vector<double> yEnter;
    vector<double> zEnter;
    vector<double> vxEnter;
    vector<double> vyEnter;
    vector<double> vzEnter;

    vector<double> tStart;
    vector<double> xStart;
    vector<double> yStart;
    vector<double> zStart;
    vector<double> vxStart;
    vector<double> vyStart;
    vector<double> vzStart;

    vector<double> tEnd;
    vector<double> xEnd;
    vector<double> yEnd;
    vector<double> zEnd;
    vector<double> vxEnd;
    vector<double> vyEnd;
    vector<double> vzEnd;

    vector<int> leftFromTheTop;

    for (t = 0; t < nTimeSteps; t++)
    {
        for (n = 1; n < 10; n++)
        {
            if (n == 4)
            {
                data >> nAtoms;
                //                 cout << "nAtoms = " << nAtoms << endl;
            }

            if (n == 2)
            {
                data >> currentTimeStep;

                cout << "currentTimeStep = " << currentTimeStep
                     << "; t = " << t << " [ "
                     << 100 * float(t + 1) / float(nTimeSteps)
                     << "% ]" << endl;
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

        // set fields
        if (t == 0)
        {
            nCollisions.resize(nAtoms, 0);
            startFromMiddle.resize(nAtoms, false);
            crossedStart.resize(nAtoms, false);

            tM.resize(nAtoms, 0);
            xM.resize(nAtoms, 0.0);
            yM.resize(nAtoms, 0.0);
            zM.resize(nAtoms, 0.0);
            vxM.resize(nAtoms, 0.0);
            vyM.resize(nAtoms, 0.0);
            vzM.resize(nAtoms, 0.0);

            tEnter.resize(nAtoms, 0);
            xEnter.resize(nAtoms, 0.0);
            yEnter.resize(nAtoms, 0.0);
            zEnter.resize(nAtoms, 0.0);
            vxEnter.resize(nAtoms, 0.0);
            vyEnter.resize(nAtoms, 0.0);
            vzEnter.resize(nAtoms, 0.0);
        }

        for (n = 0; n < nAtoms; n++)
        {
            data >> id >> typ >> x >> y >> z >> vx >> vy >> vz >> tau1 >> tau2 >> tau3 >> KE >> PE >> fx >> fy >> fz;

            if (typ == 1) // methane
            {
                if ((z >= 0) && (z <= H))
                {
                    if ((x >= xLo + xDelta + xLf) && (x <= xLo + xDelta + xLf + xL)) // within the channel
                    {
                        if (t == 0)
                        {
                            if (z >= rCut)
                            {
                                startFromMiddle[n] = true; // inside middle region
                            }
                            tM[n] = t;

                            xM[n] = x;
                            yM[n] = y;
                            zM[n] = z;

                            vxM[n] = vx * refVelocity;
                            vyM[n] = vy * refVelocity;
                            vzM[n] = vz * refVelocity;
                        }
                        else
                        {
                            if (z < rCut)
                            {
                                if (startFromMiddle[n])
                                {
                                    // the collision has started
                                    crossedStart[n] = true;

                                    tEnter[n] = tM[n];

                                    xEnter[n] = xM[n];
                                    yEnter[n] = yM[n];
                                    zEnter[n] = zM[n];

                                    vxEnter[n] = vxM[n];
                                    vyEnter[n] = vyM[n];
                                    vzEnter[n] = vzM[n];

                                    nStartedCollisions++;
                                }
                                // switch into false for re-use this molecule
                                startFromMiddle[n] = false;
                            }

                            else if (z >= rCut)
                            {
                                if (crossedStart[n])
                                {
                                    // push data when the molecule started the collision
                                    tStart.push_back(tEnter[n]);

                                    xStart.push_back(xEnter[n]);
                                    yStart.push_back(yEnter[n]);
                                    zStart.push_back(zEnter[n]);

                                    vxStart.push_back(vxEnter[n]);
                                    vyStart.push_back(vyEnter[n]);
                                    vzStart.push_back(vzEnter[n]);

                                    /* if the molecule is within the middle region and has started the collision,
                                     flag it as the END */
                                    tEnd.push_back(t);

                                    xEnd.push_back(x);
                                    yEnd.push_back(y);
                                    zEnd.push_back(z);

                                    vxEnd.push_back(vx * refVelocity);
                                    vyEnd.push_back(vy * refVelocity);
                                    vzEnd.push_back(vz * refVelocity);

                                    // switch into false for re-use this molecule
                                    crossedStart[n] = false;

                                    nCollisions[n]++;
                                    nWallCollisions++;
                                }
                                // switch into true for re-use this molecule
                                startFromMiddle[n] = true;
                            }

                            tM[n] = t;

                            xM[n] = x;
                            yM[n] = y;
                            zM[n] = z;

                            vxM[n] = vx * refVelocity;
                            vyM[n] = vy * refVelocity;
                            vzM[n] = vz * refVelocity;
                        }
                    }
                }
            }
        }

        getline(data, line);

        cout << "number of gas-wall collisions = " << nWallCollisions << endl;
        cout << endl;
    }

    cout << "Analysis Starts: " << endl;
    cout << "Number of Collisions Started: " << nStartedCollisions << endl;
    cout << "Number of Collisions Ended: " << vzEnd.size() << endl;
    cout << endl;

    vector<double> vN, vNi;
    vector<double> vT, vTi;
    vector<double> vTx, vTxi;
    vector<double> vTy, vTyi;
    double vMagT, vMagTi;

    int coTop = 0, coBottom = 0;

    for (i = 0; i < nWallCollisions; i++)
    {
        // normal
        vN.push_back(abs(vzEnd[i]));
        vNi.push_back(abs(vzStart[i]));

        // combined tangents
        vMagT = vxEnd[i] * vxEnd[i] + vyEnd[i] * vyEnd[i];
        vMagTi = vxStart[i] * vxStart[i] + vyStart[i] * vyStart[i];

        if (vMagT > 0)
        {
            vT.push_back(sqrt(vMagT));
        }
        else
        {
            vT.push_back(0.0);
        }

        if (vMagTi > 0)
        {
            vTi.push_back(sqrt(vMagTi));
        }
        else
        {
            vTi.push_back(0.0);
        }

        vTx.push_back(vxEnd[i]);
        vTxi.push_back(vxStart[i]);

        vTy.push_back(vyEnd[i]);
        vTyi.push_back(vyStart[i]);

        coBottom++;
    }

    cout << "Size of vN: " << vN.size() << endl;
    cout << "No. of Collisions at top: " << coTop << endl;
    cout << "No. of Collisions at bottom: " << coBottom << endl;

    int nPts = vN.size();
    vector<double> range = linspace(0.10, 2.0, 20);

    // Partial reflected velocity distributions - Normal component
    vector<vector<double> > NormalDis;
    vector<vector<double> > NormalPro;
    ofstream normalfile("reflected_N.txt");

    for (int n = 0; n < range.size(); n++)
    {
        vector<double> distribution(nBins, 0.0);
        vector<double> bins(nBins, 0.0);
        int count = 0;

        for (i = 0; i < nBins; i++)
        {
            bins[i] = binWidth * 0.5 + binWidth * i;
        }

        for (i = 0; i < nPts; i++)
        {
            if ((vNi[i] > vM * (range[n] - 0.1)) && (vNi[i] <= vM * (range[n] + 0.1)))
            {
                bo = floor(vN[i] / binWidth);

                if (bo >= nBins)
                {
                    bo = nBins - 1;
                }

                distribution[bo] += 1.0;
                count++;
            }
        }
        cout << "No. molecules within Normal velocity range " << (range[n] - 0.1) << " and " << (range[n] + 0.1) << " are: " << count << endl;

        vector<double> probability(nBins, 0.0);
        double areaUnderGraph = 0.0;
        for (i = 0; i < nBins - 1; i++)
        {
            areaUnderGraph += (bins[i + 1] - bins[i]) * (distribution[i] + distribution[i + 1]) * 0.5;
        }

        for (i = 0; i < nBins; i++)
        {
            probability[i] = distribution[i] / areaUnderGraph;
        }

        NormalDis.push_back(distribution);
        NormalPro.push_back(probability);

        for (int i = 0; i < nBins; i++)
        {
            normalfile << bins[i] << " " << NormalDis[n][i] << " " << NormalPro[n][i] << endl;
        }
    }

    // Partial reflected velocity distributions - Tangential component
    vector<vector<double> > TangentialDisPo;
    vector<vector<double> > TangentialDisNe;
    vector<vector<double> > TangentialProPo;
    vector<vector<double> > TangentialProNe;
    ofstream tangentialfile("reflected_Tx.txt");

    for (int n = 0; n < range.size(); n++)
    {
        vector<double> distributionPos(nBins, 0.0);
        vector<double> distributionNeg(nBins, 0.0);
        vector<double> binsPos(nBins, 0.0);
        vector<double> binsNeg(nBins, 0.0);
        int count = 0;

        for (i = 0; i < nBins; i++)
        {
            binsPos[i] = binWidth * 0.5 + binWidth * i;
            binsNeg[i] = -binsPos[i];
        }

        for (i = 0; i < nPts; i++)
        {
            if ((vTxi[i] > vM * (range[n] - 0.1)) && (vTxi[i] <= vM * (range[n] + 0.1)))
            {
                if (vTx[i] >= 0)
                {
                    bo = floor(vTx[i] / binWidth);

                    if (bo >= nBins)
                    {
                        bo = nBins - 1;
                    }

                    distributionPos[bo] += 1.0;
                }
                else
                {
                    bo = floor(abs(vTx[i]) / binWidth);

                    if (bo >= nBins)
                    {
                        bo = nBins - 1;
                    }

                    distributionNeg[bo] += 1.0;
                }
                count++;
            }
        }

        cout << "No. molecules within Tangential velocity range " << (range[n] - 0.1) << " and " << range[n] + 0.1 << " are: " << count << endl;

        vector<double> probabilityPos(nBins, 0.0);
        vector<double> probabilityNeg(nBins, 0.0);

        // scale with area
        double areaUnderGraph = 0.0;

        for (i = 0; i < nBins - 1; i++)
        {
            areaUnderGraph += (binsPos[i + 1] - binsPos[i]) * (distributionPos[i] + distributionPos[i + 1]) * 0.5;
            areaUnderGraph += abs(binsNeg[i + 1] - binsNeg[i]) * (distributionNeg[i] + distributionNeg[i + 1]) * 0.5;
        }

        for (i = 0; i < nBins; i++)
        {
            probabilityPos[i] = distributionPos[i] / areaUnderGraph;
            probabilityNeg[i] = distributionNeg[i] / areaUnderGraph;
        }

        TangentialDisNe.push_back(distributionNeg);
        TangentialDisPo.push_back(distributionPos);
        TangentialProNe.push_back(probabilityNeg);
        TangentialProPo.push_back(probabilityPos);

        for (int i = 0; i < nBins; i++)
        {
            tangentialfile << binsNeg[nBins - i - 1] << " " << TangentialDisNe[n][nBins - i - 1] << " " << TangentialProNe[n][nBins - i - 1] << endl;
        }

        for (int i = 0; i < nBins; i++)
        {
            tangentialfile << binsPos[i] << " " << TangentialDisPo[n][i] << " " << TangentialProPo[n][i] << endl;
        }
    }

    return 0;
}

// Create a vector of evenly spaced numbers.
vector<double> linspace(double min, double max, size_t N)
{
    vector<double> linspace;
    double delta = (max - min) / double(N - 1);
    for (int i = 0; i < N; i++)
    {
        linspace.push_back(min + i * delta);
    }
    return linspace;
}
