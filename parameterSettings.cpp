#include "MOPSO.h"

//  the address of the input data
const char rootAddress[] = "/home/ws/zzZyj/MOPSO/";
//const char inputAddress[]="/home/ws/zzZyj/MOPSO/data/input/test_thread";
const char *inputAddress      = catStrStr(rootAddress, "data/input/newInput/");
const char *energyFileAddress = catStrStr(rootAddress, "data/energyFile/");
const char *answerAddress     = catStrStr(rootAddress, "data/answer/newAnswer/");
const char *draftAddress      = catStrStr(rootAddress, "data/draft/");

const char *refine_1Address   = catStrStr(rootAddress, "refine1/");
const char *scoreAddress      = catStrStr(rootAddress, "rosetta_source/bin/score.linuxgccrelease");
const char *databaseAddress   = catStrStr(rootAddress, "rosetta_database/");
const char *mybinAddress      = catStrStr(rootAddress, "mybin/");

// Problem Definition
const int nVar = 183;
const int VarMin = -180;
const int VarMax = 180;
const int VarSize[] = {1, nVar};
const int VelMax = 5;             //  =20 without rama_map
const int tidSize = 10;
// MOPSO Settings
const int nPop = 5;                  // Population Size
const int nRep = 50;                // Repository Size
const int  MaxIt = 5;          // Maximum Number of Iterations

const double phi1 = 2.05;
const double phi2 = 2.05;
const double phi = phi1 + phi2;
const double chi = 2 / ( phi - 2 + sqrt( phi * phi - 4 * phi ) );    // 0.73

const double wMin = chi;                        // =chi  Inertia Weight
const double wMax = 1.2;
const double wDamp=1;                       //  Inertia Weight Damping Ratio
const double c1 = chi*phi1;                 //  Personal Learning Coefficient
const double c2 = chi*phi2;                 //  Global Learning Coefficient

const double Alpha = 0.1;       //Grid Inflation Parameter
const int nGrid = 10;               //Number of Grids per each Dimension
const int Beta = 4;                   //Leader Selection Pressure Parameter
const int Gamma = 2;             // Extra (to be deleted) Repository Member Selection Pressure

const int numObjective = 1; //  Multiple Obejectives settings

const double INF=100000000;
