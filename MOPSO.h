#ifndef MOPSO_H
#define MOPSO_H

#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <list>
#include <cstdlib>
#include <ctime>
#include <pthread.h>

using namespace std;

struct structBest;
struct Particle;
struct Rep;

typedef list<Rep> myRep;            //   a wonderful list, the data struct of the repositpry, which
typedef const char CONCH;
//store the best solutions that have been found (pateto front)


// fileDisposal
    // Input Data
char * inputSeq();
void inputParticle(Particle &particle, int k, char *seq);
    // a function that will open file and throw error if it can't open file
void openFile(const char fileDir[10000], ifstream &infile);
    // a function that will remove file and throw error if it can't remove file
void removeFile(const char * str);

    // for the object function
void printPdb(Particle * particle);
    // for debug
        // print rep in several files
void printRep(myRep & rep, int len);
        // print the particles' status in one file
void printParticleCost(Particle * particle, int iterator);
    // print answer rep
void printAnswer(myRep & rep);
    //  print run time
void printTime(int choice, time_t start, time_t finish);
    //  create a new fold for the answer
void createNewFold();

void printParticleCostForDebug(Particle * particle, int iterator);
void printAnswerForDebug(myRep & rep);
void printTimeForDebug(int choice, time_t start, time_t finish);


// initialize particle
void initializeParticle(Particle * particle);

//  MOPSOFunction
    //  apply the formulation of PSO
void PSOAdaption(Particle &particle, myRep & rep, double w,int it);

    //  decide if particle[i] is dominated, for all i
void decideDominated(Particle * particle);

    //  compare the new position with the pBest and update it
void updatePBest(Particle * particle);

    // calculate f(x)
void getAllParticleCostForTest(Particle * particle);
void getAllParticleCost(Particle * particle);

void *getAParticleCost(void* particle);

    //if the result is true, then cost1 is dominated by cost2
bool isDominated(double *cost1, double *cost2);

    //    add the new non-dominated particles into rep / eliminate some in the rep that is not so good now
void putNewParticleIntoRep(Particle * particle, myRep & rep, int iterator);

    //decide whether particle[i] can be put into rep or not, as well as eliminate the rep particles
    //which are dominated by particle[i]
    // if return true then can put
bool canPutIntoRep(double *Cost, myRep & rep);


//  MOPSOAidFunction
    //  cat c1-ss-c2  e.g. :  "mj" , 520 , "forever"  ->  mj520forever        used to conbine file name
char *catStrIntStr(const char *c1, int ss, const char *c2);
char *catStrIntStr(const char *c1, const char *c2, int ss, const char *c3);

    //  cat c1-c2
char *catStrStr(CONCH *c1, CONCH *c2);
char *catStrStr(CONCH *c1, CONCH *c2, CONCH *c3);
char *catStrStr(CONCH *c1, CONCH *c2, CONCH *c3, CONCH *c4);
char *catStrStr(CONCH *c1, CONCH *c2, CONCH *c3, CONCH *c4, CONCH *c5);
char *catStrStr(CONCH *c1, CONCH *c2, CONCH *c3, CONCH *c4, CONCH *c5, CONCH *c6);
char *catStrStr(CONCH *c1, CONCH *c2, CONCH *c3, CONCH *c4, CONCH *c5, CONCH *c6, CONCH *c7);

    //  copy p2 to p1
void cpyDoubleArray(double *& p1, double *& p2, int n);

    //  check if V is out of range
void checkV(double *& tmp, int n);

    //  check if position is out of range
void checkP(double *& tmp, double *& tmpV, int n);

    //  get the right rep_h for this update time
int rouletteWheel(myRep & rep);

    //  get the w
double getInertiaWeight(double it,double MaxIt);

    //  for test object function
double dis1(double x, double y, double z);
double dis2(double x, double y, double z);

    //  get the AA that pr represent
char * getRelativeAA(char pr);

    // get the cost for the object function 1 and 2
double getCost0(int i);
//double getCost1(int i);


//  convert rotation to coordinary
void convertRotationToCoordinary(Particle &particle);


//  Matrix Function
    //  ans = a*b
void matrixProduct(double a[4][4], double b[4][4], double ans[4][4]);

    //  set a[4][4] to be full of 0
void initialize_0(double a[4][4]);

    //  dan wei(chinese) matrix
void initialize_1(double a[4][4]);

    //  get the R in the paper
void getR(int it, double * a, double theta, double R[4][4]);

    //  get the new coordinate for itAtom
void getNewCoordinary(int itAtom, double R[4][4], double * a);

void cpyOriginTOAddO_origin(int it, double * x, double * xAddO);

void getNewCoordinaryForO(int it, double R[4][4], double * xAddO, int numAtom);

    //  newR = preR
void cpyMatrix(double preR[4][4], double newR[4][4]);

    //   dan wei hua (chinese)
void unitization(double v[3]);

    //set Q[3][3] the paper used
void getQ(double Q[4][4], double v[3], double theta);


//  const statement
extern const char rootAddress[];
extern const char *inputAddress;
extern const char *energyFileAddress;
extern const char *refine_1Address;
extern const char *answerAddress;
extern const char *draftAddress;
extern const char *scoreAddress;
extern const char *databaseAddress;
extern const char *mybinAddress;

extern const int nVar;
extern const int VarMin;
extern const int VarMax ;
extern const int VarSize[];
extern const int VelMax;             //  =20 without rama_map
extern const int tidSize;


// MOPSO Settings
extern const int nPop;                  // Population Size
extern const int nRep;                // Repository Size
extern const int  MaxIt;          // Maximum Number of Iterations

extern const double phi1;
extern const double phi2;
extern const double phi;
extern const double chi ;    // 0.73

extern const double wMin;                        //  Inertia Weight
extern const double wMax;
extern const double wDamp;                       //  Inertia Weight Damping Ratio
extern const double c1;                 //  Personal Learning Coefficient
extern const double c2;                 //  Global Learning Coefficient

extern const double Alpha;       //Grid Inflation Parameter
extern const int nGrid;               //Number of Grids per each Dimension
extern const int Beta;                   //Leader Selection Pressure Parameter
extern const int Gamma;             // Extra (to be deleted) Repository Member Selection Pressure

extern const int numObjective;      //  Multiple Objectives
extern const double INF;

//  struct definition
// pBest
struct structBest{
    double *Position, *Cost;
};

struct Particle{
    double * origin;                          //    the coordinate(x,y,z) of every atom (not include atom O),
    //  three indexes represent an Atom's position
    double * addO_origin;                   // the coordinate(x,y,z) of every atom (include the atom O)
    double * Position;                      //  the angles
    double * old_position;                  //  used in the calculation of angle change
    double * Velocity;                      //  will be used in the MOPSO
    double * Cost ;                         //  will be used int the MOPSO, store every f(x) the particle have
    int * GridIndex, * GridSubIndex;
    bool dominated;                     //   represent if the particle is dominated by some other,
    //  which means it is certainly not the best answer
    int numAA, numAtom;         //  the number of AA and the number of atoms
    int sizeOfOrigin, sizeOfAddO_origin, sizeOfPosition, sizeOfVelocity;    //  the size of arrays
    int index;
    char *seq;
    structBest Best;                    //  pBest
};

struct Rep{
    double *Cost;
    double *addO_origin;
    double *Position;
    int sizeOfAddO_origin;
    int iterator;
};


#endif // MOPSO_H
