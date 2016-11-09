#include "MOPSO.h"

int main()
{
    createNewFold();
    time_t It1, It2;
    It1 = time(NULL);
    srand(time(NULL));
    // Input Data
    Particle particle[nPop];
    char *seq;
    seq = inputSeq();                                               //input seq
    for (int i = 0; i < nPop; i++) { inputParticle(particle[i], i, seq); }  // input the particles data
    decideDominated(particle);                        //  decide if particle[i] is dominated
    printParticleCost(particle, 0);
    //printParticleCostForDebug(particle, 0);
    myRep rep;
    putNewParticleIntoRep(particle, rep, 0);          //  put particles into rep
    initializeParticle(particle);                 //  init the velocity, bestCost, bestPosition



    // MOPSO Main Loop
    for (int it = 0; it < MaxIt; it++) {
        for (int i = 0; i < nPop; i++) {
            double w = getInertiaWeight(it, MaxIt);       //  get inertia wight
            PSOAdaption(particle[i], rep, w, it);            //  apply the PSO formula
            convertRotationToCoordinary(particle[i]);    //  convert rotation to coordinary
        }
        decideDominated(particle);                     // decide new domination relationship
        printParticleCost(particle, it + 1);
        //printParticleCostForDebug(particle, it + 1);
        if (it == 0) {
            rep.clear();
        }

        putNewParticleIntoRep(particle, rep, it);        //    add the new particles into rep
        updatePBest(particle);                       //      update pBest
    }

    printAnswer(rep);
    //printAnswerForDebug(rep);
    delete (seq);
    It2 = time(NULL);
    printTime(0, It1, It2);
    //printTimeForDebug(0, It1, It2);
    return 0;
}