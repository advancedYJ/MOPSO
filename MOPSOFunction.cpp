#include "MOPSO.h"

void PSOAdaption(Particle & particle, myRep & rep, double w, int it) {
    int len = particle.sizeOfPosition;
    int h = rouletteWheel(rep);
    cpyDoubleArray(particle.old_position, particle.Position, len);
    myRep ::iterator a = rep.begin();
    for (int i=0; i<h; i++) a++;
    //  calculate the new Velocity
    for (int i=0; i<len; i++){
        double r1 = rand() / double(RAND_MAX) , r2 = rand() / double(RAND_MAX) ;
        //double r1=0.7, r2=0.7;
        if(it==0 && i%3 != 1) particle.Velocity[i] = rand() / double(RAND_MAX);
        //if(it==0 && i%3 != 1) particle.Velocity[i] = 0.7;
        particle.Velocity[i] = w * particle.Velocity[i] +
                               c1*r1* (particle.Best.Position[i] - particle.Position[i]) +
                               c2*r2 * (a->Position[i] - particle.Position[i]);
    }

    checkV(particle.Velocity, particle.sizeOfVelocity);

    //  calculate the new position
    for (int i=0; i<len; i++)   particle.Position[i] += particle.Velocity[i];

    checkP(particle.Position, particle.Velocity, particle.sizeOfPosition);
}

void getAllParticleCost(Particle * particle){
    pthread_t tids[tidSize];
    int haveRun=0;
    while (haveRun < nPop){
        for (int i=0; i<tidSize && haveRun+i<nPop; i++){
            pthread_create(&tids[i], NULL, getAParticleCost, (void *) &particle[haveRun+i]);
        }
        for (int i=0; i<tidSize && haveRun+i<nPop; i++){
            pthread_join(tids[i], NULL);
        }
        haveRun += tidSize;
    }

}

void *getAParticleCost(void *p){
    Particle *particle = (Particle *) p;
    printPdb(particle);

    particle->Cost[0] = getCost0(particle->index);

    //particle->Cost[1] = getCost1(particle->index);
}

/*void getAllParticleCostForTest(Particle * particle){                   // calculate f(x)
    for (int i=0; i<nPop; i++)
    {
        int len = particle[i].sizeOfAddO_origin;
        double ans1 = 0, ans2=0;
        for (int j=0; j<len; j+=3){
            double x= particle[i].addO_origin[j];
            double y= particle[i].addO_origin[j+1];
            double z= particle[i].addO_origin[j+2];
            ans1 += dis1(x,y,z);
            ans2 += dis2(x, y, z);
        }
        particle[i].Cost[0] = ans1;
        particle[i].Cost[1] = ans2;
    }
}
*/
bool isDominated(double *cost1, double *cost2){
    bool flag = 1;                                          //  if flag is true then  cost1 is dominated by cost2
    //  which means  cost1 >= cost2, for all the f(x)
    bool equal = 1;
    for (int i=0; i<numObjective; i++){
        if (cost1[i] < cost2[i]){   flag=0;     break;      }
        if (cost1[i] != cost2[i] )  equal = 0;
    }
    if (equal)  flag  =  0;                             //   if cost1 = cost2 ,for all the f(x), then cost1 is not
    return flag;                                          //dominated by cost2, which means 2 is not better than 1
}

void decideDominated(Particle * particle){
    getAllParticleCost(particle);
    //getAllParticleCostForTest(particle);
    for (int i=0; i<nPop; i++)
        particle[i].dominated = 0;
    for (int i=0; i<nPop; i++)
        for (int j=0; j<nPop; j++)
            if  (isDominated(particle[i].Cost, particle[j].Cost))
                particle[i].dominated = 1;
}

void putNewParticleIntoRep(Particle * particle, myRep & rep, int iterator){
    for (int i=0; i<nPop; i++)
        if (!particle[i].dominated){
            if (canPutIntoRep(particle[i].Cost, rep)){
                Rep tmp;
                tmp.addO_origin = new double [particle[i].sizeOfAddO_origin];
                tmp.Position = new double [particle[i].sizeOfPosition];
                tmp.Cost = new double [numObjective];
                tmp.sizeOfAddO_origin = particle[i].sizeOfAddO_origin;
                tmp.iterator = iterator;
                cpyDoubleArray(tmp.addO_origin, particle[i].addO_origin, particle[i].sizeOfAddO_origin);
                cpyDoubleArray(tmp.Position, particle[i].Position, particle[i].sizeOfPosition);
                cpyDoubleArray(tmp.Cost, particle[i].Cost, numObjective);
                rep.push_back(tmp);
            }
        }
}

bool canPutIntoRep(double *cost, myRep & rep){
    bool flag=1;                                                    // if return true then can put
    myRep::iterator it=rep.begin();
    myRep::iterator tmpIt;
    while (it != rep.end()){
        bool erase=0;
        if (isDominated(cost, it->Cost)){       //      there is at least one element in the rep better than
            flag = 0;                           //  the new particle, so it is also impossible that
            break;                                              //  the new one is better than the remains in the rep
        }
        if (isDominated(it->Cost, cost)){       //  the new one is better then one of the element in rep
            tmpIt = it;
            erase =1;
        }
        it++;
        if (erase) rep.erase(tmpIt);            //  if erase it ... error will occur
    }
    return flag;
}

void updatePBest(Particle * particle){
    for (int i=0; i<nPop; i++){
        int update=0;               //  update=-1  :  pre pBest is better;  update=1: new particle is better
        if (isDominated(particle[i].Cost, particle[i].Best.Cost))  update = -1;
        if (isDominated(particle[i].Best.Cost, particle[i].Cost))  update = 1;

        if (update==0)
            if (rand() % 2 ==1)
                update=1;
        if (update == 1){
            cpyDoubleArray(particle[i].Best.Position, particle[i].Position, particle[i].sizeOfPosition);
            cpyDoubleArray(particle[i].Best.Cost, particle[i].Cost, numObjective);
        }
    }
}
